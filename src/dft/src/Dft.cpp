// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2023-2025, The OpenROAD Authors

#include "dft/Dft.hh"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <fstream>
#include <memory>
#include <optional>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <variant>
#include <vector>

#include "ClockDomain.hh"
#include "DftConfig.hh"
#include "ScanArchitect.hh"
#include "ScanArchitectConfig.hh"
#include "ScanCell.hh"
#include "ScanCellFactory.hh"
#include "ScanPin.hh"
#include "ScanReplace.hh"
#include "ScanStitch.hh"
#include "ShiftRegister.hh"
#include "db_sta/dbNetwork.hh"
#include "db_sta/dbSta.hh"
#include "odb/db.h"
#include "sta/Liberty.hh"
#include "sta/PortDirection.hh"
#include "utl/Logger.h"

namespace {
constexpr char kDefaultPartition[] = "default";

sta::LibertyCell* getLibertyCell(odb::dbInst* inst, sta::dbNetwork* db_network)
{
  if (inst == nullptr || db_network == nullptr) {
    return nullptr;
  }
  sta::Cell* master_cell = db_network->dbToSta(inst->getMaster());
  return master_cell != nullptr ? db_network->libertyCell(master_cell) : nullptr;
}

bool isClockGateInstance(odb::dbInst* inst, sta::dbNetwork* db_network)
{
  sta::LibertyCell* liberty_cell = getLibertyCell(inst, db_network);
  return liberty_cell != nullptr && liberty_cell->isClockGate();
}

bool isTristateDriverInstance(odb::dbInst* inst, sta::dbNetwork* db_network)
{
  sta::LibertyCell* liberty_cell = getLibertyCell(inst, db_network);
  if (liberty_cell == nullptr) {
    return false;
  }

  sta::LibertyCellPortIterator port_it(liberty_cell);
  while (port_it.hasNext()) {
    sta::LibertyPort* port = port_it.next();
    if (port == nullptr) {
      continue;
    }
    const sta::PortDirection* dir = port->direction();
    if (dir != nullptr && dir->isAnyOutput() && port->tristateEnable() != nullptr) {
      return true;
    }
  }
  return false;
}

void collectInstancesRec(odb::dbBlock* block, std::vector<odb::dbInst*>& out)
{
  for (odb::dbInst* inst : block->getInsts()) {
    out.push_back(inst);
  }
  for (odb::dbBlock* child : block->getChildren()) {
    collectInstancesRec(child, out);
  }
}

// Warn about clock-gated scan clocks and tri-state drivers. When scan_enable
// is available (post-stitch), also check for common clock-gate test-enable pins
// connected to the scan enable net.
void warnSpecialCells(odb::dbDatabase* db,
                      sta::dbSta* sta,
                      utl::Logger* logger,
                      const std::vector<std::unique_ptr<dft::ScanChain>>& chains,
                      bool post_stitch)
{
  if (db == nullptr || sta == nullptr || logger == nullptr) {
    return;
  }
  odb::dbChip* chip = db->getChip();
  if (chip == nullptr || chip->getBlock() == nullptr) {
    return;
  }

  odb::dbBlock* block = chip->getBlock();
  sta::dbNetwork* db_network = sta->getDbNetwork();
  if (db_network == nullptr) {
    return;
  }

  // Determine scan_enable nets (if any) from stitched chains.
  std::unordered_set<odb::dbNet*> scan_enable_nets;
  if (post_stitch) {
    for (const auto& chain : chains) {
      if (!chain) {
        continue;
      }
      const std::optional<dft::ScanDriver> se_driver = chain->getScanEnable();
      if (!se_driver.has_value()) {
        continue;
      }
      odb::dbNet* net = nullptr;
      std::visit([&](auto&& term) { net = term ? term->getNet() : nullptr; },
                 se_driver->getValue());
      if (net != nullptr) {
        scan_enable_nets.insert(net);
      }
    }
  }

  // Map clock-gate instance -> number of scan flops clocked by it.
  std::unordered_map<odb::dbInst*, int> gate_to_scan_flops;

  int total_scan_flops = 0;
  for (const auto& chain : chains) {
    if (!chain) {
      continue;
    }
    for (const auto& cell : chain->getScanCells()) {
      if (!cell) {
        continue;
      }
      odb::dbInst* inst = block->findInst(std::string(cell->getName()).c_str());
      if (inst == nullptr) {
        continue;
      }
      ++total_scan_flops;

      const std::vector<odb::dbITerm*> clocks = dft::utils::GetClockPin(inst);
      if (clocks.empty() || clocks.front() == nullptr) {
        continue;
      }
      odb::dbNet* clk_net = clocks.front()->getNet();
      if (clk_net == nullptr) {
        continue;
      }

      for (odb::dbITerm* iterm : clk_net->getITerms()) {
        if (iterm == nullptr || !iterm->isOutputSignal()) {
          continue;
        }
        odb::dbInst* driver_inst = iterm->getInst();
        if (driver_inst == nullptr) {
          continue;
        }
        if (!isClockGateInstance(driver_inst, db_network)) {
          continue;
        }
        gate_to_scan_flops[driver_inst] += 1;
      }
    }
  }

  if (!gate_to_scan_flops.empty()) {
    logger->warn(
        utl::DFT,
        256,
        "Found {} clock gate instance(s) driving clock nets used by {} scan "
        "cell(s). Ensure clock gates are enabled in scan mode.",
        gate_to_scan_flops.size(),
        total_scan_flops);

    if (post_stitch) {
      int enabled_gates = 0;
      int gates_with_te_pin = 0;

      auto isTestEnablePinName = [](const char* name) -> bool {
        return name != nullptr
               && (strcasecmp(name, "SE") == 0 || strcasecmp(name, "TE") == 0
                   || strcasecmp(name, "SCAN_EN") == 0
                   || strcasecmp(name, "SCAN_ENABLE") == 0
                   || strcasecmp(name, "SCANENABLE") == 0);
      };

      for (const auto& [gate, _] : gate_to_scan_flops) {
        if (gate == nullptr) {
          continue;
        }

        bool has_te = false;
        bool te_connected = false;
        for (odb::dbITerm* iterm : gate->getITerms()) {
          if (iterm == nullptr) {
            continue;
          }
          if (!isTestEnablePinName(iterm->getMTerm()->getConstName())) {
            continue;
          }
          has_te = true;
          odb::dbNet* net = iterm->getNet();
          if (net != nullptr && scan_enable_nets.find(net) != scan_enable_nets.end()) {
            te_connected = true;
            break;
          }
        }
        if (has_te) {
          ++gates_with_te_pin;
        }
        if (te_connected) {
          ++enabled_gates;
        }
      }

      if (!scan_enable_nets.empty()) {
        logger->info(
            utl::DFT,
            257,
            "Clock gate scan enable check: {} of {} clock gate(s) have a "
            "test-enable pin (SE/TE/SCAN_ENABLE) connected to the scan enable net.",
            enabled_gates,
            gate_to_scan_flops.size());
      } else {
        logger->warn(utl::DFT,
                     258,
                     "Clock gate scan enable check skipped: scan enable net not found.");
      }

      if (gates_with_te_pin < static_cast<int>(gate_to_scan_flops.size())) {
        logger->warn(
            utl::DFT,
            259,
            "Some clock gate cells do not expose a recognized test-enable pin "
            "(SE/TE/SCAN_ENABLE). Scan-mode gate control may require synthesis/DFT insertion.",
            gates_with_te_pin);
      }
    }
  }

  if (post_stitch) {
    // Tri-state drivers are a common scan-mode hazard. Report presence so users
    // can ensure they are disabled/controlled in test mode.
    std::vector<odb::dbInst*> all_insts;
    all_insts.reserve(block->getInsts().size());
    collectInstancesRec(block, all_insts);

    int tristate_count = 0;
    for (odb::dbInst* inst : all_insts) {
      if (inst == nullptr || inst->isDoNotTouch()) {
        continue;
      }
      if (isTristateDriverInstance(inst, db_network)) {
        ++tristate_count;
      }
    }
    if (tristate_count > 0) {
      logger->warn(
          utl::DFT,
          260,
          "Detected {} tri-state driver instance(s). Ensure tri-state drivers are "
          "disabled or controlled during scan to avoid X propagation/contending.",
          tristate_count);
    }
  }
}

struct PowerDomainCrossingExample
{
  std::string chain_name;
  std::string from_inst;
  std::string to_inst;
  std::string from_domain;
  std::string to_domain;
  float from_voltage = 0.0F;
  float to_voltage = 0.0F;
};

void collectGroupInstsRec(odb::dbGroup* group,
                          odb::dbPowerDomain* domain,
                          std::unordered_map<odb::dbInst*, odb::dbPowerDomain*>& out)
{
  if (group == nullptr) {
    return;
  }
  for (odb::dbInst* inst : group->getInsts()) {
    out[inst] = domain;
  }
  for (odb::dbGroup* child : group->getGroups()) {
    collectGroupInstsRec(child, domain, out);
  }
}

void warnPowerDomainCrossings(odb::dbDatabase* db,
                              utl::Logger* logger,
                              const std::vector<std::unique_ptr<dft::ScanChain>>& chains)
{
  if (db == nullptr || logger == nullptr) {
    return;
  }
  odb::dbChip* chip = db->getChip();
  if (chip == nullptr || chip->getBlock() == nullptr) {
    return;
  }
  odb::dbBlock* block = chip->getBlock();

  odb::dbSet<odb::dbPowerDomain> power_domains = block->getPowerDomains();
  if (power_domains.empty()) {
    return;
  }

  std::unordered_map<odb::dbInst*, odb::dbPowerDomain*> inst_to_domain;
  for (odb::dbPowerDomain* pd : power_domains) {
    odb::dbGroup* group = pd->getGroup();
    if (group == nullptr) {
      continue;
    }
    collectGroupInstsRec(group, pd, inst_to_domain);
  }

  if (inst_to_domain.empty()) {
    return;
  }

  int voltage_crossings = 0;
  int switched_crossings = 0;
  int unknown_crossings = 0;
  std::vector<PowerDomainCrossingExample> examples;
  constexpr int kMaxExamples = 8;

  auto getDomainName = [](odb::dbPowerDomain* pd) -> std::string {
    return pd != nullptr ? pd->getName() : "<none>";
  };

  auto domainIsSwitched = [](odb::dbPowerDomain* pd) -> bool {
    return pd != nullptr && !pd->getPowerSwitches().empty();
  };

  for (const auto& chain : chains) {
    if (!chain) {
      continue;
    }
    const auto& cells = chain->getScanCells();
    for (std::size_t i = 1; i < cells.size(); ++i) {
      const auto& prev_cell = cells[i - 1];
      const auto& next_cell = cells[i];
      if (!prev_cell || !next_cell) {
        continue;
      }

      odb::dbInst* from_inst
          = block->findInst(std::string(prev_cell->getName()).c_str());
      odb::dbInst* to_inst
          = block->findInst(std::string(next_cell->getName()).c_str());
      if (from_inst == nullptr || to_inst == nullptr) {
        continue;
      }

      odb::dbPowerDomain* from_pd
          = inst_to_domain.count(from_inst) ? inst_to_domain[from_inst] : nullptr;
      odb::dbPowerDomain* to_pd
          = inst_to_domain.count(to_inst) ? inst_to_domain[to_inst] : nullptr;

      if (from_pd == to_pd) {
        continue;
      }

      if (from_pd == nullptr || to_pd == nullptr) {
        ++unknown_crossings;
        continue;
      }

      const float v_from = from_pd->getVoltage();
      const float v_to = to_pd->getVoltage();
      if (std::abs(v_from - v_to) > 1e-6F) {
        ++voltage_crossings;
        if (static_cast<int>(examples.size()) < kMaxExamples) {
          examples.push_back({.chain_name = chain->getName(),
                              .from_inst = from_inst->getName(),
                              .to_inst = to_inst->getName(),
                              .from_domain = getDomainName(from_pd),
                              .to_domain = getDomainName(to_pd),
                              .from_voltage = v_from,
                              .to_voltage = v_to});
        }
        continue;
      }

      if (domainIsSwitched(from_pd) || domainIsSwitched(to_pd)) {
        ++switched_crossings;
        if (static_cast<int>(examples.size()) < kMaxExamples) {
          examples.push_back({.chain_name = chain->getName(),
                              .from_inst = from_inst->getName(),
                              .to_inst = to_inst->getName(),
                              .from_domain = getDomainName(from_pd),
                              .to_domain = getDomainName(to_pd),
                              .from_voltage = v_from,
                              .to_voltage = v_to});
        }
      }
    }
  }

  if (unknown_crossings > 0) {
    logger->warn(utl::DFT,
                 261,
                 "Scan chains include {} crossing(s) between different power "
                 "domains, but power-domain assignments are incomplete; cannot "
                 "validate level shifter/isolation requirements.",
                 unknown_crossings);
  }

  if (voltage_crossings > 0) {
    logger->warn(utl::DFT,
                 262,
                 "Scan chains include {} crossing(s) between power domains at "
                 "different voltages; level shifters may be required.",
                 voltage_crossings);
  }

  if (switched_crossings > 0) {
    logger->warn(utl::DFT,
                 263,
                 "Scan chains include {} crossing(s) where at least one power "
                 "domain is switched; isolation may be required.",
                 switched_crossings);
  }

  if (!examples.empty()) {
    const auto& ex = examples.front();
    logger->report("Example power-domain crossing: chain '{}' edge '{}' -> '{}' "
                   "({} {:.3f}V -> {} {:.3f}V)",
                   ex.chain_name,
                   ex.from_inst,
                   ex.to_inst,
                   ex.from_domain,
                   ex.from_voltage,
                   ex.to_domain,
                   ex.to_voltage);
  }
}
}  // namespace

namespace dft {

Dft::Dft(odb::dbDatabase* db, sta::dbSta* sta, utl::Logger* logger)
    : db_(db),
      sta_(sta),
      logger_(logger),
      dft_config_(std::make_unique<DftConfig>())
{
}

Dft::~Dft() = default;

void Dft::reset()
{
  scan_replace_.reset();
  need_to_run_pre_dft_ = true;
  auto_exclusions_cache_valid_ = false;
  cached_auto_exclude_shift_registers_ = false;
  cached_shift_register_min_length_ = 0;
}

void Dft::pre_dft()
{
  scan_replace_
      = std::make_unique<ScanReplace>(db_,
                                      sta_,
                                      logger_,
                                      &dft_config_->getScanArchitectConfig());
  scan_replace_->collectScanCellAvailable();

  // This should always be at the end
  need_to_run_pre_dft_ = false;
}

void Dft::reportDftPlan(bool verbose)
{
  if (need_to_run_pre_dft_) {
    pre_dft();
  }

  std::vector<std::unique_ptr<ScanChain>> scan_chains = scanArchitect();

  logger_->report("***************************");
  logger_->report("Report DFT Plan");
  logger_->report("Number of chains: {:d}", scan_chains.size());
  logger_->report("Clock domain: {:s}",
                  ScanArchitectConfig::ClockMixingName(
                      dft_config_->getScanArchitectConfig().getClockMixing()));
  logger_->report("***************************\n");
  for (const auto& scan_chain : scan_chains) {
    scan_chain->report(logger_, verbose);
  }
  warnPowerDomainCrossings(db_, logger_, scan_chains);
  warnSpecialCells(db_, sta_, logger_, scan_chains, /*post_stitch=*/false);
  logger_->report("");
}

void Dft::scanReplace()
{
  if (need_to_run_pre_dft_) {
    pre_dft();
  }
  applyAutoExclusions();
  scan_replace_->scanReplace();
}

void Dft::executeDftPlan()
{
  if (need_to_run_pre_dft_) {
    pre_dft();
  }
  std::vector<std::unique_ptr<ScanChain>> scan_chains = scanArchitect();

  warnPowerDomainCrossings(db_, logger_, scan_chains);
  warnSpecialCells(db_, sta_, logger_, scan_chains, /*post_stitch=*/false);

  ScanStitch stitch(db_,
                    logger_,
                    dft_config_->getScanArchitectConfig(),
                    dft_config_->getScanStitchConfig());
  stitch.Stitch(scan_chains);

  // Post-stitch: scan_enable nets exist; re-run checks with more context.
  warnSpecialCells(db_, sta_, logger_, scan_chains, /*post_stitch=*/true);

  // Write scan chains to odb
  odb::dbBlock* db_block = db_->getChip()->getBlock();
  odb::dbDft* db_dft = db_block->getDft();

  for (const auto& chain : scan_chains) {
    odb::dbScanChain* db_sc = odb::dbScanChain::create(db_dft);
    db_sc->setName(chain->getName());
    odb::dbScanPartition* db_part = odb::dbScanPartition::create(db_sc);
    db_part->setName(kDefaultPartition);
    odb::dbScanList* db_scanlist = odb::dbScanList::create(db_part);

    for (const auto& scan_cell : chain->getScanCells()) {
      std::string inst_name(scan_cell->getName());
      odb::dbInst* db_inst = db_block->findInst(inst_name.c_str());
      odb::dbScanInst* db_scaninst = db_scanlist->add(db_inst);
      db_scaninst->setBits(scan_cell->getBits());
      ScanLoad scan_enable = scan_cell->getScanEnable();
      std::visit([&](auto&& pin) { db_scaninst->setScanEnable(pin); },
                 scan_enable.getValue());
      auto scan_in_term = scan_cell->getScanIn().getValue();
      auto scan_out_term = scan_cell->getScanOut().getValue();
      db_scaninst->setAccessPins(
          {.scan_in = scan_in_term, .scan_out = scan_out_term});

      const ClockDomain& clock_domain = scan_cell->getClockDomain();
      db_scaninst->setScanClock(clock_domain.getClockName());
      switch (clock_domain.getClockEdge()) {
        case ClockEdge::Rising:
          db_scaninst->setClockEdge(odb::dbScanInst::ClockEdge::Rising);
          break;
        case ClockEdge::Falling:
          db_scaninst->setClockEdge(odb::dbScanInst::ClockEdge::Falling);
          break;
      }
    }

    std::optional<ScanDriver> sc_enable_driver = chain->getScanEnable();
    std::optional<ScanDriver> sc_in_driver = chain->getScanIn();
    std::optional<ScanLoad> sc_out_load = chain->getScanOut();

    if (sc_enable_driver.has_value()) {
      std::visit(
          [&](auto&& sc_enable_term) { db_sc->setScanEnable(sc_enable_term); },
          sc_enable_driver.value().getValue());
    }
    if (sc_in_driver.has_value()) {
      std::visit([&](auto&& sc_in_term) { db_sc->setScanIn(sc_in_term); },
                 sc_in_driver.value().getValue());
    }
    if (sc_out_load.has_value()) {
      std::visit([&](auto&& sc_out_term) { db_sc->setScanOut(sc_out_term); },
                 sc_out_load.value().getValue());
    }
  }
}

int Dft::bufferScanEnable(const std::string& buffer_cell,
                          int max_fanout,
                          int max_levels)
{
  if (buffer_cell.empty()) {
    logger_->warn(utl::DFT,
                  280,
                  "Scan enable buffering skipped: buffer cell not specified.");
    return 0;
  }

  if (max_fanout <= 0) {
    logger_->error(
        utl::DFT, 281, "Expected max_fanout > 0 (got {}).", max_fanout);
  }
  if (max_levels <= 0) {
    logger_->error(
        utl::DFT, 282, "Expected max_levels > 0 (got {}).", max_levels);
  }

  if (db_ == nullptr) {
    logger_->error(utl::DFT, 283, "No database available.");
  }
  odb::dbChip* chip = db_->getChip();
  if (chip == nullptr || chip->getBlock() == nullptr) {
    logger_->error(utl::DFT, 284, "No design block found.");
  }
  odb::dbBlock* block = chip->getBlock();

  odb::dbNet* scan_enable_net = nullptr;
  odb::dbDft* db_dft = block->getDft();
  if (db_dft != nullptr) {
    for (odb::dbScanChain* chain : db_dft->getScanChains()) {
      if (chain == nullptr) {
        continue;
      }
      odb::dbNet* net = std::visit(
          [&](auto&& term) { return term ? term->getNet() : nullptr; },
          chain->getScanEnable());
      if (net != nullptr) {
        scan_enable_net = net;
        break;
      }
    }
  }

  if (scan_enable_net == nullptr) {
    logger_->warn(utl::DFT,
                  285,
                  "Scan enable buffering skipped: scan enable net not found "
                  "(run execute_dft_plan first).");
    return 0;
  }

  scan_enable_net->setSigType(odb::dbSigType::SCAN);

  odb::dbMaster* master = db_->findMaster(buffer_cell.c_str());
  if (master == nullptr) {
    logger_->error(utl::DFT,
                   286,
                   "Scan enable buffer cell master '{}' not found in the database.",
                   buffer_cell);
  }

  std::string in_pin;
  std::string out_pin;
  for (odb::dbMTerm* mterm : master->getMTerms()) {
    if (mterm == nullptr || mterm->getSigType() != odb::dbSigType::SIGNAL) {
      continue;
    }
    const odb::dbIoType io = mterm->getIoType();
    if (in_pin.empty() && io == odb::dbIoType::INPUT) {
      in_pin = mterm->getName();
    } else if (out_pin.empty() && io == odb::dbIoType::OUTPUT) {
      out_pin = mterm->getName();
    }
  }
  if (in_pin.empty() || out_pin.empty()) {
    logger_->error(
        utl::DFT,
        287,
        "Scan enable buffer cell '{}' must have at least one SIGNAL INPUT and "
        "one SIGNAL OUTPUT pin.",
        buffer_cell);
  }

  struct LoadLoc
  {
    odb::dbITerm* iterm = nullptr;
    int x = 0;
    int y = 0;
    std::string name;
  };

  int total_inserted = 0;

  for (int level = 0; level < max_levels; ++level) {
    std::vector<LoadLoc> loads;
    for (odb::dbITerm* iterm : scan_enable_net->getITerms()) {
      if (iterm == nullptr || iterm->getIoType() != odb::dbIoType::INPUT) {
        continue;
      }
      int x = 0;
      int y = 0;
      if (!iterm->getAvgXY(&x, &y)) {
        if (odb::dbInst* inst = iterm->getInst()) {
          inst->getLocation(x, y);
        }
      }

      std::string name;
      if (odb::dbInst* inst = iterm->getInst()) {
        name = inst->getName();
      }
      name += "/";
      if (odb::dbMTerm* mterm = iterm->getMTerm()) {
        name += mterm->getName();
      }
      loads.push_back({iterm, x, y, std::move(name)});
    }

    const int fanout = static_cast<int>(loads.size());
    if (fanout <= max_fanout) {
      if (level == 0) {
        logger_->info(utl::DFT,
                      288,
                      "Scan enable fanout={} (<= {}); no buffering needed.",
                      fanout,
                      max_fanout);
      } else {
        logger_->info(utl::DFT,
                      289,
                      "Scan enable buffering complete at level {} (fanout={}).",
                      level,
                      fanout);
      }
      break;
    }

    std::sort(loads.begin(),
              loads.end(),
              [](const LoadLoc& a, const LoadLoc& b) {
                if (a.x != b.x) {
                  return a.x < b.x;
                }
                if (a.y != b.y) {
                  return a.y < b.y;
                }
                return a.name < b.name;
              });

    for (int start = 0; start < fanout; start += max_fanout) {
      const int end = std::min(start + max_fanout, fanout);
      const int count = end - start;
      if (count <= 0) {
        continue;
      }

      long long sum_x = 0;
      long long sum_y = 0;
      for (int i = start; i < end; ++i) {
        sum_x += loads[i].x;
        sum_y += loads[i].y;
      }
      const int cx = static_cast<int>(
          std::llround(static_cast<double>(sum_x) / static_cast<double>(count)));
      const int cy = static_cast<int>(
          std::llround(static_cast<double>(sum_y) / static_cast<double>(count)));

      const int group_id = start / max_fanout;
      const std::string inst_name = "dft_scan_enable_buf_" + std::to_string(level)
                                    + "_" + std::to_string(group_id);
      odb::dbInst* buf_inst = block->findInst(inst_name.c_str());
      if (buf_inst != nullptr) {
        odb::dbInst::destroy(buf_inst);
        buf_inst = nullptr;
      }
      buf_inst = odb::dbInst::create(block, master, inst_name.c_str());
      if (buf_inst == nullptr) {
        logger_->error(utl::DFT,
                       290,
                       "Failed to create scan enable buffer instance '{}'",
                       inst_name);
      }
      buf_inst->setLocation(cx, cy);
      buf_inst->setPlacementStatus(odb::dbPlacementStatus::PLACED);

      odb::dbITerm* din = buf_inst->findITerm(in_pin.c_str());
      odb::dbITerm* dout = buf_inst->findITerm(out_pin.c_str());
      if (din == nullptr || dout == nullptr) {
        logger_->error(utl::DFT,
                       291,
                       "Scan enable buffer instance '{}' missing expected pin(s) "
                       "('{}' input, '{}' output).",
                       inst_name,
                       in_pin,
                       out_pin);
      }

      din->connect(scan_enable_net);

      const std::string net_name = "dft_scan_enable_net_" + std::to_string(level)
                                   + "_" + std::to_string(group_id);
      odb::dbNet* out_net = block->findNet(net_name.c_str());
      if (out_net == nullptr) {
        out_net = odb::dbNet::create(block, net_name.c_str());
        if (out_net == nullptr) {
          logger_->error(utl::DFT,
                         292,
                         "Failed to create scan enable buffer net '{}'",
                         net_name);
        }
      }
      out_net->setSigType(odb::dbSigType::SCAN);

      dout->connect(out_net);

      for (int i = start; i < end; ++i) {
        odb::dbITerm* iterm = loads[i].iterm;
        if (iterm == nullptr) {
          continue;
        }
        iterm->disconnect();
        iterm->connect(out_net);
      }

      ++total_inserted;
    }
  }

  if (total_inserted > 0) {
    logger_->info(utl::DFT,
                  293,
                  "Inserted {} buffer(s) for scan_enable (max_fanout={}, levels={}).",
                  total_inserted,
                  max_fanout,
                  max_levels);
  }

  return total_inserted;
}

void Dft::writeScandef(const std::string& path) const
{
  if (db_ == nullptr) {
    logger_->error(utl::DFT, 264, "No database available.");
  }

  odb::dbChip* chip = db_->getChip();
  if (chip == nullptr || chip->getBlock() == nullptr) {
    logger_->error(utl::DFT, 265, "No design block found.");
  }

  odb::dbBlock* block = chip->getBlock();
  odb::dbDft* dft = block->getDft();
  odb::dbSet<odb::dbScanChain> scan_chains = dft->getScanChains();
  if (scan_chains.empty()) {
    logger_->error(utl::DFT,
                   266,
                   "No scan chains found in the database (run execute_dft_plan first).");
  }

  std::ofstream out(path);
  if (!out) {
    logger_->error(utl::DFT, 267, "Couldn't open '{}' for writing.", path);
  }

  auto writeEndpoint = [&](const char* keyword,
                           const std::variant<odb::dbBTerm*, odb::dbITerm*>& pin,
                           bool add_semicolon) {
    out << "+ " << keyword << " ";
    std::visit(
        [&](auto&& p) {
          if (p == nullptr) {
            out << "PIN ";
            return;
          }
          if constexpr (std::is_same_v<std::decay_t<decltype(p)>, odb::dbBTerm*>) {
            out << "PIN " << p->getName();
          } else {
            out << p->getInst()->getName() << " " << p->getMTerm()->getName();
          }
        },
        pin);
    if (add_semicolon) {
      out << " ;";
    }
    out << "\n";
  };

  auto accessPinName = [](const std::variant<odb::dbBTerm*, odb::dbITerm*>& pin) {
    return std::visit(
        [](auto&& p) -> std::string {
          if (p == nullptr) {
            return "";
          }
          if constexpr (std::is_same_v<std::decay_t<decltype(p)>, odb::dbBTerm*>) {
            return std::string(p->getName());
          } else {
            return std::string(p->getMTerm()->getName());
          }
        },
        pin);
  };

  // Minimal SCANDEF/DEF-style wrapper.
  out << "VERSION 5.8 ;\n";
  out << "DIVIDERCHAR \"/\" ;\n";
  out << "BUSBITCHARS \"[]\" ;\n";
  out << "DESIGN " << block->getName() << " ;\n\n";

  out << "SCANCHAINS " << scan_chains.size() << " ;\n\n";

  for (odb::dbScanChain* scan_chain : dft->getScanChains()) {
    odb::dbSet<odb::dbScanPartition> scan_partitions
        = scan_chain->getScanPartitions();
    int chain_suffix = 0;
    for (odb::dbScanPartition* scan_partition : scan_partitions) {
      bool already_printed_floating = false;
      bool already_printed_ordered = false;

      const std::string chain_name
          = scan_partitions.size() == 1
                ? scan_chain->getName()
                : scan_chain->getName() + "_" + std::to_string(chain_suffix);

      out << "- " << chain_name << "\n";
      writeEndpoint("START", scan_chain->getScanIn(), /*add_semicolon=*/false);

      for (odb::dbScanList* scan_list : scan_partition->getScanLists()) {
        odb::dbSet<odb::dbScanInst> scan_insts = scan_list->getScanInsts();

        if (scan_insts.size() == 1 && !already_printed_floating) {
          out << "+ FLOATING\n";
          already_printed_floating = true;
          already_printed_ordered = false;
        } else if (scan_insts.size() > 1 && !already_printed_ordered) {
          out << "+ ORDERED\n";
          already_printed_floating = false;
          already_printed_ordered = true;
        }

        for (odb::dbScanInst* scan_inst : scan_insts) {
          odb::dbScanInst::AccessPins pins = scan_inst->getAccessPins();
          out << "  " << scan_inst->getInst()->getName() << " ( IN "
              << accessPinName(pins.scan_in) << " ) ( OUT "
              << accessPinName(pins.scan_out) << " )\n";
        }
      }

      out << "+ PARTITION " << scan_partition->getName() << "\n";
      writeEndpoint("STOP", scan_chain->getScanOut(), /*add_semicolon=*/true);
      out << "\n";
      ++chain_suffix;
    }
  }

  out << "END SCANCHAINS\n\n";
  out << "END DESIGN\n";

  logger_->info(utl::DFT, 268, "Wrote scan chains to '{}'", path);
}

DftConfig* Dft::getMutableDftConfig()
{
  return dft_config_.get();
}

const DftConfig& Dft::getDftConfig() const
{
  return *dft_config_;
}

void Dft::reportDftConfig() const
{
  logger_->report("DFT Config:");
  dft_config_->report(logger_);
}

std::vector<std::unique_ptr<ScanChain>> Dft::scanArchitect()
{
  applyAutoExclusions();
  std::vector<std::unique_ptr<ScanCell>> scan_cells
      = CollectScanCells(db_, sta_, dft_config_->getScanArchitectConfig(), logger_);

  // Scan Architect
  std::unique_ptr<ScanCellsBucket> scan_cells_bucket
      = std::make_unique<ScanCellsBucket>(logger_);
  scan_cells_bucket->init(dft_config_->getScanArchitectConfig(), scan_cells);

  std::unique_ptr<ScanArchitect> scan_architect
      = ScanArchitect::ConstructScanScanArchitect(
          dft_config_->getScanArchitectConfig(),
          std::move(scan_cells_bucket),
          logger_);
  scan_architect->init();
  scan_architect->architect();

  return scan_architect->getScanChains();
}

void Dft::applyAutoExclusions()
{
  ScanArchitectConfig* config
      = dft_config_->getMutableScanArchitectConfig();
  const bool enable_shift_regs = config->getAutoExcludeShiftRegisters();
  const int min_len = config->getShiftRegisterMinLength();

  const bool cache_hit
      = auto_exclusions_cache_valid_
        && cached_auto_exclude_shift_registers_ == enable_shift_regs
        && cached_shift_register_min_length_ == min_len;
  if (cache_hit) {
    return;
  }

  config->clearAutoExcludedInstances();
  if (enable_shift_regs) {
    auto detected = dft::utils::DetectShiftRegisters(db_, sta_, *config, logger_);
    config->setAutoExcludedInstances(std::move(detected.instances));
  }

  auto_exclusions_cache_valid_ = true;
  cached_auto_exclude_shift_registers_ = enable_shift_regs;
  cached_shift_register_min_length_ = min_len;
}

void Dft::scanOpt()
{
  if (need_to_run_pre_dft_) {
    pre_dft();
  }

  // Re-run scan planning using the latest placement, then re-stitch scan
  // connections. This updates scan ordering without re-running scan_replace.
  std::vector<std::unique_ptr<ScanChain>> scan_chains = scanArchitect();
  if (scan_chains.empty()) {
    logger_->warn(utl::DFT,
                  14,
                  "Scan Opt found no scan chains (run `scan_replace` first)");
    return;
  }

  ScanStitch stitch(db_,
                    logger_,
                    dft_config_->getScanArchitectConfig(),
                    dft_config_->getScanStitchConfig());
  stitch.Stitch(scan_chains);

  logger_->info(utl::DFT,
                15,
                "Scan Opt re-stitched {:d} scan chain(s)",
                scan_chains.size());
}

}  // namespace dft
