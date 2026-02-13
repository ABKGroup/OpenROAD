// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2023-2025, The OpenROAD Authors

#include "ScanStitch.hh"

#include <algorithm>
#include <cstddef>
#include <deque>
#include <functional>
#include <iterator>
#include <memory>
#include <optional>
#include <string>
#include <string_view>
#include <type_traits>
#include <unordered_set>
#include <utility>
#include <vector>

#include "ClockDomain.hh"
#include "boost/algorithm/string.hpp"
#include "odb/db.h"
#include "odb/dbTypes.h"
#include "utl/Logger.h"

namespace {
constexpr std::string_view kScanEnable = "scan-enable";
constexpr std::string_view kScanIn = "scan-in";
constexpr std::string_view kScanOut = "scan-out";
constexpr size_t kEnableNumber = 0;

constexpr std::string_view kLockupInstPrefix = "dft_lockup_";
constexpr std::string_view kLockupNetPrefix = "dft_lockup_net_";

constexpr std::string_view kScanBufferInstPrefix = "dft_scan_buf_";
constexpr std::string_view kScanBufferNetPrefix = "dft_scan_buf_net_";

odb::dbTechLayer* DefaultPortLayer(odb::dbBlock* block)
{
  if (block == nullptr) {
    return nullptr;
  }
  odb::dbTech* tech = block->getTech();
  if (tech == nullptr) {
    return nullptr;
  }
  if (odb::dbTechLayer* layer = tech->findRoutingLayer(1)) {
    return layer;
  }
  // Fall back to the first routing layer in case routing levels aren't set.
  for (odb::dbTechLayer* layer : tech->getLayers()) {
    if (layer && layer->getType() == odb::dbTechLayerType::ROUTING) {
      return layer;
    }
  }
  return nullptr;
}

void EnsureBTermHasLocation(odb::dbBlock* block,
                            odb::dbBTerm* term,
                            const odb::Point& loc,
                            utl::Logger* logger)
{
  if (block == nullptr || term == nullptr) {
    return;
  }

  int x = 0;
  int y = 0;
  if (term->getFirstPinLocation(x, y) && x == loc.x() && y == loc.y()) {
    return;
  }

  odb::dbTechLayer* layer = DefaultPortLayer(block);
  if (layer == nullptr) {
    if (logger) {
      logger->warn(utl::DFT,
                   216,
                   "Can't set pin location for port '{}' (no routing layer)",
                   term->getName());
    }
    return;
  }

  // Override any previous pin location. This is intentional when the user
  // specifies explicit (x,y) endpoints for exploration: the port's placement is
  // used as the BeginPort/EndPort location.
  std::vector<odb::dbBPin*> pins;
  for (odb::dbBPin* pin : term->getBPins()) {
    pins.push_back(pin);
  }
  for (odb::dbBPin* pin : pins) {
    odb::dbBPin::destroy(pin);
  }

  odb::dbBPin* pin = odb::dbBPin::create(term);
  pin->setPlacementStatus(odb::dbPlacementStatus::PLACED);
  odb::dbBox::create(pin, layer, loc.x(), loc.y(), loc.x(), loc.y());
}

void RemoveExistingLockups(odb::dbBlock* block)
{
  if (block == nullptr) {
    return;
  }

  std::vector<odb::dbInst*> to_remove;
  for (odb::dbInst* inst : block->getInsts()) {
    if (inst == nullptr) {
      continue;
    }
    const char* name = inst->getConstName();
    if (name == nullptr) {
      continue;
    }
    if (!boost::algorithm::starts_with(name, kLockupInstPrefix)) {
      continue;
    }
    to_remove.push_back(inst);
  }

  for (odb::dbInst* inst : to_remove) {
    odb::dbInst::destroy(inst);
  }
}

void RemoveExistingScanBuffers(odb::dbBlock* block)
{
  if (block == nullptr) {
    return;
  }

  std::vector<odb::dbInst*> to_remove;
  for (odb::dbInst* inst : block->getInsts()) {
    if (inst == nullptr) {
      continue;
    }
    const char* name = inst->getConstName();
    if (name == nullptr) {
      continue;
    }
    if (!boost::algorithm::starts_with(name, kScanBufferInstPrefix)) {
      continue;
    }
    to_remove.push_back(inst);
  }

  for (odb::dbInst* inst : to_remove) {
    odb::dbInst::destroy(inst);
  }
}

odb::dbNet* EnsureDriverNet(odb::dbBlock* block,
                            const dft::ScanDriver& driver,
                            const std::string& net_name)
{
  if (odb::dbNet* net = driver.getNet()) {
    return net;
  }
  odb::dbNet* net = block->findNet(net_name.c_str());
  if (net == nullptr) {
    net = odb::dbNet::create(block, net_name.c_str());
    if (net == nullptr) {
      return nullptr;
    }
    net->setSigType(odb::dbSigType::SCAN);
  }
  std::visit(
      [&](auto&& term) {
        if (term != nullptr) {
          term->connect(net);
        }
      },
      driver.getValue());
  return net;
}

bool IsTimingCritical(const dft::ScanCell& cell, double critical_slack)
{
  if (!cell.hasTimingSlacks()) {
    return false;
  }

  const float setup = cell.getSetupSlack();
  const float hold = cell.getHoldSlack();

  if (critical_slack <= 0.0) {
    return (setup < 0.0F) || (hold < 0.0F);
  }

  return (static_cast<double>(setup) < critical_slack)
         || (static_cast<double>(hold) < critical_slack);
}

odb::dbNet* FindClockNet(odb::dbBlock* block,
                         const dft::ScanCell& scan_cell,
                         utl::Logger* logger)
{
  odb::dbInst* inst = block->findInst(std::string(scan_cell.getName()).c_str());
  if (inst == nullptr) {
    logger->error(utl::DFT,
                  104,
                  "Scan stitch can't find instance '{}' for clock lookup",
                  scan_cell.getName());
  }
  std::vector<odb::dbITerm*> clocks = dft::utils::GetClockPin(inst);
  if (clocks.empty()) {
    logger->error(utl::DFT,
                  105,
                  "Scan stitch can't find a clock pin on instance '{}'",
                  scan_cell.getName());
  }
  odb::dbNet* net = clocks.front()->getNet();
  if (net == nullptr) {
    logger->error(utl::DFT,
                  106,
                  "Scan stitch clock pin on instance '{}' has no net",
                  scan_cell.getName());
  }
  return net;
}

void InsertTimingBufferBetween(odb::dbDatabase* db,
                               odb::dbBlock* block,
                               utl::Logger* logger,
                               const dft::ScanStitchConfig& config,
                               const dft::ScanCell& prev_cell,
                               const dft::ScanCell& next_cell,
                               size_t chain_ordinal,
                               size_t link_idx)
{
  const std::string_view buffer_cell = config.getTimingBufferCell();
  if (buffer_cell.empty()) {
    return;
  }

  odb::dbMaster* master = db->findMaster(std::string(buffer_cell).c_str());
  if (master == nullptr) {
    logger->error(utl::DFT,
                  107,
                  "Timing buffer cell master '{}' not found in the database",
                  buffer_cell);
  }

  const std::string inst_name
      = fmt::format("{}{}_{}", kScanBufferInstPrefix, chain_ordinal, link_idx);
  odb::dbInst* buf_inst = block->findInst(inst_name.c_str());
  if (buf_inst == nullptr) {
    buf_inst = odb::dbInst::create(block, master, inst_name.c_str());
    if (buf_inst == nullptr) {
      logger->error(utl::DFT,
                    108,
                    "Failed to create timing buffer instance '{}'",
                    inst_name);
    }

    if (prev_cell.isPlaced() || next_cell.isPlaced()) {
      const odb::Point a_fallback = prev_cell.getOrigin();
      const odb::Point b_fallback = next_cell.getOrigin();
      const odb::Point a = prev_cell.getScanOut().getLocation(a_fallback);
      const odb::Point b = next_cell.getScanIn().getLocation(b_fallback);
      const int x = (prev_cell.isPlaced() && next_cell.isPlaced()) ? (a.x() + b.x()) / 2
                                                                   : (next_cell.isPlaced() ? b.x() : a.x());
      const int y = (prev_cell.isPlaced() && next_cell.isPlaced()) ? (a.y() + b.y()) / 2
                                                                   : (next_cell.isPlaced() ? b.y() : a.y());
      buf_inst->setLocation(x, y);
      buf_inst->setPlacementStatus(odb::dbPlacementStatus::PLACED);
    }
  } else if (buf_inst->getMaster() != master) {
    buf_inst->swapMaster(master);
  }

  odb::dbITerm* din = buf_inst->findITerm(
      std::string(config.getTimingBufferInPin()).c_str());
  odb::dbITerm* dout = buf_inst->findITerm(
      std::string(config.getTimingBufferOutPin()).c_str());

  if (din == nullptr) {
    logger->error(utl::DFT,
                  109,
                  "Timing buffer instance '{}' missing input pin '{}'",
                  inst_name,
                  config.getTimingBufferInPin());
  }
  if (dout == nullptr) {
    logger->error(utl::DFT,
                  110,
                  "Timing buffer instance '{}' missing output pin '{}'",
                  inst_name,
                  config.getTimingBufferOutPin());
  }

  // prev.Q -> buf.A (preserve prev functional connections)
  const dft::ScanDriver prev_scan_out = prev_cell.getScanOut();
  odb::dbNet* src_net = EnsureDriverNet(
      block,
      prev_scan_out,
      fmt::format("dft_scan_src_{}_{}", chain_ordinal, link_idx));
  if (src_net == nullptr) {
    logger->error(utl::DFT,
                  111,
                  "Failed to create scan source net for timing buffer between "
                  "'{}' and '{}'",
                  prev_cell.getName(),
                  next_cell.getName());
  }
  din->connect(src_net);

  // buf.X -> next.SI (new scan net)
  const std::string net_name
      = fmt::format("{}{}_{}", kScanBufferNetPrefix, chain_ordinal, link_idx);
  odb::dbNet* out_net = block->findNet(net_name.c_str());
  if (out_net == nullptr) {
    out_net = odb::dbNet::create(block, net_name.c_str());
    if (out_net == nullptr) {
      logger->error(utl::DFT,
                    112,
                    "Failed to create timing buffer net '{}'",
                    net_name);
    }
    out_net->setSigType(odb::dbSigType::SCAN);
  } else {
    out_net->setSigType(odb::dbSigType::SCAN);
  }

  dout->connect(out_net);
  std::visit(
      [&](auto&& term) {
        if (term != nullptr) {
          term->connect(out_net);
        }
      },
      next_cell.getScanIn().getValue());
}

void InsertLockupBetween(odb::dbDatabase* db,
                         odb::dbBlock* block,
                         utl::Logger* logger,
                         const dft::ScanStitchConfig& config,
                         const dft::ScanCell& prev_cell,
                         const dft::ScanCell& next_cell,
                         size_t chain_ordinal,
                         size_t link_idx)
{
  std::string_view lockup_cell;
  std::string_view lockup_clk_pin;
  switch (next_cell.getClockDomain().getClockEdge()) {
    case dft::ClockEdge::Rising:
      lockup_cell = config.getLockupCellRising();
      lockup_clk_pin = config.getLockupClockPinRising();
      break;
    case dft::ClockEdge::Falling:
      lockup_cell = config.getLockupCellFalling();
      lockup_clk_pin = config.getLockupClockPinFalling();
      break;
  }

  if (lockup_cell.empty()) {
    logger->error(utl::DFT,
                  113,
                  "Lockup insertion requested but no lockup cell configured "
                  "for destination clock edge '{}' (use -lockup_cell_rising or "
                  "-lockup_cell_falling)",
                  next_cell.getClockDomain().getClockEdgeName());
  }
  if (lockup_clk_pin.empty()) {
    logger->error(utl::DFT,
                  114,
                  "Lockup insertion requested but no lockup clock pin "
                  "configured for destination clock edge '{}' (use "
                  "-lockup_clock_pin_rising or -lockup_clock_pin_falling)",
                  next_cell.getClockDomain().getClockEdgeName());
  }

  odb::dbMaster* master = db->findMaster(std::string(lockup_cell).c_str());
  if (master == nullptr) {
    logger->error(utl::DFT,
                  115,
                  "Lockup cell master '{}' not found in the database",
                  lockup_cell);
  }

  const std::string inst_name
      = fmt::format("{}{}_{}", kLockupInstPrefix, chain_ordinal, link_idx);
  odb::dbInst* lockup_inst = block->findInst(inst_name.c_str());
  if (lockup_inst == nullptr) {
    lockup_inst = odb::dbInst::create(block, master, inst_name.c_str());
    if (lockup_inst == nullptr) {
      logger->error(utl::DFT,
                    116,
                    "Failed to create lockup instance '{}'",
                    inst_name);
    }

    // Seed placement near the stitched cells when possible.
    if (prev_cell.isPlaced() || next_cell.isPlaced()) {
      const odb::Point a_fallback = prev_cell.getOrigin();
      const odb::Point b_fallback = next_cell.getOrigin();
      const odb::Point a = prev_cell.getScanOut().getLocation(a_fallback);
      const odb::Point b = next_cell.getScanIn().getLocation(b_fallback);
      const int x = (prev_cell.isPlaced() && next_cell.isPlaced()) ? (a.x() + b.x()) / 2
                                                                   : (next_cell.isPlaced() ? b.x() : a.x());
      const int y = (prev_cell.isPlaced() && next_cell.isPlaced()) ? (a.y() + b.y()) / 2
                                                                   : (next_cell.isPlaced() ? b.y() : a.y());
      lockup_inst->setLocation(x, y);
      lockup_inst->setPlacementStatus(odb::dbPlacementStatus::PLACED);
    }
  } else if (lockup_inst->getMaster() != master) {
    lockup_inst->swapMaster(master);
  }

  odb::dbITerm* din
      = lockup_inst->findITerm(std::string(config.getLockupInPin()).c_str());
  odb::dbITerm* dout
      = lockup_inst->findITerm(std::string(config.getLockupOutPin()).c_str());
  odb::dbITerm* clk = lockup_inst->findITerm(std::string(lockup_clk_pin).c_str());

  if (din == nullptr) {
    logger->error(utl::DFT,
                  117,
                  "Lockup instance '{}' missing data-in pin '{}'",
                  inst_name,
                  config.getLockupInPin());
  }
  if (dout == nullptr) {
    logger->error(utl::DFT,
                  118,
                  "Lockup instance '{}' missing data-out pin '{}'",
                  inst_name,
                  config.getLockupOutPin());
  }
  if (clk == nullptr) {
    logger->error(utl::DFT,
                  119,
                  "Lockup instance '{}' missing clock pin '{}'",
                  inst_name,
                  lockup_clk_pin);
  }

  // prev.Q -> lockup.D (preserve prev functional connections)
  const dft::ScanDriver prev_scan_out = prev_cell.getScanOut();
  odb::dbNet* src_net = EnsureDriverNet(
      block,
      prev_scan_out,
      fmt::format("dft_scan_src_{}_{}", chain_ordinal, link_idx));
  if (src_net == nullptr) {
    logger->error(utl::DFT,
                  120,
                  "Failed to create scan source net for lockup between '{}' "
                  "and '{}'",
                  prev_cell.getName(),
                  next_cell.getName());
  }
  din->connect(src_net);

  // lockup.CLK -> destination clock net
  odb::dbNet* clk_net = FindClockNet(block, next_cell, logger);
  clk->connect(clk_net);

  // lockup.Q -> next.SI (new scan net)
  const std::string net_name
      = fmt::format("{}{}_{}", kLockupNetPrefix, chain_ordinal, link_idx);
  odb::dbNet* out_net = block->findNet(net_name.c_str());
  if (out_net == nullptr) {
    out_net = odb::dbNet::create(block, net_name.c_str());
    if (out_net == nullptr) {
      logger->error(utl::DFT,
                    121,
                    "Failed to create lockup net '{}'",
                    net_name);
    }
    out_net->setSigType(odb::dbSigType::SCAN);
  } else {
    out_net->setSigType(odb::dbSigType::SCAN);
  }

  dout->connect(out_net);
  std::visit(
      [&](auto&& term) {
        if (term != nullptr) {
          term->connect(out_net);
        }
      },
      next_cell.getScanIn().getValue());
}
}  // namespace

namespace dft {

namespace {
std::pair<std::string, std::optional<std::string>> SplitTermIdentifier(
    std::string_view input);
}  // namespace

ScanStitch::ScanStitch(odb::dbDatabase* db,
                       utl::Logger* logger,
                       const ScanArchitectConfig& architect_config,
                       const ScanStitchConfig& config)
    : architect_config_(architect_config),
      config_(config),
      db_(db),
      logger_(logger)
{
  odb::dbChip* chip = db_->getChip();
  top_block_ = chip->getBlock();
}

void ScanStitch::Stitch(
    const std::vector<std::unique_ptr<ScanChain>>& scan_chains)
{
  if (scan_chains.empty()) {
    return;
  }

  // Defensive: ensure scan-in/out endpoint names are distinct across chains.
  // A common misconfiguration is setting scan_in/out name patterns without "{}"
  // while requesting multiple chains, which would otherwise silently create
  // multiple chains sharing the same top-level ports.
  if (scan_chains.size() > 1) {
    std::unordered_set<std::string> in_names;
    std::unordered_set<std::string> out_names;
    in_names.reserve(scan_chains.size() * 2);
    out_names.reserve(scan_chains.size() * 2);

    for (std::size_t ordinal = 0; ordinal < scan_chains.size(); ++ordinal) {
      const ScanChain& scan_chain = *scan_chains[ordinal];
      const std::optional<ScanArchitectConfig::ChainEndpoints> endpoints
          = architect_config_.getChainEndpoints(scan_chain.getName());

      std::optional<std::string_view> begin_term;
      std::optional<std::string_view> end_term;
      if (endpoints.has_value()) {
        if (endpoints->begin.has_value()
            && endpoints->begin->type
                   == ScanArchitectConfig::ChainEndpoint::Type::Term) {
          begin_term = std::string_view(endpoints->begin->term);
        }
        if (endpoints->end.has_value()
            && endpoints->end->type
                   == ScanArchitectConfig::ChainEndpoint::Type::Term) {
          end_term = std::string_view(endpoints->end->term);
        }
      }

      std::string scan_in_name;
      std::string scan_out_name;
      try {
        scan_in_name = begin_term.has_value()
                           ? std::string(begin_term.value())
                           : fmt::format(FMT_RUNTIME(config_.getInNamePattern()),
                                         ordinal);
        scan_out_name = end_term.has_value()
                            ? std::string(end_term.value())
                            : fmt::format(FMT_RUNTIME(config_.getOutNamePattern()),
                                          ordinal);
      } catch (...) {
        logger_->error(
            utl::DFT,
            247,
            "Failed to format scan-in/out port names for chain '{}' (ordinal {})",
            scan_chain.getName(),
            ordinal);
      }

      if (!in_names.insert(scan_in_name).second) {
        logger_->error(
            utl::DFT,
            248,
            "Non-unique scan-in endpoint '{}' across chains; check "
            "-scan_in_name_pattern and/or constraints file chain begin ports",
            scan_in_name);
      }
      if (!out_names.insert(scan_out_name).second) {
        logger_->error(
            utl::DFT,
            249,
            "Non-unique scan-out endpoint '{}' across chains; check "
            "-scan_out_name_pattern and/or constraints file chain end ports",
            scan_out_name);
      }
    }
  }

  RemoveExistingLockups(top_block_);
  RemoveExistingScanBuffers(top_block_);

  bool have_existing_scan_ports = false;
  try {
    const std::string in0
        = fmt::format(FMT_RUNTIME(config_.getInNamePattern()), 0);
    const std::string out0
        = fmt::format(FMT_RUNTIME(config_.getOutNamePattern()), 0);
    have_existing_scan_ports = (top_block_->findBTerm(in0.c_str()) != nullptr)
                               || (top_block_->findBTerm(out0.c_str()) != nullptr);
  } catch (...) {
    // Ignore malformed patterns here; errors will be reported when formatting
    // is attempted for a specific chain.
  }

  size_t ordinal = 0;
  for (const std::unique_ptr<ScanChain>& scan_chain : scan_chains) {
    Stitch(top_block_, *scan_chain, ordinal, have_existing_scan_ports);
    ordinal += 1;
  }
}

void ScanStitch::Stitch(odb::dbBlock* block,
                        ScanChain& scan_chain,
                        size_t ordinal,
                        bool warn_on_missing_pattern_ports)
{
  const std::optional<ScanArchitectConfig::ChainEndpoints> endpoints
      = architect_config_.getChainEndpoints(scan_chain.getName());

  auto scan_enable_name
      = fmt::format(FMT_RUNTIME(config_.getEnableNamePattern()), kEnableNumber);
  auto scan_enable_driver = FindOrCreateScanEnable(block, scan_enable_name);

  std::optional<odb::Point> begin_pt;
  std::optional<odb::Point> end_pt;
  std::optional<std::string_view> begin_term;
  std::optional<std::string_view> end_term;
  if (endpoints.has_value()) {
    if (endpoints->begin.has_value()) {
      const auto& ep = endpoints->begin.value();
      if (ep.type == ScanArchitectConfig::ChainEndpoint::Type::Point) {
        begin_pt = odb::Point(ep.point.x, ep.point.y);
      } else {
        begin_term = std::string_view(ep.term);
      }
    }
    if (endpoints->end.has_value()) {
      const auto& ep = endpoints->end.value();
      if (ep.type == ScanArchitectConfig::ChainEndpoint::Type::Point) {
        end_pt = odb::Point(ep.point.x, ep.point.y);
      } else {
        end_term = std::string_view(ep.term);
      }
    }
  }

  const std::string scan_in_name
      = begin_term.has_value()
            ? std::string(begin_term.value())
            : fmt::format(FMT_RUNTIME(config_.getInNamePattern()), ordinal);
  ScanDriver scan_in_driver = [&]() -> ScanDriver {
    if (begin_term.has_value()) {
      const auto term_info = SplitTermIdentifier(std::string_view(scan_in_name));
      if (term_info.second.has_value()) {
        return FindOrCreateScanIn(block, scan_in_name);
      }

      odb::dbBTerm* bterm = block->findBTerm(term_info.first.c_str());
      if (bterm == nullptr) {
        logger_->error(utl::DFT,
                       224,
                       "Scan chain '{}' begin port '{}' not found",
                       scan_chain.getName(),
                       term_info.first);
      }
      return ScanDriver(bterm);
    }

    const bool existed = (block->findBTerm(scan_in_name.c_str()) != nullptr);
    ScanDriver driver = FindOrCreateScanIn(block, scan_in_name);
    if (!existed && warn_on_missing_pattern_ports) {
      logger_->warn(
          utl::DFT,
          225,
          "Scan chain '{}' is missing expected scan-in port '{}' from naming "
          "pattern; creating a new top-level port (may be unplaced)",
          scan_chain.getName(),
          scan_in_name);
    }
    return driver;
  }();
  if (begin_pt.has_value()) {
    std::visit(
        [&](auto&& term) {
          if (term == nullptr) {
            return;
          }
          using T = std::decay_t<decltype(term)>;
          if constexpr (std::is_same_v<T, odb::dbBTerm*>) {
            EnsureBTermHasLocation(block, term, begin_pt.value(), logger_);
          }
        },
        scan_in_driver.getValue());
  }
  if (!begin_pt.has_value()) {
    std::visit(
        [&](auto&& term) {
          if (term == nullptr) {
            return;
          }
          using T = std::decay_t<decltype(term)>;
          if constexpr (std::is_same_v<T, odb::dbBTerm*>) {
            int x = 0;
            int y = 0;
            if (!term->getFirstPinLocation(x, y)) {
              logger_->warn(
                  utl::DFT,
                  226,
                  "Scan chain '{}' scan-in port '{}' has no pin location; "
                  "begin/end cost and plots may be misleading (place pins or "
                  "specify endpoint coordinates)",
                  scan_chain.getName(),
                  term->getName());
            }
          }
        },
        scan_in_driver.getValue());
  }

  scan_chain.setScanIn(scan_in_driver);
  scan_chain.setScanEnable(scan_enable_driver);

  const std::vector<std::unique_ptr<ScanCell>>& scan_cells
      = scan_chain.getScanCells();
  if (scan_cells.empty()) {
    return;
  }

  // All the cells in the scan chain are controlled by the same scan enable
  for (const std::unique_ptr<ScanCell>& scan_cell : scan_cells) {
    scan_cell->connectScanEnable(scan_enable_driver);
  }

  // Let's connect the first cell
  scan_cells.front()->connectScanIn(scan_in_driver);

  // Connect scan-out to scan-in along the chain. When clock mixing is enabled
  // at the planning stage, inserting lockup latches between different
  // clock/edge domains is mandatory for correctness.
  const bool insert_lockup
      = config_.getInsertLockup()
        || architect_config_.getClockMixing()
               == ScanArchitectConfig::ClockMixing::ClockMix;
  const bool insert_timing_buffers = !config_.getTimingBufferCell().empty();
  const double critical_slack = architect_config_.getTimingCriticalSlack();
  for (size_t idx = 1; idx < scan_cells.size(); idx++) {
    const ScanCell& prev_cell = *scan_cells[idx - 1];
    const ScanCell& next_cell = *scan_cells[idx];

    const bool domain_mismatch = prev_cell.getClockDomain().getClockDomainId()
                                 != next_cell.getClockDomain().getClockDomainId();
    if (insert_lockup && domain_mismatch) {
      InsertLockupBetween(
          db_, block, logger_, config_, prev_cell, next_cell, ordinal, idx);
      continue;
    }

    if (!domain_mismatch && insert_timing_buffers
        && IsTimingCritical(prev_cell, critical_slack)) {
      InsertTimingBufferBetween(
          db_, block, logger_, config_, prev_cell, next_cell, ordinal, idx);
      continue;
    }

    scan_cells[idx]->connectScanIn(scan_cells[idx - 1]->getScanOut());
  }

  // Let's connect the last cell
  const std::unique_ptr<ScanCell>& last_scan_cell = scan_cells.back();
  const std::string scan_out_name
      = end_term.has_value()
            ? std::string(end_term.value())
            : fmt::format(FMT_RUNTIME(config_.getOutNamePattern()), ordinal);
  ScanLoad scan_out_load = [&]() -> ScanLoad {
    if (end_term.has_value()) {
      const auto term_info = SplitTermIdentifier(std::string_view(scan_out_name));
      if (term_info.second.has_value()) {
        return FindOrCreateScanOut(
            block, last_scan_cell->getScanOut(), scan_out_name);
      }

      odb::dbBTerm* bterm = block->findBTerm(term_info.first.c_str());
      if (bterm == nullptr) {
        logger_->error(utl::DFT,
                       227,
                       "Scan chain '{}' end port '{}' not found",
                       scan_chain.getName(),
                       term_info.first);
      }
      if (bterm->getIoType() != odb::dbIoType::OUTPUT) {
        logger_->error(utl::DFT,
                       228,
                       "Top-level pin '{}' specified as {} is not an output port",
                       term_info.first,
                       kScanOut);
      }
      return ScanLoad(bterm);
    }

    const bool existed = (block->findBTerm(scan_out_name.c_str()) != nullptr);
    ScanLoad load = FindOrCreateScanOut(
        block, last_scan_cell->getScanOut(), scan_out_name);
    if (!existed && warn_on_missing_pattern_ports) {
      logger_->warn(
          utl::DFT,
          229,
          "Scan chain '{}' is missing expected scan-out port '{}' from naming "
          "pattern; creating a new top-level port (may be unplaced)",
          scan_chain.getName(),
          scan_out_name);
    }
    return load;
  }();
  if (end_pt.has_value()) {
    std::visit(
        [&](auto&& term) {
          if (term == nullptr) {
            return;
          }
          using T = std::decay_t<decltype(term)>;
          if constexpr (std::is_same_v<T, odb::dbBTerm*>) {
            EnsureBTermHasLocation(block, term, end_pt.value(), logger_);
          }
        },
        scan_out_load.getValue());
  }
  if (!end_pt.has_value()) {
    std::visit(
        [&](auto&& term) {
          if (term == nullptr) {
            return;
          }
          using T = std::decay_t<decltype(term)>;
          if constexpr (std::is_same_v<T, odb::dbBTerm*>) {
            int x = 0;
            int y = 0;
            if (!term->getFirstPinLocation(x, y)) {
              logger_->warn(
                  utl::DFT,
                  230,
                  "Scan chain '{}' scan-out port '{}' has no pin location; "
                  "begin/end cost and plots may be misleading (place pins or "
                  "specify endpoint coordinates)",
                  scan_chain.getName(),
                  term->getName());
            }
          }
        },
        scan_out_load.getValue());
  }
  last_scan_cell->connectScanOut(scan_out_load);
  scan_chain.setScanOut(scan_out_load);
}

namespace {
std::string UnescapeSlash(std::string_view input)
{
  std::string out;
  out.reserve(input.size());
  for (std::size_t i = 0; i < input.size(); ++i) {
    const char c = input[i];
    if (c == '\\' && i + 1 < input.size() && input[i + 1] == '/') {
      continue;  // drop the escape, keep '/'
    }
    out.push_back(c);
  }
  return out;
}

std::pair<std::string, std::optional<std::string>> SplitTermIdentifier(
    std::string_view input)
{
  size_t tracker = 0;
  size_t slash_position;
  while ((slash_position = input.find('/', tracker)) != std::string_view::npos) {
    if (slash_position != 0 && input[slash_position - 1] == '\\') {
      tracker = slash_position + 1;
      continue;
    }
    const std::string inst = UnescapeSlash(input.substr(0, slash_position));
    const std::string pin = UnescapeSlash(input.substr(slash_position + 1));
    return {inst, pin};
  }
  return {UnescapeSlash(input), std::nullopt};
}
}  // namespace

ScanDriver ScanStitch::FindOrCreateDriver(std::string_view kind,
                                          odb::dbBlock* block,
                                          const std::string& with_name)
{
  auto term_info = SplitTermIdentifier(std::string_view(with_name));

  if (term_info.second.has_value()) {  // Instance/ITerm
    auto inst = block->findInst(term_info.first.c_str());
    if (inst == nullptr) {
      logger_->error(utl::DFT,
                     34,
                     "Instance {} not found for {} port",
                     term_info.first,
                     kind);
    }
    auto iterm = inst->findITerm(term_info.second.value().c_str());
    if (iterm == nullptr) {
      logger_->error(utl::DFT,
                     35,
                     "ITerm {}/{} not found for {} port",
                     term_info.first,
                     term_info.second.value(),
                     kind);
    }
    if (iterm->getIoType() != odb::dbIoType::OUTPUT) {
      logger_->error(utl::DFT,
                     36,
                     "ITerm {}/{} for {} port is a {}",
                     term_info.first,
                     term_info.second.value(),
                     kind,
                     iterm->getIoType().getString());
    }
    return ScanDriver(iterm);
  }

  // BTerm
  auto bterm = block->findBTerm(with_name.data());
  if (bterm != nullptr) {
    // We don't actually care if it's an output, that works here too.
    return ScanDriver(bterm);
  }
  return CreateNewPort<ScanDriver>(block, with_name);
}

ScanDriver ScanStitch::FindOrCreateScanEnable(odb::dbBlock* block,
                                              const std::string& with_name)
{
  return FindOrCreateDriver(kScanEnable, block, with_name);
}

ScanDriver ScanStitch::FindOrCreateScanIn(odb::dbBlock* block,
                                          const std::string& with_name)
{
  return FindOrCreateDriver(kScanIn, block, with_name);
}

ScanLoad ScanStitch::FindOrCreateScanOut(odb::dbBlock* block,
                                         const ScanDriver& cell_scan_out,
                                         const std::string& with_name)
{
  auto term_info = SplitTermIdentifier(std::string_view(with_name));

  if (term_info.second.has_value()) {  // Instance/ITerm
    auto inst = block->findInst(term_info.first.c_str());
    if (inst == nullptr) {
      logger_->error(utl::DFT,
                     37,
                     "Instance {} not found for {} port",
                     term_info.first,
                     kScanOut);
    }
    auto iterm = inst->findITerm(term_info.second.value().c_str());
    if (iterm == nullptr) {
      logger_->error(utl::DFT,
                     38,
                     "ITerm {}/{} not found for {} port",
                     term_info.first,
                     term_info.second.value(),
                     kScanOut);
    }
    if (iterm->getIoType() != odb::dbIoType::INPUT) {
      logger_->error(utl::DFT,
                     39,
                     "ITerm {}/{} for {} port is a {}",
                     term_info.first,
                     term_info.second.value(),
                     kScanOut,
                     iterm->getIoType().getString());
    }
    return ScanLoad(iterm);
  }
  auto bterm = block->findBTerm(with_name.data());
  if (bterm != nullptr) {
    if (bterm->getIoType() != odb::dbIoType::OUTPUT) {
      logger_->error(utl::DFT,
                     40,
                     "Top-level pin '{}' specified as {} is not an output port",
                     term_info.first,
                     kScanOut);
    }
    return ScanLoad(bterm);
  }

  // Prefer creating a dedicated scan-out BTerm with the requested name.
  // If the scan-out driver is already connected to a net, attach the new BTerm
  // directly to that net to avoid creating and then orphaning a temporary net.
  odb::dbNet* scan_out_net = cell_scan_out.getNet();
  if (scan_out_net && top_block_ == scan_out_net->getBlock()) {
    return CreateNewPort<ScanLoad>(block, with_name, scan_out_net);
  }

  return CreateNewPort<ScanLoad>(block, with_name);
}

}  // namespace dft
