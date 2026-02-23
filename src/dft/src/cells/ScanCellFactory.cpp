// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2023-2025, The OpenROAD Authors

#include "ScanCellFactory.hh"

#include <algorithm>
#include <cctype>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <memory>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include "ClockDomain.hh"
#include "OneBitScanCell.hh"
#include "ScanCell.hh"
#include "OneBitScanCell.hh"
#include "ScanCell.hh"
#include "Utils.hh"
#include "db_sta/dbNetwork.hh"
#include "db_sta/dbSta.hh"
#include "odb/db.h"
#include "sta/Clock.hh"
#include "sta/FuncExpr.hh"
#include "sta/Liberty.hh"
#include "sta/LibertyClass.hh"
#include "sta/MinMax.hh"
#include "sta/NetworkClass.hh"
#include "sta/Sequential.hh"
#include "sta/Transition.hh"
#include "utl/Logger.h"

namespace dft {

namespace {

enum class TypeOfCell
{
  kOneBitCell,
  kNotSupported
};

bool EqualsIgnoreCase(const char* lhs, std::string_view rhs)
{
  if (lhs == nullptr) {
    return false;
  }
  const std::size_t n = std::strlen(lhs);
  if (n != rhs.size()) {
    return false;
  }
  for (std::size_t i = 0; i < n; ++i) {
    if (std::tolower(static_cast<unsigned char>(lhs[i]))
        != std::tolower(static_cast<unsigned char>(rhs[i]))) {
      return false;
    }
  }
  return true;
}

bool HasPrefixIgnoreCase(const char* name, std::string_view prefix)
{
  if (name == nullptr) {
    return false;
  }
  const std::size_t n = std::strlen(name);
  if (n < prefix.size()) {
    return false;
  }
  for (std::size_t i = 0; i < prefix.size(); ++i) {
    if (std::tolower(static_cast<unsigned char>(name[i]))
        != std::tolower(static_cast<unsigned char>(prefix[i]))) {
      return false;
    }
  }
  return true;
}

bool MatchesScanPrefix(const char* name, std::string_view prefix)
{
  if (!HasPrefixIgnoreCase(name, prefix)) {
    return false;
  }
  const std::size_t n = std::strlen(name);
  if (n == prefix.size()) {
    return true;
  }
  const char next = name[prefix.size()];
  if (std::isdigit(static_cast<unsigned char>(next)) || next == '[') {
    return true;
  }
  if (next == '_' && (prefix.size() + 1) < n
      && std::isdigit(static_cast<unsigned char>(name[prefix.size() + 1]))) {
    return true;
  }
  return false;
}

int ExtractPortIndex(const char* name)
{
  if (name == nullptr) {
    return -1;
  }
  const std::size_t n = std::strlen(name);
  if (n == 0) {
    return -1;
  }

  // Bracket form: FOO[12]
  if (n >= 3 && name[n - 1] == ']') {
    std::size_t l = n - 2;
    while (l > 0 && std::isdigit(static_cast<unsigned char>(name[l]))) {
      --l;
    }
    if (name[l] == '[' && (l + 1) < (n - 1)) {
      return std::atoi(name + l + 1);
    }
  }

  // Trailing digits: FOO12
  std::size_t start = n;
  while (start > 0 && std::isdigit(static_cast<unsigned char>(name[start - 1]))) {
    --start;
  }
  if (start < n) {
    return std::atoi(name + start);
  }
  return -1;
}

template <typename Pred>
std::vector<sta::LibertyPort*> CollectPortsBySignalType(
    const sta::LibertyCell* lib_cell,
    Pred predicate)
{
  std::vector<sta::LibertyPort*> out;
  sta::LibertyCellPortIterator iter(lib_cell);
  while (iter.hasNext()) {
    sta::LibertyPort* port = iter.next();
    if (port != nullptr && predicate(port->scanSignalType())) {
      out.push_back(port);
    }
  }
  return out;
}

std::vector<sta::LibertyPort*> CollectScanPortsByNameFallback(
    const sta::LibertyCell* lib_cell,
    const std::vector<std::string_view>& prefixes,
    sta::ScanSignalType set_type)
{
  std::vector<sta::LibertyPort*> out;
  sta::LibertyCellPortIterator iter(lib_cell);
  while (iter.hasNext()) {
    sta::LibertyPort* port = iter.next();
    if (port == nullptr) {
      continue;
    }
    const char* port_name = port->name();
    for (const auto& prefix : prefixes) {
      if (MatchesScanPrefix(port_name, prefix) || EqualsIgnoreCase(port_name, prefix)) {
        if (port->scanSignalType() == sta::ScanSignalType::none) {
          port->setScanSignalType(set_type);
        }
        out.push_back(port);
        break;
      }
    }
  }
  return out;
}

void SortScanPorts(std::vector<sta::LibertyPort*>& ports)
{
  std::stable_sort(ports.begin(), ports.end(), [](sta::LibertyPort* a, sta::LibertyPort* b) {
    const char* a_name = a ? a->name() : nullptr;
    const char* b_name = b ? b->name() : nullptr;
    const int a_idx = ExtractPortIndex(a_name);
    const int b_idx = ExtractPortIndex(b_name);
    if (a_idx >= 0 && b_idx >= 0 && a_idx != b_idx) {
      return a_idx < b_idx;
    }
    if (a_idx >= 0 && b_idx < 0) {
      return true;
    }
    if (a_idx < 0 && b_idx >= 0) {
      return false;
    }
    const std::string_view as = a_name ? std::string_view(a_name) : std::string_view();
    const std::string_view bs = b_name ? std::string_view(b_name) : std::string_view();
    return as < bs;
  });
}

sta::LibertyCell* GetLibertyCell(odb::dbMaster* master,
                                 sta::dbNetwork* db_network)
{
  sta::Cell* master_cell = db_network->dbToSta(master);
  return db_network->libertyCell(master_cell);
}

TypeOfCell IdentifyCell(odb::dbInst* inst, sta::dbSta* sta)
{
  sta::dbNetwork* db_network = sta->getDbNetwork();
  sta::LibertyCell* liberty_cell
      = GetLibertyCell(inst->getMaster(), db_network);
  if (liberty_cell != nullptr && liberty_cell->hasSequentials()
      && !inst->getMaster()->isBlock()) {
    // we assume that we are only dealing with one bit cells, but in the future
    // we could deal with multibit cells too
    return TypeOfCell::kOneBitCell;
  }
  return TypeOfCell::kNotSupported;
}

std::unique_ptr<ClockDomain> GetClockDomainFromClock(
    sta::LibertyCell* liberty_cell,
    sta::Clock* clock,
    odb::dbITerm* clock_pin)
{
  ClockEdge edge = ClockEdge::Rising;
  const sta::SequentialSeq& sequentials = liberty_cell->sequentials();
  // TODO: Other cells may have more than one sequential
  for (const sta::Sequential& sequential : sequentials) {
    // If the clock has a left FuncExpr, then
    sta::FuncExpr* left = sequential.clock()->left();
    sta::FuncExpr* right = sequential.clock()->right();
    if (left && !right) {
      //  When operator is NOT, left is the only operand
      edge = ClockEdge::Falling;
    } else {
      edge = ClockEdge::Rising;
    }
  }

  // TODO: Create the clock domain based on the timing instead of the name to
  // better control equivalent clocks
  return std::make_unique<ClockDomain>(clock->name(), edge);
}

std::unique_ptr<ClockDomain> FindOneBitCellClockDomain(odb::dbInst* inst,
                                                       sta::dbSta* sta,
                                                       utl::Logger* logger)
{
  std::vector<odb::dbITerm*> clock_pins = utils::GetClockPin(inst);
  if (clock_pins.empty()) {
    logger->warn(utl::DFT,
                 49,
                 "Can't find a clock pin for cell '{:s}'",
                 inst->getName());
    return nullptr;
  }
  sta::dbNetwork* db_network = sta->getDbNetwork();
  // A one bit cell should only have one clock pin

  odb::dbITerm* clock_pin = clock_pins.at(0);
  std::optional<sta::Clock*> clock = utils::GetClock(sta, clock_pin);

  if (clock.has_value()) {
    return GetClockDomainFromClock(
        GetLibertyCell(inst->getMaster(), db_network), *clock, clock_pin);
  }

  return nullptr;
}

std::vector<std::unique_ptr<ScanCell>> CreateOneBitCells(
    odb::dbInst* inst,
    sta::dbSta* sta,
    const ScanArchitectConfig& config,
    utl::Logger* logger)
{
  std::vector<std::unique_ptr<ScanCell>> out;
  sta::dbNetwork* db_network = sta->getDbNetwork();
  sta::LibertyCell* liberty_cell
      = GetLibertyCell(inst->getMaster(), db_network);
  if (liberty_cell == nullptr) {
    logger->warn(utl::DFT,
                 11,
                 "Cell master '{:s}' has no lib info. Can't find scan cell",
                 inst->getMaster()->getName());
    return out;
  }

  std::unique_ptr<ClockDomain> clock_domain
      = FindOneBitCellClockDomain(inst, sta, logger);

  if (!clock_domain) {
    logger->warn(utl::DFT,
                 5,
                 "Cell '{:s}' doesn't have a valid clock connected. Can't "
                 "create a scan cell",
                 inst->getName());
    return out;
  }

  sta::LibertyPort* scan_in_port = getLibertyScanIn(liberty_cell);
  sta::LibertyPort* scan_enable_port = getLibertyScanEnable(liberty_cell);

  if (scan_in_port == nullptr || scan_enable_port == nullptr) {
    logger->warn(
        utl::DFT,
        7,
        "Cell '{:s}' is not a scan cell. Can't use it for scan architect",
        inst->getName());
    return out;
  }

  sta::LibertyPort* scan_out_port = getLibertyScanOut(liberty_cell);
  if (scan_out_port == nullptr) {
    // Many libraries do not explicitly tag scan-out ports; fall back to the
    // functional output.
    if (config.getPreferQbarScanOut()) {
      scan_out_port = liberty_cell->findLibertyPort("QN");
      if (scan_out_port == nullptr) {
        scan_out_port = liberty_cell->findLibertyPort("Q_N");
      }
      if (scan_out_port == nullptr) {
        scan_out_port = liberty_cell->findLibertyPort("Q");
      }
    } else {
      scan_out_port = liberty_cell->findLibertyPort("Q");
      if (scan_out_port == nullptr) {
        scan_out_port = liberty_cell->findLibertyPort("QN");
      }
      if (scan_out_port == nullptr) {
        scan_out_port = liberty_cell->findLibertyPort("Q_N");
      }
    }
    if (scan_out_port == nullptr) {
      const sta::SequentialSeq& sequentials = liberty_cell->sequentials();
      if (!sequentials.empty()) {
        scan_out_port = sequentials.front()->output();
      }
    }
  }

  if (scan_out_port == nullptr) {
    logger->warn(utl::DFT,
                 50,
                 "Cell '{:s}' has no identifiable scan-out port. Can't create "
                 "a scan cell",
                 inst->getName());
    return out;
  }

  uint64_t bits_total = 0;
  const sta::SequentialSeq& sequentials = liberty_cell->sequentials();
  for (const sta::Sequential* seq : sequentials) {
    if (seq != nullptr && seq->isRegister()) {
      bits_total += 1;
    }
  }
  if (bits_total == 0) {
    bits_total = 1;
  }

  if (config.getSplitMultibitScanCells()) {
    std::vector<sta::LibertyPort*> scan_ins = CollectPortsBySignalType(
        liberty_cell,
        [](sta::ScanSignalType t) {
          return t == sta::ScanSignalType::input
                 || t == sta::ScanSignalType::input_inverted;
        });
    if (scan_ins.empty()) {
      scan_ins = CollectScanPortsByNameFallback(
          liberty_cell,
          {"SI", "SCD", "SCAN_IN", "SCANIN"},
          sta::ScanSignalType::input);
    }

    std::vector<sta::LibertyPort*> scan_outs = CollectPortsBySignalType(
        liberty_cell,
        [](sta::ScanSignalType t) {
          return t == sta::ScanSignalType::output
                 || t == sta::ScanSignalType::output_inverted;
        });
    if (scan_outs.empty()) {
      scan_outs = CollectScanPortsByNameFallback(
          liberty_cell,
          {"SO", "SCO", "SCAN_OUT", "SCANOUT"},
          sta::ScanSignalType::output);
    }

    if (scan_ins.size() > 1 && scan_outs.size() == scan_ins.size()) {
      SortScanPorts(scan_ins);
      SortScanPorts(scan_outs);

      const std::string clock_name(clock_domain->getClockName());
      const ClockEdge edge = clock_domain->getClockEdge();

      const std::size_t n = scan_ins.size();
      const uint64_t total = std::max<uint64_t>(bits_total, n);
      const uint64_t base = total / n;
      const uint64_t rem = total % n;

      out.reserve(n);
      for (std::size_t i = 0; i < n; ++i) {
        const uint64_t bits = std::max<uint64_t>(1, base + (i < rem ? 1 : 0));
        std::unique_ptr<ClockDomain> cd
            = (i == 0) ? std::move(clock_domain)
                       : std::make_unique<ClockDomain>(clock_name, edge);
        std::string name = inst->getName();
        name += "__dft_scanbit";
        name += std::to_string(i);
        out.push_back(std::make_unique<OneBitScanCell>(name,
                                                       std::move(cd),
                                                       bits,
                                                       inst,
                                                       scan_ins[i],
                                                       scan_enable_port,
                                                       scan_outs[i],
                                                       db_network,
                                                       logger));
      }
      return out;
    }
  }

  out.push_back(std::make_unique<OneBitScanCell>(inst->getName(),
                                                 std::move(clock_domain),
                                                 bits_total,
                                                 inst,
                                                 scan_in_port,
                                                 scan_enable_port,
                                                 scan_out_port,
                                                 db_network,
                                                 logger));
  return out;
}

std::vector<std::unique_ptr<ScanCell>> ScanCellFactory(
    odb::dbInst* inst,
    sta::dbSta* sta,
    const ScanArchitectConfig& config,
    utl::Logger* logger)
{
  std::vector<std::unique_ptr<ScanCell>> out;
  TypeOfCell type_of_cell = IdentifyCell(inst, sta);

  switch (type_of_cell) {
    case TypeOfCell::kOneBitCell:
      return CreateOneBitCells(inst, sta, config, logger);
    default:
      return out;
  }
}

std::optional<std::pair<float, float>> computeScanOutTimingSlacks(
    const ScanCell& cell,
    sta::dbSta* sta)
{
  sta::dbNetwork* db_network = sta->getDbNetwork();
  if (db_network == nullptr) {
    return std::nullopt;
  }

  sta::Pin* pin = nullptr;
  const ScanDriver scan_out = cell.getScanOut();
  std::visit(
      [&](auto&& term) {
        if (term != nullptr) {
          pin = db_network->dbToSta(term);
        }
      },
      scan_out.getValue());
  if (pin == nullptr) {
    return std::nullopt;
  }

  const float setup_rise
      = sta->pinSlack(pin, sta::RiseFall::rise(), sta::MinMax::max());
  const float setup_fall
      = sta->pinSlack(pin, sta::RiseFall::fall(), sta::MinMax::max());
  const float hold_rise
      = sta->pinSlack(pin, sta::RiseFall::rise(), sta::MinMax::min());
  const float hold_fall
      = sta->pinSlack(pin, sta::RiseFall::fall(), sta::MinMax::min());

  const float setup = std::min(setup_rise, setup_fall);
  const float hold = std::min(hold_rise, hold_fall);

  if (setup >= sta::INF / 2.0F && hold >= sta::INF / 2.0F) {
    return std::nullopt;
  }

  return std::make_pair(setup, hold);
}

void CollectScanCells(odb::dbBlock* block,
                      sta::dbSta* sta,
                      const ScanArchitectConfig& config,
                      bool compute_timing_slacks,
                      utl::Logger* logger,
                      std::vector<std::unique_ptr<ScanCell>>& scan_cells)
{
  for (odb::dbInst* inst : block->getInsts()) {
    if (inst->isDoNotTouch()) {
      // Do not scan
      continue;
    }
    // TODO: Support doNotScan later
    if (config.isInstanceExcluded(inst->getName(), inst->getMaster()->getName())) {
      continue;
    }

    std::vector<std::unique_ptr<ScanCell>> new_cells
        = ScanCellFactory(inst, sta, config, logger);
    for (auto& scan_cell : new_cells) {
      if (scan_cell == nullptr) {
        continue;
      }
      if (compute_timing_slacks) {
        if (auto slacks = computeScanOutTimingSlacks(*scan_cell, sta)) {
          scan_cell->setTimingSlacks(slacks->first, slacks->second);
        }
      }
      scan_cells.push_back(std::move(scan_cell));
    }
  }

  // Go inside the next blocks
  for (odb::dbBlock* next_block : block->getChildren()) {
    CollectScanCells(next_block,
                     sta,
                     config,
                     compute_timing_slacks,
                     logger,
                     scan_cells);
  }
}

}  // namespace

std::vector<std::unique_ptr<ScanCell>> CollectScanCells(odb::dbDatabase* db,
                                                        sta::dbSta* sta,
                                                        const ScanArchitectConfig& config,
                                                        utl::Logger* logger)
{
  std::vector<std::unique_ptr<ScanCell>> scan_cells;

  odb::dbChip* chip = db->getChip();
  const bool compute_timing_slacks
      = (config.getTimingWeightSetup() != 0.0
         || config.getTimingWeightHold() != 0.0);
  CollectScanCells(chip->getBlock(),
                   sta,
                   config,
                   compute_timing_slacks,
                   logger,
                   scan_cells);

  // To keep preview_dft consistent between calls and rollbacks
  std::sort(scan_cells.begin(), scan_cells.end(), [](const auto& lhs, const auto& rhs) {
    return lhs->getName() < rhs->getName();
  });

  return scan_cells;
}

}  // namespace dft
