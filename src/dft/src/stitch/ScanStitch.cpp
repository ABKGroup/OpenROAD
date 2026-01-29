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
      const odb::Point a = prev_cell.getOrigin();
      const odb::Point b = next_cell.getOrigin();
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
      const odb::Point a = prev_cell.getOrigin();
      const odb::Point b = next_cell.getOrigin();
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

  RemoveExistingLockups(top_block_);
  RemoveExistingScanBuffers(top_block_);

  size_t ordinal = 0;
  for (const std::unique_ptr<ScanChain>& scan_chain : scan_chains) {
    Stitch(top_block_, *scan_chain, ordinal);
    ordinal += 1;
  }
}

void ScanStitch::Stitch(odb::dbBlock* block,
                        ScanChain& scan_chain,
                        size_t ordinal)
{
  auto scan_enable_name
      = fmt::format(FMT_RUNTIME(config_.getEnableNamePattern()), kEnableNumber);
  auto scan_enable_driver = FindOrCreateScanEnable(block, scan_enable_name);

  auto scan_in_name
      = fmt::format(FMT_RUNTIME(config_.getInNamePattern()), ordinal);
  ScanDriver scan_in_driver = FindOrCreateScanIn(block, scan_in_name);

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
  auto scan_out_name
      = fmt::format(FMT_RUNTIME(config_.getOutNamePattern()), ordinal);
  ScanLoad scan_out_load
      = FindOrCreateScanOut(block, last_scan_cell->getScanOut(), scan_out_name);
  last_scan_cell->connectScanOut(scan_out_load);
  scan_chain.setScanOut(scan_out_load);
}

namespace {
static std::pair<std::string, std::optional<std::string>> SplitTermIdentifier(
    const std::string& input)
{
  size_t tracker = 0;
  size_t slash_position;
  while ((slash_position = input.find('/', tracker)) != std::string::npos) {
    if (slash_position != 0 && input[slash_position - 1] == '\\') {
      tracker = slash_position + 1;
      continue;
    }
    return {input.substr(0, slash_position), input.substr(slash_position + 1)};
  }
  return {input, std::nullopt};
}
}  // namespace

ScanDriver ScanStitch::FindOrCreateDriver(std::string_view kind,
                                          odb::dbBlock* block,
                                          const std::string& with_name)
{
  auto term_info = SplitTermIdentifier(with_name);

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
  auto term_info = SplitTermIdentifier(with_name);

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

  // TODO: Trace forward the scan out net so we can see if it is connected to a
  // top port or to functional logic
  odb::dbNet* scan_out_net = cell_scan_out.getNet();
  if (scan_out_net && top_block_ == scan_out_net->getBlock()) {
    // if the scan_out_net exists, and has an BTerm that is an OUTPUT, then we
    // can reuse that BTerm to act as scan_out, only if the block of the bterm
    // is the top block, otherwise we will punch a new port
    for (odb::dbBTerm* bterm : scan_out_net->getBTerms()) {
      if (bterm->getIoType() == odb::dbIoType::OUTPUT) {
        return ScanLoad(bterm);
      }
    }
  }

  return CreateNewPort<ScanLoad>(block, with_name);
}

}  // namespace dft
