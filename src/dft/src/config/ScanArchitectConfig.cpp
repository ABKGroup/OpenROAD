// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2023-2025, The OpenROAD Authors

#include "ScanArchitectConfig.hh"

#include <algorithm>
#include <cctype>
#include <cstdint>
#include <fstream>
#include <optional>
#include <sstream>
#include <string>

#include "Formatting.hh"
#include "utl/Logger.h"

namespace dft {

ScanArchitectConfig::ClockMixing ScanArchitectConfig::getClockMixing() const
{
  return clock_mixing_;
}

ScanArchitectConfig::ScanOrderMetric ScanArchitectConfig::getScanOrderMetric()
    const
{
  return scan_order_metric_;
}

ScanArchitectConfig::ScanOrderSolver ScanArchitectConfig::getScanOrderSolver()
    const
{
  return scan_order_solver_;
}

uint64_t ScanArchitectConfig::getScanOptRounds() const
{
  return scanopt_rounds_;
}

uint64_t ScanArchitectConfig::getScanOptSeed() const
{
  return scanopt_seed_;
}

double ScanArchitectConfig::getVerticalWeight() const
{
  return vertical_weight_;
}

double ScanArchitectConfig::getTimingWeightSetup() const
{
  return timing_weight_setup_;
}

double ScanArchitectConfig::getTimingWeightHold() const
{
  return timing_weight_hold_;
}

double ScanArchitectConfig::getTimingCriticalSlack() const
{
  return timing_critical_slack_;
}

const std::optional<uint64_t>& ScanArchitectConfig::getChainCount() const
{
  return chain_count_;
}

const std::optional<uint64_t>& ScanArchitectConfig::getMaxLength() const
{
  return max_length_;
}

const std::optional<uint64_t>& ScanArchitectConfig::getMaxChains() const
{
  return max_chains_;
}

void ScanArchitectConfig::setClockMixing(ClockMixing clock_mixing)
{
  clock_mixing_ = clock_mixing;
}

void ScanArchitectConfig::setScanOrderMetric(ScanOrderMetric metric)
{
  scan_order_metric_ = metric;
}

void ScanArchitectConfig::setScanOrderSolver(ScanOrderSolver solver)
{
  scan_order_solver_ = solver;
}

void ScanArchitectConfig::setScanOptRounds(uint64_t rounds)
{
  scanopt_rounds_ = rounds;
}

void ScanArchitectConfig::setScanOptSeed(uint64_t seed)
{
  scanopt_seed_ = seed;
}

void ScanArchitectConfig::setVerticalWeight(double weight)
{
  if (weight > 0.0) {
    vertical_weight_ = weight;
  }
}

void ScanArchitectConfig::setTimingWeightSetup(double weight)
{
  if (weight >= 0.0) {
    timing_weight_setup_ = weight;
  }
}

void ScanArchitectConfig::setTimingWeightHold(double weight)
{
  if (weight >= 0.0) {
    timing_weight_hold_ = weight;
  }
}

void ScanArchitectConfig::setTimingCriticalSlack(double slack)
{
  if (slack >= 0.0) {
    timing_critical_slack_ = slack;
  }
}

void ScanArchitectConfig::setChainCount(uint64_t chain_count)
{
  chain_count_ = chain_count;
}

void ScanArchitectConfig::setMaxChains(uint64_t max_chains)
{
  max_chains_ = max_chains;
}

void ScanArchitectConfig::setMaxLength(uint64_t max_length)
{
  max_length_ = max_length;
}

void ScanArchitectConfig::clearScanOrderConstraints()
{
  scan_order_groups_.clear();
  scan_order_fixed_edges_.clear();
}

void ScanArchitectConfig::setDefaultGroupPriority(int priority)
{
  default_group_priority_ = std::clamp(priority, 0, 127);
}

int ScanArchitectConfig::getDefaultGroupPriority() const
{
  return default_group_priority_;
}

const std::vector<ScanArchitectConfig::ScanOrderGroupConstraint>&
ScanArchitectConfig::getScanOrderGroups() const
{
  return scan_order_groups_;
}

const std::vector<ScanArchitectConfig::ScanOrderFixedEdgeConstraint>&
ScanArchitectConfig::getScanOrderFixedEdges() const
{
  return scan_order_fixed_edges_;
}

bool ScanArchitectConfig::loadScanOrderConstraintsFile(
    const std::string& path,
    utl::Logger* logger)
{
  std::ifstream f(path);
  if (!f) {
    if (logger) {
      logger->error(utl::DFT, 76, "Couldn't open scan constraints file '{}'", path);
    }
    return false;
  }

  clearScanOrderConstraints();

  const auto trim = [](std::string& s) {
    auto not_space = [](unsigned char c) { return !std::isspace(c); };
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), not_space));
    s.erase(std::find_if(s.rbegin(), s.rend(), not_space).base(), s.end());
  };

  std::string raw;
  int lineno = 0;
  while (std::getline(f, raw)) {
    ++lineno;
    std::string line = raw;
    if (const std::size_t hash = line.find('#'); hash != std::string::npos) {
      line = line.substr(0, hash);
    }
    trim(line);
    if (line.empty()) {
      continue;
    }

    std::istringstream iss(line);
    std::string keyword;
    iss >> keyword;
    std::transform(keyword.begin(),
                   keyword.end(),
                   keyword.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });

    if (keyword == "default_priority" || keyword == "default_group_priority") {
      int pr = 0;
      if (!(iss >> pr)) {
        if (logger) {
          logger->warn(utl::DFT,
                       77,
                       "Scan constraints parse error at {}:{}: expected integer "
                       "priority",
                       path,
                       lineno);
        }
        continue;
      }
      setDefaultGroupPriority(pr);
      continue;
    }

    if (keyword == "group") {
      int pr = default_group_priority_;
      if (!(iss >> pr)) {
        if (logger) {
          logger->warn(
              utl::DFT,
              78,
              "Scan constraints parse error at {}:{}: expected 'group <priority> "
              "<inst...>'",
              path,
              lineno);
        }
        continue;
      }
      pr = std::clamp(pr, 0, 127);
      std::vector<std::string> insts;
      std::string inst;
      while (iss >> inst) {
        insts.push_back(inst);
      }
      if (insts.empty()) {
        if (logger) {
          logger->warn(utl::DFT,
                       79,
                       "Scan constraints parse error at {}:{}: group has no "
                       "instances",
                       path,
                       lineno);
        }
        continue;
      }
      scan_order_groups_.push_back(
          {.priority = pr, .inst_names = std::move(insts)});
      continue;
    }

    if (keyword == "fixed_edge" || keyword == "fixededge" || keyword == "edge") {
      std::string from;
      std::string to;
      if (!(iss >> from >> to)) {
        if (logger) {
          logger->warn(
              utl::DFT,
              80,
              "Scan constraints parse error at {}:{}: expected 'fixed_edge <from> "
              "<to>'",
              path,
              lineno);
        }
        continue;
      }
      scan_order_fixed_edges_.push_back(
          {.from_inst = std::move(from), .to_inst = std::move(to)});
      continue;
    }

    if (logger) {
      logger->warn(utl::DFT,
                   81,
                   "Scan constraints: ignoring unknown directive '{}' at {}:{}",
                   keyword,
                   path,
                   lineno);
    }
  }

  return true;
}

void ScanArchitectConfig::report(utl::Logger* logger) const
{
  logger->report("Scan Architect Config:");
  logger->report("- Chain Count: {}", utils::FormatForReport(chain_count_));
  logger->report("- Max Length: {}", utils::FormatForReport(max_length_));
  logger->report("- Max Chains: {}", utils::FormatForReport(max_chains_));
  logger->report("- Clock Mixing: {}", ClockMixingName(clock_mixing_));
  logger->report("- Scan Order Metric: {}",
                 ScanOrderMetricName(scan_order_metric_));
  logger->report("- Scan Order Solver: {}",
                 ScanOrderSolverName(scan_order_solver_));
  logger->report("- Vertical Weight: {:.3f}", vertical_weight_);
  if (timing_weight_setup_ != 0.0 || timing_weight_hold_ != 0.0) {
    logger->report("- Timing Setup Weight: {:.3f}", timing_weight_setup_);
    logger->report("- Timing Hold Weight: {:.3f}", timing_weight_hold_);
    logger->report("- Timing Critical Slack: {:.3f}", timing_critical_slack_);
  }
  if (scan_order_solver_ == ScanOrderSolver::ScanOpt) {
    logger->report("- ScanOpt Rounds: {}", scanopt_rounds_);
    logger->report("- ScanOpt Seed: {}", scanopt_seed_);
  }
  if (!scan_order_groups_.empty() || !scan_order_fixed_edges_.empty()) {
    logger->report("- Scan Order Constraints:");
    logger->report("  - Default Priority: {}", default_group_priority_);
    logger->report("  - Groups: {}", scan_order_groups_.size());
    logger->report("  - Fixed Edges: {}", scan_order_fixed_edges_.size());
  }
}

std::string ScanArchitectConfig::ClockMixingName(
    ScanArchitectConfig::ClockMixing clock_mixing)
{
  switch (clock_mixing) {
    case ScanArchitectConfig::ClockMixing::NoMix:
      return "No Mix";
    case ScanArchitectConfig::ClockMixing::ClockMix:
      return "Clock Mix";
    default:
      return "Missing case in ClockMixingName";
  }
}

std::string ScanArchitectConfig::ScanOrderMetricName(
    ScanArchitectConfig::ScanOrderMetric metric)
{
  switch (metric) {
    case ScanArchitectConfig::ScanOrderMetric::Placement:
      return "Placement (cell-to-cell)";
    case ScanArchitectConfig::ScanOrderMetric::PinToNet:
      return "Pin-to-net (routing-aware)";
    default:
      return "Missing case in ScanOrderMetricName";
  }
}

std::string ScanArchitectConfig::ScanOrderSolverName(
    ScanArchitectConfig::ScanOrderSolver solver)
{
  switch (solver) {
    case ScanArchitectConfig::ScanOrderSolver::Heuristic:
      return "Heuristic";
    case ScanArchitectConfig::ScanOrderSolver::ScanOpt:
      return "ScanOpt";
    case ScanArchitectConfig::ScanOrderSolver::MinFeedthrough:
      return "Min-feedthrough (row-sweep)";
    default:
      return "Missing case in ScanOrderSolverName";
  }
}

}  // namespace dft
