// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2023-2025, The OpenROAD Authors

#pragma once

#include <cstdint>
#include <optional>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

#include "utl/Logger.h"

namespace dft {

class ScanArchitectConfig
{
 public:
  // TODO Add suport for mix_edges, mix_clocks, mix_clocks_not_edges
  enum class ClockMixing
  {
    NoMix,    // We create different scan chains for each clock and edge
    ClockMix  // We architect the flops of different clock and edge together
  };

  // Metric used for ordering scan cells within each scan chain.
  enum class ScanOrderMetric
  {
    Placement,  // Placement-based Manhattan (cell-to-cell)
    PinToNet    // Routing-aware incremental cost (pin-to-net)
  };

  // Solver used for ordering scan cells within each scan chain.
  enum class ScanOrderSolver
  {
    Heuristic,  // Greedy + local search heuristics (fast)
    ScanOpt,    // ScanOpt-style iterated local search (higher quality)
    MinFeedthrough  // Row-sweep DP for min-feedthrough special case
  };

  struct ScanOrderGroupConstraint
  {
    std::string name;
    int priority = 64;  // lower runs earlier
    std::vector<std::string> inst_names;
  };

  struct ScanOrderFixedEdgeConstraint
  {
    std::string from_inst;
    std::string to_inst;
  };

  struct ScanOrderBeforeConstraint
  {
    std::string before;
    std::string after;
  };

  struct Point
  {
    int x = 0;
    int y = 0;
  };

  struct ChainEndpoints
  {
    std::optional<Point> begin;
    std::optional<Point> end;
  };

  void setClockMixing(ClockMixing clock_mixing);

  // The exact number of scan chains to generate (total across the design).
  // When set, this takes priority over max_length/max_chains inference.
  void setChainCount(uint64_t chain_count);
  const std::optional<uint64_t>& getChainCount() const;

  // The max length in bits that a scan chain can have
  void setMaxLength(uint64_t max_length);
  const std::optional<uint64_t>& getMaxLength() const;

  // The max number of scan chains (total across the design) to generate
  void setMaxChains(uint64_t max_chains);
  const std::optional<uint64_t>& getMaxChains() const;

  // Max allowed chain length imbalance (percent). Constraint:
  // max(bits)/min(bits) <= 1 + max_imbalance_percent/100.
  void setMaxImbalancePercent(double percent);
  double getMaxImbalancePercent() const;

  ClockMixing getClockMixing() const;

  void setScanOrderMetric(ScanOrderMetric metric);
  ScanOrderMetric getScanOrderMetric() const;

  void setScanOrderSolver(ScanOrderSolver solver);
  ScanOrderSolver getScanOrderSolver() const;

  void setScanOptRounds(uint64_t rounds);
  uint64_t getScanOptRounds() const;

  void setScanOptSeed(uint64_t seed);
  uint64_t getScanOptSeed() const;

  // Preferred-direction tuning: vertical movement is weighted by this factor
  // relative to horizontal movement (default 1.0).
  void setVerticalWeight(double weight);
  double getVerticalWeight() const;

  // Optional timing-aware penalty (Gupta'03-style extension):
  // Edge costs are scaled by a penalty based on timing slack at the scan-out
  // driver pin of the source cell. This is a heuristic proxy for "scan
  // stitching impact on timing-critical nets".
  //
  // Criticality is derived from slack using timing_critical_slack:
  // - If timing_critical_slack == 0: only negative slack is considered critical.
  // - Else: slack below timing_critical_slack is considered critical.
  //
  // Total scale factor = 1 + setup_weight*crit(setup_slack)
  //                       + hold_weight*crit(hold_slack).
  void setTimingWeightSetup(double weight);
  double getTimingWeightSetup() const;
  void setTimingWeightHold(double weight);
  double getTimingWeightHold() const;
  void setTimingCriticalSlack(double slack);
  double getTimingCriticalSlack() const;

  // Optional scan ordering constraints (ScanOpt-style):
  // - group constraints (members must stay together in one chain)
  // - fixed edges (directed adjacency constraints)
  // - before constraints (partial order)
  // - optional chain endpoints/names and chain assignment
  //
  // Constraints are specified by scan-cell instance names (post scan_replace).
  void clearScanOrderConstraints();
  void setDefaultGroupPriority(int priority);
  int getDefaultGroupPriority() const;
  const std::vector<ScanOrderGroupConstraint>& getScanOrderGroups() const;
  const std::vector<ScanOrderFixedEdgeConstraint>& getScanOrderFixedEdges() const;
  const std::vector<ScanOrderBeforeConstraint>& getScanOrderBeforeConstraints()
      const;

  // Optional per-chain endpoint constraints (by chain name).
  std::optional<ChainEndpoints> getChainEndpoints(
      std::string_view chain_name) const;

  // Optional explicit chain names (in file order).
  const std::vector<std::string>& getChainNames() const;

  // Optional hard assignment of scan instances to a specific chain (by name).
  std::optional<std::string_view> getAssignedChainForInstance(
      std::string_view inst_name) const;

  bool loadScanOrderConstraintsFile(const std::string& path, utl::Logger* logger);

  // Prints using logger->report the config used by Scan Architect
  void report(utl::Logger* logger) const;

  static std::string ClockMixingName(ClockMixing clock_mixing);
  static std::string ScanOrderMetricName(ScanOrderMetric metric);
  static std::string ScanOrderSolverName(ScanOrderSolver solver);

 private:
  // Exact number of chains to generate.
  std::optional<uint64_t> chain_count_;

  // The max length in bits of the scan chain
  std::optional<uint64_t> max_length_;
  // The max number of chains to generate
  std::optional<uint64_t> max_chains_;

  // How we are going to mix the clocks of the scan cells
  ClockMixing clock_mixing_;

  // How we order scan cells within each chain.
  ScanOrderMetric scan_order_metric_{ScanOrderMetric::Placement};

  // Which solver to use for scan ordering.
  ScanOrderSolver scan_order_solver_{ScanOrderSolver::Heuristic};

  // ScanOpt-style solver tuning knobs.
  uint64_t scanopt_rounds_{50};
  uint64_t scanopt_seed_{1};

  // Preferred wiring direction (vertical weighting).
  double vertical_weight_{1.0};

  // Timing-aware penalty knobs.
  double timing_weight_setup_{0.0};
  double timing_weight_hold_{0.0};
  double timing_critical_slack_{0.0};

  // ScanOpt-style constraints.
  int default_group_priority_{64};
  std::vector<ScanOrderGroupConstraint> scan_order_groups_;
  std::vector<ScanOrderFixedEdgeConstraint> scan_order_fixed_edges_;
  std::vector<ScanOrderBeforeConstraint> scan_order_before_;

  // Global chain endpoint constraints (optional).
  std::vector<std::string> chain_names_;
  std::unordered_map<std::string, ChainEndpoints> chain_endpoints_by_name_;

  // Optional instance->chain assignment constraints (resolved from file).
  std::unordered_map<std::string, std::string> instance_to_chain_name_;

  // Length balance constraint (percent).
  double max_imbalance_percent_{30.0};
};

}  // namespace dft
