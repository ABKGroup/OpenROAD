// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2023-2025, The OpenROAD Authors

#pragma once

#include <cstdint>
#include <optional>
#include <string>
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
    int priority = 64;  // lower runs earlier
    std::vector<std::string> inst_names;
  };

  struct ScanOrderFixedEdgeConstraint
  {
    std::string from_inst;
    std::string to_inst;
  };

  void setClockMixing(ClockMixing clock_mixing);

  // The exact number of scan chains to generate (per clock-edge pair in NoMix,
  // total in ClockMix). When set, this takes priority over max_length/max_chains
  // inference.
  void setChainCount(uint64_t chain_count);
  const std::optional<uint64_t>& getChainCount() const;

  // The max length in bits that a scan chain can have
  void setMaxLength(uint64_t max_length);
  const std::optional<uint64_t>& getMaxLength() const;

  // The max number of scan chains (per clock in NoMixe mode) to generate
  void setMaxChains(uint64_t max_chains);
  const std::optional<uint64_t>& getMaxChains() const;

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
  //
  // Constraints are specified by scan-cell instance names (post scan_replace).
  void clearScanOrderConstraints();
  void setDefaultGroupPriority(int priority);
  int getDefaultGroupPriority() const;
  const std::vector<ScanOrderGroupConstraint>& getScanOrderGroups() const;
  const std::vector<ScanOrderFixedEdgeConstraint>& getScanOrderFixedEdges() const;
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
};

}  // namespace dft
