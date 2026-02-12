// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2023-2025, The OpenROAD Authors

#pragma once

#include <cstdint>
#include <optional>
#include <string>
#include <string_view>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "utl/Logger.h"

namespace dft {

class ScanArchitectConfig
{
 public:
  // TODO: support additional clock/edge mixing modes.
  enum class ClockMixing
  {
    NoMix,    // Separate scan chains per clock.
    ClockMix  // Mix flops of different clocks together.
  };

  // How scan cells of different edge polarity are handled within chains.
  enum class PolarityMode
  {
    Mid,     // Falling-edge cells are stitched before rising-edge cells per chain.
    Strict   // Do not mix falling-edge and rising-edge cells in the same chain.
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
    UclaScanOpt // UCLA ScanOptpack-010411 (reference implementation)
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

  // Scan-chain endpoint constraint. The point can be specified either as:
  //   - explicit coordinates in DBU, or
  //   - a design terminal name (top-level port or instance/pin).
  struct ChainEndpoint
  {
    enum class Type
    {
      Point,
      Term
    };

    Type type{Type::Point};
    Point point;
    std::string term;
  };

  struct ChainEndpoints
  {
    std::optional<ChainEndpoint> begin;
    std::optional<ChainEndpoint> end;
  };

  void setClockMixing(ClockMixing clock_mixing);
  void setPolarityMode(PolarityMode mode);

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
  PolarityMode getPolarityMode() const;

  void setScanOrderMetric(ScanOrderMetric metric);
  ScanOrderMetric getScanOrderMetric() const;

  void setScanOrderSolver(ScanOrderSolver solver);
  ScanOrderSolver getScanOrderSolver() const;

  // Default scan port naming patterns. These match the scan stitch defaults and
  // can be overridden via `set_dft_config -scan_*_name_pattern`.
  //
  // When no explicit chain endpoints are provided in the constraints file,
  // Scan Architect may use these patterns to derive implicit begin/end ports
  // for ordering cost evaluation.
  void setScanEnableNamePattern(std::string_view pattern);
  std::string_view getScanEnableNamePattern() const;
  void setScanInNamePattern(std::string_view pattern);
  std::string_view getScanInNamePattern() const;
  void setScanOutNamePattern(std::string_view pattern);
  std::string_view getScanOutNamePattern() const;

  void setScanOptRounds(uint64_t rounds);
  uint64_t getScanOptRounds() const;

  void setScanOptSeed(uint64_t seed);
  uint64_t getScanOptSeed() const;

  // Optional ScanOpt time limit in seconds (0 = unlimited).
  void setScanOptTimeLimitSeconds(double seconds);
  double getScanOptTimeLimitSeconds() const;

  // Optional temperature-controlled acceptance of non-improving moves (off by
  // default). This can help escape local minima at the expense of extra
  // randomness; use a fixed ScanOpt seed for reproducibility.
  void setScanOptTempControl(bool enable);
  bool getScanOptTempControl() const;

  // Temperature divisor when temp control is enabled. Larger values make
  // uphill acceptance rarer.
  void setScanOptTDiv(double t_div);
  double getScanOptTDiv() const;

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

  // Optional automatic exclusions.
  // When enabled, shift-register chains are detected and excluded from
  // scan_replace and scan planning.
  void setAutoExcludeShiftRegisters(bool enable);
  bool getAutoExcludeShiftRegisters() const;

  // Prefer using the complemented output (Qbar/QN) for scan-out when the
  // library does not explicitly tag a scan-out port. This can reduce added
  // load on functional Q nets for libraries whose scan architecture reuses
  // the functional output as scan-out.
  void setPreferQbarScanOut(bool enable);
  bool getPreferQbarScanOut() const;

  // Minimum number of sequential elements to consider a chain a "shift
  // register" for auto exclusion (must be >= 2).
  void setShiftRegisterMinLength(int min_length);
  int getShiftRegisterMinLength() const;

  // Replace the current auto-excluded instance set (used by automatic
  // exclusions such as shift-register detection).
  void clearAutoExcludedInstances();
  void setAutoExcludedInstances(std::unordered_set<std::string> instances);

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

  // Optional exclusion list: instances listed here are skipped by scan_replace
  // (left as functional flops) and ignored by scan planning.
  bool isInstanceExcluded(std::string_view inst_name) const;
  bool isMasterExcluded(std::string_view master_name) const;
  bool isInstanceExcluded(std::string_view inst_name,
                          std::string_view master_name) const;
  const std::unordered_set<std::string>& getExcludedInstances() const;
  const std::unordered_set<std::string>& getAutoExcludedInstances() const;
  const std::vector<std::string>& getExcludedInstancePatterns() const;
  const std::vector<std::string>& getExcludedMasterPatterns() const;

  bool loadScanOrderConstraintsFile(const std::string& path, utl::Logger* logger);

  // Prints using logger->report the config used by Scan Architect
  void report(utl::Logger* logger) const;

  static std::string ClockMixingName(ClockMixing clock_mixing);
  static std::string PolarityModeName(PolarityMode mode);
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
  ClockMixing clock_mixing_{ClockMixing::NoMix};
  PolarityMode polarity_mode_{PolarityMode::Strict};

  // How we order scan cells within each chain.
  ScanOrderMetric scan_order_metric_{ScanOrderMetric::Placement};

  // Which solver to use for scan ordering.
  ScanOrderSolver scan_order_solver_{ScanOrderSolver::ScanOpt};

  // ScanOpt-style solver tuning knobs.
  uint64_t scanopt_rounds_{500000};
  uint64_t scanopt_seed_{1};
  double scanopt_time_limit_seconds_{15.0};
  bool scanopt_temp_control_{false};
  double scanopt_t_div_{100.0};

  // Scan signal name patterns (for implicit chain endpoints).
  std::string scan_enable_name_pattern_ = "scan_enable_{}";
  std::string scan_in_name_pattern_ = "scan_in_{}";
  std::string scan_out_name_pattern_ = "scan_out_{}";

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

  // Optional exclusion list (resolved from constraints file).
  std::unordered_set<std::string> excluded_instances_;
  std::vector<std::string> excluded_instance_patterns_;
  std::vector<std::string> excluded_master_patterns_;

  // Auto exclusion list (computed from the design).
  std::unordered_set<std::string> auto_excluded_instances_;

  // Automatic shift register exclusion.
  bool auto_exclude_shift_registers_{false};
  int shift_register_min_length_{4};

  // Prefer Qbar/QN as scan-out when no scan-out is tagged in Liberty.
  bool prefer_qbar_scan_out_{false};

  // Length balance constraint (percent).
  double max_imbalance_percent_{30.0};
};

}  // namespace dft
