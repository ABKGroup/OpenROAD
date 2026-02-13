// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2023-2025, The OpenROAD Authors

%module dft

%{

#include "dft/Dft.hh"
#include "DftConfig.hh"
#include "ord/OpenRoad.hh"
#include "ScanArchitect.hh"
#include "ClockDomain.hh"

dft::Dft * getDft()
{
  return ord::OpenRoad::openRoad()->getDft();
}

utl::Logger* getLogger()
{
  return ord::OpenRoad::openRoad()->getLogger();
}

%}

%include "../../Exception.i"

// Enum: dft::ClockEdge
%typemap(typecheck) dft::ClockEdge {
  char *str = Tcl_GetStringFromObj($input, 0);
    if (strcasecmp(str, "RISING") == 0) {
    $1 = 1;
  } else if (strcasecmp(str, "FALLING") == 0) {
    $1 = 1;
  } else {
    $1 = 0;
  }
}

%typemap(in) dft::ClockEdge {
  char *str = Tcl_GetStringFromObj($input, 0);
  if (strcasecmp(str, "FALLING") == 0) {
    $1 = dft::ClockEdge::Falling;
  } else /* other values eliminated in typecheck */ {
    $1 = dft::ClockEdge::Rising;
  };
}

// Enum: dft::ScanArchitectConfig::ClockMixing
%typemap(typecheck) dft::ScanArchitectConfig::ClockMixing {
  char *str = Tcl_GetStringFromObj($input, 0);
    if (strcasecmp(str, "NO_MIX") == 0) {
    $1 = 1;
  } else if (strcasecmp(str, "CLOCK_MIX") == 0) {
    $1 = 1;
  } else {
    $1 = 0;
  }
}

%typemap(in) dft::ScanArchitectConfig::ClockMixing {
  char *str = Tcl_GetStringFromObj($input, 0);
  if (strcasecmp(str, "NO_MIX") == 0) {
    $1 = dft::ScanArchitectConfig::ClockMixing::NoMix;
  } else /* other values eliminated in typecheck */ {
    $1 = dft::ScanArchitectConfig::ClockMixing::ClockMix;
  };
}

// Enum: dft::ScanArchitectConfig::PolarityMode
%typemap(typecheck) dft::ScanArchitectConfig::PolarityMode {
  char *str = Tcl_GetStringFromObj($input, 0);
  if (strcasecmp(str, "MID") == 0) {
    $1 = 1;
  } else if (strcasecmp(str, "STRICT") == 0) {
    $1 = 1;
  } else if (strcasecmp(str, "NO_MIX") == 0) {  // alias for strict
    $1 = 1;
  } else {
    $1 = 0;
  }
}

%typemap(in) dft::ScanArchitectConfig::PolarityMode {
  char *str = Tcl_GetStringFromObj($input, 0);
  if (strcasecmp(str, "STRICT") == 0 || strcasecmp(str, "NO_MIX") == 0) {
    $1 = dft::ScanArchitectConfig::PolarityMode::Strict;
  } else /* other values eliminated in typecheck */ {
    $1 = dft::ScanArchitectConfig::PolarityMode::Mid;
  };
}

// Enum: dft::ScanArchitectConfig::ScanOrderMetric
%typemap(typecheck) dft::ScanArchitectConfig::ScanOrderMetric {
  char *str = Tcl_GetStringFromObj($input, 0);
  if (strcasecmp(str, "PLACEMENT") == 0) {
    $1 = 1;
  } else if (strcasecmp(str, "PIN_TO_NET") == 0) {
    $1 = 1;
  } else {
    $1 = 0;
  }
}

%typemap(in) dft::ScanArchitectConfig::ScanOrderMetric {
  char *str = Tcl_GetStringFromObj($input, 0);
  if (strcasecmp(str, "PIN_TO_NET") == 0) {
    $1 = dft::ScanArchitectConfig::ScanOrderMetric::PinToNet;
  } else /* other values eliminated in typecheck */ {
    $1 = dft::ScanArchitectConfig::ScanOrderMetric::Placement;
  };
}

// Enum: dft::ScanArchitectConfig::ScanOrderSolver
%typemap(typecheck) dft::ScanArchitectConfig::ScanOrderSolver {
  char *str = Tcl_GetStringFromObj($input, 0);
  if (strcasecmp(str, "HEURISTIC") == 0) {
    $1 = 1;
  } else if (strcasecmp(str, "SCANOPT") == 0) {
    $1 = 1;
  } else if (strcasecmp(str, "UCLA_SCANOPT") == 0) {
    $1 = 1;
  } else {
    $1 = 0;
  }
}

%typemap(in) dft::ScanArchitectConfig::ScanOrderSolver {
  char *str = Tcl_GetStringFromObj($input, 0);
  if (strcasecmp(str, "SCANOPT") == 0) {
    $1 = dft::ScanArchitectConfig::ScanOrderSolver::ScanOpt;
  } else if (strcasecmp(str, "UCLA_SCANOPT") == 0) {
    $1 = dft::ScanArchitectConfig::ScanOrderSolver::UclaScanOpt;
  } else /* other values eliminated in typecheck */ {
    $1 = dft::ScanArchitectConfig::ScanOrderSolver::Heuristic;
  };
}

%inline
%{

void report_dft_plan(bool verbose)
{
  getDft()->reportDftPlan(verbose);
}

void report_dft_plan_pins(bool verbose)
{
  getDft()->reportDftPlanPins(verbose);
}

void scan_replace()
{
  getDft()->scanReplace();
}

void execute_dft_plan()
{
  getDft()->executeDftPlan();
}

int buffer_scan_enable(const char* buffer_cell, int max_fanout, int max_levels)
{
  const std::string cell = (buffer_cell != nullptr) ? buffer_cell : "";
  return getDft()->bufferScanEnable(cell, max_fanout, max_levels);
}

void set_dft_config_max_length(int max_length)
{
  getDft()->getMutableDftConfig()->getMutableScanArchitectConfig()->setMaxLength(max_length);
}

void set_dft_config_chain_count(int chain_count)
{
  getDft()->getMutableDftConfig()->getMutableScanArchitectConfig()->setChainCount(chain_count);
}

void set_dft_config_max_chains(int max_chains)
{
  getDft()->getMutableDftConfig()->getMutableScanArchitectConfig()->setMaxChains(max_chains);
}

void set_dft_config_max_imbalance(double percent)
{
  getDft()
      ->getMutableDftConfig()
      ->getMutableScanArchitectConfig()
      ->setMaxImbalancePercent(percent);
}

void set_dft_config_clock_mixing(dft::ScanArchitectConfig::ClockMixing clock_mixing)
{
  getDft()->getMutableDftConfig()->getMutableScanArchitectConfig()->setClockMixing(clock_mixing);
}

void set_dft_config_polarity_mode(dft::ScanArchitectConfig::PolarityMode polarity_mode)
{
  getDft()->getMutableDftConfig()->getMutableScanArchitectConfig()->setPolarityMode(polarity_mode);
}

void set_dft_config_scan_order_metric(dft::ScanArchitectConfig::ScanOrderMetric metric)
{
  getDft()->getMutableDftConfig()->getMutableScanArchitectConfig()->setScanOrderMetric(metric);
}

void set_dft_config_scan_order_solver(dft::ScanArchitectConfig::ScanOrderSolver solver)
{
  getDft()->getMutableDftConfig()->getMutableScanArchitectConfig()->setScanOrderSolver(solver);
}

void set_dft_config_scanopt_rounds(int rounds)
{
  if (rounds > 0) {
    getDft()->getMutableDftConfig()->getMutableScanArchitectConfig()->setScanOptRounds(static_cast<uint64_t>(rounds));
  }
}

void set_dft_config_scanopt_seed(int seed)
{
  if (seed >= 0) {
    getDft()->getMutableDftConfig()->getMutableScanArchitectConfig()->setScanOptSeed(static_cast<uint64_t>(seed));
  }
}

void set_dft_config_scanopt_time_limit(double seconds)
{
  if (seconds >= 0.0) {
    getDft()->getMutableDftConfig()->getMutableScanArchitectConfig()->setScanOptTimeLimitSeconds(seconds);
  }
}

void set_dft_config_scanopt_temp_control(int enabled)
{
  getDft()
      ->getMutableDftConfig()
      ->getMutableScanArchitectConfig()
      ->setScanOptTempControl(enabled != 0);
}

void set_dft_config_scanopt_t_div(double t_div)
{
  if (t_div > 0.0) {
    getDft()->getMutableDftConfig()->getMutableScanArchitectConfig()->setScanOptTDiv(t_div);
  }
}

void set_dft_config_vertical_weight(double weight)
{
  getDft()->getMutableDftConfig()->getMutableScanArchitectConfig()->setVerticalWeight(weight);
}

void set_dft_config_timing_setup_weight(double weight)
{
  getDft()->getMutableDftConfig()->getMutableScanArchitectConfig()->setTimingWeightSetup(weight);
}

void set_dft_config_timing_hold_weight(double weight)
{
  getDft()->getMutableDftConfig()->getMutableScanArchitectConfig()->setTimingWeightHold(weight);
}

void set_dft_config_timing_critical_slack(double slack)
{
  getDft()->getMutableDftConfig()->getMutableScanArchitectConfig()->setTimingCriticalSlack(slack);
}

void set_dft_config_exclude_shift_registers(int enabled)
{
  getDft()
      ->getMutableDftConfig()
      ->getMutableScanArchitectConfig()
      ->setAutoExcludeShiftRegisters(enabled != 0);
}

void set_dft_config_prefer_qbar(int enabled)
{
  getDft()
      ->getMutableDftConfig()
      ->getMutableScanArchitectConfig()
      ->setPreferQbarScanOut(enabled != 0);
}

void set_dft_config_shift_register_min_length(int min_length)
{
  getDft()
      ->getMutableDftConfig()
      ->getMutableScanArchitectConfig()
      ->setShiftRegisterMinLength(min_length);
}

void set_dft_config_scan_order_constraints_file(const char* path_ptr)
{
  if (!path_ptr) {
    return;
  }
  std::string path(path_ptr);
  if (path.empty()) {
    return;
  }
  getDft()->getMutableDftConfig()->getMutableScanArchitectConfig()->loadScanOrderConstraintsFile(path, getLogger());
}

void set_dft_config_scan_signal_name_pattern(const char* signal_ptr, const char* pattern_ptr) {
  dft::ScanStitchConfig* config = getDft()->getMutableDftConfig()->getMutableScanStitchConfig();
  dft::ScanArchitectConfig* arch_config = getDft()->getMutableDftConfig()->getMutableScanArchitectConfig();
  std::string_view signal(signal_ptr), pattern(pattern_ptr);
  
  if (signal == "scan_in") {
    config->setInNamePattern(pattern);
    arch_config->setScanInNamePattern(pattern);
  } else if (signal == "scan_enable") {
    config->setEnableNamePattern(pattern);
    arch_config->setScanEnableNamePattern(pattern);
  } else if (signal == "scan_out") {
    config->setOutNamePattern(pattern);
    arch_config->setScanOutNamePattern(pattern);
  } else {
    getLogger()->error(utl::DFT, 6, "Internal error: unrecognized signal '{}' to set a pattern for", signal); 
  }
}

void set_dft_config_insert_lockup(int enable)
{
  getDft()
      ->getMutableDftConfig()
      ->getMutableScanStitchConfig()
      ->setInsertLockup(enable != 0);
}

void set_dft_config_lockup_cell_rising(const char* cell_ptr)
{
  if (!cell_ptr) {
    return;
  }
  getDft()
      ->getMutableDftConfig()
      ->getMutableScanStitchConfig()
      ->setLockupCellRising(cell_ptr);
}

void set_dft_config_lockup_cell_falling(const char* cell_ptr)
{
  if (!cell_ptr) {
    return;
  }
  getDft()
      ->getMutableDftConfig()
      ->getMutableScanStitchConfig()
      ->setLockupCellFalling(cell_ptr);
}

void set_dft_config_lockup_in_pin(const char* pin_ptr)
{
  if (!pin_ptr) {
    return;
  }
  getDft()
      ->getMutableDftConfig()
      ->getMutableScanStitchConfig()
      ->setLockupInPin(pin_ptr);
}

void set_dft_config_lockup_out_pin(const char* pin_ptr)
{
  if (!pin_ptr) {
    return;
  }
  getDft()
      ->getMutableDftConfig()
      ->getMutableScanStitchConfig()
      ->setLockupOutPin(pin_ptr);
}

void set_dft_config_lockup_clock_pin_rising(const char* pin_ptr)
{
  if (!pin_ptr) {
    return;
  }
  getDft()
      ->getMutableDftConfig()
      ->getMutableScanStitchConfig()
      ->setLockupClockPinRising(pin_ptr);
}

void set_dft_config_lockup_clock_pin_falling(const char* pin_ptr)
{
  if (!pin_ptr) {
    return;
  }
  getDft()
      ->getMutableDftConfig()
      ->getMutableScanStitchConfig()
      ->setLockupClockPinFalling(pin_ptr);
}

void set_dft_config_timing_buffer_cell(const char* cell_ptr)
{
  if (!cell_ptr) {
    return;
  }
  getDft()
      ->getMutableDftConfig()
      ->getMutableScanStitchConfig()
      ->setTimingBufferCell(cell_ptr);
}

void set_dft_config_timing_buffer_in_pin(const char* pin_ptr)
{
  if (!pin_ptr) {
    return;
  }
  getDft()
      ->getMutableDftConfig()
      ->getMutableScanStitchConfig()
      ->setTimingBufferInPin(pin_ptr);
}

void set_dft_config_timing_buffer_out_pin(const char* pin_ptr)
{
  if (!pin_ptr) {
    return;
  }
  getDft()
      ->getMutableDftConfig()
      ->getMutableScanStitchConfig()
      ->setTimingBufferOutPin(pin_ptr);
}

void report_dft_config() {
  getDft()->reportDftConfig();
}

void scan_opt()
{
  getDft()->scanOpt();
}

void write_scandef(const char* path_ptr)
{
  if (!path_ptr) {
    return;
  }
  std::string path(path_ptr);
  if (path.empty()) {
    return;
  }
  getDft()->writeScandef(path);
}

%}  // inline
