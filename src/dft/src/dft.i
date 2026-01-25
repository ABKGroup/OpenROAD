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
  } else if (strcasecmp(str, "MIN_FEEDTHROUGH") == 0
             || strcasecmp(str, "MIN_FEEDTHROUGH_DP") == 0) {
    $1 = 1;
  } else {
    $1 = 0;
  }
}

%typemap(in) dft::ScanArchitectConfig::ScanOrderSolver {
  char *str = Tcl_GetStringFromObj($input, 0);
  if (strcasecmp(str, "SCANOPT") == 0) {
    $1 = dft::ScanArchitectConfig::ScanOrderSolver::ScanOpt;
  } else if (strcasecmp(str, "MIN_FEEDTHROUGH") == 0
             || strcasecmp(str, "MIN_FEEDTHROUGH_DP") == 0) {
    $1 = dft::ScanArchitectConfig::ScanOrderSolver::MinFeedthrough;
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

void scan_replace()
{
  getDft()->scanReplace();
}

void execute_dft_plan()
{
  getDft()->executeDftPlan();
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

void set_dft_config_clock_mixing(dft::ScanArchitectConfig::ClockMixing clock_mixing)
{
  getDft()->getMutableDftConfig()->getMutableScanArchitectConfig()->setClockMixing(clock_mixing);
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
  std::string_view signal(signal_ptr), pattern(pattern_ptr);
  
  if (signal == "scan_in") {
    config->setInNamePattern(pattern);
  } else if (signal == "scan_enable") {
    config->setEnableNamePattern(pattern);
  } else if (signal == "scan_out") {
    config->setOutNamePattern(pattern);
  } else {
    getLogger()->error(utl::DFT, 6, "Internal error: unrecognized signal '{}' to set a pattern for", signal); 
  }
}

void report_dft_config() {
  getDft()->reportDftConfig();
}

void scan_opt()
{
  getDft()->scanOpt();
}

%}  // inline
