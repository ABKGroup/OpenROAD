// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2023-2025, The OpenROAD Authors

#pragma once

#include <string>

#include "db_sta/dbSta.hh"
#include "odb/db.h"
#include "utl/Logger.h"

namespace dft {
class ScanReplace;
class DftConfig;
class ScanChain;

// The main DFT implementation.
//
// We can split the DFT process in 3 main steps:
//
// 1) Scan Replace: Where non scan sequential cells are replaced by the scan
// equivalent.
//
// 2) Scan Architect: We create the scan chains and decide what scan cells are
// going to be in them. This is an instance of the Bin-Packing problem where we
// have elements of different size (cells of different bits) and bins (scan
// chains) where to put them.
//
// 3 Scan Stitching: Where we take the scan chains from Scan Architect and
// connect the cells together to form the scan chains in the design. This is
// done by connecting the output of each scan scell to the scan input of the
// next scan cell, based on the order defined in Scan Architect.
//
// See:
//  VLSI Test Principles and Architectures, 2006, Chapter 2.7: Scan Design Flow
//
class Dft
{
 public:
  Dft(odb::dbDatabase* db, sta::dbSta* sta, utl::Logger* logger);
  ~Dft();

  // Pre-work for insert_dft. We collect the cells that need to be
  // scan replaced. This function doesn't mutate the design.
  void pre_dft();

  // Calls pre_dft performing scan replace and scan architect.
  // A report is going to be generated and printed followed by a rollback to
  // undo the work of scan replace. Use this command/function to iterate
  // different options like max_length
  //
  // Here we do:
  //  - Scan Replace
  //  - Scan Architect
  //  - Rollback
  //  - Shows Scan Architect report
  //
  // If verbose is true, then we show all the cells that are inside the scan
  // chains
  void reportDftPlan(bool verbose);

  // Reports the scan plan in a machine-parseable format that includes scan
  // pin locations (scan-in and scan-out) for each scan cell.
  //
  // This is intended for tooling/visualization to match the scan ordering
  // objective, which is defined on scan pin locations rather than instance
  // origins.
  void reportDftPlanPins(bool verbose);

  // Inserts the scan chains into the design. For now this just replace the
  // cells in the design with scan equivalent. This functions mutates the
  // design.
  //
  // Here we do:
  //  - Scan Replace
  //
  void scanReplace();

  // Runs the complete flow for scan insertion based on the user's settings
  //
  // Here we do:
  //  - Scan Replace
  //  - Scan Architect
  //  - Scan Insertion
  //  - Store the inserted DFT (scan chains) into odb for later optimization
  void executeDftPlan();

  // Writes a SCANDEF/DEF-style SCANCHAINS section (for ATPG/external tools)
  // based on the scan chains stored in the database by execute_dft_plan.
  void writeScandef(const std::string& path) const;

  // Buffers/splits the scan enable net to reduce fanout.
  // Returns the number of inserted buffers.
  int bufferScanEnable(const std::string& buffer_cell,
                       int max_fanout,
                       int max_levels);

  // Returns a mutable version of DftConfig
  DftConfig* getMutableDftConfig();

  // Returns a const version of DftConfig
  const DftConfig& getDftConfig() const;

  // Prints to stdout
  void reportDftConfig() const;

  // Performs scan optimizations on the netlist
  void scanOpt();

 private:
  // If we need to run pre_dft to create the internal state
  bool need_to_run_pre_dft_{true};

  // Cache for automatic exclusions (e.g., shift-register detection) so we don't
  // repeatedly re-scan the netlist during a flow.
  bool auto_exclusions_cache_valid_{false};
  bool cached_auto_exclude_shift_registers_{false};
  int cached_shift_register_min_length_{0};

  // Resets the internal state
  void reset();

  // Common function to perform scan replace and scan architect. Shared between
  // report_dft_plan and execute_dft_plan
  std::vector<std::unique_ptr<ScanChain>> scanArchitect();

  // Applies any enabled automatic exclusions (e.g., shift-register detection).
  void applyAutoExclusions();

  // Global state
  odb::dbDatabase* db_;
  sta::dbSta* sta_;
  utl::Logger* logger_;

  // Internal state
  std::unique_ptr<ScanReplace> scan_replace_;
  std::unique_ptr<DftConfig> dft_config_;
};

}  // namespace dft
