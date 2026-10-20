// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2023-2025, The OpenROAD Authors

#pragma once

#include <string>
#include <unordered_set>

namespace odb {
class dbDatabase;
}

namespace sta {
class dbSta;
}

namespace utl {
class Logger;
}

namespace dft {
class ScanArchitectConfig;
}

namespace dft::utils {

struct ShiftRegisterDetectResult
{
  std::unordered_set<std::string> instances;
  int chains = 0;
};

// Detects simple "functional shift register" chains (direct Q->D connections)
// and returns the set of instance names that participate in those chains.
//
// This is used to implement the optional "exclude shift registers" behavior
// described in the requirements/comments: shift registers do not need an
// additional scan chain through them.
ShiftRegisterDetectResult DetectShiftRegisters(odb::dbDatabase* db,
                                               sta::dbSta* sta,
                                               const ScanArchitectConfig& config,
                                               utl::Logger* logger);

}  // namespace dft::utils

