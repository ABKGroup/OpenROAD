// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors

#pragma once

#include <cstddef>
#include <cstdint>
#include <string_view>
#include <utility>
#include <vector>

namespace dft {

struct UclaScanOptParams
{
  uint64_t seed = 1;

  // UCLA ScanOpt params (Optimizer1).
  uint64_t major_loops = 100;
  uint64_t n_descents = 5;
  uint64_t kick_move = 15;
  uint64_t n_near = 20;
  bool only_2opt = false;
  bool temp_control = false;
  bool zero_temp = true;
};

// Runs UCLA ScanOpt (ScanOptpack-010411) on a single scan chain with fixed
// begin/end points, returning an ordering over the input nodes.
//
// The cost model is the original UCLA ScanOpt asymmetric Manhattan:
//   dist(i->j) = |out(i) - in(j)|_1
std::vector<std::size_t> UclaScanOptOrder(
    const std::vector<std::string_view>& names,
    const std::vector<std::pair<int, int>>& scan_in_pts,
    const std::vector<std::pair<int, int>>& scan_out_pts,
    const std::pair<int, int>& begin,
    const std::pair<int, int>& end,
    const UclaScanOptParams& params);

}  // namespace dft

