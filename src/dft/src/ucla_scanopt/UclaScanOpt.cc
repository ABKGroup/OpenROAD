// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors

#include "UclaScanOpt.hh"

#include <limits>
#include <stdexcept>
#include <string>

#include <ScanOpt/optimizer1.h>

namespace dft {
namespace {

class FlatScanChain : public abkscanopt::ScanChain
{
 public:
  struct CellData
  {
    int in_x = 0;
    int in_y = 0;
    int out_x = 0;
    int out_y = 0;
    std::string name;
  };

  FlatScanChain(const std::vector<CellData>& cells,
                unsigned begin_idx,
                unsigned end_idx)
      : abkscanopt::ScanChain()
  {
    for (const auto& c : cells) {
      _addCell(c.in_x, c.in_y, c.out_x, c.out_y, c.name);
    }

    _normalizePO();

    // Fix endpoints via legalins/legalouts; Optimizer1 keeps path endpoints
    // fixed to the initial path endpoints.
    _addLegalin(static_cast<int>(begin_idx));
    _addLegalout(static_cast<int>(end_idx));

    _checkIndexRanges();
    _checkInOutLegality();

    _path.clear();
    _path.reserve(cells.size());
    for (unsigned i = 0; i < cells.size(); ++i) {
      _path.push_back(i);
    }
    _pathValid = true;
    computeInverse();
  }
};

unsigned clampToUnsigned(uint64_t v)
{
  return static_cast<unsigned>(
      std::min<uint64_t>(v, std::numeric_limits<unsigned>::max()));
}

}  // namespace

std::vector<std::size_t> UclaScanOptOrder(
    const std::vector<std::string_view>& names,
    const std::vector<std::pair<int, int>>& scan_in_pts,
    const std::vector<std::pair<int, int>>& scan_out_pts,
    const std::pair<int, int>& begin,
    const std::pair<int, int>& end,
    const UclaScanOptParams& params)
{
  const std::size_t n = names.size();
  if (scan_in_pts.size() != n || scan_out_pts.size() != n) {
    throw std::runtime_error("UclaScanOptOrder: size mismatch");
  }
  if (n <= 1) {
    return n == 1 ? std::vector<std::size_t>{0} : std::vector<std::size_t>{};
  }

  std::vector<FlatScanChain::CellData> cells;
  cells.reserve(n + 2);

  // Begin dummy: its out point becomes the distinguished zero OUT.
  {
    FlatScanChain::CellData c;
    c.in_x = begin.first;
    c.in_y = begin.second;
    c.out_x = begin.first;
    c.out_y = begin.second;
    c.name = "__begin__";
    cells.push_back(std::move(c));
  }

  // Real scan cells.
  for (std::size_t i = 0; i < n; ++i) {
    FlatScanChain::CellData c;
    c.in_x = scan_in_pts[i].first;
    c.in_y = scan_in_pts[i].second;
    c.out_x = scan_out_pts[i].first;
    c.out_y = scan_out_pts[i].second;
    c.name = std::string(names[i]);
    cells.push_back(std::move(c));
  }

  // End dummy: its in point becomes the distinguished zero IN.
  {
    FlatScanChain::CellData c;
    c.in_x = end.first;
    c.in_y = end.second;
    c.out_x = end.first;
    c.out_y = end.second;
    c.name = "__end__";
    cells.push_back(std::move(c));
  }

  const unsigned begin_idx = 0;
  const unsigned end_idx = static_cast<unsigned>(cells.size() - 1);

  FlatScanChain chain(cells, begin_idx, end_idx);

  RandomRawUnsigned randuns(clampToUnsigned(params.seed));

  abkscanopt::Optimizer::Params p;
  p.majorLoops = clampToUnsigned(params.major_loops);
  p.nDescents = clampToUnsigned(params.n_descents);
  p.kickMove = clampToUnsigned(params.kick_move);
  p.nnear = clampToUnsigned(params.n_near);
  p.only2Opt = params.only_2opt;
  p.temp_control = params.temp_control;
  p.zeroTemp = params.zero_temp;

  // Runs optimization in ctor and updates chain path.
  abkscanopt::Optimizer1 opt(chain, randuns, p);
  (void) opt;

  const auto& path = chain.getPath();
  if (path.size() != n + 2) {
    throw std::runtime_error("UclaScanOptOrder: unexpected path length");
  }
  if (path.front() != begin_idx || path.back() != end_idx) {
    throw std::runtime_error("UclaScanOptOrder: endpoints moved unexpectedly");
  }

  std::vector<std::size_t> order;
  order.reserve(n);
  for (std::size_t k = 1; k + 1 < path.size(); ++k) {
    const unsigned idx = path[k];
    if (idx == begin_idx || idx == end_idx) {
      throw std::runtime_error("UclaScanOptOrder: internal endpoint");
    }
    // Real scan cells are offset by +1 due to the begin dummy.
    order.push_back(static_cast<std::size_t>(idx - 1));
  }

  if (order.size() != n) {
    throw std::runtime_error("UclaScanOptOrder: ordering size mismatch");
  }

  return order;
}

}  // namespace dft
