// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors

#include "UclaScanOpt.hh"

#include <algorithm>
#include <limits>
#include <mutex>
#include <stdexcept>
#include <string>

#include <ABKCommon/abkseed.h>
#include <ScanOpt/optimizer.h>
#include <ScanOpt/scanTourDZ.h>

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

  void restorePath(const std::vector<unsigned>& path)
  {
    _path = path;
    _pathValid = true;
    computeInverse();
  }

  std::vector<unsigned>& mutablePath() { return _path; }
};

unsigned clampToUnsigned(uint64_t v)
{
  return static_cast<unsigned>(
      std::min<uint64_t>(v, std::numeric_limits<unsigned>::max()));
}

unsigned clampToDeterministicSeed(uint64_t v)
{
  const unsigned s = clampToUnsigned(v);
  // ABKCommon uses UINT_MAX as the sentinel for "no explicit seed".
  if (s == std::numeric_limits<unsigned>::max()) {
    return std::numeric_limits<unsigned>::max() - 1;
  }
  return s;
}

void configureSeedHandlerOnce(uint64_t seed)
{
  static std::once_flag seed_once;
  std::call_once(seed_once, [seed]() {
    // UCLApack's ABKCommon SeedHandler writes a `seeds.out` lock/log file in the
    // CWD by default. Disable that so embedded use (and parallel runs) don't
    // create/contend on a global file.
    SeedHandler::turnOffLogging();

    // Force a deterministic external seed so any ABKCommon RNGs created without
    // an explicit seed (or using multipartite locIdent seeds) remain
    // deterministic across runs.
    SeedHandler::overrideExternalSeed(clampToDeterministicSeed(seed));
  });
}

double optimizeLevelSafe(FlatScanChain& chain,
                         RandomRawUnsigned& randuns,
                         const abkscanopt::Optimizer::Params& params)
{
  abkscanopt::ScanTourDZ::OptParams opt_params;
  opt_params.nDescents = params.nDescents;
  opt_params.kickMove = params.kickMove;
  opt_params.only2Opt = params.only2Opt;
  opt_params.temp_control = params.temp_control;
  opt_params.bZeroFix = true;

  const unsigned collapse_size = chain.getNumCells() - 1;
  const unsigned in_idx = chain.getPath()[0];
  const unsigned out_idx = chain.getPath()[collapse_size];

  unsigned nnear = params.nnear;
  if (nnear > collapse_size - 1) {
    nnear = collapse_size - 1;
  }

  abkscanopt::ScanTourDZ tour(chain, nnear, randuns, opt_params);
  if (!chain.isSubpath()) {
    tour.scoOptAll();

    unsigned lesser = in_idx;
    unsigned greater = out_idx;
    if (lesser > greater) {
      std::swap(lesser, greater);
    }

    std::vector<unsigned> collapsed_path = tour.getPath();
    auto zero_it = std::find(collapsed_path.begin(), collapsed_path.end(), 0U);
    if (zero_it == collapsed_path.end()) {
      throw std::runtime_error("UclaScanOptOrder: missing distinguished zero");
    }

    // UCLApack's Optimizer1 expects the distinguished-zero cell to remain in
    // position 0 of the tour. Some cases (notably when there are no partial
    // order constraints) can return a rotated tour where the zero cell is not
    // first, which later triggers an ABKCommon fatal error while "uncollapsing"
    // the endpoints. Canonicalize by rotating the cycle so the zero cell is
    // always first.
    if (zero_it != collapsed_path.begin()) {
      std::rotate(collapsed_path.begin(), zero_it, collapsed_path.end());
    }

    std::vector<unsigned>& path = chain.mutablePath();
    for (unsigned i = 1; i < collapse_size; ++i) {
      const unsigned idx = collapsed_path[i];
      if (idx == 0U) {
        throw std::runtime_error(
            "UclaScanOptOrder: internal distinguished zero after rotation");
      }
      path[i] = abkscanopt::ScanCells::uncollapseIndex(idx, lesser, greater);
    }

    chain.computeInverse();
  }

  return tour.scoTourCost();
}

void optimizeChainSafe(FlatScanChain& chain,
                       RandomRawUnsigned& randuns,
                       const abkscanopt::Optimizer::Params& params)
{
  const unsigned major_loops = std::max(1U, params.majorLoops);

  std::vector<unsigned> best_path = chain.getPath();
  double best_cost = std::numeric_limits<double>::infinity();

  std::vector<unsigned> last_path;
  double last_cost = std::numeric_limits<double>::infinity();
  bool have_last = false;

  for (unsigned loop = 0; loop < major_loops; ++loop) {
    const double cost = optimizeLevelSafe(chain, randuns, params);

    if (cost < best_cost) {
      best_cost = cost;
      best_path = chain.getPath();
    }

    if (params.zeroTemp) {
      if (!have_last || cost < last_cost) {
        have_last = true;
        last_cost = cost;
        last_path = chain.getPath();
      } else {
        chain.restorePath(last_path);
      }
    }
  }

  if (!best_path.empty()) {
    chain.restorePath(best_path);
  }
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

  configureSeedHandlerOnce(params.seed);
  RandomRawUnsigned randuns(clampToDeterministicSeed(params.seed));

  abkscanopt::Optimizer::Params p;
  p.majorLoops = clampToUnsigned(params.major_loops);
  p.nDescents = clampToUnsigned(params.n_descents);
  p.kickMove = clampToUnsigned(params.kick_move);
  p.nnear = clampToUnsigned(params.n_near);
  p.only2Opt = params.only_2opt;
  p.temp_control = params.temp_control;
  p.zeroTemp = params.zero_temp;

  optimizeChainSafe(chain, randuns, p);

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
