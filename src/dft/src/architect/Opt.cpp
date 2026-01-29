// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2024-2025, The OpenROAD Authors

#include "Opt.hh"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <limits>
#include <map>
#include <memory>
#include <optional>
#include <random>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <vector>

#include "ScanCell.hh"
#include "boost/geometry/core/access.hpp"
#include "boost/geometry/core/cs.hpp"
#include "boost/geometry/geometries/point.hpp"
#include "boost/geometry/geometry.hpp"  // NOLINT(misc-include-cleaner)
#include "boost/geometry/index/parameters.hpp"
#include "boost/geometry/index/predicates.hpp"
#include "boost/geometry/index/rtree.hpp"
#include "odb/dbShape.h"
#include "odb/geom.h"
#include "utl/Logger.h"

namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

namespace dft {

namespace {
constexpr int64_t kInfDistance = std::numeric_limits<int64_t>::max() / 8;
constexpr int32_t kScanOptLargeCost = std::numeric_limits<int32_t>::max() / 16;

// Manhattan nearest-neighbor scan (O(n^2)) is exact but quadratic; keep it for
// moderately sized chains.
constexpr std::size_t kQuadraticHeuristicMaxCells = 10000;

// Euclidean rtree nearest query candidate count (used as a fallback for large
// chains where quadratic heuristics are too expensive).
constexpr std::size_t kNearestCandidateCount = 256;

// Bound the runtime of 2-opt (O(n^2)). The quadratic local-search is only used
// on smaller chains.
constexpr std::size_t kTwoOptMaxCellsFor1Pass = 12000;
constexpr std::size_t kTwoOptMaxCellsFor2Passes = 8000;
constexpr std::size_t kTwoOptMaxCellsFor3Passes = 4000;

constexpr std::size_t kFarthestInsertionMaxCells = kQuadraticHeuristicMaxCells;

// Bound memory used by ScanOpt-style cost matrices.
constexpr std::size_t kScanOptMaxMatrixCells = 6000;

// Only run the more expensive local-search operators on smaller chains.
constexpr std::size_t kScanOptTwoOptMaxCells = 2500;
constexpr std::size_t kScanOptSwapMaxCells = 4000;

int64_t manhattanDist(const odb::Point& a,
                      const odb::Point& b,
                      double vertical_weight);

double timingCriticalityFromSlack(float slack, double critical_slack)
{
  if (critical_slack <= 0.0) {
    return slack < 0.0F ? 1.0 : 0.0;
  }
  const double ratio
      = (critical_slack - static_cast<double>(slack)) / critical_slack;
  return std::clamp(ratio, 0.0, 1.0);
}

double timingMultiplierForCell(const ScanArchitectConfig& config,
                               const ScanCell& cell)
{
  const double setup_w = config.getTimingWeightSetup();
  const double hold_w = config.getTimingWeightHold();
  if (setup_w == 0.0 && hold_w == 0.0) {
    return 1.0;
  }
  if (!cell.hasTimingSlacks()) {
    return 1.0;
  }

  const double critical_slack = config.getTimingCriticalSlack();
  const double setup_c
      = timingCriticalityFromSlack(cell.getSetupSlack(), critical_slack);
  const double hold_c
      = timingCriticalityFromSlack(cell.getHoldSlack(), critical_slack);
  const double factor = 1.0 + setup_w * setup_c + hold_w * hold_c;
  return std::clamp(factor, 0.0, 1.0e6);
}

int64_t scaleEdgeCost(int64_t base_cost, double factor)
{
  if (factor == 1.0) {
    return base_cost;
  }
  if (base_cost <= 0) {
    return base_cost;
  }
  const double scaled = static_cast<double>(base_cost) * factor;
  if (scaled >= static_cast<double>(kInfDistance)) {
    return kInfDistance;
  }
  return static_cast<int64_t>(std::llround(scaled));
}

std::vector<std::size_t> minFeedthroughRowSweepOrder(
    const std::vector<odb::Point>& origins,
    const std::vector<std::string_view>& names,
    double vertical_weight)
{
  const std::size_t n = origins.size();
  std::map<int, std::vector<std::size_t>> rows_by_y;
  for (std::size_t i = 0; i < n; ++i) {
    rows_by_y[origins[i].y()].push_back(i);
  }

  struct Row
  {
    int y = 0;
    std::vector<std::size_t> idxs;
    std::size_t left = 0;
    std::size_t right = 0;
    int64_t span = 0;
  };

  std::vector<Row> rows;
  rows.reserve(rows_by_y.size());
  for (auto& [y, idxs] : rows_by_y) {
    std::sort(idxs.begin(), idxs.end(), [&](std::size_t a, std::size_t b) {
      const int ax = origins[a].x();
      const int bx = origins[b].x();
      if (ax != bx) {
        return ax < bx;
      }
      return names[a] < names[b];
    });

    Row row;
    row.y = y;
    row.idxs = std::move(idxs);
    row.left = row.idxs.front();
    row.right = row.idxs.back();
    row.span = std::abs(static_cast<int64_t>(origins[row.right].x())
                        - static_cast<int64_t>(origins[row.left].x()));
    rows.push_back(std::move(row));
  }

  if (rows.empty()) {
    return {};
  }

  const std::size_t m = rows.size();
  // Direction 0: traverse left->right (entry=left, exit=right)
  // Direction 1: traverse right->left (entry=right, exit=left)
  const auto entry_idx = [&](std::size_t row_i, int dir) -> std::size_t {
    return dir == 0 ? rows[row_i].left : rows[row_i].right;
  };
  const auto exit_idx = [&](std::size_t row_i, int dir) -> std::size_t {
    return dir == 0 ? rows[row_i].right : rows[row_i].left;
  };

  std::vector<std::array<int64_t, 2>> dp(m, {kInfDistance, kInfDistance});
  std::vector<std::array<int, 2>> parent(m, {-1, -1});

  for (int dir = 0; dir < 2; ++dir) {
    dp[0][dir] = rows[0].span;
  }

  for (std::size_t i = 1; i < m; ++i) {
    for (int dir = 0; dir < 2; ++dir) {
      const std::size_t cur_entry = entry_idx(i, dir);
      int64_t best = kInfDistance;
      int best_prev = -1;
      for (int prev_dir = 0; prev_dir < 2; ++prev_dir) {
        const std::size_t prev_exit = exit_idx(i - 1, prev_dir);
        const int64_t link
            = manhattanDist(origins[prev_exit], origins[cur_entry], vertical_weight);
        const int64_t cand = dp[i - 1][prev_dir] + link + rows[i].span;
        if (cand < best || (cand == best && prev_dir < best_prev)) {
          best = cand;
          best_prev = prev_dir;
        }
      }
      dp[i][dir] = best;
      parent[i][dir] = best_prev;
    }
  }

  int best_dir = 0;
  if (dp[m - 1][1] < dp[m - 1][0]) {
    best_dir = 1;
  }

  std::vector<int> dirs(m, 0);
  int dir = best_dir;
  for (std::size_t i = m; i-- > 0;) {
    dirs[i] = dir;
    if (i == 0) {
      break;
    }
    dir = parent[i][dir] >= 0 ? parent[i][dir] : 0;
  }

  std::vector<std::size_t> order;
  order.reserve(n);
  for (std::size_t i = 0; i < m; ++i) {
    if (dirs[i] == 0) {
      order.insert(order.end(), rows[i].idxs.begin(), rows[i].idxs.end());
    } else {
      order.insert(order.end(), rows[i].idxs.rbegin(), rows[i].idxs.rend());
    }
  }

  return order;
}

int64_t manhattanDist(const odb::Point& a,
                      const odb::Point& b,
                      double vertical_weight)
{
  const int64_t dx = std::abs(static_cast<int64_t>(a.x())
                              - static_cast<int64_t>(b.x()));
  const int64_t dy = std::abs(static_cast<int64_t>(a.y())
                              - static_cast<int64_t>(b.y()));
  const int64_t wy = static_cast<int64_t>(std::llround(vertical_weight * dy));
  return dx + wy;
}

int64_t manhattanPointToRectDist(const odb::Point& p,
                                 const odb::Rect& r,
                                 double vertical_weight)
{
  int64_t dx = 0;
  if (p.x() < r.xMin()) {
    dx = static_cast<int64_t>(r.xMin()) - p.x();
  } else if (p.x() > r.xMax()) {
    dx = static_cast<int64_t>(p.x()) - r.xMax();
  }

  int64_t dy = 0;
  if (p.y() < r.yMin()) {
    dy = static_cast<int64_t>(r.yMin()) - p.y();
  } else if (p.y() > r.yMax()) {
    dy = static_cast<int64_t>(p.y()) - r.yMax();
  }
  const int64_t wy = static_cast<int64_t>(std::llround(vertical_weight * dy));
  return dx + wy;
}

odb::Point scanPinLocation(const ScanPin& pin, const odb::Point& fallback)
{
  return std::visit(
      overloaded{[&](odb::dbITerm* iterm) -> odb::Point {
                   if (iterm == nullptr) {
                     return fallback;
                   }
                   int x = 0;
                   int y = 0;
                   if (iterm->getAvgXY(&x, &y)) {
                     return odb::Point(x, y);
                   }
                   odb::dbInst* inst = iterm->getInst();
                   return inst ? inst->getLocation() : fallback;
                 },
                 [&](odb::dbBTerm* bterm) -> odb::Point {
                   if (bterm == nullptr) {
                     return fallback;
                   }
                   int x = 0;
                   int y = 0;
                   if (bterm->getFirstPinLocation(x, y)) {
                     return odb::Point(x, y);
                   }
                   return fallback;
                 }},
      pin.getValue());
}

struct NetAccessGeometry
{
  std::vector<odb::Rect> boxes;
  std::vector<odb::Point> terminals;
};

NetAccessGeometry buildNetAccessGeometry(odb::dbNet* net)
{
  NetAccessGeometry geom;
  if (net == nullptr) {
    return geom;
  }

  for (odb::dbGuide* guide : net->getGuides()) {
    geom.boxes.push_back(guide->getBox());
  }

  if (geom.boxes.empty()) {
    if (odb::dbWire* wire = net->getWire()) {
      odb::dbWireShapeItr itr;
      odb::dbShape shape;
      for (itr.begin(wire); itr.next(shape);) {
        geom.boxes.push_back(shape.getBox());
      }
    }
  }

  if (geom.boxes.empty()) {
    for (odb::dbITerm* iterm : net->getITerms()) {
      int x = 0;
      int y = 0;
      if (iterm->getAvgXY(&x, &y)) {
        geom.terminals.emplace_back(x, y);
      } else if (odb::dbInst* inst = iterm->getInst()) {
        geom.terminals.emplace_back(inst->getLocation());
      }
    }
    for (odb::dbBTerm* bterm : net->getBTerms()) {
      int x = 0;
      int y = 0;
      if (bterm->getFirstPinLocation(x, y)) {
        geom.terminals.emplace_back(x, y);
      }
    }
  }

  return geom;
}

int64_t pinToNetDistance(const odb::Point& pin,
                         const NetAccessGeometry& geom,
                         double vertical_weight)
{
  int64_t best = kInfDistance;
  for (const odb::Rect& box : geom.boxes) {
    best = std::min(best, manhattanPointToRectDist(pin, box, vertical_weight));
    if (best == 0) {
      return 0;
    }
  }
  for (const odb::Point& term : geom.terminals) {
    best = std::min(best, manhattanDist(pin, term, vertical_weight));
    if (best == 0) {
      return 0;
    }
  }
  return best == kInfDistance ? 0 : best;
}

struct ScanOptMatrix
{
  std::size_t n = 0;
  std::vector<int32_t> costs;

  int32_t get(std::size_t src, std::size_t dst) const
  {
    return costs[src * n + dst];
  }
};

int64_t scanOptPathCost(const ScanOptMatrix& m,
                        const std::vector<std::size_t>& order)
{
  int64_t total = 0;
  for (std::size_t i = 1; i < order.size(); ++i) {
    total += static_cast<int64_t>(m.get(order[i - 1], order[i]));
  }
  return total;
}

std::vector<std::size_t> scanOptGreedyOrder(
    const ScanOptMatrix& m,
    std::size_t start,
    const std::vector<std::string_view>& names,
    utl::Logger* logger)
{
  const std::size_t n = m.n;
  std::vector<std::size_t> order;
  order.reserve(n);
  std::vector<bool> used(n, false);

  std::size_t cur = start;
  used[cur] = true;
  order.push_back(cur);

  while (order.size() < n) {
    bool found = false;
    std::size_t best = 0;
    int32_t best_cost = kScanOptLargeCost;

    for (std::size_t i = 0; i < n; ++i) {
      if (used[i]) {
        continue;
      }
      const int32_t cost = m.get(cur, i);
      if (!found || cost < best_cost
          || (cost == best_cost && names[i] < names[best])) {
        found = true;
        best = i;
        best_cost = cost;
      }
    }

    if (!found) {
      logger->error(utl::DFT, 75, "Couldn't find next scan cell to order");
    }

    used[best] = true;
    order.push_back(best);
    cur = best;
  }

  return order;
}

bool scanOptAdjacentSwapImprove(const ScanOptMatrix& m,
                                const std::vector<std::string_view>& names,
                                std::vector<std::size_t>& order)
{
  bool improved = false;
  for (std::size_t pos = 1; pos + 1 < order.size(); ++pos) {
    const std::size_t prev = order[pos - 1];
    const std::size_t a = order[pos];
    const std::size_t b = order[pos + 1];
    const bool has_next = (pos + 2 < order.size());
    const std::size_t next = has_next ? order[pos + 2] : 0;

    const int64_t before = static_cast<int64_t>(m.get(prev, a))
                           + static_cast<int64_t>(m.get(a, b))
                           + (has_next ? static_cast<int64_t>(m.get(b, next))
                                       : 0);
    const int64_t after = static_cast<int64_t>(m.get(prev, b))
                          + static_cast<int64_t>(m.get(b, a))
                          + (has_next ? static_cast<int64_t>(m.get(a, next))
                                      : 0);

    if (after < before || (after == before && names[b] < names[a])) {
      std::swap(order[pos], order[pos + 1]);
      improved = true;
    }
  }
  return improved;
}

bool scanOptBestRelocateMove(const ScanOptMatrix& m,
                             const std::vector<std::string_view>& names,
                             std::vector<std::size_t>& order)
{
  const std::size_t n = order.size();
  if (n < 3) {
    return false;
  }

  bool found = false;
  int64_t best_delta = 0;
  std::size_t best_from = 0;
  std::size_t best_after = 0;

  // Keep the start fixed for determinism (consistent with existing heuristics).
  for (std::size_t from = 1; from < n; ++from) {
    const std::size_t node = order[from];

    const bool has_prev = (from > 0);
    const bool has_next = (from + 1 < n);
    const std::size_t prev = has_prev ? order[from - 1] : 0;
    const std::size_t next = has_next ? order[from + 1] : 0;

    int64_t remove_delta = 0;
    if (has_prev && has_next) {
      remove_delta = static_cast<int64_t>(m.get(prev, next))
                     - static_cast<int64_t>(m.get(prev, node))
                     - static_cast<int64_t>(m.get(node, next));
    } else if (has_prev) {
      remove_delta = -static_cast<int64_t>(m.get(prev, node));
    } else if (has_next) {
      remove_delta = -static_cast<int64_t>(m.get(node, next));
    }

    for (std::size_t after = 0; after < n; ++after) {
      if (after == from || after + 1 == from) {
        continue;  // no-op or invalid
      }
      if (after >= n - 1 && after != n - 1) {
        continue;
      }

      const std::size_t ins_prev = order[after];
      const bool ins_has_next = (after + 1 < n);
      const std::size_t ins_next = ins_has_next ? order[after + 1] : 0;

      int64_t insert_delta = 0;
      if (ins_has_next) {
        insert_delta = static_cast<int64_t>(m.get(ins_prev, node))
                       + static_cast<int64_t>(m.get(node, ins_next))
                       - static_cast<int64_t>(m.get(ins_prev, ins_next));
      } else {
        insert_delta = static_cast<int64_t>(m.get(ins_prev, node));
      }

      const int64_t delta = remove_delta + insert_delta;
      if (delta < best_delta
          || (delta == best_delta && found
              && names[node] < names[order[best_from]])) {
        found = true;
        best_delta = delta;
        best_from = from;
        best_after = after;
      }
    }
  }

  if (!found || best_delta >= 0) {
    return false;
  }

  const std::size_t node = order[best_from];
  order.erase(order.begin() + static_cast<std::ptrdiff_t>(best_from));
  std::size_t insert_index = 0;
  if (best_after < best_from) {
    insert_index = best_after + 1;
  } else {
    insert_index = best_after;
  }
  order.insert(order.begin() + static_cast<std::ptrdiff_t>(insert_index), node);

  return true;
}

bool scanOptBestSwapMove(const ScanOptMatrix& m,
                         const std::vector<std::string_view>& names,
                         std::vector<std::size_t>& order)
{
  const std::size_t n = order.size();
  if (n < 4) {
    return false;
  }

  bool found = false;
  int64_t best_delta = 0;
  std::size_t best_i = 0;
  std::size_t best_j = 0;

  // Keep the start fixed (i starts at 1).
  for (std::size_t i = 1; i + 1 < n; ++i) {
    for (std::size_t j = i + 1; j < n; ++j) {
      if (i == j) {
        continue;
      }

      const std::size_t ai = order[i];
      const std::size_t aj = order[j];

      const bool has_pi = (i > 0);
      const bool has_ni = (i + 1 < n);
      const bool has_pj = (j > 0);
      const bool has_nj = (j + 1 < n);

      const std::size_t pi = has_pi ? order[i - 1] : 0;
      const std::size_t ni = has_ni ? order[i + 1] : 0;
      const std::size_t pj = has_pj ? order[j - 1] : 0;
      const std::size_t nj = has_nj ? order[j + 1] : 0;

      int64_t before = 0;
      int64_t after = 0;

      if (has_pi) {
        before += m.get(pi, ai);
        after += m.get(pi, aj);
      }
      if (has_ni) {
        if (i + 1 == j) {
          before += m.get(ai, aj);
          after += m.get(aj, ai);
        } else {
          before += m.get(ai, ni);
          after += m.get(aj, ni);
        }
      }

      if (has_pj) {
        if (j - 1 != i) {
          before += m.get(pj, aj);
          after += m.get(pj, ai);
        }
      }
      if (has_nj) {
        before += m.get(aj, nj);
        after += m.get(ai, nj);
      }

      const int64_t delta = after - before;
      if (delta < best_delta
          || (delta == best_delta && found
              && (names[aj] < names[order[best_j]]
                  || (names[aj] == names[order[best_j]]
                      && names[ai] < names[order[best_i]])))) {
        found = true;
        best_delta = delta;
        best_i = i;
        best_j = j;
      }
    }
  }

  if (!found || best_delta >= 0) {
    return false;
  }

  std::swap(order[best_i], order[best_j]);
  return true;
}

bool scanOptFirstTwoOptMove(const ScanOptMatrix& m,
                            std::vector<std::size_t>& order)
{
  const std::size_t n = order.size();
  if (n < 4) {
    return false;
  }

  // Prefix sums for forward edges and reversed-adjacent edges.
  std::vector<int64_t> forward_prefix(n, 0);
  std::vector<int64_t> rev_prefix(n, 0);
  for (std::size_t i = 1; i < n; ++i) {
    forward_prefix[i] = forward_prefix[i - 1] + m.get(order[i - 1], order[i]);
    rev_prefix[i]
        = rev_prefix[i - 1] + m.get(order[i], order[i - 1]);  // reversed edge
  }

  const auto segment_forward = [&](std::size_t from,
                                   std::size_t to) -> int64_t {
    return forward_prefix[to] - forward_prefix[from];
  };
  const auto segment_reversed = [&](std::size_t from,
                                    std::size_t to) -> int64_t {
    return rev_prefix[to] - rev_prefix[from];
  };

  for (std::size_t i = 0; i + 2 < n; ++i) {
    for (std::size_t j = i + 1; j < n; ++j) {
      if (j <= i + 1) {
        continue;
      }

      const std::size_t a = order[i];
      const std::size_t b = order[i + 1];
      const std::size_t c = order[j];
      const bool has_d = (j + 1 < n);
      const std::size_t d = has_d ? order[j + 1] : 0;

      const int64_t old_inside = segment_forward(i + 1, j);
      const int64_t new_inside = segment_reversed(i + 1, j);

      const int64_t old_edges
          = static_cast<int64_t>(m.get(a, b))
            + (has_d ? static_cast<int64_t>(m.get(c, d)) : 0);
      const int64_t new_edges
          = static_cast<int64_t>(m.get(a, c))
            + (has_d ? static_cast<int64_t>(m.get(b, d)) : 0);

      const int64_t before = old_edges + old_inside;
      const int64_t after = new_edges + new_inside;

      if (after < before) {
        std::reverse(order.begin() + static_cast<std::ptrdiff_t>(i + 1),
                     order.begin() + static_cast<std::ptrdiff_t>(j + 1));
        return true;
      }
    }
  }
  return false;
}

void scanOptDescent(const ScanOptMatrix& m,
                    const std::vector<std::string_view>& names,
                    std::vector<std::size_t>& order)
{
  // A small number of local-search passes tends to give most of the benefit.
  for (int pass = 0; pass < 4; ++pass) {
    bool improved = scanOptAdjacentSwapImprove(m, names, order);
    improved |= scanOptBestRelocateMove(m, names, order);
    if (order.size() <= kScanOptSwapMaxCells) {
      improved |= scanOptBestSwapMove(m, names, order);
    }
    if (order.size() <= kScanOptTwoOptMaxCells) {
      improved |= scanOptFirstTwoOptMove(m, order);
    }
    if (!improved) {
      break;
    }
  }
}

void scanOptDoubleBridgeKick(std::vector<std::size_t>& order,
                             std::mt19937_64& rng)
{
  const std::size_t n = order.size();
  if (n < 8) {
    std::shuffle(order.begin() + 1, order.end(), rng);
    return;
  }

  std::uniform_int_distribution<std::size_t> dist(1, n - 2);
  std::size_t a = dist(rng);
  std::size_t b = dist(rng);
  std::size_t c = dist(rng);
  std::size_t d = dist(rng);

  std::array<std::size_t, 4> cuts{a, b, c, d};
  std::sort(cuts.begin(), cuts.end());
  a = cuts[0];
  b = cuts[1];
  c = cuts[2];
  d = cuts[3];

  if (a == b || b == c || c == d) {
    return;
  }

  std::vector<std::size_t> kicked;
  kicked.reserve(n);
  kicked.insert(kicked.end(), order.begin(), order.begin() + a);
  kicked.insert(kicked.end(), order.begin() + b, order.begin() + c);
  kicked.insert(kicked.end(), order.begin() + a, order.begin() + b);
  kicked.insert(kicked.end(), order.begin() + c, order.begin() + d);
  kicked.insert(kicked.end(), order.begin() + d, order.end());

  order.swap(kicked);
}

struct UnionFind
{
  explicit UnionFind(std::size_t n) : parent(n), rank(n, 0)
  {
    for (std::size_t i = 0; i < n; ++i) {
      parent[i] = i;
    }
  }

  std::size_t find(std::size_t x)
  {
    if (parent[x] != x) {
      parent[x] = find(parent[x]);
    }
    return parent[x];
  }

  void unite(std::size_t a, std::size_t b)
  {
    a = find(a);
    b = find(b);
    if (a == b) {
      return;
    }
    if (rank[a] < rank[b]) {
      parent[a] = b;
    } else if (rank[a] > rank[b]) {
      parent[b] = a;
    } else {
      parent[b] = a;
      ++rank[a];
    }
  }

  std::vector<std::size_t> parent;
  std::vector<std::size_t> rank;
};

bool hasScanOrderConstraints(const ScanArchitectConfig& config)
{
  return !config.getScanOrderGroups().empty()
         || !config.getScanOrderFixedEdges().empty()
         || !config.getScanOrderBeforeConstraints().empty();
}

std::vector<std::size_t> orderNodesByCost(
    std::size_t n,
    std::size_t start,
    const std::vector<std::string_view>& names,
    const ScanArchitectConfig& config,
    utl::Logger* logger,
    const std::function<int64_t(std::size_t, std::size_t)>& cost_fn,
    const std::function<int64_t(std::size_t)>& terminal_cost)
{
  if (n <= 1) {
    return n == 1 ? std::vector<std::size_t>{0} : std::vector<std::size_t>{};
  }

  const bool use_scanopt
      = (config.getScanOrderSolver()
             == ScanArchitectConfig::ScanOrderSolver::ScanOpt
         && n <= kScanOptMaxMatrixCells && !terminal_cost);

  if (use_scanopt) {
    ScanOptMatrix m;
    m.n = n;
    m.costs.resize(n * n);
    for (std::size_t i = 0; i < n; ++i) {
      for (std::size_t j = 0; j < n; ++j) {
        if (i == j) {
          m.costs[i * n + j] = kScanOptLargeCost;
          continue;
        }
        const int64_t cost = cost_fn(i, j);
        m.costs[i * n + j]
            = static_cast<int32_t>(std::min<int64_t>(cost, kScanOptLargeCost));
      }
    }

    std::vector<std::size_t> best_order
        = scanOptGreedyOrder(m, start, names, logger);
    scanOptDescent(m, names, best_order);
    int64_t best_cost = scanOptPathCost(m, best_order);

    std::mt19937_64 rng(config.getScanOptSeed());
    const uint64_t rounds = config.getScanOptRounds();
    for (uint64_t r = 0; r < rounds; ++r) {
      std::vector<std::size_t> cand = best_order;
      scanOptDoubleBridgeKick(cand, rng);
      scanOptDescent(m, names, cand);
      const int64_t cand_cost = scanOptPathCost(m, cand);
      if (cand_cost < best_cost) {
        best_cost = cand_cost;
        best_order.swap(cand);
      }
    }

    return best_order;
  }

  // Greedy nearest-neighbor ordering (deterministic).
  std::vector<std::size_t> order;
  order.reserve(n);
  std::vector<bool> used(n, false);

  std::size_t cur = start;
  used[cur] = true;
  order.push_back(cur);

  while (order.size() < n) {
    bool found = false;
    std::size_t best = 0;
    int64_t best_cost = kInfDistance;
    const bool last_step = (order.size() + 1 == n);
    for (std::size_t i = 0; i < n; ++i) {
      if (used[i]) {
        continue;
      }
      int64_t cost = cost_fn(cur, i);
      if (terminal_cost && last_step) {
        cost += terminal_cost(i);
      }
      if (!found || cost < best_cost
          || (cost == best_cost && names[i] < names[best])) {
        found = true;
        best = i;
        best_cost = cost;
      }
    }
    if (!found) {
      logger->error(utl::DFT, 180, "Couldn't find next node to order");
    }
    used[best] = true;
    order.push_back(best);
    cur = best;
  }

  // Simple local improvement: repeated adjacent swaps.
  for (int pass = 0; pass < 2; ++pass) {
    bool improved = false;
    for (std::size_t pos = 1; pos + 1 < order.size(); ++pos) {
      const std::size_t prev = order[pos - 1];
      const std::size_t a = order[pos];
      const std::size_t b = order[pos + 1];
      const bool has_next = (pos + 2 < order.size());
      const std::size_t next = has_next ? order[pos + 2] : 0;

      const int64_t before
          = cost_fn(prev, a) + cost_fn(a, b) + (has_next ? cost_fn(b, next) : 0)
            + (!has_next && terminal_cost ? terminal_cost(b) : 0);
      const int64_t after
          = cost_fn(prev, b) + cost_fn(b, a) + (has_next ? cost_fn(a, next) : 0)
            + (!has_next && terminal_cost ? terminal_cost(a) : 0);
      if (after < before) {
        std::swap(order[pos], order[pos + 1]);
        improved = true;
      }
    }
    if (!improved) {
      break;
    }
  }

  return order;
}

void optimizeScanWirelengthWithConstraints(
    std::vector<std::unique_ptr<ScanCell>>& cells,
    const ScanArchitectConfig& config,
    const std::vector<odb::Point>& origins,
    const std::vector<odb::Point>& scan_in_pts,
    const std::vector<odb::Point>& scan_out_pts,
    const std::vector<std::string_view>& names,
    utl::Logger* logger,
    const std::function<int64_t(std::size_t, std::size_t)>& edge_cost,
    const std::optional<odb::Point>& begin,
    const std::optional<odb::Point>& end_pt,
    double vertical_weight)
{
  const std::size_t n = cells.size();
  if (n < 2) {
    return;
  }
  if (!hasScanOrderConstraints(config)) {
    return;
  }

  std::unordered_map<std::string_view, std::size_t> name_to_idx;
  name_to_idx.reserve(n * 2);
  for (std::size_t i = 0; i < n; ++i) {
    name_to_idx.emplace(names[i], i);
  }

  UnionFind uf(n);

  for (const auto& group : config.getScanOrderGroups()) {
    std::optional<std::size_t> first;
    for (const std::string& inst : group.inst_names) {
      const auto it = name_to_idx.find(std::string_view(inst));
      if (it == name_to_idx.end()) {
        continue;
      }
      if (!first.has_value()) {
        first = it->second;
      } else {
        uf.unite(first.value(), it->second);
      }
    }
  }

  std::vector<std::pair<std::size_t, std::size_t>> fixed_edges;
  fixed_edges.reserve(config.getScanOrderFixedEdges().size());
  for (const auto& edge : config.getScanOrderFixedEdges()) {
    const auto it_from = name_to_idx.find(std::string_view(edge.from_inst));
    const auto it_to = name_to_idx.find(std::string_view(edge.to_inst));
    if (it_from == name_to_idx.end() || it_to == name_to_idx.end()) {
      continue;
    }
    uf.unite(it_from->second, it_to->second);
    fixed_edges.emplace_back(it_from->second, it_to->second);
  }

  std::unordered_map<std::size_t, std::vector<std::size_t>> comps;
  comps.reserve(n);
  for (std::size_t i = 0; i < n; ++i) {
    comps[uf.find(i)].push_back(i);
  }

  const int default_pr = config.getDefaultGroupPriority();
  std::unordered_map<std::size_t, int> comp_priority;
  comp_priority.reserve(comps.size());
  for (const auto& [root, _] : comps) {
    comp_priority[root] = default_pr;
  }

  // Assign priorities to components; if multiple groups map into the same
  // component, pick the smallest priority and warn on conflicts.
  for (const auto& group : config.getScanOrderGroups()) {
    std::optional<std::size_t> any_idx;
    for (const std::string& inst : group.inst_names) {
      const auto it = name_to_idx.find(std::string_view(inst));
      if (it != name_to_idx.end()) {
        any_idx = it->second;
        break;
      }
    }
    if (!any_idx.has_value()) {
      continue;
    }
    const std::size_t root = uf.find(any_idx.value());
    const int prev = comp_priority[root];
    const int pr = std::clamp(group.priority, 0, 127);
    if (prev != default_pr && prev != pr && logger) {
      logger->warn(utl::DFT,
                   181,
                   "Scan constraints: merged groups with priorities {} and {} "
                   "into one component (using {}).",
                   prev,
                   pr,
                   std::min(prev, pr));
    }
    comp_priority[root] = std::min(prev, pr);
  }

  struct Component
  {
    int priority = 0;
    std::string_view rep;
    std::vector<std::size_t> members;
    std::vector<std::size_t> ordered;
    std::size_t entry = 0;
    std::size_t exit = 0;
  };

  std::vector<Component> components;
  components.reserve(comps.size());
  for (auto& [root, members] : comps) {
    std::stable_sort(members.begin(),
                     members.end(),
                     [&](std::size_t a, std::size_t b) {
      return names[a] < names[b];
    });
    Component comp;
    comp.priority = comp_priority[root];
    comp.rep = names[members.front()];
    comp.members = std::move(members);
    components.push_back(std::move(comp));
  }

  std::stable_sort(components.begin(),
                   components.end(),
                   [&](const Component& a, const Component& b) {
                     if (a.priority != b.priority) {
                       return a.priority < b.priority;
                     }
                     return a.rep < b.rep;
                   });

  // Precompute component-internal ordering (fixed edges enforced).
  const auto build_internal_order = [&](Component& comp) {
    if (comp.members.size() == 1) {
      comp.ordered = comp.members;
      comp.entry = comp.members.front();
      comp.exit = comp.members.front();
      return;
    }

    std::vector<char> in_comp(n, 0);
    for (const std::size_t m : comp.members) {
      in_comp[m] = 1;
    }

    std::unordered_map<std::size_t, std::size_t> succ;
    std::unordered_map<std::size_t, std::size_t> pred;
    succ.reserve(comp.members.size());
    pred.reserve(comp.members.size());

    for (const auto& [from, to] : fixed_edges) {
      if (!in_comp[from] || !in_comp[to]) {
        continue;
      }
      const auto succ_it = succ.find(from);
      if (succ_it != succ.end() && succ_it->second != to) {
        if (logger) {
          logger->error(utl::DFT,
                        182,
                        "Scan constraints infeasible: fixed-edge conflict on '{}' "
                        "(multiple successors: '{}' and '{}').",
                        names[from],
                        names[succ_it->second],
                        names[to]);
        }
      }
      const auto pred_it = pred.find(to);
      if (pred_it != pred.end() && pred_it->second != from) {
        if (logger) {
          logger->error(utl::DFT,
                        183,
                        "Scan constraints infeasible: fixed-edge conflict on '{}' "
                        "(multiple predecessors: '{}' and '{}').",
                        names[to],
                        names[pred_it->second],
                        names[from]);
        }
      }
      succ[from] = to;
      pred[to] = from;
    }

    // Detect cycles in strict fixed-edge constraints (not representable in a
    // directed Hamiltonian path).
    std::vector<char> color(n, 0);  // 0=unseen, 1=visiting, 2=done
    for (const std::size_t start : comp.members) {
      if (color[start] != 0) {
        continue;
      }
      std::vector<std::size_t> stack;
      std::size_t cur = start;
      while (true) {
        if (color[cur] == 0) {
          color[cur] = 1;
          stack.push_back(cur);
          const auto it = succ.find(cur);
          if (it == succ.end()) {
            break;
          }
          cur = it->second;
          continue;
        }
        if (color[cur] == 1) {
          if (logger) {
            logger->error(
                utl::DFT,
                194,
                "Scan constraints infeasible: fixed-edge cycle detected within "
                "a constrained component (cycle includes '{}').",
                names[cur]);
          }
        }
        break;  // color==2
      }
      for (const std::size_t v : stack) {
        color[v] = 2;
      }
    }

    std::vector<char> visited(n, 0);
    std::vector<std::vector<std::size_t>> segments;

    const auto emit_segment_from = [&](std::size_t start) {
      std::vector<std::size_t> seg;
      std::size_t cur = start;
      while (in_comp[cur] && !visited[cur]) {
        visited[cur] = 1;
        seg.push_back(cur);
        const auto it = succ.find(cur);
        if (it == succ.end()) {
          break;
        }
        cur = it->second;
      }
      if (!seg.empty()) {
        segments.push_back(std::move(seg));
      }
    };

    // Start with nodes that have no predecessor to form maximal paths.
    for (const std::size_t m : comp.members) {
      if (pred.find(m) == pred.end()) {
        emit_segment_from(m);
      }
    }
    // Break remaining cycles arbitrarily (best effort).
    for (const std::size_t m : comp.members) {
      if (!visited[m]) {
        emit_segment_from(m);
      }
    }

    if (segments.empty()) {
      comp.ordered = comp.members;
      comp.entry = comp.members.front();
      comp.exit = comp.members.back();
      return;
    }

    if (segments.size() == 1) {
      comp.ordered = segments.front();
      comp.entry = comp.ordered.front();
      comp.exit = comp.ordered.back();
      return;
    }

    const std::size_t sn = segments.size();
    std::vector<std::string_view> seg_names;
    seg_names.reserve(sn);
    std::vector<std::size_t> seg_entry(sn, 0);
    std::vector<std::size_t> seg_exit(sn, 0);
    for (std::size_t si = 0; si < sn; ++si) {
      seg_entry[si] = segments[si].front();
      seg_exit[si] = segments[si].back();
      seg_names.push_back(names[seg_entry[si]]);
    }

    std::size_t seg_start = 0;
    int64_t best_key = std::numeric_limits<int64_t>::max();
    for (std::size_t si = 0; si < sn; ++si) {
      const odb::Point& p = origins[seg_entry[si]];
      const int64_t key = static_cast<int64_t>(p.x()) + static_cast<int64_t>(p.y());
      if (key < best_key
          || (key == best_key && seg_names[si] < seg_names[seg_start])) {
        best_key = key;
        seg_start = si;
      }
    }

    const auto seg_cost = [&](std::size_t a, std::size_t b) -> int64_t {
      return edge_cost(seg_exit[a], seg_entry[b]);
    };
    const std::vector<std::size_t> seg_order
        = orderNodesByCost(sn, seg_start, seg_names, config, logger, seg_cost, {});

    comp.ordered.clear();
    for (const std::size_t si : seg_order) {
      comp.ordered.insert(
          comp.ordered.end(), segments[si].begin(), segments[si].end());
    }
    comp.entry = comp.ordered.front();
    comp.exit = comp.ordered.back();
  };

  for (auto& comp : components) {
    build_internal_order(comp);
  }

  std::vector<std::size_t> final_order;
  final_order.reserve(n);

  if (!config.getScanOrderBeforeConstraints().empty()) {
    // Build component-level DAG from before constraints.
    const std::size_t cn = components.size();
    std::unordered_map<std::size_t, std::size_t> root_to_comp;
    root_to_comp.reserve(cn * 2);
    for (std::size_t ci = 0; ci < cn; ++ci) {
      root_to_comp.emplace(uf.find(components[ci].members.front()), ci);
    }

    std::unordered_map<std::string_view, std::size_t> group_to_comp;
    group_to_comp.reserve(config.getScanOrderGroups().size() * 2);
    for (const auto& group : config.getScanOrderGroups()) {
      if (group.name.empty()) {
        continue;
      }
      std::optional<std::size_t> any_idx;
      for (const std::string& inst : group.inst_names) {
        const auto it = name_to_idx.find(std::string_view(inst));
        if (it != name_to_idx.end()) {
          any_idx = it->second;
          break;
        }
      }
      if (!any_idx.has_value()) {
        continue;
      }
      const std::size_t root = uf.find(any_idx.value());
      const auto rit = root_to_comp.find(root);
      if (rit != root_to_comp.end()) {
        group_to_comp.emplace(group.name, rit->second);
      }
    }

    const auto resolve_token
        = [&](std::string_view tok) -> std::optional<std::size_t> {
      const auto it = name_to_idx.find(tok);
      if (it != name_to_idx.end()) {
        const auto rit = root_to_comp.find(uf.find(it->second));
        if (rit != root_to_comp.end()) {
          return rit->second;
        }
      }
      const auto git = group_to_comp.find(tok);
      if (git != group_to_comp.end()) {
        return git->second;
      }
      return std::nullopt;
    };

    std::vector<std::vector<std::size_t>> succs(cn);
    for (const auto& bc : config.getScanOrderBeforeConstraints()) {
      const auto a = resolve_token(bc.before);
      const auto b = resolve_token(bc.after);
      if (!a.has_value() || !b.has_value()) {
        continue;
      }
      if (a.value() == b.value()) {
        continue;
      }
      succs[a.value()].push_back(b.value());
    }

    // Deduplicate edges and compute indegrees.
    std::vector<int> indeg(cn, 0);
    for (auto& vs : succs) {
      std::sort(vs.begin(), vs.end());
      vs.erase(std::unique(vs.begin(), vs.end()), vs.end());
      for (const std::size_t to : vs) {
        indeg[to]++;
      }
    }

    std::vector<char> used_comp(cn, 0);
    std::optional<std::size_t> prev_exit;
    for (std::size_t step = 0; step < cn; ++step) {
      const std::size_t remaining = cn - step;
      bool found = false;
      std::size_t best = 0;
      int64_t best_cost = kInfDistance;

      for (std::size_t ci = 0; ci < cn; ++ci) {
        if (used_comp[ci] || indeg[ci] != 0) {
          continue;
        }

        int64_t cost = 0;
        if (prev_exit.has_value()) {
          cost = edge_cost(prev_exit.value(), components[ci].entry);
        } else if (begin.has_value()) {
          cost = manhattanDist(begin.value(),
                               scan_in_pts[components[ci].entry],
                               vertical_weight);
        } else {
          const odb::Point& p = origins[components[ci].entry];
          cost = static_cast<int64_t>(p.x()) + static_cast<int64_t>(p.y());
        }

        if (remaining == 1 && end_pt.has_value()) {
          cost += manhattanDist(scan_out_pts[components[ci].exit],
                                end_pt.value(),
                                vertical_weight);
        }

        if (!found || cost < best_cost
            || (cost == best_cost
                && (components[ci].priority < components[best].priority
                    || (components[ci].priority == components[best].priority
                        && components[ci].rep < components[best].rep)))) {
          best = ci;
          best_cost = cost;
          found = true;
        }
      }

      if (!found) {
        if (logger) {
          logger->error(
              utl::DFT,
              190,
              "Scan constraints: before-constraint cycle detected while "
              "ordering components.");
        }
        return;
      }

      used_comp[best] = 1;
      final_order.insert(final_order.end(),
                         components[best].ordered.begin(),
                         components[best].ordered.end());
      prev_exit = components[best].exit;
      for (const std::size_t to : succs[best]) {
        indeg[to]--;
      }
    }
  } else {
    // Order components bucketed by priority (ScanOpt-style).
    std::optional<std::size_t> prev_exit;
    std::size_t idx = 0;
    while (idx < components.size()) {
      const int pr = components[idx].priority;
      std::size_t bucket_end = idx;
      while (bucket_end < components.size() && components[bucket_end].priority == pr) {
        ++bucket_end;
      }

      const std::size_t cn = bucket_end - idx;
      if (cn == 1) {
        final_order.insert(final_order.end(),
                           components[idx].ordered.begin(),
                           components[idx].ordered.end());
        prev_exit = components[idx].exit;
        idx = bucket_end;
        continue;
      }

      std::vector<std::string_view> comp_names;
      comp_names.reserve(cn);
      for (std::size_t ci = idx; ci < bucket_end; ++ci) {
        comp_names.push_back(components[ci].rep);
      }

      std::size_t start_comp = 0;
      if (prev_exit.has_value()) {
        bool found = false;
        int64_t best = kInfDistance;
        for (std::size_t off = 0; off < cn; ++off) {
          const Component& c = components[idx + off];
          const int64_t cost = edge_cost(prev_exit.value(), c.entry);
          if (!found || cost < best
              || (cost == best && c.rep < components[idx + start_comp].rep)) {
            best = cost;
            start_comp = off;
            found = true;
          }
        }
      } else {
        int64_t best_key = std::numeric_limits<int64_t>::max();
        for (std::size_t off = 0; off < cn; ++off) {
          const Component& c = components[idx + off];
          int64_t key = 0;
          if (begin.has_value()) {
            key = manhattanDist(begin.value(),
                                scan_in_pts[c.entry],
                                vertical_weight);
          } else {
            const odb::Point& p = origins[c.entry];
            key = static_cast<int64_t>(p.x()) + static_cast<int64_t>(p.y());
          }
          if (key < best_key
              || (key == best_key && c.rep < components[idx + start_comp].rep)) {
            best_key = key;
            start_comp = off;
          }
        }
      }

      const auto comp_cost = [&](std::size_t a, std::size_t b) -> int64_t {
        const Component& ca = components[idx + a];
        const Component& cb = components[idx + b];
        return edge_cost(ca.exit, cb.entry);
      };

      std::function<int64_t(std::size_t)> terminal_cost;
      const bool last_bucket = (bucket_end == components.size());
      if (last_bucket && end_pt.has_value()) {
        terminal_cost = [&](std::size_t a) -> int64_t {
          const Component& ca = components[idx + a];
          return manhattanDist(scan_out_pts[ca.exit],
                               end_pt.value(),
                               vertical_weight);
        };
      }

      const std::vector<std::size_t> comp_order = orderNodesByCost(
          cn, start_comp, comp_names, config, logger, comp_cost, terminal_cost);

      for (const std::size_t off : comp_order) {
        Component& c = components[idx + off];
        final_order.insert(final_order.end(), c.ordered.begin(), c.ordered.end());
        prev_exit = c.exit;
      }

      idx = bucket_end;
    }
  }

  if (final_order.size() != n) {
    logger->warn(utl::DFT,
                 184,
                 "Scan constraints: internal error while ordering; expected {} "
                 "cells but got {}. Leaving order unchanged.",
                 n,
                 final_order.size());
    return;
  }

  std::vector<std::unique_ptr<ScanCell>> ordered;
  ordered.reserve(n);
  for (const std::size_t i : final_order) {
    ordered.emplace_back(std::move(cells[i]));
  }
  cells.swap(ordered);
}

void OptimizeScanWirelengthPinToNet(std::vector<std::unique_ptr<ScanCell>>& cells,
                                    const ScanArchitectConfig& config,
                                    utl::Logger* logger,
                                    const std::optional<ScanArchitectConfig::ChainEndpoints>&
                                        endpoints)
{
  const std::size_t n = cells.size();
  if (n < 2) {
    return;
  }

  std::optional<odb::Point> begin;
  std::optional<odb::Point> end;
  if (endpoints.has_value()) {
    if (endpoints->begin.has_value()) {
      begin = odb::Point(endpoints->begin->x, endpoints->begin->y);
    }
    if (endpoints->end.has_value()) {
      end = odb::Point(endpoints->end->x, endpoints->end->y);
    }
  }

  std::vector<odb::Point> origins;
  origins.reserve(n);
  std::vector<odb::Point> scan_in_pts;
  scan_in_pts.reserve(n);
  std::vector<odb::Point> scan_out_pts;
  scan_out_pts.reserve(n);
  std::vector<std::string_view> names;
  names.reserve(n);
  std::vector<NetAccessGeometry> net_geoms;
  net_geoms.reserve(n);
  std::vector<double> timing_mul;
  timing_mul.reserve(n);

  for (const auto& cell : cells) {
    const odb::Point origin = cell->getOrigin();
    origins.emplace_back(origin);
    names.emplace_back(cell->getName());

    scan_in_pts.emplace_back(scanPinLocation(cell->getScanIn(), origin));
    scan_out_pts.emplace_back(scanPinLocation(cell->getScanOut(), origin));
    net_geoms.emplace_back(buildNetAccessGeometry(cell->getScanOut().getNet()));
    timing_mul.emplace_back(timingMultiplierForCell(config, *cell));
  }

  const double vertical_weight = config.getVerticalWeight();

  std::size_t start_index = 0;
  int64_t lowest = std::numeric_limits<int64_t>::max();
  if (begin.has_value()) {
    for (std::size_t i = 0; i < n; ++i) {
      const int64_t score
          = manhattanDist(begin.value(), scan_in_pts[i], vertical_weight);
      if (score < lowest
          || (score == lowest && names[i] < names[start_index])) {
        start_index = i;
        lowest = score;
      }
    }
  } else {
    for (std::size_t i = 0; i < n; ++i) {
      const odb::Point& p = origins[i];
      const int64_t score
          = static_cast<int64_t>(p.x()) + static_cast<int64_t>(p.y());
      if (score < lowest
          || (score == lowest && names[i] < names[start_index])) {
        start_index = i;
        lowest = score;
      }
    }
  }

  const auto edge_cost = [&](std::size_t src, std::size_t dst) -> int64_t {
    const int64_t dist
        = pinToNetDistance(scan_in_pts[dst], net_geoms[src], vertical_weight);
    if (dist != 0 || !net_geoms[src].boxes.empty()
        || !net_geoms[src].terminals.empty()) {
      return scaleEdgeCost(dist, timing_mul[src]);
    }
    // No routing/pin geometry available; fall back to pin-based Manhattan.
    return scaleEdgeCost(
        manhattanDist(scan_out_pts[src], scan_in_pts[dst], vertical_weight),
        timing_mul[src]);
  };

  if (hasScanOrderConstraints(config)) {
    optimizeScanWirelengthWithConstraints(
        cells,
        config,
        origins,
        scan_in_pts,
        scan_out_pts,
        names,
        logger,
        edge_cost,
        begin,
        end,
        vertical_weight);
    return;
  }

  if (config.getScanOrderSolver()
      == ScanArchitectConfig::ScanOrderSolver::ScanOpt) {
    if (end.has_value()) {
      logger->warn(
          utl::DFT,
          186,
          "ScanOpt ordering does not currently support EndPort costs; falling "
          "back to heuristic ordering.");
    } else if (n > kScanOptMaxMatrixCells) {
      logger->warn(
          utl::DFT,
          72,
          "ScanOpt ordering requested for {} cells, which exceeds the current "
          "matrix limit {}. Falling back to heuristic ordering.",
          n,
          kScanOptMaxMatrixCells);
    } else {
      ScanOptMatrix m;
      m.n = n;
      m.costs.resize(n * n);
      for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
          if (i == j) {
            m.costs[i * n + j] = kScanOptLargeCost;
            continue;
          }
          const int64_t cost = edge_cost(i, j);
          m.costs[i * n + j]
              = static_cast<int32_t>(std::min<int64_t>(cost, kScanOptLargeCost));
        }
      }

      std::vector<std::size_t> best_order
          = scanOptGreedyOrder(m, start_index, names, logger);
      scanOptDescent(m, names, best_order);
      int64_t best_cost = scanOptPathCost(m, best_order);

      std::mt19937_64 rng(config.getScanOptSeed());
      const uint64_t rounds = config.getScanOptRounds();
      for (uint64_t r = 0; r < rounds; ++r) {
        std::vector<std::size_t> cand = best_order;
        scanOptDoubleBridgeKick(cand, rng);
        scanOptDescent(m, names, cand);
        const int64_t cand_cost = scanOptPathCost(m, cand);
        if (cand_cost < best_cost) {
          best_cost = cand_cost;
          best_order.swap(cand);
        }
      }

      std::vector<std::unique_ptr<ScanCell>> ordered;
      ordered.reserve(n);
      for (const std::size_t idx : best_order) {
        ordered.emplace_back(std::move(cells[idx]));
      }
      std::swap(cells, ordered);
      return;
    }
  }

  std::vector<std::size_t> order;
  order.reserve(n);
  std::vector<bool> used(n, false);

  std::size_t cur = start_index;
  used[cur] = true;
  order.push_back(cur);

  while (order.size() < n) {
    bool found = false;
    std::size_t best = 0;
    int64_t best_cost = kInfDistance;
    const bool last_step = (order.size() + 1 == n);

    for (std::size_t i = 0; i < n; ++i) {
      if (used[i]) {
        continue;
      }
      int64_t cost = edge_cost(cur, i);
      if (end.has_value() && last_step) {
        cost += manhattanDist(scan_out_pts[i], end.value(), vertical_weight);
      }
      if (!found || cost < best_cost
          || (cost == best_cost && names[i] < names[best])) {
        found = true;
        best = i;
        best_cost = cost;
      }
    }

    if (!found) {
      logger->error(utl::DFT, 17, "Couldn't find next scan cell to order");
    }

    used[best] = true;
    order.push_back(best);
    cur = best;
  }

  // Local improvement: repeated adjacent swaps (directional 2-opt-lite).
  for (int pass = 0; pass < 3; ++pass) {
    bool improved = false;
    for (std::size_t pos = 1; pos + 1 < order.size(); ++pos) {
      const std::size_t prev = order[pos - 1];
      const std::size_t a = order[pos];
      const std::size_t b = order[pos + 1];
      const bool has_next = (pos + 2 < order.size());
      const std::size_t next = has_next ? order[pos + 2] : 0;

      const int64_t before
          = edge_cost(prev, a) + edge_cost(a, b)
            + (has_next ? edge_cost(b, next) : 0)
            + (!has_next && end.has_value()
                   ? manhattanDist(scan_out_pts[b], end.value(), vertical_weight)
                   : 0);
      const int64_t after
          = edge_cost(prev, b) + edge_cost(b, a)
            + (has_next ? edge_cost(a, next) : 0)
            + (!has_next && end.has_value()
                   ? manhattanDist(scan_out_pts[a], end.value(), vertical_weight)
                   : 0);

      if (after < before) {
        std::swap(order[pos], order[pos + 1]);
        improved = true;
      }
    }
    if (!improved) {
      break;
    }
  }

  std::vector<std::unique_ptr<ScanCell>> ordered;
  ordered.reserve(n);
  for (const std::size_t idx : order) {
    ordered.emplace_back(std::move(cells[idx]));
  }
  std::swap(cells, ordered);
}
}  // namespace

void OptimizeScanWirelength(std::vector<std::unique_ptr<ScanCell>>& cells,
                            const ScanArchitectConfig& config,
                            utl::Logger* logger)
{
  OptimizeScanWirelength(cells, config, logger, std::nullopt);
}

void OptimizeScanWirelength(
    std::vector<std::unique_ptr<ScanCell>>& cells,
    const ScanArchitectConfig& config,
    utl::Logger* logger,
    const std::optional<ScanArchitectConfig::ChainEndpoints>& endpoints)
{
  // Nothing to order
  if (cells.empty()) {
    return;
  }
  // No point running this if the cells aren't placed yet
  for (const auto& cell : cells) {
    if (!cell->isPlaced()) {
      return;
    }
  }

  if (config.getScanOrderMetric()
      == ScanArchitectConfig::ScanOrderMetric::PinToNet) {
    OptimizeScanWirelengthPinToNet(cells, config, logger, endpoints);
    return;
  }

  const double vertical_weight = config.getVerticalWeight();
  const auto manhattan = [&](const odb::Point& a, const odb::Point& b) -> int64_t {
    return manhattanDist(a, b, vertical_weight);
  };

  const auto path_len = [&](const std::vector<std::size_t>& order,
                            const std::vector<odb::Point>& scan_in_pts,
                            const std::vector<odb::Point>& scan_out_pts) -> int64_t {
    int64_t total = 0;
    for (std::size_t i = 1; i < order.size(); ++i) {
      total += manhattan(scan_out_pts[order[i - 1]], scan_in_pts[order[i]]);
    }
    return total;
  };

  const auto no_next_scan_cell = [&]() -> void {
    logger->error(utl::DFT, 16, "Couldn't find next scan cell to order");
  };

  const std::size_t n = cells.size();
  std::vector<odb::Point> origins;
  origins.reserve(n);
  std::vector<odb::Point> scan_in_pts;
  scan_in_pts.reserve(n);
  std::vector<odb::Point> scan_out_pts;
  scan_out_pts.reserve(n);
  std::vector<std::string_view> names;
  names.reserve(n);
  std::vector<double> timing_mul;
  timing_mul.reserve(n);
  for (const auto& cell : cells) {
    origins.emplace_back(cell->getOrigin());
    scan_in_pts.emplace_back(scanPinLocation(cell->getScanIn(), origins.back()));
    scan_out_pts.emplace_back(
        scanPinLocation(cell->getScanOut(), origins.back()));
    names.emplace_back(cell->getName());
    timing_mul.emplace_back(timingMultiplierForCell(config, *cell));
  }

  std::optional<odb::Point> begin;
  std::optional<odb::Point> end;
  if (endpoints.has_value()) {
    if (endpoints->begin.has_value()) {
      begin = odb::Point(endpoints->begin->x, endpoints->begin->y);
    }
    if (endpoints->end.has_value()) {
      end = odb::Point(endpoints->end->x, endpoints->end->y);
    }
  }

  const auto edge_cost = [&](std::size_t src, std::size_t dst) -> int64_t {
    return scaleEdgeCost(manhattan(scan_out_pts[src], scan_in_pts[dst]),
                         timing_mul[src]);
  };

  bool directional_pins = false;
  for (std::size_t i = 0; i < n; ++i) {
    if (scan_in_pts[i].x() != scan_out_pts[i].x()
        || scan_in_pts[i].y() != scan_out_pts[i].y()) {
      directional_pins = true;
      break;
    }
  }

  if (hasScanOrderConstraints(config)) {
    optimizeScanWirelengthWithConstraints(
        cells,
        config,
        origins,
        scan_in_pts,
        scan_out_pts,
        names,
        logger,
        edge_cost,
        begin,
        end,
        vertical_weight);
    return;
  }

  // Define the starting node as the lower leftmost, so we don't accidentally
  // start somewhere in the middle
  size_t start_index = 0;
  int64_t lowest_dist = std::numeric_limits<int64_t>::max();

  if (begin.has_value()) {
    for (size_t i = 0; i < n; i++) {
      const int64_t dist = manhattan(begin.value(), scan_in_pts[i]);
      if (dist < lowest_dist
          || (dist == lowest_dist && names[i] < names[start_index])) {
        start_index = i;
        lowest_dist = dist;
      }
    }
  } else {
    for (size_t i = 0; i < n; i++) {
      // Find the lower leftmost cell by looking for the cell with the lowest
      // manhattan distance to the origin.
      const auto& origin = scan_in_pts[i];
      const int64_t dist
          = static_cast<int64_t>(origin.x()) + static_cast<int64_t>(origin.y());
      if (dist < lowest_dist
          || (dist == lowest_dist && names[i] < names[start_index])) {
        start_index = i;
        lowest_dist = dist;
      }
    }
  }

  if (config.getScanOrderSolver()
      == ScanArchitectConfig::ScanOrderSolver::ScanOpt) {
    if (end.has_value()) {
      logger->warn(
          utl::DFT,
          187,
          "ScanOpt ordering does not currently support EndPort costs; falling "
          "back to heuristic ordering.");
    } else if (n > kScanOptMaxMatrixCells) {
      logger->warn(
          utl::DFT,
          185,
          "ScanOpt ordering requested for {} cells, which exceeds the current "
          "matrix limit {}. Falling back to heuristic ordering.",
          n,
          kScanOptMaxMatrixCells);
    } else {
      ScanOptMatrix m;
      m.n = n;
      m.costs.resize(n * n);
      for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
          if (i == j) {
            m.costs[i * n + j] = kScanOptLargeCost;
            continue;
          }
          const int64_t cost = edge_cost(i, j);
          m.costs[i * n + j]
              = static_cast<int32_t>(std::min<int64_t>(cost, kScanOptLargeCost));
        }
      }

      std::vector<std::size_t> best_order
          = scanOptGreedyOrder(m, start_index, names, logger);
      scanOptDescent(m, names, best_order);
      int64_t best_cost = scanOptPathCost(m, best_order);

      std::mt19937_64 rng(config.getScanOptSeed());
      const uint64_t rounds = config.getScanOptRounds();
      for (uint64_t r = 0; r < rounds; ++r) {
        std::vector<std::size_t> cand = best_order;
        scanOptDoubleBridgeKick(cand, rng);
        scanOptDescent(m, names, cand);
        const int64_t cand_cost = scanOptPathCost(m, cand);
        if (cand_cost < best_cost) {
          best_cost = cand_cost;
          best_order.swap(cand);
        }
      }

      std::vector<std::unique_ptr<ScanCell>> ordered;
      ordered.reserve(n);
      for (const std::size_t idx : best_order) {
        ordered.emplace_back(std::move(cells[idx]));
      }
      std::swap(cells, ordered);
      return;
    }
  }

  if (config.getScanOrderSolver()
      == ScanArchitectConfig::ScanOrderSolver::MinFeedthrough) {
    if (config.getTimingWeightSetup() != 0.0
        || config.getTimingWeightHold() != 0.0) {
      logger->warn(
          utl::DFT,
          188,
          "Min-feedthrough ordering does not currently support timing-aware "
          "penalties; falling back to heuristic ordering.");
    } else {
      const std::vector<std::size_t> order
          = minFeedthroughRowSweepOrder(origins, names, vertical_weight);
      if (order.size() == n) {
        std::vector<std::unique_ptr<ScanCell>> ordered;
        ordered.reserve(n);
        for (const std::size_t idx : order) {
          ordered.emplace_back(std::move(cells[idx]));
        }
        std::swap(cells, ordered);
        return;
      }
      logger->warn(utl::DFT,
                   189,
                   "Min-feedthrough ordering failed (expected {} cells but got "
                   "{}); falling back to heuristic ordering.",
                   n,
                   order.size());
    }
  }

  // Timing-aware ordering makes the objective asymmetric (penalizes outgoing
  // edges from timing-critical sources). EndPort costs also make the objective
  // directional. When scan-in/out pins differ, scan_out->scan_in costs are also
  // directional. Use a directed greedy heuristic with a small local-improvement
  // pass rather than symmetric 2-opt/farthest-insertion moves.
  if (config.getTimingWeightSetup() != 0.0 || config.getTimingWeightHold() != 0.0
      || end.has_value() || directional_pins) {
    std::vector<std::size_t> order;
    order.reserve(n);
    std::vector<bool> used(n, false);

    std::size_t cur = start_index;
    used[cur] = true;
    order.push_back(cur);

    while (order.size() < n) {
      bool found = false;
      std::size_t best = 0;
      int64_t best_cost = kInfDistance;
      const bool last_step = (order.size() + 1 == n);

      for (std::size_t i = 0; i < n; ++i) {
        if (used[i]) {
          continue;
        }
        int64_t cost = edge_cost(cur, i);
        if (end.has_value() && last_step) {
          cost += manhattan(scan_out_pts[i], end.value());
        }
        if (!found || cost < best_cost
            || (cost == best_cost && names[i] < names[best])) {
          found = true;
          best = i;
          best_cost = cost;
        }
      }

      if (!found) {
        no_next_scan_cell();
      }

      used[best] = true;
      order.push_back(best);
      cur = best;
    }

    // Local improvement: repeated adjacent swaps (directional 2-opt-lite).
    for (int pass = 0; pass < 3; ++pass) {
      bool improved = false;
      for (std::size_t pos = 1; pos + 1 < order.size(); ++pos) {
        const std::size_t prev = order[pos - 1];
        const std::size_t a = order[pos];
        const std::size_t b = order[pos + 1];
        const bool has_next = (pos + 2 < order.size());
        const std::size_t next = has_next ? order[pos + 2] : 0;

        const int64_t before
            = edge_cost(prev, a) + edge_cost(a, b)
              + (has_next ? edge_cost(b, next) : 0)
              + (!has_next && end.has_value() ? manhattan(scan_out_pts[b], end.value())
                                              : 0);
        const int64_t after
            = edge_cost(prev, b) + edge_cost(b, a)
              + (has_next ? edge_cost(a, next) : 0)
              + (!has_next && end.has_value() ? manhattan(scan_out_pts[a], end.value())
                                              : 0);

        if (after < before) {
          std::swap(order[pos], order[pos + 1]);
          improved = true;
        }
      }
      if (!improved) {
        break;
      }
    }

    std::vector<std::unique_ptr<ScanCell>> ordered;
    ordered.reserve(n);
    for (const std::size_t idx : order) {
      ordered.emplace_back(std::move(cells[idx]));
    }
    std::swap(cells, ordered);
    return;
  }

  const auto greedy_nn_order
      = [&](std::size_t start) -> std::vector<std::size_t> {
    std::vector<std::size_t> order;
    order.reserve(n);
    std::vector<bool> used(n, false);

    std::size_t cur = start;
    used[cur] = true;
    order.push_back(cur);

    while (order.size() < n) {
      bool found = false;
      std::size_t best = cur;
      int64_t best_dist = std::numeric_limits<int64_t>::max();
      for (std::size_t i = 0; i < n; ++i) {
        if (used[i]) {
          continue;
        }
        const int64_t dist = manhattan(scan_in_pts[cur], scan_in_pts[i]);
        if (!found || dist < best_dist
            || (dist == best_dist && names[i] < names[best])) {
          best = i;
          best_dist = dist;
          found = true;
        }
      }
      if (!found) {
        no_next_scan_cell();
      }
      used[best] = true;
      cur = best;
      order.push_back(cur);
    }
    return order;
  };

  const auto farthest_insertion_order
      = [&](std::size_t start) -> std::vector<std::size_t> {
    std::vector<std::size_t> order;
    order.reserve(n);
    std::vector<bool> in_path(n, false);

    order.push_back(start);
    in_path[start] = true;

    if (n == 1) {
      return order;
    }

    // Seed the path with the cell farthest from the start so the initial
    // segment spans the placement.
    std::size_t farthest = start;
    int64_t farthest_dist = -1;
    for (std::size_t i = 0; i < n; ++i) {
      if (i == start) {
        continue;
      }
      const int64_t dist = manhattan(scan_in_pts[start], scan_in_pts[i]);
      if (dist > farthest_dist
          || (dist == farthest_dist && names[i] < names[farthest])) {
        farthest = i;
        farthest_dist = dist;
      }
    }

    order.push_back(farthest);
    in_path[farthest] = true;

    std::vector<int64_t> nearest_dist(n, std::numeric_limits<int64_t>::max());
    for (std::size_t i = 0; i < n; ++i) {
      if (in_path[i]) {
        nearest_dist[i] = 0;
        continue;
      }
      nearest_dist[i] = std::min(manhattan(scan_in_pts[i], scan_in_pts[start]),
                                 manhattan(scan_in_pts[i], scan_in_pts[farthest]));
    }

    while (order.size() < n) {
      // Pick the cell farthest from the current path (farthest insertion).
      std::size_t next = start;
      bool found = false;
      int64_t best_score = -1;
      for (std::size_t i = 0; i < n; ++i) {
        if (in_path[i]) {
          continue;
        }
        const int64_t score = nearest_dist[i];
        if (!found || score > best_score
            || (score == best_score && names[i] < names[next])) {
          next = i;
          best_score = score;
          found = true;
        }
      }
      if (!found) {
        no_next_scan_cell();
      }

      // Insert it at the position that minimises the path length increase.
      std::size_t best_pos = order.size();  // append by default
      int64_t best_delta = manhattan(scan_in_pts[order.back()], scan_in_pts[next]);

      for (std::size_t pos = 0; pos + 1 < order.size(); ++pos) {
        const auto a = order[pos];
        const auto b = order[pos + 1];
        const int64_t delta = manhattan(scan_in_pts[a], scan_in_pts[next])
                              + manhattan(scan_in_pts[next], scan_in_pts[b])
                              - manhattan(scan_in_pts[a], scan_in_pts[b]);
        if (delta < best_delta || (delta == best_delta && pos + 1 < best_pos)) {
          best_delta = delta;
          best_pos = pos + 1;
        }
      }

      order.insert(order.begin() + static_cast<std::ptrdiff_t>(best_pos), next);
      in_path[next] = true;

      // Update the nearest distance cache for remaining nodes.
      for (std::size_t i = 0; i < n; ++i) {
        if (in_path[i]) {
          continue;
        }
        const int64_t dist = manhattan(scan_in_pts[i], scan_in_pts[next]);
        nearest_dist[i] = std::min(nearest_dist[i], dist);
      }
    }

    // Keep the start fixed as the first element.
    if (!order.empty() && order.front() != start) {
      const auto it = std::find(order.begin(), order.end(), start);
      if (it != order.end()) {
        std::rotate(order.begin(), it, order.end());
      }
    }

    return order;
  };

  const auto two_opt_improve
      = [&](std::vector<std::size_t>& order, int max_passes) -> void {
    if (max_passes <= 0 || order.size() < 4) {
      return;
    }

    const int64_t before = path_len(order, scan_in_pts, scan_out_pts);
    bool improved_any = false;

    for (int pass = 0; pass < max_passes; ++pass) {
      bool pass_improved = false;
      for (std::size_t i = 0; i + 3 < order.size(); ++i) {
        for (std::size_t j = i + 2; j + 1 < order.size(); ++j) {
          const auto a = order[i];
          const auto b = order[i + 1];
          const auto c = order[j];
          const auto d = order[j + 1];
          const int64_t old_cost = manhattan(scan_in_pts[a], scan_in_pts[b])
                                   + manhattan(scan_in_pts[c], scan_in_pts[d]);
          const int64_t new_cost = manhattan(scan_in_pts[a], scan_in_pts[c])
                                   + manhattan(scan_in_pts[b], scan_in_pts[d]);
          if (new_cost < old_cost) {
            std::reverse(order.begin() + static_cast<std::ptrdiff_t>(i + 1),
                         order.begin() + static_cast<std::ptrdiff_t>(j + 1));
            pass_improved = true;
            improved_any = true;
          }
        }
      }
      if (!pass_improved) {
        break;
      }
    }

    if (improved_any) {
      const int64_t after = path_len(order, scan_in_pts, scan_out_pts);
      debugPrint(logger,
                 utl::DFT,
                 "scan_chain_opt",
                 1,
                 "OptimizeScanWirelength: 2-opt improved path length {} -> {} "
                 "({} cells)",
                 before,
                 after,
                 order.size());
    }
  };

  int two_opt_passes = 0;
  if (n <= kTwoOptMaxCellsFor3Passes) {
    two_opt_passes = 3;
  } else if (n <= kTwoOptMaxCellsFor2Passes) {
    two_opt_passes = 2;
  } else if (n <= kTwoOptMaxCellsFor1Pass) {
    two_opt_passes = 1;
  }

  // Quadratic heuristics for small/medium chains.
  if (n <= kQuadraticHeuristicMaxCells) {
    std::vector<std::size_t> best_order = greedy_nn_order(start_index);
    two_opt_improve(best_order, two_opt_passes);
    int64_t best_cost = path_len(best_order, scan_in_pts, scan_out_pts);

    if (n <= kFarthestInsertionMaxCells) {
      std::vector<std::size_t> fi_order = farthest_insertion_order(start_index);
      two_opt_improve(fi_order, two_opt_passes);
      const int64_t fi_cost = path_len(fi_order, scan_in_pts, scan_out_pts);
      if (fi_cost < best_cost) {
        best_order = std::move(fi_order);
      }
    }

    std::vector<std::unique_ptr<ScanCell>> ordered;
    ordered.reserve(n);
    for (const std::size_t idx : best_order) {
      ordered.emplace_back(std::move(cells[idx]));
    }
    std::swap(cells, ordered);
    return;
  }

  // Fallback for large chains: rtree greedy ordering with a Manhattan-aware
  // choice among the nearest Euclidean candidates.
  // Get points in a form ready to insert into index
  using Point = bg::model::point<int, 2, bg::cs::cartesian>;
  std::vector<std::pair<Point, size_t>> transformed;

  for (size_t i = 0; i < n; i++) {
    const auto& p = scan_in_pts[i];
    transformed.emplace_back(Point(p.x(), p.y()), i);
  }
  // Update the index
  bgi::rtree<std::pair<Point, size_t>, bgi::rstar<4>> rtree(transformed);
  auto cursor = transformed[start_index];

  // Search nearest neighbours. The rtree nearest search is Euclidean, so we
  // evaluate a handful of nearest Euclidean candidates and pick the best by
  // Manhattan distance (the scan-chain cost proxy we care about).
  std::vector<std::unique_ptr<ScanCell>> ordered;
  ordered.reserve(cells.size());

  ordered.emplace_back(std::move(cells[cursor.second]));
  rtree.remove(cursor);

  while (ordered.size() < cells.size()) {
    bool found = false;
    std::pair<Point, size_t> best = cursor;
    int64_t best_dist = std::numeric_limits<int64_t>::max();
    odb::Point cursor_pt(bg::get<0>(cursor.first), bg::get<1>(cursor.first));

    for (auto it
         = rtree.qbegin(bgi::nearest(cursor.first, kNearestCandidateCount));
         it != rtree.qend();
         ++it) {
      const auto cand = *it;
      odb::Point cand_pt(bg::get<0>(cand.first), bg::get<1>(cand.first));
      const int64_t dist = manhattan(cursor_pt, cand_pt);
      if (!found || dist < best_dist
          || (dist == best_dist && cand.second < best.second)) {
        best = cand;
        best_dist = dist;
        found = true;
      }
    }

    if (!found) {
      no_next_scan_cell();
    }

    cursor = best;
    ordered.emplace_back(std::move(cells[cursor.second]));
    rtree.remove(cursor);
  }

  // Replace with ordered vector.
  std::swap(cells, ordered);
}

}  // namespace dft
