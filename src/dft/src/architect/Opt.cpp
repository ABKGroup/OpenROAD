// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2024-2025, The OpenROAD Authors

#include "Opt.hh"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <limits>
#include <memory>
#include <string_view>
#include <utility>
#include <vector>

#include "ClockDomain.hh"
#include "ScanCell.hh"
#include "boost/geometry/geometries/register/point.hpp"
#include "boost/geometry/geometry.hpp"
#include "boost/geometry/index/rtree.hpp"
#include "odb/dbShape.h"
#include "odb/geom.h"
#include "utl/Logger.h"

namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

namespace dft {

namespace {
constexpr int64_t kInfDistance = std::numeric_limits<int64_t>::max() / 8;

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

int64_t manhattanDist(const odb::Point& a, const odb::Point& b)
{
  const int64_t dx = static_cast<int64_t>(a.x()) - static_cast<int64_t>(b.x());
  const int64_t dy = static_cast<int64_t>(a.y()) - static_cast<int64_t>(b.y());
  return std::abs(dx) + std::abs(dy);
}

int64_t manhattanPointToRectDist(const odb::Point& p, const odb::Rect& r)
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
  return dx + dy;
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

int64_t pinToNetDistance(const odb::Point& pin, const NetAccessGeometry& geom)
{
  int64_t best = kInfDistance;
  for (const odb::Rect& box : geom.boxes) {
    best = std::min(best, manhattanPointToRectDist(pin, box));
    if (best == 0) {
      return 0;
    }
  }
  for (const odb::Point& term : geom.terminals) {
    best = std::min(best, manhattanDist(pin, term));
    if (best == 0) {
      return 0;
    }
  }
  return best == kInfDistance ? 0 : best;
}

void OptimizeScanWirelengthPinToNet(std::vector<std::unique_ptr<ScanCell>>& cells,
                                    const ScanArchitectConfig& config,
                                    utl::Logger* logger)
{
  (void) config;
  const std::size_t n = cells.size();
  if (n < 2) {
    return;
  }

  std::vector<odb::Point> origins;
  origins.reserve(n);
  std::vector<odb::Point> scan_in_pts;
  scan_in_pts.reserve(n);
  std::vector<std::string_view> names;
  names.reserve(n);
  std::vector<NetAccessGeometry> net_geoms;
  net_geoms.reserve(n);

  for (const auto& cell : cells) {
    const odb::Point origin = cell->getOrigin();
    origins.emplace_back(origin);
    names.emplace_back(cell->getName());

    scan_in_pts.emplace_back(scanPinLocation(cell->getScanIn(), origin));
    net_geoms.emplace_back(buildNetAccessGeometry(cell->getScanOut().getNet()));
  }

  std::size_t start_index = 0;
  int64_t lowest = std::numeric_limits<int64_t>::max();
  for (std::size_t i = 0; i < n; ++i) {
    const odb::Point& p = origins[i];
    const int64_t score
        = static_cast<int64_t>(p.x()) + static_cast<int64_t>(p.y());
    if (score < lowest || (score == lowest && names[i] < names[start_index])) {
      start_index = i;
      lowest = score;
    }
  }

  const auto edge_cost = [&](std::size_t src, std::size_t dst) -> int64_t {
    const int64_t dist = pinToNetDistance(scan_in_pts[dst], net_geoms[src]);
    if (dist != 0 || !net_geoms[src].boxes.empty()
        || !net_geoms[src].terminals.empty()) {
      return dist;
    }
    // No routing/pin geometry available; fall back to placement-based distance.
    return manhattanDist(origins[src], origins[dst]);
  };

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

    for (std::size_t i = 0; i < n; ++i) {
      if (used[i]) {
        continue;
      }
      const int64_t cost = edge_cost(cur, i);
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
            + (has_next ? edge_cost(b, next) : 0);
      const int64_t after
          = edge_cost(prev, b) + edge_cost(b, a)
            + (has_next ? edge_cost(a, next) : 0);

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
    OptimizeScanWirelengthPinToNet(cells, config, logger);
    return;
  }

  const auto manhattan = [](const odb::Point& a, const odb::Point& b) -> int64_t {
    const int64_t dx
        = static_cast<int64_t>(a.x()) - static_cast<int64_t>(b.x());
    const int64_t dy
        = static_cast<int64_t>(a.y()) - static_cast<int64_t>(b.y());
    return std::abs(dx) + std::abs(dy);
  };

  const auto path_len = [&](const std::vector<std::size_t>& order,
                            const std::vector<odb::Point>& origins) -> int64_t {
    int64_t total = 0;
    for (std::size_t i = 1; i < order.size(); ++i) {
      total += manhattan(origins[order[i - 1]], origins[order[i]]);
    }
    return total;
  };

  const auto no_next_scan_cell = [&]() -> void {
    logger->error(utl::DFT, 16, "Couldn't find next scan cell to order");
  };

  const std::size_t n = cells.size();
  std::vector<odb::Point> origins;
  origins.reserve(n);
  std::vector<std::string_view> names;
  names.reserve(n);
  for (const auto& cell : cells) {
    origins.emplace_back(cell->getOrigin());
    names.emplace_back(cell->getName());
  }

  // Define the starting node as the lower leftmost, so we don't accidentally
  // start somewhere in the middle
  size_t start_index = 0;
  int64_t lowest_dist = std::numeric_limits<int64_t>::max();

  for (size_t i = 0; i < n; i++) {
    // Find the lower leftmost cell by looking for the cell with the lowest
    // manhattan distance to the origin.
    const auto& origin = origins[i];
    const int64_t dist
        = static_cast<int64_t>(origin.x()) + static_cast<int64_t>(origin.y());
    if (dist < lowest_dist
        || (dist == lowest_dist && names[i] < names[start_index])) {
      start_index = i;
      lowest_dist = dist;
    }
  }

  const auto greedy_nn_order = [&](std::size_t start) -> std::vector<std::size_t> {
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
        const int64_t dist = manhattan(origins[cur], origins[i]);
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

  const auto farthest_insertion_order =
      [&](std::size_t start) -> std::vector<std::size_t> {
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
      const int64_t dist = manhattan(origins[start], origins[i]);
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
      nearest_dist[i] = std::min(manhattan(origins[i], origins[start]),
                                 manhattan(origins[i], origins[farthest]));
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
      int64_t best_delta = manhattan(origins[order.back()], origins[next]);

      for (std::size_t pos = 0; pos + 1 < order.size(); ++pos) {
        const auto a = order[pos];
        const auto b = order[pos + 1];
        const int64_t delta = manhattan(origins[a], origins[next])
                              + manhattan(origins[next], origins[b])
                              - manhattan(origins[a], origins[b]);
        if (delta < best_delta
            || (delta == best_delta && pos + 1 < best_pos)) {
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
        const int64_t dist = manhattan(origins[i], origins[next]);
        if (dist < nearest_dist[i]) {
          nearest_dist[i] = dist;
        }
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

  const auto two_opt_improve = [&](std::vector<std::size_t>& order,
                                  int max_passes) -> void {
    if (max_passes <= 0 || order.size() < 4) {
      return;
    }

    const int64_t before = path_len(order, origins);
    bool improved_any = false;

    for (int pass = 0; pass < max_passes; ++pass) {
      bool pass_improved = false;
      for (std::size_t i = 0; i + 3 < order.size(); ++i) {
        for (std::size_t j = i + 2; j + 1 < order.size(); ++j) {
          const auto a = order[i];
          const auto b = order[i + 1];
          const auto c = order[j];
          const auto d = order[j + 1];
          const int64_t old_cost = manhattan(origins[a], origins[b])
                                  + manhattan(origins[c], origins[d]);
          const int64_t new_cost = manhattan(origins[a], origins[c])
                                  + manhattan(origins[b], origins[d]);
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
      const int64_t after = path_len(order, origins);
      debugPrint(logger,
                 utl::DFT,
                 "scan_chain_opt",
                 1,
                 "OptimizeScanWirelength: 2-opt improved path length {} -> {} ({} cells)",
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
    int64_t best_cost = path_len(best_order, origins);
    two_opt_improve(best_order, two_opt_passes);
    best_cost = path_len(best_order, origins);

    if (n <= kFarthestInsertionMaxCells) {
      std::vector<std::size_t> fi_order = farthest_insertion_order(start_index);
      two_opt_improve(fi_order, two_opt_passes);
      const int64_t fi_cost = path_len(fi_order, origins);
      if (fi_cost < best_cost) {
        best_cost = fi_cost;
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
    const auto& origin = origins[i];
    transformed.emplace_back(Point(origin.x(), origin.y()), i);
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

    for (auto it = rtree.qbegin(bgi::nearest(cursor.first, kNearestCandidateCount));
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
