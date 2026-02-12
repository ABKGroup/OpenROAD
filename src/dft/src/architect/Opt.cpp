// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2024-2025, The OpenROAD Authors

#include "Opt.hh"

#include <algorithm>
#include <array>
#include <chrono>
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
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <vector>

#include "ScanCell.hh"
#include "UclaScanOpt.hh"
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

// ScanOpt local-search operators. These are O(n^2), but ScanOpt is only enabled
// when we can afford building the full O(n^2) cost matrix (kScanOptMaxMatrixCells),
// so allow them up to that same size to avoid leaving obvious long-edge/crossing
// artifacts ("jumps") on larger single-chain designs.
constexpr std::size_t kScanOptTwoOptMaxCells = kScanOptMaxMatrixCells;
constexpr std::size_t kScanOptSwapMaxCells = kScanOptMaxMatrixCells;

int64_t manhattanDist(const odb::Point& a,
                      const odb::Point& b,
                      double vertical_weight);

std::optional<std::chrono::steady_clock::time_point> scanOptDeadline(
    const ScanArchitectConfig& config)
{
  const double limit_s = config.getScanOptTimeLimitSeconds();
  if (limit_s <= 0.0) {
    return std::nullopt;
  }
  const auto dur = std::chrono::duration_cast<std::chrono::steady_clock::duration>(
      std::chrono::duration<double>(limit_s));
  return std::chrono::steady_clock::now() + dur;
}

bool scanOptTimeExpired(
    const std::optional<std::chrono::steady_clock::time_point>& deadline)
{
  return deadline.has_value()
         && std::chrono::steady_clock::now() >= deadline.value();
}

int64_t estimateLocalManhattanScale(const std::vector<odb::Point>& pts,
                                    double vertical_weight)
{
  const std::size_t n = pts.size();
  if (n < 2) {
    return 1;
  }

  using Pt = bg::model::point<int, 2, bg::cs::cartesian>;
  std::vector<std::pair<Pt, std::size_t>> data;
  data.reserve(n);
  for (std::size_t i = 0; i < n; ++i) {
    data.emplace_back(Pt(pts[i].x(), pts[i].y()), i);
  }

  bgi::rtree<std::pair<Pt, std::size_t>, bgi::rstar<4>> rtree(data);

  static constexpr std::size_t kNearest = 8;
  std::vector<int64_t> nearest;
  nearest.reserve(n);

  for (std::size_t i = 0; i < n; ++i) {
    const Pt q(pts[i].x(), pts[i].y());
    int64_t best = kInfDistance;

    for (auto it = rtree.qbegin(bgi::nearest(q, kNearest + 1));
         it != rtree.qend();
         ++it) {
      const auto cand = *it;
      if (cand.second == i) {
        continue;
      }
      best = std::min(best,
                      manhattanDist(pts[i], pts[cand.second], vertical_weight));
    }
    if (best < kInfDistance) {
      nearest.push_back(best);
    }
  }

  if (nearest.empty()) {
    return 1;
  }

  const std::size_t mid = nearest.size() / 2;
  std::nth_element(nearest.begin(), nearest.begin() + mid, nearest.end());
  return std::max<int64_t>(nearest[mid], 1);
}

int64_t jumpPenaltyFromManhattan(int64_t manhattan, int64_t local_scale)
{
  if (manhattan <= 0 || local_scale <= 0) {
    return 0;
  }

  // Penalize long edges superlinearly to avoid visually/physically "jumpy"
  // chains. This acts as a soft proxy for minimizing maximum hop length while
  // remaining a simple additive edge objective.
  // Smaller divisor => stronger penalty. Favor suppressing long "jump" edges
  // even at the expense of some total wirelength, since large hops are both
  // visually obvious and often physically undesirable.
  static constexpr int64_t kPenaltyDivisor = 1;
  static constexpr __int128 kPenaltyMultiplier = 16;
  const __int128 num
      = static_cast<__int128>(manhattan) * manhattan * kPenaltyMultiplier;
  const __int128 den
      = static_cast<__int128>(local_scale) * kPenaltyDivisor;
  if (den <= 0) {
    return 0;
  }

  const __int128 q = num / den;
  if (q >= static_cast<__int128>(kInfDistance)) {
    return kInfDistance;
  }
  return static_cast<int64_t>(q);
}

// If chain endpoints are not fixed, the chain "break" can be chosen freely.
// Rotating the ordering to drop the most expensive edge (i -> i+1) can improve
// both total wirelength proxy and worst-edge proxy, without changing relative
// order (only the start/end).
template <typename EdgeCost>
void rotateOrderToDropWorstEdge(std::vector<std::size_t>& order,
                                const std::vector<std::string_view>& names,
                                EdgeCost edge_cost)
{
  const std::size_t n = order.size();
  if (n < 2) {
    return;
  }

  // Consider the "closure" edge from end->start. The original ordering is
  // equivalent to choosing this closure as the break (not present in the path).
  const int64_t closure_cost = edge_cost(order.back(), order.front());

  int64_t max_internal_cost = -1;
  std::vector<std::size_t> candidates;
  candidates.reserve(4);
  for (std::size_t i = 0; i + 1 < n; ++i) {
    const int64_t cost = edge_cost(order[i], order[i + 1]);
    if (cost > max_internal_cost) {
      max_internal_cost = cost;
      candidates.clear();
      candidates.push_back(i);
    } else if (cost == max_internal_cost) {
      candidates.push_back(i);
    }
  }

  // Only rotate if it strictly improves the objective (removing an internal
  // edge larger than the would-be closure).
  if (max_internal_cost <= closure_cost || candidates.empty()) {
    return;
  }

  // Deterministic tie-break: pick the cut that yields the lexicographically
  // smallest new start cell name.
  std::size_t best = candidates.front();
  std::string_view best_start = names[order[best + 1]];
  for (std::size_t k = 1; k < candidates.size(); ++k) {
    const std::size_t i = candidates[k];
    const std::string_view start = names[order[i + 1]];
    if (start < best_start) {
      best_start = start;
      best = i;
    }
  }

  std::rotate(order.begin(),
              order.begin() + static_cast<std::ptrdiff_t>(best + 1),
              order.end());
}

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
                   const odb::Rect bbox = iterm->getBBox();
                   return odb::Point(bbox.xMin(), bbox.yMin());
                 },
                 [&](odb::dbBTerm* bterm) -> odb::Point {
                   if (bterm == nullptr) {
                     return fallback;
                   }
                   const odb::Rect bbox = bterm->getBBox();
                   return odb::Point(bbox.xMin(), bbox.yMin());
                 }},
      pin.getValue());
}

namespace {
std::string UnescapeSlash(std::string_view input)
{
  std::string out;
  out.reserve(input.size());
  for (std::size_t i = 0; i < input.size(); ++i) {
    const char c = input[i];
    if (c == '\\' && i + 1 < input.size() && input[i + 1] == '/') {
      continue;  // drop the escape, keep '/'
    }
    out.push_back(c);
  }
  return out;
}

std::pair<std::string, std::optional<std::string>> SplitTermIdentifier(
    std::string_view input)
{
  std::size_t tracker = 0;
  std::size_t slash_position;
  while ((slash_position = input.find('/', tracker)) != std::string_view::npos) {
    if (slash_position != 0 && input[slash_position - 1] == '\\') {
      tracker = slash_position + 1;
      continue;
    }
    const std::string inst = UnescapeSlash(input.substr(0, slash_position));
    const std::string pin = UnescapeSlash(input.substr(slash_position + 1));
    return {inst, pin};
  }
  return {UnescapeSlash(input), std::nullopt};
}

odb::dbBlock* InferBlockFromPlacedCells(
    const std::vector<std::unique_ptr<ScanCell>>& cells)
{
  if (cells.empty()) {
    return nullptr;
  }

  // Any scan pin on a scan cell should be an ITerm, and can be used to recover
  // the owning block. This is only used to resolve endpoint terminal names.
  odb::dbBlock* block = nullptr;
  const ScanLoad scan_in = cells.front()->getScanIn();
  std::visit(
      [&](auto&& term) {
        if (term == nullptr) {
          return;
        }
        using T = std::decay_t<decltype(term)>;
        if constexpr (std::is_same_v<T, odb::dbITerm*>) {
          if (odb::dbInst* inst = term->getInst()) {
            block = inst->getBlock();
          }
        } else if constexpr (std::is_same_v<T, odb::dbBTerm*>) {
          block = term->getBlock();
        }
      },
      scan_in.getValue());

  return block;
}

std::optional<odb::Point> ResolveEndpointTerm(odb::dbBlock* block,
                                              std::string_view term,
                                              utl::Logger* logger)
{
  if (block == nullptr) {
    return std::nullopt;
  }

  const auto term_info = SplitTermIdentifier(term);
  if (term_info.second.has_value()) {
    odb::dbInst* inst = block->findInst(term_info.first.c_str());
    if (inst == nullptr) {
      if (logger) {
        logger->warn(utl::DFT,
                     210,
                     "Scan constraints: endpoint instance '{}' not found "
                     "(term '{}'); ignoring endpoint.",
                     term_info.first,
                     term);
      }
      return std::nullopt;
    }
    odb::dbITerm* iterm = inst->findITerm(term_info.second->c_str());
    if (iterm == nullptr) {
      if (logger) {
        logger->warn(
            utl::DFT,
            211,
            "Scan constraints: endpoint iterm '{}/{}' not found; ignoring "
            "endpoint.",
            term_info.first,
            *term_info.second);
      }
      return std::nullopt;
    }
    const odb::Rect bbox = iterm->getBBox();
    return odb::Point(bbox.xMin(), bbox.yMin());
  }

  odb::dbBTerm* bterm = block->findBTerm(term_info.first.c_str());
  if (bterm == nullptr) {
    if (logger) {
      logger->warn(utl::DFT,
                   212,
                   "Scan constraints: endpoint port '{}' not found; ignoring "
                   "endpoint.",
                   term_info.first);
    }
    return std::nullopt;
  }
  const odb::Rect bbox = bterm->getBBox();
  return odb::Point(bbox.xMin(), bbox.yMin());
}

std::optional<odb::Point> EndpointPoint(
    const ScanArchitectConfig::ChainEndpoint& endpoint,
    odb::dbBlock* block,
    utl::Logger* logger)
{
  using Type = ScanArchitectConfig::ChainEndpoint::Type;
  switch (endpoint.type) {
    case Type::Point:
      return odb::Point(endpoint.point.x, endpoint.point.y);
    case Type::Term:
      return ResolveEndpointTerm(block, endpoint.term, logger);
    default:
      return std::nullopt;
  }
}
}  // namespace

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
      const odb::Rect bbox = iterm->getBBox();
      geom.terminals.emplace_back(bbox.xMin(), bbox.yMin());
    }
    for (odb::dbBTerm* bterm : net->getBTerms()) {
      const odb::Rect bbox = bterm->getBBox();
      geom.terminals.emplace_back(bbox.xMin(), bbox.yMin());
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

std::vector<std::size_t> scanOptGreedyOrderEndFixed(
    const ScanOptMatrix& m,
    std::size_t start,
    std::size_t fixed_end,
    const std::vector<std::string_view>& names,
    utl::Logger* logger)
{
  const std::size_t n = m.n;
  if (n < 2) {
    return {};
  }
  if (fixed_end >= n) {
    logger->error(utl::DFT,
                  206,
                  "Internal error: fixed_end index {} out of range for n={}",
                  fixed_end,
                  n);
  }
  if (start == fixed_end) {
    logger->error(utl::DFT,
                  207,
                  "Internal error: start index equals fixed_end index {}",
                  start);
  }

  std::vector<std::size_t> order;
  order.reserve(n);
  std::vector<bool> used(n, false);

  used[fixed_end] = true;

  std::size_t cur = start;
  used[cur] = true;
  order.push_back(cur);

  while (order.size() + 1 < n) {  // leave room for fixed_end
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
      logger->error(utl::DFT, 208, "Couldn't find next scan cell to order");
    }

    used[best] = true;
    order.push_back(best);
    cur = best;
  }

  order.push_back(fixed_end);
  return order;
}

bool scanOptAdjacentSwapImprove(const ScanOptMatrix& m,
                                const std::vector<std::string_view>& names,
                                std::vector<std::size_t>& order,
                                bool end_fixed)
{
  bool improved = false;
  const std::size_t n = order.size();
  if (n < 3) {
    return false;
  }
  const std::size_t limit = end_fixed ? (n - 1) : n;
  for (std::size_t pos = 1; pos + 1 < limit; ++pos) {
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
                             std::vector<std::size_t>& order,
                             bool end_fixed)
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
  const std::size_t last_movable = end_fixed ? (n - 2) : (n - 1);
  for (std::size_t from = 1; from <= last_movable; ++from) {
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

    for (std::size_t after = 0; after <= last_movable; ++after) {
      if (after == from || after + 1 == from) {
        continue;  // no-op or invalid
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

bool scanOptBestSegmentRelocateMove(const ScanOptMatrix& m,
                                    const std::vector<std::string_view>& names,
                                    std::vector<std::size_t>& order,
                                    bool end_fixed,
                                    std::size_t seg_len)
{
  const std::size_t n = order.size();
  if (seg_len < 2 || n < seg_len + 2) {
    return false;
  }

  bool found = false;
  int64_t best_delta = 0;
  std::size_t best_from = 0;
  std::size_t best_after = 0;
  std::string_view best_name;

  const std::size_t last_movable = end_fixed ? (n - 2) : (n - 1);
  if (last_movable < seg_len) {
    return false;
  }

  // Keep the start fixed (from starts at 1).
  for (std::size_t from = 1; from + seg_len - 1 <= last_movable; ++from) {
    const std::size_t to = from + seg_len - 1;
    const std::size_t prev = order[from - 1];
    const std::size_t first = order[from];
    const std::size_t last = order[to];
    const bool has_next = (to + 1 < n);
    const std::size_t next = has_next ? order[to + 1] : 0;

    const int64_t remove_before
        = static_cast<int64_t>(m.get(prev, first))
          + (has_next ? static_cast<int64_t>(m.get(last, next)) : 0);
    const int64_t remove_after
        = has_next ? static_cast<int64_t>(m.get(prev, next)) : 0;
    const int64_t remove_delta = remove_after - remove_before;

    for (std::size_t after = 0; after <= last_movable; ++after) {
      // Skip no-op and illegal insertions (within or immediately adjacent to
      // the segment).
      if (after + 1 >= from && after <= to) {
        continue;
      }

      const std::size_t ins_prev = order[after];
      const bool has_ins_next = (after + 1 < n);
      const std::size_t ins_next = has_ins_next ? order[after + 1] : 0;

      const int64_t insert_before
          = has_ins_next ? static_cast<int64_t>(m.get(ins_prev, ins_next)) : 0;
      const int64_t insert_after
          = static_cast<int64_t>(m.get(ins_prev, first))
            + (has_ins_next ? static_cast<int64_t>(m.get(last, ins_next)) : 0);
      const int64_t insert_delta = insert_after - insert_before;

      const int64_t delta = remove_delta + insert_delta;
      if (delta < best_delta
          || (delta == best_delta && found
              && (names[first] < best_name
                  || (names[first] == best_name
                      && (from < best_from
                          || (from == best_from && after < best_after)))))) {
        found = true;
        best_delta = delta;
        best_from = from;
        best_after = after;
        best_name = names[first];
      }
    }
  }

  if (!found || best_delta >= 0) {
    return false;
  }

  const std::size_t best_to = best_from + seg_len - 1;
  std::vector<std::size_t> segment(
      order.begin() + static_cast<std::ptrdiff_t>(best_from),
      order.begin() + static_cast<std::ptrdiff_t>(best_to + 1));
  order.erase(order.begin() + static_cast<std::ptrdiff_t>(best_from),
              order.begin() + static_cast<std::ptrdiff_t>(best_to + 1));

  std::size_t after = best_after;
  if (after > best_to) {
    after -= seg_len;
  }
  const std::size_t insert_index = after + 1;
  order.insert(order.begin() + static_cast<std::ptrdiff_t>(insert_index),
               segment.begin(),
               segment.end());

  return true;
}

bool scanOptBestSwapMove(const ScanOptMatrix& m,
                         const std::vector<std::string_view>& names,
                         std::vector<std::size_t>& order,
                         bool end_fixed)
{
  const std::size_t n = order.size();
  if (n < 4) {
    return false;
  }

  bool found = false;
  int64_t best_delta = 0;
  std::size_t best_i = 0;
  std::size_t best_j = 0;
  const std::size_t j_limit = end_fixed ? (n - 1) : n;

  // Keep the start fixed (i starts at 1).
  for (std::size_t i = 1; i + 1 < n; ++i) {
    for (std::size_t j = i + 1; j < j_limit; ++j) {
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

bool scanOptBestTwoOptMove(
    const ScanOptMatrix& m,
    const std::vector<std::string_view>& names,
    std::vector<std::size_t>& order,
    bool end_fixed,
    const std::optional<std::chrono::steady_clock::time_point>& deadline)
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

  bool found = false;
  int64_t best_delta = 0;
  std::size_t best_i = 0;
  std::size_t best_j = 0;

  const std::size_t j_limit = end_fixed ? (n - 1) : n;
  for (std::size_t i = 0; i + 2 < j_limit; ++i) {
    // Poll the deadline periodically; this loop can be O(n^2).
    if ((i & 63U) == 0U && scanOptTimeExpired(deadline)) {
      break;
    }
    for (std::size_t j = i + 2; j < j_limit; ++j) {
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
      const int64_t delta = after - before;

      if (delta < best_delta
          || (delta == best_delta && found
              && (names[c] < names[order[best_j]]
                  || (names[c] == names[order[best_j]]
                      && names[a] < names[order[best_i]])))) {
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

  std::reverse(order.begin() + static_cast<std::ptrdiff_t>(best_i + 1),
               order.begin() + static_cast<std::ptrdiff_t>(best_j + 1));
  return true;
}

// Targeted 2-opt that tries to break the current worst edge (by ScanOpt cost),
// rather than doing a full O(n^2) best-improvement search. This is a cheap way
// to reduce visually-obvious long hops on large chains under a time limit.
bool scanOptWorstEdgeTwoOptImprove(
    const ScanOptMatrix& m,
    const std::vector<std::string_view>& names,
    std::vector<std::size_t>& order,
    bool end_fixed,
    const std::optional<std::chrono::steady_clock::time_point>& deadline)
{
  const std::size_t n = order.size();
  if (n < 4) {
    return false;
  }

  // j indexes the start of the second edge (c->d), so j+1 must be valid.
  // Exclude the last node in the open-path case to avoid out-of-bounds.
  const std::size_t j_limit = n - 1;

  std::size_t worst_i = 0;
  int64_t worst_cost = -1;
  bool have_worst = false;
  for (std::size_t i = 0; i + 2 < j_limit; ++i) {
    const int64_t cost = m.get(order[i], order[i + 1]);
    if (!have_worst || cost > worst_cost) {
      have_worst = true;
      worst_cost = cost;
      worst_i = i;
    }
  }
  if (!have_worst) {
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

  const std::size_t i = worst_i;
  bool found = false;
  int64_t best_delta = 0;
  std::size_t best_j = 0;

  for (std::size_t j = i + 2; j < j_limit; ++j) {
    if ((j & 255U) == 0U && scanOptTimeExpired(deadline)) {
      break;
    }

    const std::size_t a = order[i];
    const std::size_t b = order[i + 1];
    const std::size_t c = order[j];
    const std::size_t d = order[j + 1];

    const int64_t old_inside = segment_forward(i + 1, j);
    const int64_t new_inside = segment_reversed(i + 1, j);

    const int64_t old_edges
        = static_cast<int64_t>(m.get(a, b)) + static_cast<int64_t>(m.get(c, d));
    const int64_t new_edges
        = static_cast<int64_t>(m.get(a, c)) + static_cast<int64_t>(m.get(b, d));

    const int64_t before = old_edges + old_inside;
    const int64_t after = new_edges + new_inside;
    const int64_t delta = after - before;

    if (delta < best_delta
        || (delta == best_delta && found
            && (names[c] < names[order[best_j]]
                || (names[c] == names[order[best_j]]
                    && names[a] < names[order[i]])))) {
      found = true;
      best_delta = delta;
      best_j = j;
    }
  }

  if (!found || best_delta >= 0) {
    return false;
  }

  std::reverse(order.begin() + static_cast<std::ptrdiff_t>(i + 1),
               order.begin() + static_cast<std::ptrdiff_t>(best_j + 1));
  return true;
}

// Targeted Or-opt style segment relocation focused on the current worst edge.
// This is a direction-preserving move (no internal reversal) and can use larger
// segments than scanOptBestSegmentRelocateMove without O(n^2) work.
bool scanOptWorstEdgeSegmentRelocateImprove(
    const ScanOptMatrix& m,
    const std::vector<std::string_view>& names,
    std::vector<std::size_t>& order,
    bool end_fixed,
    const std::optional<std::chrono::steady_clock::time_point>& deadline,
    std::size_t max_seg_len)
{
  const std::size_t n = order.size();
  if (n < 4 || max_seg_len < 2) {
    return false;
  }

  const std::size_t j_limit = end_fixed ? (n - 1) : n;
  if (j_limit < 4) {
    return false;
  }

  std::size_t worst_i = 0;
  int64_t worst_cost = -1;
  bool have_worst = false;
  for (std::size_t i = 0; i + 2 < j_limit; ++i) {
    const int64_t cost = m.get(order[i], order[i + 1]);
    if (!have_worst || cost > worst_cost) {
      have_worst = true;
      worst_cost = cost;
      worst_i = i;
    }
  }
  if (!have_worst) {
    return false;
  }

  // Try relocating a segment that starts at the head of the worst edge. This
  // guarantees we attempt to remove that edge from the path.
  const std::size_t from = worst_i + 1;
  const std::size_t last_movable = end_fixed ? (n - 2) : (n - 1);
  if (from == 0 || from > last_movable) {
    return false;
  }

  bool found = false;
  int64_t best_delta = 0;
  std::size_t best_to = 0;
  std::size_t best_after = 0;
  std::string_view best_name;

  const std::size_t seg_limit = std::min(max_seg_len, last_movable - from + 1);
  for (std::size_t seg_len = 2; seg_len <= seg_limit; ++seg_len) {
    if (scanOptTimeExpired(deadline)) {
      break;
    }

    const std::size_t to = from + seg_len - 1;
    const std::size_t prev = order[from - 1];
    const std::size_t first = order[from];
    const std::size_t last = order[to];
    const bool has_next = (to + 1 < n);
    const std::size_t next = has_next ? order[to + 1] : 0;

    const int64_t remove_before = static_cast<int64_t>(m.get(prev, first));
    const int64_t remove_after
        = has_next ? static_cast<int64_t>(m.get(last, next)) : 0;
    const int64_t remove_join
        = has_next ? static_cast<int64_t>(m.get(prev, next)) : 0;
    const int64_t remove_delta = remove_join - remove_before - remove_after;

    for (std::size_t after = 0; after <= last_movable; ++after) {
      if ((after & 255U) == 0U && scanOptTimeExpired(deadline)) {
        break;
      }

      // Disallow inserting inside the segment or immediately before it (no-op).
      if (after >= from - 1 && after <= to) {
        continue;
      }

      const std::size_t ins_prev = order[after];
      const bool has_ins_next = (after + 1 < n);
      const std::size_t ins_next = has_ins_next ? order[after + 1] : 0;

      const int64_t insert_before
          = has_ins_next ? static_cast<int64_t>(m.get(ins_prev, ins_next)) : 0;
      const int64_t insert_after
          = static_cast<int64_t>(m.get(ins_prev, first))
            + (has_ins_next ? static_cast<int64_t>(m.get(last, ins_next)) : 0);
      const int64_t insert_delta = insert_after - insert_before;

      const int64_t delta = remove_delta + insert_delta;
      if (delta < best_delta
          || (delta == best_delta && found
              && (names[first] < best_name
                  || (names[first] == best_name
                      && (to < best_to
                          || (to == best_to && after < best_after)))))) {
        found = true;
        best_delta = delta;
        best_to = to;
        best_after = after;
        best_name = names[first];
      }
    }
  }

  if (!found || best_delta >= 0) {
    return false;
  }

  std::vector<std::size_t> segment(
      order.begin() + static_cast<std::ptrdiff_t>(from),
      order.begin() + static_cast<std::ptrdiff_t>(best_to + 1));
  order.erase(order.begin() + static_cast<std::ptrdiff_t>(from),
              order.begin() + static_cast<std::ptrdiff_t>(best_to + 1));

  std::size_t after = best_after;
  if (after > best_to) {
    after -= segment.size();
  }
  const std::size_t insert_index = after + 1;
  order.insert(order.begin() + static_cast<std::ptrdiff_t>(insert_index),
               segment.begin(),
               segment.end());
  return true;
}

// Direction-preserving 3-opt "subtour swap" focused on the current worst edge.
// This swaps two adjacent segments [a..b-1] and [b..c] into [b..c][a..b-1],
// changing only the boundary edges (no internal reversal).
//
// This move corresponds to the direction-preserving 3-opt reconnection from
// Boese/Kahng/Tsayy (1994) and can be effective at removing large "jump" edges
// that remain after 2-opt / Or-opt style relocations.
bool scanOptWorstEdgeSegmentSwapImprove(
    const ScanOptMatrix& m,
    const std::vector<std::string_view>& names,
    std::vector<std::size_t>& order,
    bool end_fixed,
    const std::optional<std::chrono::steady_clock::time_point>& deadline,
    std::size_t max_seg_len)
{
  const std::size_t n = order.size();
  if (n < 4 || max_seg_len == 0) {
    return false;
  }

  const std::size_t last_movable = end_fixed ? (n - 2) : (n - 1);
  if (last_movable < 2) {
    return false;
  }

  // Choose the worst "swappable" edge. We keep order[0] fixed, so the boundary
  // index b must be >= 2 (so segment1 starts at a >= 1). Also require segment2
  // to start at b <= last_movable.
  std::size_t worst_i = 0;
  int64_t worst_cost = -1;
  bool have_worst = false;
  for (std::size_t i = 1; i + 1 <= last_movable; ++i) {
    const int64_t cost = m.get(order[i], order[i + 1]);
    if (!have_worst || cost > worst_cost) {
      have_worst = true;
      worst_cost = cost;
      worst_i = i;
    }
  }
  if (!have_worst) {
    return false;
  }

  // Boundary between the two segments (edge X->B).
  const std::size_t b = worst_i + 1;
  if (b < 2 || b > last_movable) {
    return false;
  }

  const std::size_t max_left = std::min(max_seg_len, b - 1);  // len1 <= b-1
  const std::size_t max_right
      = std::min(max_seg_len, last_movable - b + 1);  // len2 <= ...
  if (max_left == 0 || max_right == 0) {
    return false;
  }

  bool found = false;
  int64_t best_delta = 0;
  std::size_t best_len1 = 0;
  std::size_t best_len2 = 0;
  std::string_view best_key;

  const std::size_t x = order[b - 1];
  const std::size_t bnode = order[b];

  for (std::size_t len1 = 1; len1 <= max_left; ++len1) {
    if (scanOptTimeExpired(deadline)) {
      break;
    }
    const std::size_t a = b - len1;
    if (a < 1) {
      continue;  // keep start fixed
    }
    const std::size_t p = order[a - 1];
    const std::size_t anode = order[a];

    for (std::size_t len2 = 1; len2 <= max_right; ++len2) {
      if ((len2 & 1023U) == 0U && scanOptTimeExpired(deadline)) {
        break;
      }

      const std::size_t c = b + len2 - 1;
      if (c < b || c > last_movable) {
        continue;
      }

      const std::size_t cnode = order[c];
      const bool has_n = (c + 1 < n);
      const std::size_t nnode = has_n ? order[c + 1] : 0;

      const int64_t before = static_cast<int64_t>(m.get(p, anode))
                             + static_cast<int64_t>(m.get(x, bnode))
                             + (has_n ? static_cast<int64_t>(m.get(cnode, nnode))
                                      : 0);
      const int64_t after = static_cast<int64_t>(m.get(p, bnode))
                            + static_cast<int64_t>(m.get(cnode, anode))
                            + (has_n ? static_cast<int64_t>(m.get(x, nnode))
                                     : 0);
      const int64_t delta = after - before;

      // Deterministic tie-break for reproducibility.
      const std::string_view key = names[bnode];
      if (delta < best_delta
          || (delta == best_delta && found
              && (key < best_key
                  || (key == best_key
                      && (len1 < best_len1
                          || (len1 == best_len1 && len2 < best_len2)))))) {
        found = true;
        best_delta = delta;
        best_len1 = len1;
        best_len2 = len2;
        best_key = key;
      }
    }
  }

  if (!found || best_delta >= 0) {
    return false;
  }

  const std::size_t a = b - best_len1;
  const std::size_t c = b + best_len2 - 1;
  std::rotate(order.begin() + static_cast<std::ptrdiff_t>(a),
              order.begin() + static_cast<std::ptrdiff_t>(b),
              order.begin() + static_cast<std::ptrdiff_t>(c + 1));
  return true;
}

void scanOptDescent(const ScanOptMatrix& m,
                    const std::vector<std::string_view>& names,
                    std::vector<std::size_t>& order,
                    bool end_fixed,
                    const std::optional<std::chrono::steady_clock::time_point>&
                        deadline)
{
  // Iterate local-search passes until convergence or time limit.
  for (int pass = 0; pass < 500; ++pass) {
    if (scanOptTimeExpired(deadline)) {
      break;
    }
    bool improved = scanOptAdjacentSwapImprove(m, names, order, end_fixed);
    if (scanOptTimeExpired(deadline)) {
      break;
    }
    improved |= scanOptBestRelocateMove(m, names, order, end_fixed);
    if (scanOptTimeExpired(deadline)) {
      break;
    }
    improved |= scanOptWorstEdgeTwoOptImprove(m, names, order, end_fixed, deadline);
    if (scanOptTimeExpired(deadline)) {
      break;
    }
    improved |= scanOptWorstEdgeSegmentRelocateImprove(
        m, names, order, end_fixed, deadline, 256);
    if (scanOptTimeExpired(deadline)) {
      break;
    }
    improved |= scanOptWorstEdgeSegmentSwapImprove(
        m, names, order, end_fixed, deadline, 256);
    if (scanOptTimeExpired(deadline)) {
      break;
    }
    // Or-opt: relocate short segments to improve locality and reduce long
    // "jump" edges (direction-preserving 3-opt neighborhood).
    improved |= scanOptBestSegmentRelocateMove(m, names, order, end_fixed, 10);
    if (scanOptTimeExpired(deadline)) {
      break;
    }
    improved |= scanOptBestSegmentRelocateMove(m, names, order, end_fixed, 8);
    if (scanOptTimeExpired(deadline)) {
      break;
    }
    improved |= scanOptBestSegmentRelocateMove(m, names, order, end_fixed, 6);
    if (scanOptTimeExpired(deadline)) {
      break;
    }
    improved |= scanOptBestSegmentRelocateMove(m, names, order, end_fixed, 5);
    if (scanOptTimeExpired(deadline)) {
      break;
    }
    improved |= scanOptBestSegmentRelocateMove(m, names, order, end_fixed, 4);
    if (scanOptTimeExpired(deadline)) {
      break;
    }
    improved |= scanOptBestSegmentRelocateMove(m, names, order, end_fixed, 3);
    if (scanOptTimeExpired(deadline)) {
      break;
    }
    improved |= scanOptBestSegmentRelocateMove(m, names, order, end_fixed, 2);
    if (scanOptTimeExpired(deadline)) {
      break;
    }
    if (order.size() <= kScanOptSwapMaxCells) {
      improved |= scanOptBestSwapMove(m, names, order, end_fixed);
    }
    if (scanOptTimeExpired(deadline)) {
      break;
    }
    if (order.size() <= kScanOptTwoOptMaxCells) {
      improved |= scanOptBestTwoOptMove(m, names, order, end_fixed, deadline);
    }
    if (!improved) {
      break;
    }
  }
}

void scanOptDoubleBridgeKick(std::vector<std::size_t>& order,
                             std::mt19937_64& rng,
                             bool end_fixed)
{
  const std::size_t n = order.size();
  if (n < 8) {
    if (end_fixed && n > 2) {
      std::shuffle(order.begin() + 1, order.end() - 1, rng);
    } else {
      std::shuffle(order.begin() + 1, order.end(), rng);
    }
    return;
  }

  if (!end_fixed) {
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
    return;
  }

  // End-fixed kick: preserve order[0] (start) and order[n-1] (fixed end).
  if (n < 9) {
    std::shuffle(order.begin() + 1, order.end() - 1, rng);
    return;
  }

  std::vector<std::size_t> core(order.begin() + 1, order.end() - 1);
  const std::size_t m = core.size();
  if (m < 8) {
    std::shuffle(core.begin(), core.end(), rng);
    std::vector<std::size_t> kicked;
    kicked.reserve(n);
    kicked.push_back(order.front());
    kicked.insert(kicked.end(), core.begin(), core.end());
    kicked.push_back(order.back());
    order.swap(kicked);
    return;
  }

  std::uniform_int_distribution<std::size_t> dist(1, m - 1);
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
  kicked.push_back(order.front());
  kicked.insert(kicked.end(), core.begin(), core.begin() + a);
  kicked.insert(kicked.end(), core.begin() + b, core.begin() + c);
  kicked.insert(kicked.end(), core.begin() + a, core.begin() + b);
  kicked.insert(kicked.end(), core.begin() + c, core.begin() + d);
  kicked.insert(kicked.end(), core.begin() + d, core.end());
  kicked.push_back(order.back());

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

  const bool end_fixed = static_cast<bool>(terminal_cost);
  const std::size_t scanopt_n = end_fixed ? (n + 1) : n;
  const bool use_scanopt
      = (config.getScanOrderSolver()
             == ScanArchitectConfig::ScanOrderSolver::ScanOpt
         && scanopt_n <= kScanOptMaxMatrixCells);

  if (use_scanopt) {
    const auto deadline = scanOptDeadline(config);

    ScanOptMatrix m;
    m.n = scanopt_n;
    m.costs.resize(scanopt_n * scanopt_n);
    for (std::size_t i = 0; i < scanopt_n; ++i) {
      for (std::size_t j = 0; j < scanopt_n; ++j) {
        if (i == j) {
          m.costs[i * scanopt_n + j] = kScanOptLargeCost;
          continue;
        }
        int64_t cost = kInfDistance;
        if (end_fixed && i == n) {
          cost = kInfDistance;
        } else if (end_fixed && j == n) {
          cost = terminal_cost(i);
        } else {
          cost = cost_fn(i, j);
        }
        m.costs[i * scanopt_n + j]
            = static_cast<int32_t>(std::min<int64_t>(cost, kScanOptLargeCost));
      }
    }

    std::vector<std::string_view> scanopt_names = names;
    static constexpr std::string_view kEndName = "__end__";
    if (end_fixed) {
      scanopt_names.push_back(kEndName);
    }

    std::vector<std::size_t> best_order
        = end_fixed ? scanOptGreedyOrderEndFixed(m,
                                                 start,
                                                 n,
                                                 scanopt_names,
                                                 logger)
                    : scanOptGreedyOrder(m, start, scanopt_names, logger);
    scanOptDescent(m, scanopt_names, best_order, end_fixed, deadline);
    int64_t best_cost = scanOptPathCost(m, best_order);

    std::mt19937_64 rng(config.getScanOptSeed());
    const uint64_t rounds = config.getScanOptRounds();
    if (!config.getScanOptTempControl()) {
      for (uint64_t r = 0; r < rounds; ++r) {
        if (scanOptTimeExpired(deadline)) {
          break;
        }
        std::vector<std::size_t> cand = best_order;
        scanOptDoubleBridgeKick(cand, rng, end_fixed);
        scanOptDescent(m, scanopt_names, cand, end_fixed, deadline);
        const int64_t cand_cost = scanOptPathCost(m, cand);
        if (cand_cost < best_cost) {
          best_cost = cand_cost;
          best_order.swap(cand);
        }
      }
    } else {
      constexpr int kNoImproveBeforeTemp = 3;
      constexpr int kTempSteps = 3;

      std::vector<std::size_t> cur_order = best_order;
      int64_t cur_cost = best_cost;

      double temperature = 0.0;
      int no_improve = 0;
      int temp_steps_left = 0;

      const double t_div = config.getScanOptTDiv();
      std::uniform_real_distribution<double> u01(0.0, 1.0);

      for (uint64_t r = 0; r < rounds; ++r) {
        if (scanOptTimeExpired(deadline)) {
          break;
        }
        std::vector<std::size_t> cand = cur_order;
        scanOptDoubleBridgeKick(cand, rng, end_fixed);
        scanOptDescent(m, scanopt_names, cand, end_fixed, deadline);
        const int64_t cand_cost = scanOptPathCost(m, cand);

        const int64_t delta = cand_cost - cur_cost;
        if (delta < 0) {
          cur_order.swap(cand);
          cur_cost = cand_cost;
          temperature = 0.0;
          no_improve = 0;
          temp_steps_left = 0;
        } else {
          no_improve++;
          if (temperature == 0.0 && no_improve >= kNoImproveBeforeTemp) {
            temperature = static_cast<double>(cur_cost) / t_div;
            temp_steps_left = kTempSteps;
            no_improve = 0;
          }

          bool accept = false;
          if (temperature > 0.0) {
            const double prob
                = std::exp(-static_cast<double>(delta) / temperature);
            accept = u01(rng) < prob;
            if (--temp_steps_left <= 0) {
              temperature = 0.0;
              no_improve = 0;
            }
          }
          if (accept) {
            cur_order.swap(cand);
            cur_cost = cand_cost;
          }
        }

        if (cur_cost < best_cost) {
          best_cost = cur_cost;
          best_order = cur_order;
        }
      }
    }

    if (end_fixed) {
      if (!best_order.empty() && best_order.back() == n) {
        best_order.pop_back();
      } else {
        logger->warn(utl::DFT,
                     209,
                     "Internal error: expected fixed-end node at end of "
                     "ScanOpt order; leaving order unchanged.");
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
      const odb::Point& p = scan_in_pts[seg_entry[si]];
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
          const odb::Point& p = scan_in_pts[components[ci].entry];
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
            const odb::Point& p = scan_in_pts[c.entry];
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
  odb::dbBlock* endpoint_block = nullptr;
  if (endpoints.has_value()) {
    if ((endpoints->begin.has_value()
         && endpoints->begin->type
                == ScanArchitectConfig::ChainEndpoint::Type::Term)
        || (endpoints->end.has_value()
            && endpoints->end->type
                   == ScanArchitectConfig::ChainEndpoint::Type::Term)) {
      endpoint_block = InferBlockFromPlacedCells(cells);
    }
  }
  if (endpoints.has_value()) {
    if (endpoints->begin.has_value()) {
      begin = EndpointPoint(*endpoints->begin, endpoint_block, logger);
    }
    if (endpoints->end.has_value()) {
      end = EndpointPoint(*endpoints->end, endpoint_block, logger);
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
  const int64_t local_scale
      = estimateLocalManhattanScale(scan_in_pts, vertical_weight);

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
      const odb::Point& p = scan_in_pts[i];
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
    // Routing-aware + placement-aware hybrid cost:
    // - pin-to-net distance approximates incremental branch length to attach a
    //   new scan-in load onto the already-routed scan-out (often Q) net, and
    // - direct pin-to-pin Manhattan encourages spatial locality to avoid
    //   visually/physically "jumpy" chains when many candidates tie at 0 in the
    //   pin-to-net metric (e.g., long functional nets spanning the core).
    const int64_t p2n
        = pinToNetDistance(scan_in_pts[dst], net_geoms[src], vertical_weight);
    const int64_t man
        = manhattanDist(scan_out_pts[src], scan_in_pts[dst], vertical_weight);
    const int64_t jump_pen = jumpPenaltyFromManhattan(man, local_scale);
    if (p2n != 0 || !net_geoms[src].boxes.empty()
        || !net_geoms[src].terminals.empty()) {
      return scaleEdgeCost(p2n + man + jump_pen, timing_mul[src]);
    }
    // No routing/pin geometry available; fall back to pin-based Manhattan.
    return scaleEdgeCost(man + jump_pen, timing_mul[src]);
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
    const auto deadline = scanOptDeadline(config);

    const bool have_begin = begin.has_value();
    const bool have_end = end.has_value();
    const bool end_fixed = have_end;
    const std::size_t scanopt_n
        = n + (have_begin ? 1 : 0) + (have_end ? 1 : 0);
    if (scanopt_n > kScanOptMaxMatrixCells) {
      logger->warn(
          utl::DFT,
          72,
          "ScanOpt ordering requested for {} cells, which exceeds the current "
          "matrix limit {}. Falling back to heuristic ordering.",
          scanopt_n,
          kScanOptMaxMatrixCells);
    } else {
      const std::size_t begin_node = have_begin ? n : 0;
      const std::size_t end_node = have_end ? (n + (have_begin ? 1 : 0)) : 0;
      const std::size_t start_node = have_begin ? begin_node : start_index;

      ScanOptMatrix m;
      m.n = scanopt_n;
      m.costs.resize(scanopt_n * scanopt_n);
      for (std::size_t i = 0; i < scanopt_n; ++i) {
        for (std::size_t j = 0; j < scanopt_n; ++j) {
          if (i == j) {
            m.costs[i * scanopt_n + j] = kScanOptLargeCost;
            continue;
          }
          int64_t cost = kInfDistance;
          if (have_begin && i == begin_node) {
            if (j < n) {
              const int64_t man
                  = manhattanDist(begin.value(), scan_in_pts[j], vertical_weight);
              const int64_t jump_pen = jumpPenaltyFromManhattan(man, local_scale);
              cost = man + jump_pen;
            } else {
              cost = kInfDistance;
            }
          } else if (have_begin && j == begin_node) {
            cost = kInfDistance;  // never enter begin node
          } else if (end_fixed && i == end_node) {
            cost = kInfDistance;  // don't leave the fixed end node
          } else if (end_fixed && j == end_node) {
            if (i < n) {
              const int64_t man
                  = manhattanDist(scan_out_pts[i], end.value(), vertical_weight);
              const int64_t jump_pen = jumpPenaltyFromManhattan(man, local_scale);
              cost = scaleEdgeCost(man + jump_pen, timing_mul[i]);
            } else {
              cost = kInfDistance;
            }
          } else {
            cost = edge_cost(i, j);
          }
          m.costs[i * scanopt_n + j]
              = static_cast<int32_t>(std::min<int64_t>(cost, kScanOptLargeCost));
        }
      }

      std::vector<std::string_view> scanopt_names = names;
      static constexpr std::string_view kBeginName = "__begin__";
      static constexpr std::string_view kEndName = "__end__";
      if (have_begin) {
        scanopt_names.push_back(kBeginName);
      }
      if (end_fixed) {
        scanopt_names.push_back(kEndName);
      }

      std::vector<std::size_t> best_order
          = end_fixed ? scanOptGreedyOrderEndFixed(m,
                                                   start_node,
                                                   end_node,
                                                   scanopt_names,
                                                   logger)
                      : scanOptGreedyOrder(m,
                                           start_node,
                                           scanopt_names,
                                           logger);
      scanOptDescent(m, scanopt_names, best_order, end_fixed, deadline);
      int64_t best_cost = scanOptPathCost(m, best_order);

      std::mt19937_64 rng(config.getScanOptSeed());
      const uint64_t rounds = config.getScanOptRounds();
      if (!config.getScanOptTempControl()) {
        for (uint64_t r = 0; r < rounds; ++r) {
          if (scanOptTimeExpired(deadline)) {
            break;
          }
          std::vector<std::size_t> cand = best_order;
          scanOptDoubleBridgeKick(cand, rng, end_fixed);
          scanOptDescent(m, scanopt_names, cand, end_fixed, deadline);
          const int64_t cand_cost = scanOptPathCost(m, cand);
          if (cand_cost < best_cost) {
            best_cost = cand_cost;
            best_order.swap(cand);
          }
        }
      } else {
        constexpr int kNoImproveBeforeTemp = 3;
        constexpr int kTempSteps = 3;

        std::vector<std::size_t> cur_order = best_order;
        int64_t cur_cost = best_cost;

        double temperature = 0.0;
        int no_improve = 0;
        int temp_steps_left = 0;

        const double t_div = config.getScanOptTDiv();
        std::uniform_real_distribution<double> u01(0.0, 1.0);

        for (uint64_t r = 0; r < rounds; ++r) {
          if (scanOptTimeExpired(deadline)) {
            break;
          }
          std::vector<std::size_t> cand = cur_order;
          scanOptDoubleBridgeKick(cand, rng, end_fixed);
          scanOptDescent(m, scanopt_names, cand, end_fixed, deadline);
          const int64_t cand_cost = scanOptPathCost(m, cand);

          const int64_t delta = cand_cost - cur_cost;
          if (delta < 0) {
            cur_order.swap(cand);
            cur_cost = cand_cost;
            temperature = 0.0;
            no_improve = 0;
            temp_steps_left = 0;
          } else {
            no_improve++;
            if (temperature == 0.0 && no_improve >= kNoImproveBeforeTemp) {
              temperature = static_cast<double>(cur_cost) / t_div;
              temp_steps_left = kTempSteps;
              no_improve = 0;
            }

            bool accept = false;
            if (temperature > 0.0) {
              const double prob
                  = std::exp(-static_cast<double>(delta) / temperature);
              accept = u01(rng) < prob;
              if (--temp_steps_left <= 0) {
                temperature = 0.0;
                no_improve = 0;
              }
            }
            if (accept) {
              cur_order.swap(cand);
              cur_cost = cand_cost;
            }
          }

          if (cur_cost < best_cost) {
            best_cost = cur_cost;
            best_order = cur_order;
          }
        }
      }

      std::vector<std::unique_ptr<ScanCell>> ordered;
      ordered.reserve(n);
      if (!begin.has_value() && !end.has_value()) {
        rotateOrderToDropWorstEdge(best_order, names, edge_cost);
      }
      for (const std::size_t idx : best_order) {
        if (idx >= n) {
          continue;
        }
        ordered.emplace_back(std::move(cells[idx]));
      }
      std::swap(cells, ordered);
      return;
    }
  }

  const bool end_fixed = end.has_value();
  const auto terminal_cost = [&](std::size_t idx) -> int64_t {
    if (!end_fixed) {
      return 0;
    }
    return manhattanDist(scan_out_pts[idx], end.value(), vertical_weight);
  };
  const auto no_next_scan_cell = [&]() -> void {
    logger->error(utl::DFT, 17, "Couldn't find next scan cell to order");
  };

  const auto path_cost = [&](const std::vector<std::size_t>& order) -> int64_t {
    int64_t total = 0;
    for (std::size_t i = 1; i < order.size(); ++i) {
      total += edge_cost(order[i - 1], order[i]);
    }
    if (end_fixed && !order.empty()) {
      total += terminal_cost(order.back());
    }
    return total;
  };

  const auto greedy_nn_order = [&]() -> std::vector<std::size_t> {
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
        if (end_fixed && last_step) {
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
        no_next_scan_cell();
      }

      used[best] = true;
      order.push_back(best);
      cur = best;
    }

    return order;
  };

  const auto farthest_insertion_order = [&]() -> std::vector<std::size_t> {
    std::vector<std::size_t> order;
    order.reserve(n);
    std::vector<bool> in_path(n, false);

    order.push_back(start_index);
    in_path[start_index] = true;

    if (n == 1) {
      return order;
    }

    std::size_t farthest = start_index;
    int64_t farthest_dist = -1;
    for (std::size_t i = 0; i < n; ++i) {
      if (i == start_index) {
        continue;
      }
      const int64_t dist = manhattanDist(scan_in_pts[start_index],
                                         scan_in_pts[i],
                                         vertical_weight);
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
      nearest_dist[i]
          = std::min(manhattanDist(scan_in_pts[i],
                                   scan_in_pts[start_index],
                                   vertical_weight),
                     manhattanDist(scan_in_pts[i],
                                   scan_in_pts[farthest],
                                   vertical_weight));
    }

    while (order.size() < n) {
      std::size_t next = start_index;
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

      // Insert at the position that minimises the path length increase.
      std::size_t best_pos = order.size();  // append by default
      int64_t best_delta = edge_cost(order.back(), next);
      if (end_fixed) {
        best_delta += terminal_cost(next) - terminal_cost(order.back());
      }

      for (std::size_t pos = 0; pos + 1 < order.size(); ++pos) {
        const auto a = order[pos];
        const auto b = order[pos + 1];
        const int64_t delta
            = edge_cost(a, next) + edge_cost(next, b) - edge_cost(a, b);
        if (delta < best_delta || (delta == best_delta && pos + 1 < best_pos)) {
          best_delta = delta;
          best_pos = pos + 1;
        }
      }

      order.insert(order.begin() + static_cast<std::ptrdiff_t>(best_pos), next);
      in_path[next] = true;

      for (std::size_t i = 0; i < n; ++i) {
        if (in_path[i]) {
          continue;
        }
        const int64_t dist
            = manhattanDist(scan_in_pts[i], scan_in_pts[next], vertical_weight);
        nearest_dist[i] = std::min(nearest_dist[i], dist);
      }
    }

    // Keep the start fixed as the first element.
    if (!order.empty() && order.front() != start_index) {
      const auto it = std::find(order.begin(), order.end(), start_index);
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

    const int64_t before = path_cost(order);
    bool improved_any = false;

    const auto two_opt_first_improve
        = [&](std::vector<std::size_t>& ord) -> bool {
      const std::size_t nn = ord.size();
      if (nn < 4) {
        return false;
      }

      std::vector<int64_t> forward_prefix(nn, 0);
      std::vector<int64_t> rev_prefix(nn, 0);
      for (std::size_t i = 1; i < nn; ++i) {
        forward_prefix[i]
            = forward_prefix[i - 1] + edge_cost(ord[i - 1], ord[i]);
        rev_prefix[i] = rev_prefix[i - 1] + edge_cost(ord[i], ord[i - 1]);
      }

      const auto segment_forward = [&](std::size_t from,
                                       std::size_t to) -> int64_t {
        return forward_prefix[to] - forward_prefix[from];
      };
      const auto segment_reversed = [&](std::size_t from,
                                        std::size_t to) -> int64_t {
        return rev_prefix[to] - rev_prefix[from];
      };

      for (std::size_t i = 0; i + 2 < nn; ++i) {
        for (std::size_t j = i + 2; j < nn; ++j) {
          const std::size_t a = ord[i];
          const std::size_t b = ord[i + 1];
          const std::size_t c = ord[j];
          const bool has_d = (j + 1 < nn);
          const std::size_t d = has_d ? ord[j + 1] : 0;

          const int64_t old_inside = segment_forward(i + 1, j);
          const int64_t new_inside = segment_reversed(i + 1, j);

          const int64_t old_edges
              = edge_cost(a, b) + (has_d ? edge_cost(c, d) : 0)
                + (!has_d && end_fixed ? terminal_cost(c) : 0);
          const int64_t new_edges
              = edge_cost(a, c) + (has_d ? edge_cost(b, d) : 0)
                + (!has_d && end_fixed ? terminal_cost(b) : 0);

          if (new_edges + new_inside < old_edges + old_inside) {
            std::reverse(ord.begin() + static_cast<std::ptrdiff_t>(i + 1),
                         ord.begin() + static_cast<std::ptrdiff_t>(j + 1));
            return true;
          }
        }
      }
      return false;
    };

    for (int pass = 0; pass < max_passes; ++pass) {
      if (!two_opt_first_improve(order)) {
        break;
      }
      improved_any = true;
    }

    if (improved_any) {
      const int64_t after = path_cost(order);
      debugPrint(logger,
                 utl::DFT,
                 "scan_chain_opt",
                 1,
                 "OptimizeScanWirelengthPinToNet: 2-opt improved path length "
                 "{} -> {} ({} cells)",
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

  std::vector<std::size_t> best_order;
  if (n <= kQuadraticHeuristicMaxCells) {
    best_order = greedy_nn_order();
    two_opt_improve(best_order, two_opt_passes);
    int64_t best_cost = path_cost(best_order);

    if (n <= kFarthestInsertionMaxCells) {
      std::vector<std::size_t> fi_order = farthest_insertion_order();
      two_opt_improve(fi_order, two_opt_passes);
      const int64_t fi_cost = path_cost(fi_order);
      if (fi_cost < best_cost) {
        best_cost = fi_cost;
        best_order = std::move(fi_order);
      }
    }
  } else {
    // Large-chain fallback: rtree greedy ordering with a Manhattan-aware choice
    // among the nearest Euclidean candidates.
    using Point = bg::model::point<int, 2, bg::cs::cartesian>;
    std::vector<std::pair<Point, std::size_t>> transformed;
    transformed.reserve(n);
    for (std::size_t i = 0; i < n; ++i) {
      const auto& p = scan_in_pts[i];
      transformed.emplace_back(Point(p.x(), p.y()), i);
    }
    bgi::rtree<std::pair<Point, std::size_t>, bgi::rstar<4>> rtree(transformed);
    auto cursor = transformed[start_index];
    rtree.remove(cursor);

    best_order.reserve(n);
    best_order.push_back(cursor.second);

    while (best_order.size() < n) {
      bool found = false;
      std::pair<Point, std::size_t> best = cursor;
      int64_t best_cost = kInfDistance;
      const std::size_t cursor_idx = cursor.second;
      const odb::Point cursor_out = scan_out_pts[cursor_idx];
      const Point query_pt(cursor_out.x(), cursor_out.y());
      const bool last_step = (best_order.size() + 1 == n);

      for (auto it
           = rtree.qbegin(bgi::nearest(query_pt, kNearestCandidateCount));
           it != rtree.qend();
           ++it) {
        const auto cand = *it;
        int64_t cost = edge_cost(cursor_idx, cand.second);
        if (end_fixed && last_step) {
          cost += terminal_cost(cand.second);
        }
        if (!found || cost < best_cost
            || (cost == best_cost && names[cand.second] < names[best.second])) {
          best = cand;
          best_cost = cost;
          found = true;
        }
      }

      if (!found) {
        no_next_scan_cell();
      }

      cursor = best;
      rtree.remove(cursor);
      best_order.push_back(cursor.second);
    }
  }

  std::vector<std::unique_ptr<ScanCell>> ordered;
  ordered.reserve(n);
  if (!begin.has_value() && !end.has_value()) {
    rotateOrderToDropWorstEdge(best_order, names, edge_cost);
  }
  for (const std::size_t idx : best_order) {
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
  // If UCLA ScanOpt is selected but we can't use it for the chosen metric,
  // map it to the in-tree ScanOpt solver.
  std::optional<ScanArchitectConfig> cfg_override;
  const ScanArchitectConfig* cfg = &config;

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
    if (config.getScanOrderSolver()
        == ScanArchitectConfig::ScanOrderSolver::UclaScanOpt) {
      cfg_override = config;
      cfg_override->setScanOrderSolver(
          ScanArchitectConfig::ScanOrderSolver::ScanOpt);
      cfg = &cfg_override.value();
      logger->warn(
          utl::DFT,
          223,
          "UCLA_SCANOPT is only supported for placement metric; falling back "
          "to SCANOPT for PIN_TO_NET.");
    }
    OptimizeScanWirelengthPinToNet(cells, *cfg, logger, endpoints);
    return;
  }

  if (config.getScanOrderSolver()
      == ScanArchitectConfig::ScanOrderSolver::UclaScanOpt) {
    cfg_override = config;
    cfg_override->setScanOrderSolver(ScanArchitectConfig::ScanOrderSolver::ScanOpt);
    cfg = &cfg_override.value();
  }

  const double vertical_weight = cfg->getVerticalWeight();
  const auto manhattan = [&](const odb::Point& a, const odb::Point& b) -> int64_t {
    return manhattanDist(a, b, vertical_weight);
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
    timing_mul.emplace_back(timingMultiplierForCell(*cfg, *cell));
  }
  const int64_t local_scale
      = estimateLocalManhattanScale(scan_in_pts, vertical_weight);

  std::optional<odb::Point> begin;
  std::optional<odb::Point> end;
  odb::dbBlock* endpoint_block = nullptr;
  if (endpoints.has_value()) {
    if ((endpoints->begin.has_value()
         && endpoints->begin->type
                == ScanArchitectConfig::ChainEndpoint::Type::Term)
        || (endpoints->end.has_value()
            && endpoints->end->type
                   == ScanArchitectConfig::ChainEndpoint::Type::Term)) {
      endpoint_block = InferBlockFromPlacedCells(cells);
    }
  }
  if (endpoints.has_value()) {
    if (endpoints->begin.has_value()) {
      begin = EndpointPoint(*endpoints->begin, endpoint_block, logger);
    }
    if (endpoints->end.has_value()) {
      end = EndpointPoint(*endpoints->end, endpoint_block, logger);
    }
  }

  const auto edge_cost = [&](std::size_t src, std::size_t dst) -> int64_t {
    const int64_t man = manhattan(scan_out_pts[src], scan_in_pts[dst]);
    const int64_t jump_pen = jumpPenaltyFromManhattan(man, local_scale);
    return scaleEdgeCost(man + jump_pen, timing_mul[src]);
  };

  // UCLA ScanOpt only supports unconstrained placement ordering with fixed
  // begin/end points.
  if (config.getScanOrderSolver()
      == ScanArchitectConfig::ScanOrderSolver::UclaScanOpt) {
    if (hasScanOrderConstraints(config)) {
      logger->warn(
          utl::DFT,
          187,
          "UCLA_SCANOPT does not support scan order constraints; falling back "
          "to SCANOPT.");
    } else if (begin.has_value() && end.has_value()) {
      UclaScanOptParams p;
      p.seed = cfg->getScanOptSeed();
      // Map OpenROAD ScanOpt rounds (default 500k) to UCLA majorLoops (default 100).
      p.major_loops = std::max<uint64_t>(1, cfg->getScanOptRounds() / 5000);
      p.n_descents = 5;
      p.kick_move = 15;
      p.n_near = 20;
      p.only_2opt = false;
      p.temp_control = cfg->getScanOptTempControl();
      p.zero_temp = !cfg->getScanOptTempControl();

      std::vector<std::pair<int, int>> in_pts;
      std::vector<std::pair<int, int>> out_pts;
      in_pts.reserve(n);
      out_pts.reserve(n);
      for (std::size_t i = 0; i < n; ++i) {
        in_pts.emplace_back(scan_in_pts[i].x(), scan_in_pts[i].y());
        out_pts.emplace_back(scan_out_pts[i].x(), scan_out_pts[i].y());
      }

      std::vector<std::size_t> order;
      try {
        order = UclaScanOptOrder(
            names,
            in_pts,
            out_pts,
            std::make_pair(begin->x(), begin->y()),
            std::make_pair(end->x(), end->y()),
            p);
      } catch (const std::exception& e) {
        logger->warn(utl::DFT,
                     188,
                     "UCLA_SCANOPT failed ({}); falling back to SCANOPT.",
                     e.what());
        order.clear();
      }

      if (order.size() == n) {
        std::vector<std::unique_ptr<ScanCell>> ordered;
        ordered.reserve(n);
        for (const std::size_t i : order) {
          ordered.emplace_back(std::move(cells[i]));
        }
        cells.swap(ordered);
        return;
      }
    } else {
      logger->warn(
          utl::DFT,
          189,
          "UCLA_SCANOPT requires fixed begin/end points; falling back to "
          "SCANOPT.");
    }
  }

  if (hasScanOrderConstraints(*cfg)) {
    optimizeScanWirelengthWithConstraints(
        cells,
        *cfg,
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
    const auto deadline = scanOptDeadline(config);

    const bool have_begin = begin.has_value();
    const bool have_end = end.has_value();
    const bool end_fixed = have_end;
    const std::size_t scanopt_n
        = n + (have_begin ? 1 : 0) + (have_end ? 1 : 0);
    if (scanopt_n > kScanOptMaxMatrixCells) {
      logger->warn(
          utl::DFT,
          185,
          "ScanOpt ordering requested for {} cells, which exceeds the current "
          "matrix limit {}. Falling back to heuristic ordering.",
          scanopt_n,
          kScanOptMaxMatrixCells);
    } else {
      const std::size_t begin_node = have_begin ? n : 0;
      const std::size_t end_node = have_end ? (n + (have_begin ? 1 : 0)) : 0;
      const std::size_t start_node = have_begin ? begin_node : start_index;

      ScanOptMatrix m;
      m.n = scanopt_n;
      m.costs.resize(scanopt_n * scanopt_n);
      for (std::size_t i = 0; i < scanopt_n; ++i) {
        for (std::size_t j = 0; j < scanopt_n; ++j) {
          if (i == j) {
            m.costs[i * scanopt_n + j] = kScanOptLargeCost;
            continue;
          }
          int64_t cost = kInfDistance;
          if (have_begin && i == begin_node) {
            if (j < n) {
              const int64_t man
                  = manhattanDist(begin.value(), scan_in_pts[j], vertical_weight);
              const int64_t jump_pen = jumpPenaltyFromManhattan(man, local_scale);
              cost = man + jump_pen;
            } else {
              cost = kInfDistance;
            }
          } else if (have_begin && j == begin_node) {
            cost = kInfDistance;  // never enter begin node
          } else if (end_fixed && i == end_node) {
            cost = kInfDistance;  // don't leave the fixed end node
          } else if (end_fixed && j == end_node) {
            if (i < n) {
              const int64_t man
                  = manhattanDist(scan_out_pts[i], end.value(), vertical_weight);
              const int64_t jump_pen = jumpPenaltyFromManhattan(man, local_scale);
              cost = scaleEdgeCost(man + jump_pen, timing_mul[i]);
            } else {
              cost = kInfDistance;
            }
          } else {
            cost = edge_cost(i, j);
          }
          m.costs[i * scanopt_n + j]
              = static_cast<int32_t>(std::min<int64_t>(cost, kScanOptLargeCost));
        }
      }

      std::vector<std::string_view> scanopt_names = names;
      static constexpr std::string_view kBeginName = "__begin__";
      static constexpr std::string_view kEndName = "__end__";
      if (have_begin) {
        scanopt_names.push_back(kBeginName);
      }
      if (end_fixed) {
        scanopt_names.push_back(kEndName);
      }

      std::vector<std::size_t> best_order
          = end_fixed ? scanOptGreedyOrderEndFixed(m,
                                                   start_node,
                                                   end_node,
                                                   scanopt_names,
                                                   logger)
                      : scanOptGreedyOrder(m,
                                           start_node,
                                           scanopt_names,
                                           logger);
      scanOptDescent(m, scanopt_names, best_order, end_fixed, deadline);
      int64_t best_cost = scanOptPathCost(m, best_order);

      std::mt19937_64 rng(config.getScanOptSeed());
      const uint64_t rounds = config.getScanOptRounds();
      if (!config.getScanOptTempControl()) {
        for (uint64_t r = 0; r < rounds; ++r) {
          if (scanOptTimeExpired(deadline)) {
            break;
          }
          std::vector<std::size_t> cand = best_order;
          scanOptDoubleBridgeKick(cand, rng, end_fixed);
          scanOptDescent(m, scanopt_names, cand, end_fixed, deadline);
          const int64_t cand_cost = scanOptPathCost(m, cand);
          if (cand_cost < best_cost) {
            best_cost = cand_cost;
            best_order.swap(cand);
          }
        }
      } else {
        constexpr int kNoImproveBeforeTemp = 3;
        constexpr int kTempSteps = 3;

        std::vector<std::size_t> cur_order = best_order;
        int64_t cur_cost = best_cost;

        double temperature = 0.0;
        int no_improve = 0;
        int temp_steps_left = 0;

        const double t_div = config.getScanOptTDiv();
        std::uniform_real_distribution<double> u01(0.0, 1.0);

        for (uint64_t r = 0; r < rounds; ++r) {
          if (scanOptTimeExpired(deadline)) {
            break;
          }
          std::vector<std::size_t> cand = cur_order;
          scanOptDoubleBridgeKick(cand, rng, end_fixed);
          scanOptDescent(m, scanopt_names, cand, end_fixed, deadline);
          const int64_t cand_cost = scanOptPathCost(m, cand);

          const int64_t delta = cand_cost - cur_cost;
          if (delta < 0) {
            cur_order.swap(cand);
            cur_cost = cand_cost;
            temperature = 0.0;
            no_improve = 0;
            temp_steps_left = 0;
          } else {
            no_improve++;
            if (temperature == 0.0 && no_improve >= kNoImproveBeforeTemp) {
              temperature = static_cast<double>(cur_cost) / t_div;
              temp_steps_left = kTempSteps;
              no_improve = 0;
            }

            bool accept = false;
            if (temperature > 0.0) {
              const double prob
                  = std::exp(-static_cast<double>(delta) / temperature);
              accept = u01(rng) < prob;
              if (--temp_steps_left <= 0) {
                temperature = 0.0;
                no_improve = 0;
              }
            }
            if (accept) {
              cur_order.swap(cand);
              cur_cost = cand_cost;
            }
          }

          if (cur_cost < best_cost) {
            best_cost = cur_cost;
            best_order = cur_order;
          }
        }
      }

      std::vector<std::unique_ptr<ScanCell>> ordered;
      ordered.reserve(n);
      if (!begin.has_value() && !end.has_value()) {
        rotateOrderToDropWorstEdge(best_order, names, edge_cost);
      }
      for (const std::size_t idx : best_order) {
        if (idx >= n) {
          continue;
        }
        ordered.emplace_back(std::move(cells[idx]));
      }
      std::swap(cells, ordered);
      return;
    }
  }

  // Timing-aware ordering makes the objective asymmetric (penalizes outgoing
  // edges from timing-critical sources). EndPort costs also make the objective
  // directional. Use a directed greedy heuristic with a small local-improvement
  // pass rather than symmetric 2-opt/farthest-insertion moves.
  if (config.getTimingWeightSetup() != 0.0 || config.getTimingWeightHold() != 0.0
      || end.has_value()) {
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
    if (!begin.has_value() && !end.has_value()) {
      rotateOrderToDropWorstEdge(order, names, edge_cost);
    }
    for (const std::size_t idx : order) {
      ordered.emplace_back(std::move(cells[idx]));
    }
    std::swap(cells, ordered);
	    return;
	  }

	  const bool end_fixed = end.has_value();
	  const auto terminal_cost = [&](std::size_t idx) -> int64_t {
	    if (!end_fixed) {
	      return 0;
	    }
	    return manhattan(scan_out_pts[idx], end.value());
	  };

	  const auto path_cost = [&](const std::vector<std::size_t>& order) -> int64_t {
	    if (order.empty()) {
	      return 0;
	    }
	    int64_t total = 0;
	    if (begin.has_value()) {
	      total += manhattan(begin.value(), scan_in_pts[order.front()]);
	    }
	    for (std::size_t i = 1; i < order.size(); ++i) {
	      total += edge_cost(order[i - 1], order[i]);
	    }
	    if (end_fixed) {
	      total += terminal_cost(order.back());
	    }
	    return total;
	  };

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
		      int64_t best_cost = std::numeric_limits<int64_t>::max();
		      const bool last_step = (order.size() + 1 == n);
		      for (std::size_t i = 0; i < n; ++i) {
		        if (used[i]) {
		          continue;
		        }
		        int64_t cost = edge_cost(cur, i);
		        if (end_fixed && last_step) {
		          cost += terminal_cost(i);
		        }
		        if (!found || cost < best_cost
		            || (cost == best_cost && names[i] < names[best])) {
		          best = i;
		          best_cost = cost;
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
      int64_t best_delta = edge_cost(order.back(), next);
      if (end_fixed) {
        best_delta += terminal_cost(next) - terminal_cost(order.back());
      }

      for (std::size_t pos = 0; pos + 1 < order.size(); ++pos) {
        const auto a = order[pos];
        const auto b = order[pos + 1];
        const int64_t delta
            = edge_cost(a, next) + edge_cost(next, b) - edge_cost(a, b);
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

		    const int64_t before = path_cost(order);
		    bool improved_any = false;

	    const auto two_opt_first_improve
	        = [&](std::vector<std::size_t>& ord) -> bool {
	      const std::size_t nn = ord.size();
	      if (nn < 4) {
	        return false;
	      }
	      // Prefix sums for forward edges and reversed-adjacent edges.
	      std::vector<int64_t> forward_prefix(nn, 0);
	      std::vector<int64_t> rev_prefix(nn, 0);
	      for (std::size_t i = 1; i < nn; ++i) {
	        forward_prefix[i]
	            = forward_prefix[i - 1] + edge_cost(ord[i - 1], ord[i]);
	        rev_prefix[i]
	            = rev_prefix[i - 1] + edge_cost(ord[i], ord[i - 1]);
	      }

	      const auto segment_forward = [&](std::size_t from,
	                                       std::size_t to) -> int64_t {
	        return forward_prefix[to] - forward_prefix[from];
	      };
	      const auto segment_reversed = [&](std::size_t from,
	                                        std::size_t to) -> int64_t {
	        return rev_prefix[to] - rev_prefix[from];
	      };

	      for (std::size_t i = 0; i + 2 < nn; ++i) {
	        for (std::size_t j = i + 2; j < nn; ++j) {
	          const std::size_t a = ord[i];
	          const std::size_t b = ord[i + 1];
		          const std::size_t c = ord[j];
		          const bool has_d = (j + 1 < nn);
		          const std::size_t d = has_d ? ord[j + 1] : 0;

	          const int64_t old_inside = segment_forward(i + 1, j);
	          const int64_t new_inside = segment_reversed(i + 1, j);

		          const int64_t old_edges
		              = edge_cost(a, b) + (has_d ? edge_cost(c, d) : 0)
		                + (!has_d && end_fixed ? terminal_cost(c) : 0);
		          const int64_t new_edges
		              = edge_cost(a, c) + (has_d ? edge_cost(b, d) : 0)
		                + (!has_d && end_fixed ? terminal_cost(b) : 0);

		          if (new_edges + new_inside < old_edges + old_inside) {
		            std::reverse(ord.begin() + static_cast<std::ptrdiff_t>(i + 1),
		                         ord.begin() + static_cast<std::ptrdiff_t>(j + 1));
	            return true;
	          }
	        }
	      }
	      return false;
	    };

	    for (int pass = 0; pass < max_passes; ++pass) {
	      if (!two_opt_first_improve(order)) {
	        break;
	      }
	      improved_any = true;
	    }

		    if (improved_any) {
		      const int64_t after = path_cost(order);
		      debugPrint(logger,
		                 utl::DFT,
	                 "scan_chain_opt",
	                 1,
	                 "OptimizeScanWirelength: 2-opt improved path cost {} -> {} "
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
	    int64_t best_cost = path_cost(best_order);

	    if (n <= kFarthestInsertionMaxCells) {
	      std::vector<std::size_t> fi_order = farthest_insertion_order(start_index);
	      two_opt_improve(fi_order, two_opt_passes);
	      const int64_t fi_cost = path_cost(fi_order);
	      if (fi_cost < best_cost) {
	        best_order = std::move(fi_order);
	      }
	    }

    std::vector<std::unique_ptr<ScanCell>> ordered;
    ordered.reserve(n);
    if (!begin.has_value() && !end.has_value()) {
      rotateOrderToDropWorstEdge(best_order, names, edge_cost);
    }
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
		    const std::size_t cursor_idx = cursor.second;
		    const odb::Point cursor_out = scan_out_pts[cursor_idx];
		    const Point query_pt(cursor_out.x(), cursor_out.y());
		    const bool last_step = (ordered.size() + 1 == cells.size());

		    for (auto it
		         = rtree.qbegin(bgi::nearest(query_pt, kNearestCandidateCount));
		         it != rtree.qend();
		         ++it) {
		      const auto cand = *it;
		      int64_t dist = edge_cost(cursor_idx, cand.second);
		      if (end_fixed && last_step) {
		        dist += terminal_cost(cand.second);
		      }
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
