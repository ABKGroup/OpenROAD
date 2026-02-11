// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2023-2025, The OpenROAD Authors

#include "ScanArchitectHeuristic.hh"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <iterator>
#include <limits>
#include <memory>
#include <numeric>
#include <optional>
#include <string>
#include <string_view>
#include <type_traits>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "Opt.hh"
#include "ScanArchitect.hh"
#include "ScanArchitectConfig.hh"
#include "ScanCell.hh"
#include "odb/geom.h"
#include "spdlog/fmt/fmt.h"
#include "utl/Logger.h"

namespace dft {

namespace {

struct PlacedScanCell
{
  std::unique_ptr<ScanCell> cell;
  odb::Point origin;
  uint64_t bits = 0;
  std::string_view name;
};

struct PlacedBundle
{
  std::vector<std::unique_ptr<ScanCell>> cells;
  odb::Point origin;
  uint64_t bits = 0;
  std::string name;
  bool all_placed = true;
  std::optional<std::size_t> fixed_chain;
};

int64_t manhattanDist(const odb::Point& a, const odb::Point& b)
{
  const int64_t dx = static_cast<int64_t>(a.x()) - static_cast<int64_t>(b.x());
  const int64_t dy = static_cast<int64_t>(a.y()) - static_cast<int64_t>(b.y());
  return std::abs(dx) + std::abs(dy);
}

struct UnionFind
{
  std::vector<std::size_t> parent;
  std::vector<std::size_t> rank;

  explicit UnionFind(std::size_t n) : parent(n), rank(n, 0)
  {
    std::iota(parent.begin(), parent.end(), 0);
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
      rank[a]++;
    }
  }
};

struct ChainRef
{
  std::size_t hash_domain = 0;
  std::size_t local_index = 0;
};

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

odb::dbBlock* InferBlockFromBundles(const std::vector<PlacedBundle>& bundles)
{
  for (const PlacedBundle& bundle : bundles) {
    if (bundle.cells.empty() || bundle.cells.front() == nullptr) {
      continue;
    }
    const ScanLoad scan_in = bundle.cells.front()->getScanIn();
    odb::dbBlock* block = nullptr;
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
    if (block != nullptr) {
      return block;
    }
  }
  return nullptr;
}

std::optional<odb::Point> ResolveEndpointTerm(odb::dbBlock* block,
                                              std::string_view term)
{
  if (block == nullptr) {
    return std::nullopt;
  }

  const auto term_info = SplitTermIdentifier(term);
  if (term_info.second.has_value()) {
    odb::dbInst* inst = block->findInst(term_info.first.c_str());
    if (inst == nullptr) {
      return std::nullopt;
    }
    odb::dbITerm* iterm = inst->findITerm(term_info.second->c_str());
    if (iterm == nullptr) {
      return std::nullopt;
    }
    const odb::Rect bbox = iterm->getBBox();
    return odb::Point(bbox.xMin(), bbox.yMin());
  }

  odb::dbBTerm* bterm = block->findBTerm(term_info.first.c_str());
  if (bterm == nullptr) {
    return std::nullopt;
  }
  const odb::Rect bbox = bterm->getBBox();
  return odb::Point(bbox.xMin(), bbox.yMin());
}

struct ChainEndpointPoints
{
  std::optional<odb::Point> begin;
  std::optional<odb::Point> end;
};

enum class AnchorAxis
{
  None,
  X,
  Y
};

struct AnchorOrder
{
  AnchorAxis axis{AnchorAxis::None};
  std::vector<std::size_t> order;
  std::vector<std::size_t> pos;
};

AnchorOrder ComputeAnchorOrderForAxis(
    const std::vector<std::optional<odb::Point>>& chain_anchor,
    AnchorAxis axis)
{
  AnchorOrder out;
  const std::size_t chain_count = chain_anchor.size();
  if (chain_count < 2) {
    return out;
  }

  std::size_t have = 0;
  bool first = true;
  int min_x = 0;
  int max_x = 0;
  int min_y = 0;
  int max_y = 0;
  for (const auto& a : chain_anchor) {
    if (!a.has_value()) {
      continue;
    }
    const odb::Point& p = a.value();
    if (first) {
      min_x = max_x = p.x();
      min_y = max_y = p.y();
      first = false;
    } else {
      min_x = std::min(min_x, p.x());
      max_x = std::max(max_x, p.x());
      min_y = std::min(min_y, p.y());
      max_y = std::max(max_y, p.y());
    }
    ++have;
  }

  if (have < 2) {
    return out;
  }

  const int range = (axis == AnchorAxis::X) ? (max_x - min_x) : (max_y - min_y);
  if (range == 0) {
    return out;
  }

  out.axis = axis;
  out.order.resize(chain_count);
  std::iota(out.order.begin(), out.order.end(), 0);
  std::stable_sort(out.order.begin(),
                   out.order.end(),
                   [&](std::size_t a, std::size_t b) {
                     const auto& pa = chain_anchor[a];
                     const auto& pb = chain_anchor[b];
                     if (pa.has_value() && pb.has_value()) {
                       const int va
                           = (axis == AnchorAxis::Y) ? pa->y() : pa->x();
                       const int vb
                           = (axis == AnchorAxis::Y) ? pb->y() : pb->x();
                       if (va != vb) {
                         return va < vb;
                       }
                       const int sa
                           = (axis == AnchorAxis::Y) ? pa->x() : pa->y();
                       const int sb
                           = (axis == AnchorAxis::Y) ? pb->x() : pb->y();
                       if (sa != sb) {
                         return sa < sb;
                       }
                     }
                     return a < b;
                   });

  out.pos.assign(chain_count, 0);
  for (std::size_t i = 0; i < out.order.size(); ++i) {
    out.pos[out.order[i]] = i;
  }
  return out;
}

AnchorOrder ComputeAnchorOrder(
    const std::vector<std::optional<odb::Point>>& chain_anchor)
{
  AnchorOrder out;
  const std::size_t chain_count = chain_anchor.size();
  if (chain_count < 2) {
    return out;
  }

  std::size_t have = 0;
  bool first = true;
  int min_x = 0;
  int max_x = 0;
  int min_y = 0;
  int max_y = 0;
  for (const auto& a : chain_anchor) {
    if (!a.has_value()) {
      continue;
    }
    const odb::Point& p = a.value();
    if (first) {
      min_x = max_x = p.x();
      min_y = max_y = p.y();
      first = false;
    } else {
      min_x = std::min(min_x, p.x());
      max_x = std::max(max_x, p.x());
      min_y = std::min(min_y, p.y());
      max_y = std::max(max_y, p.y());
    }
    ++have;
  }

  if (have < 2) {
    return out;
  }

  const int range_x = max_x - min_x;
  const int range_y = max_y - min_y;
  if (range_x == 0 && range_y == 0) {
    return out;
  }

  // Prefer Y on ties since standard-cell placements are typically row-based
  // and scan ordering is particularly sensitive to large vertical "jumps".
  const AnchorAxis axis = (range_y >= range_x) ? AnchorAxis::Y : AnchorAxis::X;
  return ComputeAnchorOrderForAxis(chain_anchor, axis);
}

std::vector<ChainEndpointPoints> ResolveChainEndpointsForDomain(
    const std::vector<PlacedBundle>& bundles,
    std::size_t hash_domain,
    const ScanArchitectConfig& config,
    const std::unordered_map<std::string_view, ChainRef>& chain_name_to_ref,
    const std::unordered_map<std::string_view, std::size_t>&
        chain_name_to_ordinal,
    std::size_t chain_count)
{
  std::vector<std::string_view> chain_names(chain_count);
  for (const auto& [name, ref] : chain_name_to_ref) {
    if (ref.hash_domain != hash_domain) {
      continue;
    }
    if (ref.local_index < chain_count) {
      chain_names[ref.local_index] = name;
    }
  }

  odb::dbBlock* block = InferBlockFromBundles(bundles);
  std::vector<ChainEndpointPoints> out(chain_count);
  for (std::size_t c = 0; c < chain_count; ++c) {
    if (chain_names[c].empty()) {
      continue;
    }

    std::optional<ScanArchitectConfig::ChainEndpoints> endpoints
        = config.getChainEndpoints(chain_names[c]);
    // If no explicit endpoints were provided for this chain, treat the scan
    // in/out ports (from name patterns) as implicit endpoints so clustering can
    // avoid large begin/end IO "jumps" by default.
    const auto ord_it = chain_name_to_ordinal.find(chain_names[c]);
    if (ord_it != chain_name_to_ordinal.end()) {
      ScanArchitectConfig::ChainEndpoints implicit;
      if (endpoints.has_value()) {
        implicit = endpoints.value();
      }
      bool changed = false;
      if (!implicit.begin.has_value()
          && !config.getScanInNamePattern().empty()) {
        ScanArchitectConfig::ChainEndpoint begin;
        begin.type = ScanArchitectConfig::ChainEndpoint::Type::Term;
        begin.term = fmt::format(FMT_RUNTIME(config.getScanInNamePattern()),
                                 ord_it->second);
        implicit.begin = std::move(begin);
        changed = true;
      }
      if (!implicit.end.has_value()
          && !config.getScanOutNamePattern().empty()) {
        ScanArchitectConfig::ChainEndpoint end;
        end.type = ScanArchitectConfig::ChainEndpoint::Type::Term;
        end.term = fmt::format(FMT_RUNTIME(config.getScanOutNamePattern()),
                               ord_it->second);
        implicit.end = std::move(end);
        changed = true;
      }
      if (changed) {
        endpoints = std::move(implicit);
      }
    }

    if (!endpoints.has_value()) {
      continue;
    }

    if (endpoints->begin.has_value()) {
      const auto& ep = endpoints->begin.value();
      if (ep.type == ScanArchitectConfig::ChainEndpoint::Type::Point) {
        out[c].begin = odb::Point(ep.point.x, ep.point.y);
      } else {
        out[c].begin = ResolveEndpointTerm(block, ep.term);
      }
    }
    if (endpoints->end.has_value()) {
      const auto& ep = endpoints->end.value();
      if (ep.type == ScanArchitectConfig::ChainEndpoint::Type::Point) {
        out[c].end = odb::Point(ep.point.x, ep.point.y);
      } else {
        out[c].end = ResolveEndpointTerm(block, ep.term);
      }
    }
  }
  return out;
}

void enforceMaxImbalanceOrDie(const std::vector<PlacedBundle>& bundles,
                              std::vector<std::size_t>& assignment,
                              std::size_t chain_count,
                              uint64_t max_length,
                              bool hard_max_length,
                              const ScanArchitectConfig& config,
                              utl::Logger* logger,
                              std::size_t hash_domain,
                              const std::vector<std::optional<odb::Point>>*
                                  chain_anchor)
{
  if (chain_count <= 1 || bundles.empty()) {
    return;
  }

  const double allowed_ratio
      = 1.0 + (std::max(0.0, config.getMaxImbalancePercent()) / 100.0);

  std::optional<AnchorOrder> anchor_order;
  if (chain_anchor != nullptr) {
    AnchorOrder o = ComputeAnchorOrder(*chain_anchor);
    if (o.axis != AnchorAxis::None && o.order.size() == chain_count
        && o.pos.size() == chain_count) {
      anchor_order = std::move(o);
    }
  }

  auto compute_chain_stats = [&]() {
    std::vector<uint64_t> used(chain_count, 0);
    std::vector<int64_t> sum_x(chain_count, 0);
    std::vector<int64_t> sum_y(chain_count, 0);
    std::vector<uint64_t> sum_bits(chain_count, 0);
    for (std::size_t i = 0; i < bundles.size(); ++i) {
      const std::size_t c = assignment[i];
      if (c >= chain_count) {
        if (logger) {
          logger->error(utl::DFT,
                        123,
                        "Scan constraints internal error: bundle assignment "
                        "index {} is out of range for chain_count {} (hash "
                        "domain {}).",
                        c,
                        chain_count,
                        hash_domain);
        }
      }
      used[c] += bundles[i].bits;
      sum_x[c] += static_cast<int64_t>(bundles[i].origin.x())
                  * static_cast<int64_t>(bundles[i].bits);
      sum_y[c] += static_cast<int64_t>(bundles[i].origin.y())
                  * static_cast<int64_t>(bundles[i].bits);
      sum_bits[c] += bundles[i].bits;
    }
    return std::make_tuple(std::move(used),
                           std::move(sum_x),
                           std::move(sum_y),
                           std::move(sum_bits));
  };

  auto get_min_max = [&](const std::vector<uint64_t>& used,
                         std::size_t& min_c,
                         std::size_t& max_c,
                         uint64_t& min_bits,
                         uint64_t& max_bits) -> void {
    min_c = 0;
    max_c = 0;
    min_bits = std::numeric_limits<uint64_t>::max();
    max_bits = 0;
    for (std::size_t c = 0; c < used.size(); ++c) {
      const uint64_t v = used[c];
      if (v < min_bits) {
        min_bits = v;
        min_c = c;
      }
      if (v > max_bits) {
        max_bits = v;
        max_c = c;
      }
    }
  };

  auto violates = [&](const std::vector<uint64_t>& used) -> bool {
    uint64_t min_bits = std::numeric_limits<uint64_t>::max();
    uint64_t max_bits = 0;
    for (const uint64_t v : used) {
      min_bits = std::min(min_bits, v);
      max_bits = std::max(max_bits, v);
    }
    if (max_bits == 0) {
      return false;  // empty domain
    }
    if (min_bits == 0) {
      return true;
    }
    const double ratio = static_cast<double>(max_bits) / static_cast<double>(min_bits);
    return ratio > allowed_ratio + 1e-12;
  };

  auto [used_bits, sum_x, sum_y, sum_bits] = compute_chain_stats();
  if (!violates(used_bits)) {
    return;
  }

  // Build per-chain bundle lists for move operations.
  std::vector<std::vector<std::size_t>> chain_to_bundles(chain_count);
  chain_to_bundles.reserve(chain_count);
  for (std::size_t i = 0; i < bundles.size(); ++i) {
    chain_to_bundles[assignment[i]].push_back(i);
  }

  const auto is_movable = [&](std::size_t bundle_idx) -> bool {
    return !bundles[bundle_idx].fixed_chain.has_value();
  };

  const auto centroid = [&](std::size_t c) -> odb::Point {
    if (sum_bits[c] == 0) {
      return odb::Point(0, 0);
    }
    const int64_t cx = sum_x[c] / static_cast<int64_t>(sum_bits[c]);
    const int64_t cy = sum_y[c] / static_cast<int64_t>(sum_bits[c]);
    return odb::Point(static_cast<int>(cx), static_cast<int>(cy));
  };

  const auto chain_point = [&](std::size_t c) -> odb::Point {
    if (chain_anchor != nullptr && c < chain_anchor->size()
        && (*chain_anchor)[c].has_value()) {
      return (*chain_anchor)[c].value();
    }
    return centroid(c);
  };

  const auto move_bundle = [&](std::size_t bundle_idx,
                               std::size_t from_chain,
                               std::size_t to_chain) -> void {
    assignment[bundle_idx] = to_chain;
    used_bits[from_chain] -= bundles[bundle_idx].bits;
    used_bits[to_chain] += bundles[bundle_idx].bits;
    sum_x[from_chain] -= static_cast<int64_t>(bundles[bundle_idx].origin.x())
                         * static_cast<int64_t>(bundles[bundle_idx].bits);
    sum_y[from_chain] -= static_cast<int64_t>(bundles[bundle_idx].origin.y())
                         * static_cast<int64_t>(bundles[bundle_idx].bits);
    sum_bits[from_chain] -= bundles[bundle_idx].bits;
    sum_x[to_chain] += static_cast<int64_t>(bundles[bundle_idx].origin.x())
                       * static_cast<int64_t>(bundles[bundle_idx].bits);
    sum_y[to_chain] += static_cast<int64_t>(bundles[bundle_idx].origin.y())
                       * static_cast<int64_t>(bundles[bundle_idx].bits);
    sum_bits[to_chain] += bundles[bundle_idx].bits;

    auto& from = chain_to_bundles[from_chain];
    auto it = std::find(from.begin(), from.end(), bundle_idx);
    if (it != from.end()) {
      from.erase(it);
    }
    chain_to_bundles[to_chain].push_back(bundle_idx);
  };

  // First, ensure there are no empty chains when there are bits to pack.
  for (;;) {
    std::size_t empty_chain = chain_count;
    for (std::size_t c = 0; c < chain_count; ++c) {
      if (used_bits[c] == 0) {
        empty_chain = c;
        break;
      }
    }
    if (empty_chain == chain_count) {
      break;
    }

    std::size_t max_c = 0;
    uint64_t max_bits = 0;
    for (std::size_t c = 0; c < chain_count; ++c) {
      if (used_bits[c] > max_bits) {
        max_bits = used_bits[c];
        max_c = c;
      }
    }

    bool moved = false;
    std::size_t best_bundle = 0;
    uint64_t best_bits = 0;
    for (const std::size_t bi : chain_to_bundles[max_c]) {
      if (!is_movable(bi)) {
        continue;
      }
      if (chain_to_bundles[max_c].size() <= 1) {
        continue;  // don't empty the donor
      }
      const uint64_t bbits = bundles[bi].bits;
      if (hard_max_length && bbits > max_length) {
        continue;
      }
      if (bbits > best_bits) {
        best_bits = bbits;
        best_bundle = bi;
        moved = true;
      }
    }

    if (!moved) {
      if (logger) {
        logger->error(
            utl::DFT,
            124,
            "Scan architect constraints infeasible in hash domain {}: cannot "
            "avoid empty scan chains while satisfying constraints.",
            hash_domain);
      }
    }
    move_bundle(best_bundle, max_c, empty_chain);
  }

  // Rebalance by moving a bundle from the longest to the shortest chain.
  constexpr std::size_t kMaxMoves = 100000;
  std::size_t moves = 0;
  while (violates(used_bits)) {
    if (moves++ > kMaxMoves) {
      if (logger) {
        logger->error(utl::DFT,
                      125,
                      "Scan architect internal error: max_imbalance rebalancing "
                      "did not converge (hash domain {}).",
                      hash_domain);
      }
    }

    std::size_t min_c = 0;
    std::size_t max_c = 0;
    uint64_t min_bits = 0;
    uint64_t max_bits = 0;
    get_min_max(used_bits, min_c, max_c, min_bits, max_bits);
    if (min_bits == 0) {
      // Should have been handled above, but keep a safety check.
      if (logger) {
        logger->error(
            utl::DFT,
            126,
            "Scan architect constraints infeasible in hash domain {}: cannot "
            "satisfy max_imbalance due to empty scan chains.",
            hash_domain);
      }
    }

    const std::size_t receiver = [&]() -> std::size_t {
      if (!anchor_order.has_value()) {
        return min_c;
      }
      const std::size_t pos_max = anchor_order->pos[max_c];
      const std::size_t pos_min = anchor_order->pos[min_c];
      if (pos_max < pos_min && pos_max + 1 < anchor_order->order.size()) {
        return anchor_order->order[pos_max + 1];
      }
      if (pos_max > pos_min && pos_max > 0) {
        return anchor_order->order[pos_max - 1];
      }
      return min_c;
    }();

    bool found = false;
    std::size_t best_bundle = 0;
    double best_ratio = std::numeric_limits<double>::infinity();
    int64_t best_spatial_delta = std::numeric_limits<int64_t>::max();
    std::string_view best_name;

    const odb::Point recv_pt = chain_point(receiver);
    const odb::Point donor_pt = chain_point(max_c);

    for (const std::size_t bi : chain_to_bundles[max_c]) {
      if (!is_movable(bi)) {
        continue;
      }
      if (chain_to_bundles[max_c].size() <= 1) {
        continue;
      }
      const uint64_t bbits = bundles[bi].bits;
      if (hard_max_length && used_bits[receiver] + bbits > max_length) {
        continue;
      }

      const uint64_t new_donor = max_bits - bbits;
      const uint64_t new_recv = used_bits[receiver] + bbits;
      if (new_donor == 0) {
        continue;
      }

      uint64_t new_min = std::numeric_limits<uint64_t>::max();
      uint64_t new_max = 0;
      for (std::size_t c = 0; c < chain_count; ++c) {
        uint64_t v = used_bits[c];
        if (c == max_c) {
          v = new_donor;
        } else if (c == receiver) {
          v = new_recv;
        }
        new_min = std::min(new_min, v);
        new_max = std::max(new_max, v);
      }
      if (new_min == 0) {
        continue;
      }
      const double ratio
          = static_cast<double>(new_max) / static_cast<double>(new_min);

      // Spatial tie-break when multiple moves yield the same imbalance ratio:
      // prefer moves that keep chains geographically contiguous.
      const int64_t dist_to_recv = manhattanDist(bundles[bi].origin, recv_pt);
      const int64_t dist_to_donor = manhattanDist(bundles[bi].origin, donor_pt);
      const int64_t spatial_delta = dist_to_recv - dist_to_donor;
      const std::string_view name = bundles[bi].name;

      if (!found || ratio < best_ratio - 1e-12
          || (std::abs(ratio - best_ratio) <= 1e-12
              && (spatial_delta < best_spatial_delta
                  || (spatial_delta == best_spatial_delta && name < best_name)))) {
        found = true;
        best_ratio = ratio;
        best_bundle = bi;
        best_spatial_delta = spatial_delta;
        best_name = name;
      }
    }

    if (!found) {
      if (logger) {
        logger->error(
            utl::DFT,
            127,
            "Scan architect constraints infeasible in hash domain {}: cannot "
            "satisfy max_imbalance={:.1f}% with the current grouping/assignment "
            "constraints.",
            hash_domain,
            config.getMaxImbalancePercent());
      }
    }

    move_bundle(best_bundle, max_c, receiver);
  }
}

// Multi-chain partitioning heuristic guardrail:
// If a chain contains a small "outlier" cluster separated by a very large
// coordinate gap, move the outlier cluster into a neighbouring chain (based on
// begin/end anchor ordering) to avoid long hop artifacts ("big jumps") that
// cannot be removed by intra-chain ordering alone.
void mitigateOutlierYGaps(
    const std::vector<PlacedBundle>& bundles,
    std::vector<std::size_t>& assignment,
    std::size_t chain_count,
    uint64_t max_length,
    bool hard_max_length,
    utl::Logger* logger,
    std::size_t hash_domain,
    const std::vector<std::optional<odb::Point>>& chain_anchor)
{
  if (chain_count <= 1 || bundles.empty()) {
    return;
  }

  const AnchorOrder y_order
      = ComputeAnchorOrderForAxis(chain_anchor, AnchorAxis::Y);
  const AnchorOrder x_order
      = ComputeAnchorOrderForAxis(chain_anchor, AnchorAxis::X);
  const auto axis_order = [&](AnchorAxis axis) -> const AnchorOrder* {
    const AnchorOrder& o = (axis == AnchorAxis::X) ? x_order : y_order;
    if (o.axis == AnchorAxis::None || o.order.size() != chain_count
        || o.pos.size() != chain_count) {
      return nullptr;
    }
    return &o;
  };

  constexpr int kMaxPasses = 8;
  for (int pass = 0; pass < kMaxPasses; ++pass) {
    std::vector<uint64_t> used_bits(chain_count, 0);
    std::vector<int64_t> sum_x(chain_count, 0);
    std::vector<int64_t> sum_y(chain_count, 0);
    std::vector<uint64_t> sum_bits(chain_count, 0);
    std::vector<std::vector<std::size_t>> chain_to_bundles(chain_count);
    chain_to_bundles.reserve(chain_count);
    for (std::size_t i = 0; i < bundles.size(); ++i) {
      const std::size_t c = assignment[i];
      if (c >= chain_count) {
        continue;
      }
      used_bits[c] += bundles[i].bits;
      sum_x[c] += static_cast<int64_t>(bundles[i].origin.x())
                  * static_cast<int64_t>(bundles[i].bits);
      sum_y[c] += static_cast<int64_t>(bundles[i].origin.y())
                  * static_cast<int64_t>(bundles[i].bits);
      sum_bits[c] += bundles[i].bits;
      chain_to_bundles[c].push_back(i);
    }

    std::vector<std::optional<odb::Point>> chain_centroid(chain_count);
    for (std::size_t c = 0; c < chain_count; ++c) {
      if (sum_bits[c] == 0) {
        continue;
      }
      const int64_t cx = sum_x[c] / static_cast<int64_t>(sum_bits[c]);
      const int64_t cy = sum_y[c] / static_cast<int64_t>(sum_bits[c]);
      chain_centroid[c] = odb::Point(static_cast<int>(cx), static_cast<int>(cy));
    }

    struct GapCandidate
    {
      AnchorAxis axis{AnchorAxis::None};
      std::size_t chain = 0;
      int max_gap = 0;
      std::size_t split_idx = 0;  // gap between split_idx and split_idx+1
      std::size_t outlier_size = 0;
      bool low_outliers = true;  // outliers are on the low side of the gap
    };

    std::optional<GapCandidate> best;
    for (const AnchorAxis axis : {AnchorAxis::Y, AnchorAxis::X}) {
      const auto bundle_axis = [&](std::size_t idx) -> int {
        if (idx >= bundles.size()) {
          return 0;
        }
        return (axis == AnchorAxis::X) ? bundles[idx].origin.x()
                                       : bundles[idx].origin.y();
      };

      for (std::size_t c = 0; c < chain_count; ++c) {
        const auto& idxs = chain_to_bundles[c];
        const std::size_t m = idxs.size();
        if (m < 32) {
          continue;
        }

        std::vector<std::pair<int, std::size_t>> coords;
        coords.reserve(m);
        for (const std::size_t bi : idxs) {
          coords.emplace_back(bundle_axis(bi), bi);
        }
        std::sort(coords.begin(),
                  coords.end(),
                  [&](const auto& a, const auto& b) {
                    if (a.first != b.first) {
                      return a.first < b.first;
                    }
                    return bundles[a.second].name < bundles[b.second].name;
                  });

        std::vector<int> gaps;
        gaps.reserve(m);
        int max_gap = 0;
        std::size_t split_idx = 0;
        for (std::size_t i = 0; i + 1 < coords.size(); ++i) {
          const int gap = coords[i + 1].first - coords[i].first;
          gaps.push_back(gap);
          if (gap > max_gap) {
            max_gap = gap;
            split_idx = i;
          }
        }
        if (max_gap <= 0 || gaps.empty()) {
          continue;
        }

        const std::size_t left = split_idx + 1;
        const std::size_t right = m - left;
        const std::size_t outlier_size = std::min(left, right);
        const std::size_t max_outlier = std::max<std::size_t>(8, m / 20);  // <=5%
        if (outlier_size == 0 || outlier_size > max_outlier) {
          continue;
        }

        // Require the max gap to be much larger than typical within-chain spacing.
        const std::size_t mid = gaps.size() / 2;
        std::nth_element(gaps.begin(), gaps.begin() + mid, gaps.end());
        const int median_gap = gaps[mid];
        if (median_gap > 0 && max_gap < 10 * median_gap) {
          continue;
        }

        GapCandidate cand;
        cand.axis = axis;
        cand.chain = c;
        cand.max_gap = max_gap;
        cand.split_idx = split_idx;
        cand.outlier_size = outlier_size;
        cand.low_outliers = (left <= right);

        if (!best.has_value() || cand.max_gap > best->max_gap
            || (cand.max_gap == best->max_gap
                && (cand.outlier_size < best->outlier_size
                    || (cand.outlier_size == best->outlier_size
                        && cand.chain < best->chain)))) {
          best = cand;
        }
      }
    }

    if (!best.has_value()) {
      break;
    }

    const std::size_t donor = best->chain;
    const AnchorAxis axis = best->axis;
    const auto bundle_axis = [&](std::size_t idx) -> int {
      if (idx >= bundles.size()) {
        return 0;
      }
      return (axis == AnchorAxis::X) ? bundles[idx].origin.x()
                                     : bundles[idx].origin.y();
    };

    const auto& donor_idxs = chain_to_bundles[donor];
    std::vector<std::pair<int, std::size_t>> coords;
    coords.reserve(donor_idxs.size());
    for (const std::size_t bi : donor_idxs) {
      coords.emplace_back(bundle_axis(bi), bi);
    }
    std::sort(coords.begin(),
              coords.end(),
              [&](const auto& a, const auto& b) {
                if (a.first != b.first) {
                  return a.first < b.first;
                }
                return bundles[a.second].name < bundles[b.second].name;
              });

    const std::size_t left = best->split_idx + 1;
    const std::size_t right = coords.size() - left;
    const bool move_left = best->low_outliers;
    const std::size_t outlier_count = move_left ? left : right;

    uint64_t outlier_bits = 0;
    std::vector<std::size_t> outliers;
    outliers.reserve(outlier_count);
    if (move_left) {
      for (std::size_t i = 0; i < left; ++i) {
        outliers.push_back(coords[i].second);
        outlier_bits += bundles[coords[i].second].bits;
      }
    } else {
      for (std::size_t i = left; i < coords.size(); ++i) {
        outliers.push_back(coords[i].second);
        outlier_bits += bundles[coords[i].second].bits;
      }
    }

    const auto outlier_center = [&]() -> std::optional<odb::Point> {
      if (outliers.empty() || outlier_bits == 0) {
        return std::nullopt;
      }
      int64_t sx = 0;
      int64_t sy = 0;
      uint64_t sb = 0;
      for (const std::size_t bi : outliers) {
        const uint64_t bits = bundles[bi].bits;
        sx += static_cast<int64_t>(bundles[bi].origin.x())
              * static_cast<int64_t>(bits);
        sy += static_cast<int64_t>(bundles[bi].origin.y())
              * static_cast<int64_t>(bits);
        sb += bits;
      }
      if (sb == 0) {
        return std::nullopt;
      }
      return odb::Point(static_cast<int>(sx / static_cast<int64_t>(sb)),
                        static_cast<int>(sy / static_cast<int64_t>(sb)));
    }();

    const auto chain_point = [&](std::size_t c) -> std::optional<odb::Point> {
      if (c >= chain_count) {
        return std::nullopt;
      }
      if (chain_centroid[c].has_value()) {
        return chain_centroid[c];
      }
      if (c < chain_anchor.size() && chain_anchor[c].has_value()) {
        return chain_anchor[c];
      }
      return std::nullopt;
    };

    const auto feasible_target = [&](std::size_t c) -> bool {
      if (c >= chain_count || c == donor) {
        return false;
      }
      return !(hard_max_length && used_bits[c] + outlier_bits > max_length);
    };

    const std::size_t target = [&]() -> std::size_t {
      const std::optional<odb::Point> center = outlier_center;

      // Prefer moving to an adjacent chain in anchor order when available for
      // this axis; otherwise fall back to nearest centroid/anchor chain.
      if (const AnchorOrder* ord = axis_order(axis)) {
        const std::size_t pos = ord->pos[donor];
        std::vector<std::size_t> candidates;
        candidates.reserve(2);
        if (pos > 0) {
          candidates.push_back(ord->order[pos - 1]);
        }
        if (pos + 1 < ord->order.size()) {
          candidates.push_back(ord->order[pos + 1]);
        }
        // Filter infeasible.
        candidates.erase(
            std::remove_if(candidates.begin(),
                           candidates.end(),
                           [&](std::size_t c) { return !feasible_target(c); }),
            candidates.end());
        if (!candidates.empty()) {
          // If we have a center, choose the adjacent chain closest to it.
          if (center.has_value()) {
            std::size_t best_c = candidates.front();
            int64_t best_d = std::numeric_limits<int64_t>::max();
            for (const std::size_t c : candidates) {
              const auto pt = chain_point(c);
              if (!pt.has_value()) {
                continue;
              }
              const int64_t d = manhattanDist(center.value(), pt.value());
              if (d < best_d || (d == best_d && c < best_c)) {
                best_d = d;
                best_c = c;
              }
            }
            return best_c;
          }
          // Otherwise, follow the low/high hint to preserve ordering.
          if (best->low_outliers && pos > 0
              && feasible_target(ord->order[pos - 1])) {
            return ord->order[pos - 1];
          }
          if (!best->low_outliers && pos + 1 < ord->order.size()
              && feasible_target(ord->order[pos + 1])) {
            return ord->order[pos + 1];
          }
          return candidates.front();
        }
      }

      if (!center.has_value()) {
        return donor;
      }
      std::size_t best_c = donor;
      int64_t best_d = std::numeric_limits<int64_t>::max();
      for (std::size_t c = 0; c < chain_count; ++c) {
        if (!feasible_target(c)) {
          continue;
        }
        const auto pt = chain_point(c);
        if (!pt.has_value()) {
          continue;
        }
        const int64_t d = manhattanDist(center.value(), pt.value());
        if (d < best_d || (d == best_d && c < best_c)) {
          best_d = d;
          best_c = c;
        }
      }
      return best_c;
    }();

    if (target == donor) {
      break;
    }

    bool moved_any = false;
    for (const std::size_t bi : outliers) {
      if (bundles[bi].fixed_chain.has_value()) {
        continue;
      }
      assignment[bi] = target;
      moved_any = true;
    }

    if (moved_any && logger) {
      debugPrint(logger,
                 utl::DFT,
                 "partition",
                 1,
                 "Moved {} outlier bundle(s) ({} bits) from chain {} to chain {} "
                 "to mitigate a large {}-axis gap ({} DBU) in hash domain {}.",
                 outliers.size(),
                 outlier_bits,
                 donor,
                 target,
                 (axis == AnchorAxis::X ? "X" : "Y"),
                 best->max_gap,
                 hash_domain);
    }
    if (!moved_any) {
      break;
    }
  }
}

int computeWorstAxisGap(const std::vector<PlacedBundle>& bundles,
                        const std::vector<std::size_t>& assignment,
                        std::size_t chain_count,
                        AnchorAxis axis)
{
  if (chain_count == 0 || bundles.empty()) {
    return 0;
  }

  std::vector<std::vector<int>> coords(chain_count);
  for (std::size_t i = 0; i < bundles.size(); ++i) {
    const std::size_t c = assignment[i];
    if (c >= chain_count) {
      continue;
    }
    const int v = (axis == AnchorAxis::X) ? bundles[i].origin.x()
                                          : bundles[i].origin.y();
    coords[c].push_back(v);
  }

  int worst = 0;
  for (std::size_t c = 0; c < chain_count; ++c) {
    auto& v = coords[c];
    if (v.size() < 2) {
      continue;
    }
    std::sort(v.begin(), v.end());
    for (std::size_t i = 0; i + 1 < v.size(); ++i) {
      worst = std::max(worst, v[i + 1] - v[i]);
    }
  }
  return worst;
}

int computeWorstGap(const std::vector<PlacedBundle>& bundles,
                    const std::vector<std::size_t>& assignment,
                    std::size_t chain_count)
{
  return std::max(computeWorstAxisGap(bundles, assignment, chain_count, AnchorAxis::X),
                  computeWorstAxisGap(bundles, assignment, chain_count, AnchorAxis::Y));
}

int64_t computeWorstManhattanDiameter(const std::vector<PlacedBundle>& bundles,
                                      const std::vector<std::size_t>& assignment,
                                      std::size_t chain_count)
{
  if (chain_count == 0 || bundles.empty()) {
    return 0;
  }

  std::vector<int> min_x(chain_count, std::numeric_limits<int>::max());
  std::vector<int> max_x(chain_count, std::numeric_limits<int>::min());
  std::vector<int> min_y(chain_count, std::numeric_limits<int>::max());
  std::vector<int> max_y(chain_count, std::numeric_limits<int>::min());
  std::vector<bool> seen(chain_count, false);

  for (std::size_t i = 0; i < bundles.size(); ++i) {
    const std::size_t c = assignment[i];
    if (c >= chain_count) {
      continue;
    }
    const odb::Point& p = bundles[i].origin;
    min_x[c] = std::min(min_x[c], p.x());
    max_x[c] = std::max(max_x[c], p.x());
    min_y[c] = std::min(min_y[c], p.y());
    max_y[c] = std::max(max_y[c], p.y());
    seen[c] = true;
  }

  int64_t worst = 0;
  for (std::size_t c = 0; c < chain_count; ++c) {
    if (!seen[c]) {
      continue;
    }
    const int64_t dx
        = static_cast<int64_t>(max_x[c]) - static_cast<int64_t>(min_x[c]);
    const int64_t dy
        = static_cast<int64_t>(max_y[c]) - static_cast<int64_t>(min_y[c]);
    worst = std::max(worst, dx + dy);
  }
  return worst;
}

uint32_t hilbertXYToIndex(uint32_t n, uint32_t x, uint32_t y)
{
  uint64_t d = 0;
  for (uint32_t s = n / 2; s > 0; s /= 2) {
    const uint32_t rx = (x & s) ? 1U : 0U;
    const uint32_t ry = (y & s) ? 1U : 0U;
    d += static_cast<uint64_t>(s) * static_cast<uint64_t>(s)
         * static_cast<uint64_t>((3U * rx) ^ ry);
    if (ry == 0U) {
      if (rx == 1U) {
        x = (n - 1U) - x;
        y = (n - 1U) - y;
      }
      std::swap(x, y);
    }
  }
  return static_cast<uint32_t>(d);
}

std::vector<std::size_t> hilbertSweepAssignBundles(
    const std::vector<PlacedBundle>& bundles,
    std::size_t chain_count,
    uint64_t max_length,
    bool hard_max_length,
    std::size_t hash_domain,
    const std::vector<std::optional<odb::Point>>& chain_anchor,
    utl::Logger* logger)
{
  std::vector<std::size_t> assignment(bundles.size(), 0);
  if (chain_count == 0 || bundles.empty()) {
    return assignment;
  }

  std::vector<uint64_t> used_bits(chain_count, 0);
  std::vector<std::size_t> used_items(chain_count, 0);
  uint64_t total_bits = 0;

  for (std::size_t i = 0; i < bundles.size(); ++i) {
    total_bits += bundles[i].bits;
    if (!bundles[i].fixed_chain.has_value()) {
      continue;
    }
    const std::size_t c = bundles[i].fixed_chain.value();
    if (c >= chain_count) {
      if (logger) {
        logger->error(
            utl::DFT,
            244,
            "Scan constraints internal error: bundle fixed_chain index {} is "
            "out of range for chain_count {} (hash domain {}).",
            c,
            chain_count,
            hash_domain);
      }
      continue;
    }
    assignment[i] = c;
    used_bits[c] += bundles[i].bits;
    used_items[c] += 1;
    if (hard_max_length && used_bits[c] > max_length) {
      if (logger) {
        logger->error(
            utl::DFT,
            245,
            "Scan architect constraints infeasible in hash domain {}: chain {} "
            "overflows max_length={} due to fixed assignments (used {}).",
            hash_domain,
            c,
            max_length,
            used_bits[c]);
      }
    }
  }

  std::vector<std::size_t> chain_order(chain_count);
  std::iota(chain_order.begin(), chain_order.end(), 0);
  const AnchorOrder ord = ComputeAnchorOrderForAxis(chain_anchor, AnchorAxis::Y);
  if (ord.axis != AnchorAxis::None && ord.order.size() == chain_count) {
    chain_order = ord.order;
  }

  int min_x = std::numeric_limits<int>::max();
  int max_x = std::numeric_limits<int>::min();
  int min_y = std::numeric_limits<int>::max();
  int max_y = std::numeric_limits<int>::min();
  for (const auto& b : bundles) {
    min_x = std::min(min_x, b.origin.x());
    max_x = std::max(max_x, b.origin.x());
    min_y = std::min(min_y, b.origin.y());
    max_y = std::max(max_y, b.origin.y());
  }

  constexpr int kHilbertBits = 16;
  static constexpr uint32_t kHilbertN = 1U << kHilbertBits;
  static constexpr uint32_t kHilbertMax = kHilbertN - 1U;

  const auto scale = [&](int v, int v_min, int v_max) -> uint32_t {
    if (v_max <= v_min) {
      return 0;
    }
    const uint64_t num
        = static_cast<uint64_t>(static_cast<int64_t>(v) - v_min) * kHilbertMax;
    const uint64_t den
        = static_cast<uint64_t>(static_cast<int64_t>(v_max) - v_min);
    return static_cast<uint32_t>(num / den);
  };

  std::vector<std::size_t> free_bundles;
  free_bundles.reserve(bundles.size());
  for (std::size_t i = 0; i < bundles.size(); ++i) {
    if (!bundles[i].fixed_chain.has_value()) {
      free_bundles.push_back(i);
    }
  }

  std::sort(free_bundles.begin(),
            free_bundles.end(),
            [&](std::size_t a, std::size_t b) {
              const uint32_t ax = scale(bundles[a].origin.x(), min_x, max_x);
              const uint32_t ay = scale(bundles[a].origin.y(), min_y, max_y);
              const uint32_t bx = scale(bundles[b].origin.x(), min_x, max_x);
              const uint32_t by = scale(bundles[b].origin.y(), min_y, max_y);
              const uint32_t ha = hilbertXYToIndex(kHilbertN, ax, ay);
              const uint32_t hb = hilbertXYToIndex(kHilbertN, bx, by);
              if (ha != hb) {
                return ha < hb;
              }
              return bundles[a].name < bundles[b].name;
            });

  const uint64_t target_bits
      = (chain_count == 0) ? 0 : ((total_bits + chain_count - 1) / chain_count);
  std::size_t cur_pos = 0;

  for (std::size_t k = 0; k < free_bundles.size(); ++k) {
    const std::size_t idx = free_bundles[k];
    const uint64_t bits = bundles[idx].bits;

    const auto remaining_free = free_bundles.size() - k;
    while (cur_pos + 1 < chain_order.size()) {
      const std::size_t c = chain_order[cur_pos];
      if (hard_max_length && used_bits[c] + bits > max_length) {
        ++cur_pos;
        continue;
      }
      const std::size_t remaining_chains = chain_order.size() - cur_pos - 1;
      if (used_bits[c] >= target_bits && remaining_free > remaining_chains) {
        ++cur_pos;
        continue;
      }
      break;
    }

    std::size_t c = chain_order[std::min(cur_pos, chain_order.size() - 1)];
    if (hard_max_length && used_bits[c] + bits > max_length) {
      bool found = false;
      for (std::size_t p = cur_pos; p < chain_order.size(); ++p) {
        const std::size_t cand = chain_order[p];
        if (used_bits[cand] + bits <= max_length) {
          c = cand;
          cur_pos = p;
          found = true;
          break;
        }
      }
      if (!found) {
        if (logger) {
          logger->error(
              utl::DFT,
              246,
              "Scan architect constraints infeasible in hash domain {}: cannot "
              "pack constrained scan items into {} chains with max_length={}.",
              hash_domain,
              chain_count,
              max_length);
        }
      }
    }

    assignment[idx] = c;
    used_bits[c] += bits;
    used_items[c] += 1;
  }

  // Ensure every chain has at least one item when possible. If any chain is
  // empty, move a single non-fixed bundle from the heaviest chain.
  if (chain_count > 1) {
    std::optional<std::size_t> max_chain;
    uint64_t max_used = 0;
    for (std::size_t c = 0; c < chain_count; ++c) {
      if (!max_chain.has_value() || used_bits[c] > max_used) {
        max_chain = c;
        max_used = used_bits[c];
      }
    }
    if (max_chain.has_value()) {
      for (std::size_t c = 0; c < chain_count; ++c) {
        if (used_items[c] != 0) {
          continue;
        }
        std::optional<std::size_t> move_idx;
        for (std::size_t i = 0; i < bundles.size(); ++i) {
          if (assignment[i] != max_chain.value()) {
            continue;
          }
          if (bundles[i].fixed_chain.has_value()) {
            continue;
          }
          move_idx = i;
          break;
        }
        if (!move_idx.has_value()) {
          continue;
        }
        assignment[move_idx.value()] = c;
        used_items[c] = 1;
        used_bits[c] += bundles[move_idx.value()].bits;
        used_bits[max_chain.value()] -= bundles[move_idx.value()].bits;
      }
    }
  }

  return assignment;
}

std::vector<std::size_t> axisSweepAssignBundles(
    const std::vector<PlacedBundle>& bundles,
    std::size_t chain_count,
    uint64_t max_length,
    bool hard_max_length,
    std::size_t hash_domain,
    const std::vector<std::optional<odb::Point>>& chain_anchor,
    AnchorAxis axis,
    utl::Logger* logger)
{
  std::vector<std::size_t> assignment(bundles.size(), 0);
  if (chain_count == 0 || bundles.empty()) {
    return assignment;
  }

  std::vector<uint64_t> used_bits(chain_count, 0);
  std::vector<std::size_t> used_items(chain_count, 0);
  uint64_t total_bits = 0;

  for (std::size_t i = 0; i < bundles.size(); ++i) {
    total_bits += bundles[i].bits;
    if (!bundles[i].fixed_chain.has_value()) {
      continue;
    }
    const std::size_t c = bundles[i].fixed_chain.value();
    if (c >= chain_count) {
      if (logger) {
        logger->error(
            utl::DFT,
            231,
            "Scan constraints internal error: bundle fixed_chain index {} is "
            "out of range for chain_count {} (hash domain {}).",
            c,
            chain_count,
            hash_domain);
      }
    }
    assignment[i] = c;
    used_bits[c] += bundles[i].bits;
    used_items[c] += 1;
    if (hard_max_length && used_bits[c] > max_length) {
      if (logger) {
        logger->error(
            utl::DFT,
            232,
            "Scan architect constraints infeasible in hash domain {}: chain {} "
            "overflows max_length={} due to fixed assignments (used {}).",
            hash_domain,
            c,
            max_length,
            used_bits[c]);
      }
    }
  }

  std::vector<std::size_t> chain_order(chain_count);
  std::iota(chain_order.begin(), chain_order.end(), 0);
  const AnchorOrder ord = ComputeAnchorOrderForAxis(chain_anchor, axis);
  if (ord.axis != AnchorAxis::None && ord.order.size() == chain_count) {
    chain_order = ord.order;
  }

  std::vector<std::size_t> free_bundles;
  free_bundles.reserve(bundles.size());
  for (std::size_t i = 0; i < bundles.size(); ++i) {
    if (!bundles[i].fixed_chain.has_value()) {
      free_bundles.push_back(i);
    }
  }
  std::sort(free_bundles.begin(),
            free_bundles.end(),
            [&](std::size_t a, std::size_t b) {
              const int ax = (axis == AnchorAxis::X) ? bundles[a].origin.x()
                                                     : bundles[a].origin.y();
              const int bx = (axis == AnchorAxis::X) ? bundles[b].origin.x()
                                                     : bundles[b].origin.y();
              if (ax != bx) {
                return ax < bx;
              }
              return bundles[a].name < bundles[b].name;
            });

  const uint64_t target_bits
      = (chain_count == 0) ? 0 : ((total_bits + chain_count - 1) / chain_count);
  std::size_t cur_pos = 0;

  for (std::size_t k = 0; k < free_bundles.size(); ++k) {
    const std::size_t idx = free_bundles[k];
    const uint64_t bits = bundles[idx].bits;

    const auto remaining_free = free_bundles.size() - k;
    while (cur_pos + 1 < chain_order.size()) {
      const std::size_t c = chain_order[cur_pos];
      if (hard_max_length && used_bits[c] + bits > max_length) {
        ++cur_pos;
        continue;
      }
      const std::size_t remaining_chains = chain_order.size() - cur_pos - 1;
      if (used_bits[c] >= target_bits && remaining_free > remaining_chains) {
        ++cur_pos;
        continue;
      }
      break;
    }

    std::size_t c = chain_order[std::min(cur_pos, chain_order.size() - 1)];
    if (hard_max_length && used_bits[c] + bits > max_length) {
      bool found = false;
      for (std::size_t p = cur_pos; p < chain_order.size(); ++p) {
        const std::size_t cand = chain_order[p];
        if (used_bits[cand] + bits <= max_length) {
          c = cand;
          cur_pos = p;
          found = true;
          break;
        }
      }
      if (!found) {
        if (logger) {
          logger->error(
              utl::DFT,
              233,
              "Scan architect constraints infeasible in hash domain {}: cannot "
              "pack constrained scan items into {} chains with max_length={}.",
              hash_domain,
              chain_count,
              max_length);
        }
      }
    }

    assignment[idx] = c;
    used_bits[c] += bits;
    used_items[c] += 1;
  }

  // Ensure every chain has at least one item when possible. If any chain is
  // empty, move a single non-fixed bundle from the heaviest chain.
  if (chain_count > 1) {
    std::optional<std::size_t> max_chain;
    uint64_t max_used = 0;
    for (std::size_t c = 0; c < chain_count; ++c) {
      if (!max_chain.has_value() || used_bits[c] > max_used) {
        max_chain = c;
        max_used = used_bits[c];
      }
    }
    if (max_chain.has_value()) {
      for (std::size_t c = 0; c < chain_count; ++c) {
        if (used_items[c] != 0) {
          continue;
        }
        // Find a movable bundle currently in max_chain.
        std::optional<std::size_t> move_idx;
        for (std::size_t i = 0; i < bundles.size(); ++i) {
          if (assignment[i] != max_chain.value()) {
            continue;
          }
          if (bundles[i].fixed_chain.has_value()) {
            continue;
          }
          move_idx = i;
          break;
        }
        if (!move_idx.has_value()) {
          continue;
        }
        assignment[move_idx.value()] = c;
        used_items[c] = 1;
        used_bits[c] += bundles[move_idx.value()].bits;
        used_bits[max_chain.value()] -= bundles[move_idx.value()].bits;
      }
    }
  }

  return assignment;
}

std::vector<PlacedBundle> buildConstraintBundles(
    std::vector<std::unique_ptr<ScanCell>> domain_cells,
    std::size_t hash_domain,
    const ScanArchitectConfig& config,
    const std::unordered_map<std::string_view, ChainRef>& chain_name_to_ref,
    utl::Logger* logger)
{
  std::vector<PlacedBundle> bundles;
  if (domain_cells.empty()) {
    return bundles;
  }

  const std::size_t n = domain_cells.size();
  std::vector<std::string> names;
  names.reserve(n);
  std::unordered_map<std::string, std::size_t> name_to_idx;
  name_to_idx.reserve(n * 2);

  for (std::size_t i = 0; i < n; ++i) {
    std::string name(domain_cells[i]->getName());
    names.push_back(name);
    name_to_idx.emplace(std::move(name), i);
  }

  UnionFind uf(n);

  // Union group members.
  for (const auto& group : config.getScanOrderGroups()) {
    std::optional<std::size_t> first;
    for (const std::string& inst : group.inst_names) {
      auto it = name_to_idx.find(inst);
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

  // Union fixed-edge endpoints (must be in same chain).
  for (const auto& edge : config.getScanOrderFixedEdges()) {
    auto it_from = name_to_idx.find(edge.from_inst);
    auto it_to = name_to_idx.find(edge.to_inst);
    if (it_from == name_to_idx.end() || it_to == name_to_idx.end()) {
      continue;
    }
    uf.unite(it_from->second, it_to->second);
  }

  std::unordered_map<std::size_t, std::vector<std::size_t>> comps;
  comps.reserve(n);
  for (std::size_t i = 0; i < n; ++i) {
    comps[uf.find(i)].push_back(i);
  }

  // Deterministic component order by smallest member name.
  std::vector<std::pair<std::string, std::size_t>> comp_order;
  comp_order.reserve(comps.size());
  for (const auto& [root, idxs] : comps) {
    std::string best = names[idxs.front()];
    for (const std::size_t idx : idxs) {
      if (names[idx] < best) {
        best = names[idx];
      }
    }
    comp_order.emplace_back(best, root);
  }
  std::sort(comp_order.begin(), comp_order.end(), [](const auto& a, const auto& b) {
    return a.first < b.first;
  });

  std::vector<bool> moved(n, false);
  bundles.reserve(comp_order.size());
  for (const auto& [rep_name, root] : comp_order) {
    const auto& idxs = comps.at(root);
    PlacedBundle bundle;
    bundle.name = rep_name;

    int64_t sum_x = 0;
    int64_t sum_y = 0;
    uint64_t sum_bits = 0;
    bool all_placed = true;

    for (const std::size_t idx : idxs) {
      if (moved[idx]) {
        continue;
      }
      moved[idx] = true;
      std::unique_ptr<ScanCell> cell = std::move(domain_cells[idx]);
      if (!cell) {
        continue;
      }
      const uint64_t bits = cell->getBits();
      sum_bits += bits;
      const odb::Point origin = cell->getOrigin();
      sum_x += static_cast<int64_t>(origin.x()) * static_cast<int64_t>(bits);
      sum_y += static_cast<int64_t>(origin.y()) * static_cast<int64_t>(bits);
      all_placed &= cell->isPlaced();
      bundle.cells.push_back(std::move(cell));
    }

    bundle.bits = sum_bits;
    bundle.all_placed = all_placed;
    if (sum_bits != 0) {
      const int64_t cx = sum_x / static_cast<int64_t>(sum_bits);
      const int64_t cy = sum_y / static_cast<int64_t>(sum_bits);
      bundle.origin = odb::Point(static_cast<int>(cx), static_cast<int>(cy));
    } else if (!bundle.cells.empty()) {
      bundle.origin = bundle.cells.front()->getOrigin();
    }

    if (bundle.bits == 0) {
      continue;
    }

    // Sort cells within the bundle by name for determinism before later
    // ordering re-sorts them.
    std::stable_sort(
        bundle.cells.begin(), bundle.cells.end(), [](const auto& a, const auto& b) {
          return a->getName() < b->getName();
        });

    // Resolve fixed chain assignment for this bundle (if any). A conflict here
    // means constraints force these instances into the same chain, but they
    // were assigned to different chains.
    std::optional<std::string_view> fixed_chain_name;
    for (const auto& cell : bundle.cells) {
      const std::optional<std::string_view> asn
          = config.getAssignedChainForInstance(cell->getName());
      if (!asn.has_value()) {
        continue;
      }
      const auto it = chain_name_to_ref.find(asn.value());
      if (it == chain_name_to_ref.end()) {
        if (logger) {
          logger->error(utl::DFT,
                        161,
                        "Scan constraints: instance '{}' is assigned to unknown "
                        "chain '{}' (in hash domain {}).",
                        cell->getName(),
                        asn.value(),
                        hash_domain);
        }
      }
      if (config.getClockMixing() == ScanArchitectConfig::ClockMixing::NoMix
          && it->second.hash_domain != hash_domain) {
        if (logger) {
          logger->error(
              utl::DFT,
              162,
              "Scan constraints: instance '{}' in hash domain {} is assigned to "
              "chain '{}' in a different hash domain {}.",
              cell->getName(),
              hash_domain,
              asn.value(),
              it->second.hash_domain);
        }
      }
      if (!fixed_chain_name.has_value()) {
        fixed_chain_name = asn.value();
        bundle.fixed_chain = it->second.local_index;
      } else if (fixed_chain_name.value() != asn.value()) {
        if (logger) {
          logger->error(
              utl::DFT,
              163,
              "Scan constraints: instances forced into one chain cannot be "
              "assigned to different chains ('{}' vs '{}').",
              fixed_chain_name.value(),
              asn.value());
        }
      }
    }

    bundles.push_back(std::move(bundle));
  }

  return bundles;
}

std::vector<std::vector<std::unique_ptr<ScanCell>>> clusterPlacedScanCells(
    std::vector<std::unique_ptr<ScanCell>> domain_cells,
    std::size_t hash_domain,
    const ScanArchitectConfig& config,
    const std::unordered_map<std::string_view, ChainRef>& chain_name_to_ref,
    const std::unordered_map<std::string_view, std::size_t>&
        chain_name_to_ordinal,
    std::size_t chain_count,
    uint64_t max_length,
    bool hard_max_length,
    utl::Logger* logger)
{
  // Backward-compatible wrapper: treat every cell as its own bundle.
  std::vector<PlacedBundle> bundles
      = buildConstraintBundles(
          std::move(domain_cells), hash_domain, config, chain_name_to_ref, logger);

  std::vector<std::vector<std::unique_ptr<ScanCell>>> clustered(chain_count);
  if (chain_count == 0 || bundles.empty()) {
    return clustered;
  }

  const std::vector<ChainEndpointPoints> chain_endpoints
      = ResolveChainEndpointsForDomain(
          bundles,
          hash_domain,
          config,
          chain_name_to_ref,
          chain_name_to_ordinal,
          chain_count);
  const auto endpoint_bias = [&](std::size_t bundle_idx,
                                 std::size_t chain_idx) -> int64_t {
    if (bundle_idx >= bundles.size() || chain_idx >= chain_endpoints.size()) {
      return 0;
    }
    int count = 0;
    int64_t sum = 0;
    if (chain_endpoints[chain_idx].begin.has_value()) {
      sum += manhattanDist(bundles[bundle_idx].origin,
                           chain_endpoints[chain_idx].begin.value());
      ++count;
    }
    if (chain_endpoints[chain_idx].end.has_value()) {
      sum += manhattanDist(bundles[bundle_idx].origin,
                           chain_endpoints[chain_idx].end.value());
      ++count;
    }
    if (count == 0) {
      return 0;
    }
    // Endpoint locations are a proxy for pin access; bias clustering toward
    // chains that are geometrically compatible with their begin/end ports.
    return (sum / count) / 4;  // weight 0.25 vs centroid distance
  };

  std::vector<std::optional<odb::Point>> chain_anchor(chain_count);
  for (std::size_t c = 0; c < chain_count; ++c) {
    const auto& ep = chain_endpoints[c];
    if (ep.begin.has_value() && ep.end.has_value()) {
      const int64_t ax = (static_cast<int64_t>(ep.begin->x())
                          + static_cast<int64_t>(ep.end->x()))
                         / 2;
      const int64_t ay = (static_cast<int64_t>(ep.begin->y())
                          + static_cast<int64_t>(ep.end->y()))
                         / 2;
      chain_anchor[c] = odb::Point(static_cast<int>(ax), static_cast<int>(ay));
    } else if (ep.begin.has_value()) {
      chain_anchor[c] = ep.begin;
    } else if (ep.end.has_value()) {
      chain_anchor[c] = ep.end;
    }
  }

  uint64_t total_bits = 0;
  uint64_t max_item_bits = 0;
  for (const auto& b : bundles) {
    total_bits += b.bits;
    max_item_bits = std::max(max_item_bits, b.bits);
  }

  if (max_length == 0) {
    // Defensive: fall back to a single chain assignment.
    for (auto& b : bundles) {
      for (auto& cell : b.cells) {
        clustered.front().push_back(std::move(cell));
      }
    }
    return clustered;
  }

  if (max_item_bits > max_length) {
    if (hard_max_length) {
      logger->error(utl::DFT,
                    164,
                    "Scan architect constraints infeasible in hash domain {}: "
                    "max_length={} is smaller than the largest constrained scan "
                    "item ({} bits).",
                    hash_domain,
                    max_length,
                    max_item_bits);
    }
    // Soft max_length (inferred): bump capacity to fit the largest item so we
    // never need to overflow during packing.
    max_length = max_item_bits;
  }
  if (hard_max_length && total_bits > chain_count * max_length) {
    logger->error(utl::DFT,
                  165,
                  "Scan architect constraints infeasible in hash domain {}: "
                  "{} bits over {} chains with max_length={} (capacity {}).",
                  hash_domain,
                  total_bits,
                  chain_count,
                  max_length,
                  chain_count * max_length);
  }

  const std::size_t n = bundles.size();
  if (config.getChainCount().has_value() && total_bits != 0
      && chain_count > n) {
    logger->error(utl::DFT,
                  166,
                  "Scan architect constraints infeasible in hash domain {}: "
                  "chain_count={} exceeds the number of constrained scan items "
                  "({}); cannot create non-empty chains.",
                  hash_domain,
                  chain_count,
                  n);
  }

  // Seed each chain either from a fixed assignment or from farthest-first
  // selection among free bundles.
  std::vector<std::vector<std::size_t>> fixed_in_chain(chain_count);
  std::vector<std::size_t> free_bundles;
  free_bundles.reserve(n);
  for (std::size_t i = 0; i < n; ++i) {
    if (bundles[i].fixed_chain.has_value()) {
      const std::size_t c = bundles[i].fixed_chain.value();
      if (c >= chain_count) {
        logger->error(
            utl::DFT,
            167,
            "Scan constraints internal error: bundle fixed_chain index {} is "
            "out of range for chain_count {} (hash domain {}).",
            c,
            chain_count,
            hash_domain);
      }
      fixed_in_chain[c].push_back(i);
    } else {
      free_bundles.push_back(i);
    }
  }

  const std::size_t k = std::min(chain_count, n);

  std::vector<bool> chain_has_fixed(chain_count, false);
  for (std::size_t c = 0; c < chain_count; ++c) {
    chain_has_fixed[c] = !fixed_in_chain[c].empty();
  }

  std::vector<std::size_t> unseeded_chains;
  unseeded_chains.reserve(chain_count);
  for (std::size_t c = 0; c < chain_count; ++c) {
    if (!chain_has_fixed[c]) {
      unseeded_chains.push_back(c);
    }
  }

  // Seed bundles to avoid empty-chain local minima in K-means. Favor chains
  // with begin/end endpoints by picking initial seeds closest to their anchors.
  std::vector<std::optional<std::size_t>> seed_chain_for_bundle(n);
  std::vector<bool> bundle_seeded(n, false);
  std::vector<std::size_t> remaining_chains;
  remaining_chains.reserve(unseeded_chains.size());

  for (const std::size_t c : unseeded_chains) {
    if (!chain_anchor[c].has_value() || free_bundles.empty()) {
      remaining_chains.push_back(c);
      continue;
    }

    std::size_t best = 0;
    bool found = false;
    int64_t best_dist = std::numeric_limits<int64_t>::max();
    for (const std::size_t idx : free_bundles) {
      if (bundle_seeded[idx]) {
        continue;
      }
      const int64_t dist
          = manhattanDist(bundles[idx].origin, chain_anchor[c].value());
      if (!found || dist < best_dist
          || (dist == best_dist && bundles[idx].name < bundles[best].name)) {
        found = true;
        best = idx;
        best_dist = dist;
      }
    }
    if (!found) {
      remaining_chains.push_back(c);
      continue;
    }
    seed_chain_for_bundle[best] = c;
    bundle_seeded[best] = true;
  }

  std::vector<std::size_t> seed_indices;
  seed_indices.reserve(std::min(unseeded_chains.size(), free_bundles.size()));
  for (std::size_t i = 0; i < n; ++i) {
    if (seed_chain_for_bundle[i].has_value()) {
      seed_indices.push_back(i);
    }
  }

  const std::size_t target_seed_count = std::min(unseeded_chains.size(), k);
  const std::size_t have_seed_count = seed_indices.size();
  const std::size_t need
      = (have_seed_count < target_seed_count) ? (target_seed_count - have_seed_count)
                                              : 0;

  std::vector<std::size_t> extra_seeds;
  extra_seeds.reserve(need);
  if (need > 0 && !free_bundles.empty()) {
    // Farthest-first among remaining bundles.
    if (seed_indices.empty()) {
      std::size_t first = free_bundles.front();
      int64_t lowest = std::numeric_limits<int64_t>::max();
      for (const std::size_t idx : free_bundles) {
        if (bundle_seeded[idx]) {
          continue;
        }
        const auto& p = bundles[idx].origin;
        const int64_t score
            = static_cast<int64_t>(p.x()) + static_cast<int64_t>(p.y());
        if (score < lowest
            || (score == lowest && bundles[idx].name < bundles[first].name)) {
          first = idx;
          lowest = score;
        }
      }
      seed_indices.push_back(first);
      extra_seeds.push_back(first);
      bundle_seeded[first] = true;
    }

    while (extra_seeds.size() < need) {
      std::size_t best = 0;
      bool found = false;
      int64_t best_score = -1;
      for (const std::size_t idx : free_bundles) {
        if (bundle_seeded[idx]) {
          continue;
        }
        int64_t nearest = std::numeric_limits<int64_t>::max();
        for (const std::size_t s : seed_indices) {
          nearest = std::min(nearest,
                             manhattanDist(bundles[idx].origin, bundles[s].origin));
        }
        if (!found || nearest > best_score
            || (nearest == best_score && bundles[idx].name < bundles[best].name)) {
          best = idx;
          best_score = nearest;
          found = true;
        }
      }
      if (!found) {
        break;
      }
      seed_indices.push_back(best);
      extra_seeds.push_back(best);
      bundle_seeded[best] = true;
    }
  }

  // Assign extra seeds to remaining chains deterministically.
  for (std::size_t i = 0; i < extra_seeds.size() && i < remaining_chains.size();
       ++i) {
    seed_chain_for_bundle[extra_seeds[i]] = remaining_chains[i];
  }

  std::vector<odb::Point> centroids(chain_count, bundles.front().origin);
  // Initialize centroids from fixed bundles when present.
  for (std::size_t c = 0; c < chain_count; ++c) {
    if (fixed_in_chain[c].empty()) {
      continue;
    }
    int64_t sum_x = 0;
    int64_t sum_y = 0;
    uint64_t sum_bits = 0;
    for (const std::size_t idx : fixed_in_chain[c]) {
      sum_x += static_cast<int64_t>(bundles[idx].origin.x())
               * static_cast<int64_t>(bundles[idx].bits);
      sum_y += static_cast<int64_t>(bundles[idx].origin.y())
               * static_cast<int64_t>(bundles[idx].bits);
      sum_bits += bundles[idx].bits;
    }
    if (sum_bits != 0) {
      centroids[c]
          = odb::Point(static_cast<int>(sum_x / static_cast<int64_t>(sum_bits)),
                       static_cast<int>(sum_y / static_cast<int64_t>(sum_bits)));
    } else {
      centroids[c] = bundles[fixed_in_chain[c].front()].origin;
    }
  }
  // Initialize remaining centroids from seeds/anchors when available.
  std::vector<bool> chain_has_seed(chain_count, false);
  for (std::size_t i = 0; i < n; ++i) {
    if (!seed_chain_for_bundle[i].has_value()) {
      continue;
    }
    const std::size_t c = seed_chain_for_bundle[i].value();
    chain_has_seed[c] = true;
    if (fixed_in_chain[c].empty()) {
      centroids[c] = bundles[i].origin;
    }
  }
  for (std::size_t c = 0; c < chain_count; ++c) {
    if (!fixed_in_chain[c].empty()) {
      continue;
    }
    if (!chain_has_seed[c] && chain_anchor[c].has_value()) {
      // If we didn't get a seed bundle (e.g., due to constraints), at least
      // start K-means near the desired begin/end location.
      centroids[c] = chain_anchor[c].value();
    }
  }

  std::vector<std::size_t> order(n);
  std::iota(order.begin(), order.end(), 0);
  std::stable_sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) {
      if (bundles[a].bits != bundles[b].bits) {
        return bundles[a].bits > bundles[b].bits;  // pack large items first
      }
      return bundles[a].name < bundles[b].name;
    });

  std::vector<std::size_t> assignment(n, 0);
  std::vector<std::size_t> prev_assignment(n, 0);

  constexpr int kMaxIters = 5;
  for (int iter = 0; iter < kMaxIters; ++iter) {
    prev_assignment = assignment;
    std::vector<uint64_t> used_bits(chain_count, 0);
    std::vector<bool> is_assigned(n, false);

    // Apply fixed chain assignments.
    for (std::size_t i = 0; i < n; ++i) {
      if (!bundles[i].fixed_chain.has_value()) {
        continue;
      }
      const std::size_t c = bundles[i].fixed_chain.value();
      assignment[i] = c;
      is_assigned[i] = true;
      used_bits[c] += bundles[i].bits;
      if (hard_max_length && used_bits[c] > max_length) {
        logger->error(
            utl::DFT,
            168,
            "Scan architect constraints infeasible in hash domain {}: chain {} "
            "overflows max_length={} due to fixed assignments (used {}).",
            hash_domain,
            c,
            max_length,
            used_bits[c]);
      }
    }

    // Pin centroid seed bundles (only for chains without fixed bundles).
    for (std::size_t i = 0; i < n; ++i) {
      if (!seed_chain_for_bundle[i].has_value()) {
        continue;
      }
      if (is_assigned[i]) {
        continue;
      }
      const std::size_t c = seed_chain_for_bundle[i].value();
      assignment[i] = c;
      is_assigned[i] = true;
      used_bits[c] += bundles[i].bits;
      if (hard_max_length && used_bits[c] > max_length) {
        logger->error(
            utl::DFT,
            169,
            "Scan architect constraints infeasible in hash domain {}: seed "
            "assignment overflows max_length={} for chain {} (used {}).",
            hash_domain,
            max_length,
            c,
            used_bits[c]);
      }
    }

    for (const std::size_t idx : order) {
      if (is_assigned[idx]) {
        continue;
      }

      bool found = false;
      std::size_t best_chain = 0;
      int64_t best_dist = std::numeric_limits<int64_t>::max();
      uint64_t best_used = std::numeric_limits<uint64_t>::max();

      for (std::size_t c = 0; c < chain_count; ++c) {
        if (hard_max_length && used_bits[c] + bundles[idx].bits > max_length) {
          continue;
        }
        const int64_t dist = manhattanDist(bundles[idx].origin, centroids[c])
                             + endpoint_bias(idx, c);
        if (!found || dist < best_dist
            || (dist == best_dist
                && (used_bits[c] < best_used
                    || (used_bits[c] == best_used && c < best_chain)))) {
          found = true;
          best_chain = c;
          best_dist = dist;
          best_used = used_bits[c];
        }
      }

      if (!found) {
        if (hard_max_length) {
        logger->error(
            utl::DFT,
            170,
            "Scan architect constraints infeasible in hash domain {}: cannot "
            "pack constrained scan items into {} chains with max_length={}.",
            hash_domain,
            chain_count,
              max_length);
        }
        // Soft max_length: pick the least-overflowing chain.
        best_chain = 0;
        uint64_t best_overflow = std::numeric_limits<uint64_t>::max();
        for (std::size_t c = 0; c < chain_count; ++c) {
          const uint64_t new_bits = used_bits[c] + bundles[idx].bits;
          const uint64_t overflow
              = (new_bits > max_length) ? (new_bits - max_length) : 0;
          if (overflow < best_overflow
              || (overflow == best_overflow && c < best_chain)) {
            best_overflow = overflow;
            best_chain = c;
          }
        }
      }

      assignment[idx] = best_chain;
      is_assigned[idx] = true;
      used_bits[best_chain] += bundles[idx].bits;
    }

    // Recompute centroids (weighted by bits). Keep centroid for empty chains.
    std::vector<int64_t> sum_x(chain_count, 0);
    std::vector<int64_t> sum_y(chain_count, 0);
    std::vector<uint64_t> sum_bits(chain_count, 0);
    for (std::size_t i = 0; i < n; ++i) {
      const std::size_t c = assignment[i];
      sum_x[c] += static_cast<int64_t>(bundles[i].origin.x())
                  * static_cast<int64_t>(bundles[i].bits);
      sum_y[c] += static_cast<int64_t>(bundles[i].origin.y())
                  * static_cast<int64_t>(bundles[i].bits);
      sum_bits[c] += bundles[i].bits;
    }
    for (std::size_t c = 0; c < chain_count; ++c) {
      if (sum_bits[c] == 0) {
        continue;
      }
      const int64_t cx = sum_x[c] / static_cast<int64_t>(sum_bits[c]);
      const int64_t cy = sum_y[c] / static_cast<int64_t>(sum_bits[c]);
      centroids[c] = odb::Point(static_cast<int>(cx), static_cast<int>(cy));
    }

    const bool changed = (assignment != prev_assignment);
    if (!changed) {
      break;
    }
  }

  // If the placement-aware K-means style assignment produces chains with very
  // large internal coordinate gaps, try a simple axis-sweep partitioning and
  // keep whichever assignment has the smaller worst (X/Y) gap. This is a cheap
  // guardrail against highly non-local multi-chain solutions that lead to
  // obvious long-hop artifacts even after intra-chain ordering.
  const int kmeans_gap = computeWorstGap(bundles, assignment, chain_count);
  const int64_t kmeans_diam
      = computeWorstManhattanDiameter(bundles, assignment, chain_count);
  std::vector<std::size_t> sweep_y = axisSweepAssignBundles(bundles,
                                                            chain_count,
                                                            max_length,
                                                            hard_max_length,
                                                            hash_domain,
                                                            chain_anchor,
                                                            AnchorAxis::Y,
                                                            logger);
  const int sweep_y_gap = computeWorstGap(bundles, sweep_y, chain_count);
  const int64_t sweep_y_diam
      = computeWorstManhattanDiameter(bundles, sweep_y, chain_count);
  std::vector<std::size_t> sweep_x = axisSweepAssignBundles(bundles,
                                                            chain_count,
                                                            max_length,
                                                            hard_max_length,
                                                            hash_domain,
                                                            chain_anchor,
                                                            AnchorAxis::X,
                                                            logger);
  const int sweep_x_gap = computeWorstGap(bundles, sweep_x, chain_count);
  const int64_t sweep_x_diam
      = computeWorstManhattanDiameter(bundles, sweep_x, chain_count);

  std::vector<std::size_t> sweep_h = hilbertSweepAssignBundles(bundles,
                                                               chain_count,
                                                               max_length,
                                                               hard_max_length,
                                                               hash_domain,
                                                               chain_anchor,
                                                               logger);
  const int sweep_h_gap = computeWorstGap(bundles, sweep_h, chain_count);
  const int64_t sweep_h_diam
      = computeWorstManhattanDiameter(bundles, sweep_h, chain_count);

  int best_gap = kmeans_gap;
  int64_t best_diam = kmeans_diam;
  std::vector<std::size_t> best_assignment = assignment;
  const char* best_name = "kmeans";
  if (sweep_y_diam < best_diam
      || (sweep_y_diam == best_diam && sweep_y_gap < best_gap)) {
    best_diam = sweep_y_diam;
    best_gap = sweep_y_gap;
    best_assignment = std::move(sweep_y);
    best_name = "sweep_y";
  }
  if (sweep_x_diam < best_diam
      || (sweep_x_diam == best_diam && sweep_x_gap < best_gap)) {
    best_diam = sweep_x_diam;
    best_gap = sweep_x_gap;
    best_assignment = std::move(sweep_x);
    best_name = "sweep_x";
  }
  if (sweep_h_diam < best_diam
      || (sweep_h_diam == best_diam && sweep_h_gap < best_gap)) {
    best_diam = sweep_h_diam;
    best_gap = sweep_h_gap;
    best_assignment = std::move(sweep_h);
    best_name = "hilbert";
  }
  if ((best_diam < kmeans_diam || best_gap < kmeans_gap) && logger) {
    debugPrint(logger,
               utl::DFT,
               "partition",
               1,
               "Using {} partitioning in hash domain {} (worst_diam {} DBU, "
               "worst_gap {} DBU vs kmeans diam {}, gap {}).",
               best_name,
               hash_domain,
               best_diam,
               best_gap,
               kmeans_diam,
               kmeans_gap);
  }
  assignment = std::move(best_assignment);

  mitigateOutlierYGaps(bundles,
                       assignment,
                       chain_count,
                       max_length,
                       hard_max_length,
                       logger,
                       hash_domain,
                       chain_anchor);

  enforceMaxImbalanceOrDie(bundles,
                           assignment,
                           chain_count,
                           max_length,
                           hard_max_length,
                           config,
                           logger,
                           hash_domain,
                           &chain_anchor);

  for (std::size_t i = 0; i < n; ++i) {
    PlacedBundle& b = bundles[i];
    for (auto& cell : b.cells) {
      clustered[assignment[i]].push_back(std::move(cell));
    }
  }

  return clustered;
}

std::vector<std::vector<std::unique_ptr<ScanCell>>> packScanCellsDeterministic(
    std::vector<std::unique_ptr<ScanCell>> domain_cells,
    std::size_t hash_domain,
    const ScanArchitectConfig& config,
    const std::unordered_map<std::string_view, ChainRef>& chain_name_to_ref,
    const std::unordered_map<std::string_view, std::size_t>&
        chain_name_to_ordinal,
    std::size_t chain_count,
    uint64_t max_length,
    bool hard_max_length,
    utl::Logger* logger)
{
  std::vector<PlacedBundle> bundles = buildConstraintBundles(
      std::move(domain_cells), hash_domain, config, chain_name_to_ref, logger);
  std::vector<std::vector<std::unique_ptr<ScanCell>>> chain_cells(chain_count);
  if (chain_count == 0 || bundles.empty()) {
    return chain_cells;
  }

  const std::vector<ChainEndpointPoints> chain_endpoints
      = ResolveChainEndpointsForDomain(
          bundles,
          hash_domain,
          config,
          chain_name_to_ref,
          chain_name_to_ordinal,
          chain_count);
  const auto endpoint_bias = [&](std::size_t bundle_idx,
                                 std::size_t chain_idx) -> int64_t {
    if (bundle_idx >= bundles.size() || chain_idx >= chain_endpoints.size()) {
      return 0;
    }
    int count = 0;
    int64_t sum = 0;
    if (chain_endpoints[chain_idx].begin.has_value()) {
      sum += manhattanDist(bundles[bundle_idx].origin,
                           chain_endpoints[chain_idx].begin.value());
      ++count;
    }
    if (chain_endpoints[chain_idx].end.has_value()) {
      sum += manhattanDist(bundles[bundle_idx].origin,
                           chain_endpoints[chain_idx].end.value());
      ++count;
    }
    if (count == 0) {
      return 0;
    }
    return (sum / count) / 4;  // weight 0.25 vs centroid distance
  };

  std::vector<std::optional<odb::Point>> chain_anchor(chain_count);
  for (std::size_t c = 0; c < chain_count; ++c) {
    const auto& ep = chain_endpoints[c];
    if (ep.begin.has_value() && ep.end.has_value()) {
      const int64_t ax = (static_cast<int64_t>(ep.begin->x())
                          + static_cast<int64_t>(ep.end->x()))
                         / 2;
      const int64_t ay = (static_cast<int64_t>(ep.begin->y())
                          + static_cast<int64_t>(ep.end->y()))
                         / 2;
      chain_anchor[c] = odb::Point(static_cast<int>(ax), static_cast<int>(ay));
    } else if (ep.begin.has_value()) {
      chain_anchor[c] = ep.begin;
    } else if (ep.end.has_value()) {
      chain_anchor[c] = ep.end;
    }
  }

  uint64_t total_bits = 0;
  uint64_t max_item_bits = 0;
  for (const auto& b : bundles) {
    total_bits += b.bits;
    max_item_bits = std::max(max_item_bits, b.bits);
  }

  if (max_length == 0) {
    for (auto& b : bundles) {
      for (auto& cell : b.cells) {
        chain_cells.front().push_back(std::move(cell));
      }
    }
    return chain_cells;
  }

  if (max_item_bits > max_length) {
    if (hard_max_length) {
      logger->error(utl::DFT,
                    171,
                    "Scan architect constraints infeasible in hash domain {}: "
                    "max_length={} is smaller than the largest constrained scan "
                    "item ({} bits).",
                    hash_domain,
                    max_length,
                    max_item_bits);
    }
    max_length = max_item_bits;
  }

  if (config.getChainCount().has_value() && total_bits != 0
      && chain_count > bundles.size()) {
    logger->error(utl::DFT,
                  172,
                  "Scan architect constraints infeasible in hash domain {}: "
                  "chain_count={} exceeds the number of constrained scan items "
                  "({}); cannot create non-empty chains.",
                  hash_domain,
                  chain_count,
                  bundles.size());
  }

  std::vector<uint64_t> used_bits(chain_count, 0);
  std::vector<int64_t> sum_x(chain_count, 0);
  std::vector<int64_t> sum_y(chain_count, 0);
  std::vector<uint64_t> sum_bits(chain_count, 0);
  std::vector<std::size_t> assignment(bundles.size(), 0);

  auto update_centroid = [&](std::size_t c) -> odb::Point {
    if (sum_bits[c] == 0) {
      return odb::Point(0, 0);
    }
    const int64_t cx = sum_x[c] / static_cast<int64_t>(sum_bits[c]);
    const int64_t cy = sum_y[c] / static_cast<int64_t>(sum_bits[c]);
    return odb::Point(static_cast<int>(cx), static_cast<int>(cy));
  };

  std::vector<odb::Point> centroids(chain_count, odb::Point(0, 0));

  // Place fixed bundles first.
  for (std::size_t i = 0; i < bundles.size(); ++i) {
    if (!bundles[i].fixed_chain.has_value()) {
      continue;
    }
    const std::size_t c = bundles[i].fixed_chain.value();
    if (c >= chain_count) {
      logger->error(utl::DFT,
                    173,
                    "Scan constraints internal error: bundle fixed_chain index "
                    "{} is out of range for chain_count {} (hash domain {}).",
                    c,
                    chain_count,
                    hash_domain);
    }
    if (hard_max_length && used_bits[c] + bundles[i].bits > max_length) {
      logger->error(utl::DFT,
                    174,
                    "Scan architect constraints infeasible in hash domain {}: "
                    "fixed assignments overflow max_length={} for chain {}.",
                    hash_domain,
                    max_length,
                    c);
    }
    used_bits[c] += bundles[i].bits;
    assignment[i] = c;
    sum_x[c] += static_cast<int64_t>(bundles[i].origin.x())
                * static_cast<int64_t>(bundles[i].bits);
    sum_y[c] += static_cast<int64_t>(bundles[i].origin.y())
                * static_cast<int64_t>(bundles[i].bits);
    sum_bits[c] += bundles[i].bits;
  }

  for (std::size_t c = 0; c < chain_count; ++c) {
    if (sum_bits[c] == 0) {
      if (chain_anchor[c].has_value()) {
        centroids[c] = chain_anchor[c].value();
      } else {
        centroids[c] = bundles.front().origin;
      }
    } else {
      centroids[c] = update_centroid(c);
    }
  }

  // Pack remaining bundles (best-fit by centroid distance, then used bits).
  std::vector<std::size_t> order(bundles.size());
  std::iota(order.begin(), order.end(), 0);
  std::stable_sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) {
    if (bundles[a].bits != bundles[b].bits) {
      return bundles[a].bits > bundles[b].bits;
    }
    return bundles[a].name < bundles[b].name;
  });

  for (const std::size_t i : order) {
    if (bundles[i].fixed_chain.has_value()) {
      continue;
    }

    bool found = false;
    std::size_t best_chain = 0;
    int64_t best_dist = std::numeric_limits<int64_t>::max();
    uint64_t best_used = std::numeric_limits<uint64_t>::max();

    for (std::size_t c = 0; c < chain_count; ++c) {
      if (hard_max_length && used_bits[c] + bundles[i].bits > max_length) {
        continue;
      }
      const int64_t dist
          = manhattanDist(bundles[i].origin, centroids[c]) + endpoint_bias(i, c);
      if (!found || dist < best_dist
          || (dist == best_dist
              && (used_bits[c] < best_used
                  || (used_bits[c] == best_used && c < best_chain)))) {
        found = true;
        best_chain = c;
        best_dist = dist;
        best_used = used_bits[c];
      }
    }

    if (!found) {
      if (hard_max_length) {
        logger->error(
            utl::DFT,
            175,
            "Scan architect constraints infeasible in hash domain {}: cannot "
            "pack constrained scan items into {} chains with max_length={}.",
            hash_domain,
            chain_count,
            max_length);
      }
      // Soft max_length: fall back to least-overflowing chain.
      best_chain = 0;
      uint64_t best_overflow = std::numeric_limits<uint64_t>::max();
      for (std::size_t c = 0; c < chain_count; ++c) {
        const uint64_t new_bits = used_bits[c] + bundles[i].bits;
        const uint64_t overflow
            = (new_bits > max_length) ? (new_bits - max_length) : 0;
        if (overflow < best_overflow
            || (overflow == best_overflow && c < best_chain)) {
          best_overflow = overflow;
          best_chain = c;
        }
      }
    }

    assignment[i] = best_chain;
    used_bits[best_chain] += bundles[i].bits;
    sum_x[best_chain] += static_cast<int64_t>(bundles[i].origin.x())
                         * static_cast<int64_t>(bundles[i].bits);
    sum_y[best_chain] += static_cast<int64_t>(bundles[i].origin.y())
                         * static_cast<int64_t>(bundles[i].bits);
    sum_bits[best_chain] += bundles[i].bits;
    centroids[best_chain] = update_centroid(best_chain);
  }

  mitigateOutlierYGaps(bundles,
                       assignment,
                       chain_count,
                       max_length,
                       hard_max_length,
                       logger,
                       hash_domain,
                       chain_anchor);

  enforceMaxImbalanceOrDie(bundles,
                           assignment,
                           chain_count,
                           max_length,
                           hard_max_length,
                           config,
                           logger,
                           hash_domain,
                           &chain_anchor);

  for (std::size_t i = 0; i < bundles.size(); ++i) {
    const std::size_t c = assignment[i];
    for (auto& cell : bundles[i].cells) {
      chain_cells[c].push_back(std::move(cell));
    }
  }

  return chain_cells;
}

}  // namespace

ScanArchitectHeuristic::ScanArchitectHeuristic(
    const ScanArchitectConfig& config,
    std::unique_ptr<ScanCellsBucket> scan_cells_bucket,
    utl::Logger* logger)
    : ScanArchitect(config, std::move(scan_cells_bucket)), logger_(logger)
{
}

void ScanArchitectHeuristic::architect()
{
  const bool hard_max_length = config_.getMaxLength().has_value();

  // Build chain name -> (hash_domain,index) map for validating assignments.
  std::unordered_map<std::string_view, ChainRef> chain_name_to_ref;
  std::size_t total_chains = 0;
  for (const auto& [hash_domain, scan_chains] : hash_domain_scan_chains_) {
    (void) hash_domain;
    total_chains += scan_chains.size();
  }
  chain_name_to_ref.reserve(total_chains * 2);
  for (const auto& [hash_domain, scan_chains] : hash_domain_scan_chains_) {
    for (std::size_t i = 0; i < scan_chains.size(); ++i) {
      const std::string_view name = scan_chains[i]->getName();
      chain_name_to_ref.emplace(name, ChainRef{hash_domain, i});
    }
  }

  // Map chain names to the global ordinal used by scan signal name patterns
  // (scan_in_{i}/scan_out_{i}). This matches the ScanStitch ordering when
  // `ScanArchitect::getScanChains()` returns chains in creation order.
  std::unordered_map<std::string_view, std::size_t> chain_name_to_ordinal;
  chain_name_to_ordinal.reserve(total_chains * 2);
  std::size_t ordinal = 0;
  for (const auto& [hash_domain, _] : hash_domain_to_limits_) {
    (void) _;
    const auto it = hash_domain_scan_chains_.find(hash_domain);
    if (it == hash_domain_scan_chains_.end()) {
      continue;
    }
    for (const auto& chain : it->second) {
      if (chain == nullptr) {
        continue;
      }
      chain_name_to_ordinal.emplace(chain->getName(), ordinal);
      ++ordinal;
    }
  }

  // Interpret scanopt_time_limit as a total time budget across all chains
  // (rather than per-chain). Split it evenly here by default to keep runtime
  // roughly stable as the number of chains increases.
  const double scanopt_time_limit_total_s
      = config_.getScanOptTimeLimitSeconds();
  const bool split_scanopt_time_limit
      = config_.getScanOrderSolver()
            == ScanArchitectConfig::ScanOrderSolver::ScanOpt
        && scanopt_time_limit_total_s > 0.0 && total_chains > 1;
  const double scanopt_time_limit_per_chain_s
      = split_scanopt_time_limit
            ? (scanopt_time_limit_total_s / static_cast<double>(total_chains))
            : scanopt_time_limit_total_s;
  if (split_scanopt_time_limit) {
    logger_->info(utl::DFT,
                  202,
                  "ScanOpt time limit {:.1f}s split across {} chains => {:.1f}s "
                  "per chain",
                  scanopt_time_limit_total_s,
                  total_chains,
                  scanopt_time_limit_per_chain_s);
  }

  // Pop all scan cells up-front so we can validate cross-domain constraints.
  std::unordered_map<std::size_t, std::vector<std::unique_ptr<ScanCell>>>
      cells_by_domain;
  cells_by_domain.reserve(hash_domain_scan_chains_.size() * 2);
  std::unordered_map<std::string_view, std::size_t> inst_to_domain;
  inst_to_domain.reserve(1024);

	  for (auto& [hash_domain, scan_chains] : hash_domain_scan_chains_) {
    (void) scan_chains;
    std::vector<std::unique_ptr<ScanCell>> domain_cells;
    domain_cells.reserve(scan_cells_bucket_->numberOfCells(hash_domain));
    while (scan_cells_bucket_->numberOfCells(hash_domain)) {
      std::unique_ptr<ScanCell> cell = scan_cells_bucket_->pop(hash_domain);
      inst_to_domain.emplace(cell->getName(), hash_domain);
      domain_cells.push_back(std::move(cell));
    }
    cells_by_domain.emplace(hash_domain, std::move(domain_cells));
  }

  // Validate assignment and "must be in same chain" constraints early.
  {
    const std::size_t n = inst_to_domain.size();
    std::unordered_map<std::string_view, std::size_t> name_to_idx;
    name_to_idx.reserve(n * 2);
    std::vector<std::string_view> names;
    names.reserve(n);
    for (const auto& [name, _] : inst_to_domain) {
      name_to_idx.emplace(name, names.size());
      names.push_back(name);
    }

    UnionFind uf(n);

    for (const auto& group : config_.getScanOrderGroups()) {
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

    for (const auto& edge : config_.getScanOrderFixedEdges()) {
      const auto it_from = name_to_idx.find(std::string_view(edge.from_inst));
      const auto it_to = name_to_idx.find(std::string_view(edge.to_inst));
      if (it_from == name_to_idx.end() || it_to == name_to_idx.end()) {
        continue;
      }
      uf.unite(it_from->second, it_to->second);
    }

    std::vector<std::optional<std::size_t>> root_domain(n);
    std::vector<std::optional<std::string_view>> root_chain(n);

    for (const auto& [inst, domain] : inst_to_domain) {
      const std::size_t idx = name_to_idx.find(inst)->second;
      const std::size_t root = uf.find(idx);

      // Validate that constrained components do not span hash domains when
      // clock mixing is disabled.
      if (config_.getClockMixing() == ScanArchitectConfig::ClockMixing::NoMix) {
        if (!root_domain[root].has_value()) {
          root_domain[root] = domain;
        } else if (root_domain[root].value() != domain) {
          logger_->error(
              utl::DFT,
              176,
              "Scan constraints infeasible: instances forced into one chain "
              "span multiple hash domains ({} and {}) while clock_mixing=NoMix.",
              root_domain[root].value(),
              domain);
        }
      }

      const std::optional<std::string_view> asn
          = config_.getAssignedChainForInstance(inst);
      if (!asn.has_value()) {
        continue;
      }

      const auto cit = chain_name_to_ref.find(asn.value());
      if (cit == chain_name_to_ref.end()) {
        logger_->error(
            utl::DFT,
            177,
            "Scan constraints: instance '{}' is assigned to unknown chain '{}'.",
            inst,
            asn.value());
      }

      if (config_.getClockMixing() == ScanArchitectConfig::ClockMixing::NoMix
          && cit->second.hash_domain != domain) {
        logger_->error(
            utl::DFT,
            178,
            "Scan constraints infeasible: instance '{}' in hash domain {} is "
            "assigned to chain '{}' in hash domain {} while clock_mixing=NoMix.",
            inst,
            domain,
            asn.value(),
            cit->second.hash_domain);
      }

      if (!root_chain[root].has_value()) {
        root_chain[root] = asn.value();
      } else if (root_chain[root].value() != asn.value()) {
        logger_->error(
            utl::DFT,
            179,
            "Scan constraints infeasible: instances forced into one chain are "
            "assigned to different chains ('{}' vs '{}').",
            root_chain[root].value(),
            asn.value());
      }
    }
  }

  // For each hash_domain, lets distribute the scan cells over the scan chains
  for (auto& [hash_domain, scan_chains] : hash_domain_scan_chains_) {
    const uint64_t max_length
        = hash_domain_to_limits_.find(hash_domain)->second.max_length;

    // Collect all cells for this hash domain so we can do placement-aware
    // partitioning when multiple chains are requested.
    std::vector<std::unique_ptr<ScanCell>> domain_cells
        = std::move(cells_by_domain[hash_domain]);

    const bool do_spatial_partition
        = (scan_chains.size() > 1 && !domain_cells.empty()
           && std::all_of(domain_cells.begin(),
                          domain_cells.end(),
                          [](const std::unique_ptr<ScanCell>& c) {
                            return c->isPlaced();
                          }));

    std::vector<std::vector<std::unique_ptr<ScanCell>>> chain_cells;
    if (do_spatial_partition) {
      // Placement-aware clustering with iterative reassignment (swap/move)
      // before running the intra-chain TSP heuristic.
      chain_cells = clusterPlacedScanCells(
          std::move(domain_cells),
          hash_domain,
          config_,
          chain_name_to_ref,
          chain_name_to_ordinal,
          scan_chains.size(),
          max_length,
          hard_max_length,
          logger_);
    } else {
      chain_cells = packScanCellsDeterministic(std::move(domain_cells),
                                               hash_domain,
                                               config_,
                                               chain_name_to_ref,
                                               chain_name_to_ordinal,
                                               scan_chains.size(),
                                               max_length,
                                               hard_max_length,
                                               logger_);
    }

    for (std::size_t i = 0; i < scan_chains.size(); ++i) {
      for (auto& cell : chain_cells[i]) {
        scan_chains[i]->add(std::move(cell));
      }
    }

	    for (auto& current_chain : scan_chains) {
	      std::optional<ScanArchitectConfig::ChainEndpoints> endpoints
	          = config_.getChainEndpoints(current_chain->getName());
	      // If the user did not provide per-chain endpoints, treat the scan
	      // in/out ports (from name patterns) as implicit endpoints so the
	      // ordering objective includes begin->first and last->end costs.
	      const auto ord_it = chain_name_to_ordinal.find(current_chain->getName());
	      if (ord_it != chain_name_to_ordinal.end()) {
	        ScanArchitectConfig::ChainEndpoints implicit;
	        if (endpoints.has_value()) {
	          implicit = endpoints.value();
	        }
	        bool changed = false;
	        if (!implicit.begin.has_value()
	            && !config_.getScanInNamePattern().empty()) {
	          ScanArchitectConfig::ChainEndpoint begin;
	          begin.type = ScanArchitectConfig::ChainEndpoint::Type::Term;
	          begin.term = fmt::format(
	              FMT_RUNTIME(config_.getScanInNamePattern()), ord_it->second);
	          implicit.begin = std::move(begin);
	          changed = true;
	        }
	        if (!implicit.end.has_value()
	            && !config_.getScanOutNamePattern().empty()) {
	          ScanArchitectConfig::ChainEndpoint end;
	          end.type = ScanArchitectConfig::ChainEndpoint::Type::Term;
	          end.term = fmt::format(
	              FMT_RUNTIME(config_.getScanOutNamePattern()), ord_it->second);
	          implicit.end = std::move(end);
	          changed = true;
	        }
	        if (changed) {
	          endpoints = std::move(implicit);
	        }
	      }
	      const std::string chain_name = std::string(current_chain->getName());
	      current_chain->sortScanCells(
	          [this,
	           endpoints,
	           chain_name,
	           split_scanopt_time_limit,
	           scanopt_time_limit_per_chain_s](
	              std::vector<std::unique_ptr<ScanCell>>& falling,
	              std::vector<std::unique_ptr<ScanCell>>& rising,
	              std::vector<std::unique_ptr<ScanCell>>& sorted) {
	            if (!falling.empty() && !rising.empty()) {
	              std::unordered_set<std::string_view> falling_names;
	              std::unordered_set<std::string_view> rising_names;
	              falling_names.reserve(falling.size() * 2);
	              rising_names.reserve(rising.size() * 2);
              for (const auto& c : falling) {
                falling_names.insert(c->getName());
              }
              for (const auto& c : rising) {
                rising_names.insert(c->getName());
              }

              auto token_has = [&](std::string_view tok) -> std::pair<bool, bool> {
                bool has_falling = falling_names.find(tok) != falling_names.end();
                bool has_rising = rising_names.find(tok) != rising_names.end();
                for (const auto& group : config_.getScanOrderGroups()) {
                  if (group.name.empty() || group.name != tok) {
                    continue;
                  }
                  for (const std::string& inst : group.inst_names) {
                    const std::string_view v(inst);
                    has_falling |= falling_names.find(v) != falling_names.end();
                    has_rising |= rising_names.find(v) != rising_names.end();
                    if (has_falling && has_rising) {
                      break;
                    }
                  }
                  break;
                }
                return {has_falling, has_rising};
              };

              // Fixed-edge constraints crossing polarity boundary are not supported
              // with the current falling-then-rising chain structure.
              for (const auto& edge : config_.getScanOrderFixedEdges()) {
                const std::string_view from(edge.from_inst);
                const std::string_view to(edge.to_inst);
                const bool from_falling = falling_names.find(from) != falling_names.end();
                const bool from_rising = rising_names.find(from) != rising_names.end();
                const bool to_falling = falling_names.find(to) != falling_names.end();
                const bool to_rising = rising_names.find(to) != rising_names.end();
                if (!(from_falling || from_rising) || !(to_falling || to_rising)) {
                  continue;
                }
                if (from_rising && to_falling) {
                  logger_->error(
                      utl::DFT,
                      191,
                      "Scan constraints infeasible in chain '{}': fixed_edge '{}' -> '{}' "
                      "violates polarity (rising before falling).",
                      chain_name,
                      from,
                      to);
                }
                if (from_falling && to_rising) {
                  logger_->error(
                      utl::DFT,
                      192,
                      "Scan constraints infeasible in chain '{}': fixed_edge '{}' -> '{}' "
                      "spans falling->rising boundary, which is not supported with polarity "
                      "ordering.",
                      chain_name,
                      from,
                      to);
                }
              }

              // Polarity hard-constraint check for partial ordering ("before").
              // If any rising instance/group is constrained to appear before any
              // falling instance/group within the same chain, the plan is infeasible.
              for (const auto& bc : config_.getScanOrderBeforeConstraints()) {
                const auto [a_falling, a_rising] = token_has(bc.before);
                const auto [b_falling, b_rising] = token_has(bc.after);
                if (!(a_falling || a_rising) || !(b_falling || b_rising)) {
                  continue;
                }
                if (a_rising && b_falling) {
                  logger_->error(
                      utl::DFT,
                      193,
                      "Scan constraints infeasible in chain '{}': before constraint '{}' -> "
                      "'{}' violates polarity (rising before falling).",
                      chain_name,
                      bc.before,
                      bc.after);
                }
              }
            }

            sorted.reserve(falling.size() + rising.size());

            std::optional<ScanArchitectConfig::ChainEndpoints> falling_eps;
            std::optional<ScanArchitectConfig::ChainEndpoints> rising_eps;
            if (endpoints.has_value()) {
              if (!falling.empty() && !rising.empty()) {
                if (endpoints->begin.has_value()) {
                  falling_eps = ScanArchitectConfig::ChainEndpoints{
                      .begin = endpoints->begin, .end = std::nullopt};
                }
                if (endpoints->end.has_value()) {
                  rising_eps = ScanArchitectConfig::ChainEndpoints{
                      .begin = std::nullopt, .end = endpoints->end};
                }
              } else if (!falling.empty()) {
                falling_eps = endpoints;
              } else if (!rising.empty()) {
                rising_eps = endpoints;
	            }
	            }

	            // Sort to reduce wire length
	            if (split_scanopt_time_limit) {
	              const auto sum_bits
	                  = [](const std::vector<std::unique_ptr<ScanCell>>& v)
	                  -> uint64_t {
	                uint64_t bits = 0;
	                for (const auto& c : v) {
	                  bits += c->getBits();
	                }
	                return bits;
	              };
	              const uint64_t falling_bits = sum_bits(falling);
	              const uint64_t rising_bits = sum_bits(rising);
	              const uint64_t total_bits = falling_bits + rising_bits;

	              double falling_s = scanopt_time_limit_per_chain_s;
	              double rising_s = scanopt_time_limit_per_chain_s;
	              if (total_bits != 0 && !falling.empty() && !rising.empty()) {
	                falling_s
	                    = scanopt_time_limit_per_chain_s
	                      * (static_cast<double>(falling_bits)
	                         / static_cast<double>(total_bits));
	                rising_s = scanopt_time_limit_per_chain_s - falling_s;
	              }

	              ScanArchitectConfig falling_cfg = config_;
	              ScanArchitectConfig rising_cfg = config_;
	              falling_cfg.setScanOptTimeLimitSeconds(falling_s);
	              rising_cfg.setScanOptTimeLimitSeconds(rising_s);

	              OptimizeScanWirelength(
	                  falling, falling_cfg, logger_, falling_eps);
	              OptimizeScanWirelength(rising, rising_cfg, logger_, rising_eps);
	            } else {
	              OptimizeScanWirelength(falling, config_, logger_, falling_eps);
	              OptimizeScanWirelength(rising, config_, logger_, rising_eps);
	            }
	            // Falling edge first
	            std::move(falling.begin(),
	                      falling.end(),
	                      std::back_inserter(sorted));
            std::move(rising.begin(), rising.end(), std::back_inserter(sorted));
          });
	    }
	  }

	  // Global max_imbalance check across all scan chains.
	  {
	    const double allowed_ratio
	        = 1.0
	          + (std::max(0.0, config_.getMaxImbalancePercent()) / 100.0);
	    uint64_t min_bits = std::numeric_limits<uint64_t>::max();
	    uint64_t max_bits = 0;
	    for (const auto& [hash_domain, scan_chains] : hash_domain_scan_chains_) {
	      (void) hash_domain;
	      for (const auto& chain : scan_chains) {
	        const uint64_t bits = chain->getBits();
	        if (bits == 0) {
	          continue;
	        }
	        min_bits = std::min(min_bits, bits);
	        max_bits = std::max(max_bits, bits);
	      }
	    }
	    if (min_bits != std::numeric_limits<uint64_t>::max()) {
	      const double ratio
	          = static_cast<double>(max_bits) / static_cast<double>(min_bits);
	      if (ratio > allowed_ratio + 1e-12) {
	        logger_->error(
	            utl::DFT,
	            197,
	            "Scan architect constraints infeasible: max_imbalance={:.1f}% "
	            "violated across final chains (min_bits={}, max_bits={}, ratio={:.3f}).",
	            config_.getMaxImbalancePercent(),
	            min_bits,
	            max_bits,
	            ratio);
	      }
	    }
	  }
	}

}  // namespace dft
