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
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "Opt.hh"
#include "ScanArchitect.hh"
#include "ScanArchitectConfig.hh"
#include "ScanCell.hh"
#include "odb/geom.h"
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

void enforceMaxImbalanceOrDie(const std::vector<PlacedBundle>& bundles,
                              std::vector<std::size_t>& assignment,
                              std::size_t chain_count,
                              uint64_t max_length,
                              bool hard_max_length,
                              const ScanArchitectConfig& config,
                              utl::Logger* logger,
                              std::size_t hash_domain)
{
  if (chain_count <= 1 || bundles.empty()) {
    return;
  }

  const double allowed_ratio
      = 1.0 + (std::max(0.0, config.getMaxImbalancePercent()) / 100.0);

  auto compute_used_bits = [&]() -> std::vector<uint64_t> {
    std::vector<uint64_t> used(chain_count, 0);
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
    }
    return used;
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

  std::vector<uint64_t> used_bits = compute_used_bits();
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

  const auto move_bundle = [&](std::size_t bundle_idx,
                               std::size_t from_chain,
                               std::size_t to_chain) -> void {
    assignment[bundle_idx] = to_chain;
    used_bits[from_chain] -= bundles[bundle_idx].bits;
    used_bits[to_chain] += bundles[bundle_idx].bits;

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

    bool found = false;
    std::size_t best_bundle = 0;
    double best_ratio = std::numeric_limits<double>::infinity();

    for (const std::size_t bi : chain_to_bundles[max_c]) {
      if (!is_movable(bi)) {
        continue;
      }
      if (chain_to_bundles[max_c].size() <= 1) {
        continue;
      }
      const uint64_t bbits = bundles[bi].bits;
      if (hard_max_length && used_bits[min_c] + bbits > max_length) {
        continue;
      }

      const uint64_t new_max = max_bits - bbits;
      const uint64_t new_min = min_bits + bbits;
      if (new_max == 0) {
        continue;
      }
      const double ratio
          = static_cast<double>(std::max(new_max, new_min))
            / static_cast<double>(std::min(new_max, new_min));
      if (!found || ratio < best_ratio) {
        found = true;
        best_ratio = ratio;
        best_bundle = bi;
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

    move_bundle(best_bundle, max_c, min_c);
  }
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

  std::vector<std::size_t> seed_indices;
  seed_indices.reserve(std::min(unseeded_chains.size(), free_bundles.size()));

  if (!unseeded_chains.empty() && !free_bundles.empty()) {
    // Deterministic centroid seeds among free bundles: pick lower-leftmost,
    // then farthest-from-nearest.
    std::size_t first = free_bundles.front();
    int64_t lowest = std::numeric_limits<int64_t>::max();
    for (const std::size_t idx : free_bundles) {
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

    const std::size_t need = std::min(unseeded_chains.size(), k);
    while (seed_indices.size() < need) {
      std::size_t best = 0;
      bool found = false;
      int64_t best_score = -1;
      for (const std::size_t idx : free_bundles) {
        const bool is_seed = std::find(seed_indices.begin(),
                                       seed_indices.end(),
                                       idx)
                             != seed_indices.end();
        if (is_seed) {
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
    }
  }

  std::vector<std::optional<std::size_t>> seed_chain_for_bundle(n);
  for (std::size_t i = 0; i < seed_indices.size(); ++i) {
    if (i >= unseeded_chains.size()) {
      break;
    }
    seed_chain_for_bundle[seed_indices[i]] = unseeded_chains[i];
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
  // Initialize remaining centroids from seeds when available.
  for (std::size_t i = 0; i < seed_indices.size(); ++i) {
    const std::size_t bidx = seed_indices[i];
    const std::size_t c = seed_chain_for_bundle[bidx].value();
    centroids[c] = bundles[bidx].origin;
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
        const int64_t dist = manhattanDist(bundles[idx].origin, centroids[c]);
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

  enforceMaxImbalanceOrDie(bundles,
                           assignment,
                           chain_count,
                           max_length,
                           hard_max_length,
                           config,
                           logger,
                           hash_domain);

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
    centroids[c] = update_centroid(c);
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
      const int64_t dist = manhattanDist(bundles[i].origin, centroids[c]);
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

  enforceMaxImbalanceOrDie(bundles,
                           assignment,
                           chain_count,
                           max_length,
                           hard_max_length,
                           config,
                           logger,
                           hash_domain);

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
          scan_chains.size(),
          max_length,
          hard_max_length,
          logger_);
    } else {
      chain_cells = packScanCellsDeterministic(std::move(domain_cells),
                                               hash_domain,
                                               config_,
                                               chain_name_to_ref,
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
      const std::optional<ScanArchitectConfig::ChainEndpoints> endpoints
          = config_.getChainEndpoints(current_chain->getName());
      const std::string chain_name = std::string(current_chain->getName());
      current_chain->sortScanCells(
          [this, endpoints, chain_name](std::vector<std::unique_ptr<ScanCell>>& falling,
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
            OptimizeScanWirelength(falling, config_, logger_, falling_eps);
            OptimizeScanWirelength(rising, config_, logger_, rising_eps);
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
