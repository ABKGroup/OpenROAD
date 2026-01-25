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

std::vector<PlacedBundle> buildConstraintBundles(
    std::vector<std::unique_ptr<ScanCell>> domain_cells,
    const ScanArchitectConfig& config,
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

    bundles.push_back(std::move(bundle));
  }

  return bundles;
}

std::vector<std::vector<std::unique_ptr<ScanCell>>> clusterPlacedScanCells(
    std::vector<std::unique_ptr<ScanCell>> domain_cells,
    const ScanArchitectConfig& config,
    std::size_t chain_count,
    uint64_t max_length,
    utl::Logger* logger)
{
  // Backward-compatible wrapper: treat every cell as its own bundle.
  std::vector<PlacedBundle> bundles
      = buildConstraintBundles(std::move(domain_cells), config, logger);

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
    logger->warn(utl::DFT,
                 69,
                 "Scan architect max_length={} is smaller than the largest "
                 "scan item ({} bits); packing will overflow.",
                 max_length,
                 max_item_bits);
  }
  if (total_bits > chain_count * max_length) {
    logger->warn(utl::DFT,
                 70,
                 "Scan architect chain budget is infeasible ({} bits over "
                 "{} chains with max_length={} -> capacity {}); packing will "
                 "overflow.",
                 total_bits,
                 chain_count,
                 max_length,
                 chain_count * max_length);
  }

  const std::size_t n = bundles.size();
  const std::size_t k = std::min(chain_count, n);

  // Deterministic centroid seeds: pick lower-leftmost, then farthest-from-nearest.
  std::vector<std::size_t> seed_indices;
  seed_indices.reserve(k);

  std::size_t first = 0;
  int64_t lowest = std::numeric_limits<int64_t>::max();
  for (std::size_t i = 0; i < n; ++i) {
    const auto& p = bundles[i].origin;
    const int64_t score
        = static_cast<int64_t>(p.x()) + static_cast<int64_t>(p.y());
    if (score < lowest || (score == lowest && bundles[i].name < bundles[first].name)) {
      first = i;
      lowest = score;
    }
  }
  seed_indices.push_back(first);

  while (seed_indices.size() < k) {
    std::size_t best = 0;
    bool found = false;
    int64_t best_score = -1;
    for (std::size_t i = 0; i < n; ++i) {
      const bool is_seed
          = std::find(seed_indices.begin(), seed_indices.end(), i) != seed_indices.end();
      if (is_seed) {
        continue;
      }
      int64_t nearest = std::numeric_limits<int64_t>::max();
      for (const std::size_t s : seed_indices) {
        nearest = std::min(nearest,
                           manhattanDist(bundles[i].origin, bundles[s].origin));
      }
      if (!found || nearest > best_score
          || (nearest == best_score && bundles[i].name < bundles[best].name)) {
        best = i;
        best_score = nearest;
        found = true;
      }
    }
    if (!found) {
      break;
    }
    seed_indices.push_back(best);
  }

  std::vector<odb::Point> centroids(chain_count, bundles[first].origin);
  for (std::size_t i = 0; i < k; ++i) {
    centroids[i] = bundles[seed_indices[i]].origin;
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
    bool overflowed = false;

    // Pin the centroid seed cells to ensure we use all chains when possible.
    for (std::size_t c = 0; c < k; ++c) {
      const std::size_t idx = seed_indices[c];
      assignment[idx] = c;
      is_assigned[idx] = true;
      used_bits[c] += bundles[idx].bits;
      if (used_bits[c] > max_length) {
        overflowed = true;
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
        if (used_bits[c] + bundles[idx].bits > max_length) {
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
        // No chain has remaining capacity; pick the least-overflowing chain.
        overflowed = true;
        best_chain = 0;
        uint64_t best_overflow = std::numeric_limits<uint64_t>::max();
        for (std::size_t c = 0; c < chain_count; ++c) {
          const uint64_t new_bits = used_bits[c] + bundles[idx].bits;
          const uint64_t overflow = (new_bits > max_length) ? (new_bits - max_length) : 0;
          if (overflow < best_overflow || (overflow == best_overflow && c < best_chain)) {
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
      if (overflowed) {
        logger->warn(utl::DFT,
                     71,
                     "Scan architect chain packing overflowed max_length={} "
                     "during clustering.",
                     max_length);
      }
      break;
    }
  }

  for (std::size_t i = 0; i < n; ++i) {
    PlacedBundle& b = bundles[i];
    for (auto& cell : b.cells) {
      clustered[assignment[i]].push_back(std::move(cell));
    }
  }

  return clustered;
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
  // For each hash_domain, lets distribute the scan cells over the scan chains
  for (auto& [hash_domain, scan_chains] : hash_domain_scan_chains_) {
    const uint64_t max_length
        = hash_domain_to_limits_.find(hash_domain)->second.max_length;

    // Collect all cells for this hash domain so we can do placement-aware
    // partitioning when multiple chains are requested.
    std::vector<std::unique_ptr<ScanCell>> domain_cells;
    domain_cells.reserve(scan_cells_bucket_->numberOfCells(hash_domain));
    while (scan_cells_bucket_->numberOfCells(hash_domain)) {
      domain_cells.push_back(scan_cells_bucket_->pop(hash_domain));
    }

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
          std::move(domain_cells), config_, scan_chains.size(), max_length, logger_);
    } else {
      // Fallback: deterministic first-fit packing.
      chain_cells.resize(scan_chains.size());

      std::vector<PlacedBundle> bundles
          = buildConstraintBundles(std::move(domain_cells), config_, logger_);
      std::size_t cursor = 0;
      for (std::size_t c = 0; c < scan_chains.size(); ++c) {
        uint64_t used = 0;
        while (cursor < bundles.size()
               && used + bundles[cursor].bits <= max_length) {
          used += bundles[cursor].bits;
          for (auto& cell : bundles[cursor].cells) {
            chain_cells[c].push_back(std::move(cell));
          }
          ++cursor;
        }
      }

      if (cursor < bundles.size()) {
        logger_->warn(utl::DFT,
                      68,
                      "Scan architect could not pack all scan cells into the "
                      "requested chain budget ({} remaining).",
                      bundles.size() - cursor);
        for (; cursor < bundles.size(); ++cursor) {
          for (auto& cell : bundles[cursor].cells) {
            chain_cells.back().push_back(std::move(cell));
          }
        }
      }
    }

    for (std::size_t i = 0; i < scan_chains.size(); ++i) {
      for (auto& cell : chain_cells[i]) {
        scan_chains[i]->add(std::move(cell));
      }
    }

    for (auto& current_chain : scan_chains) {
      current_chain->sortScanCells(
          [this](std::vector<std::unique_ptr<ScanCell>>& falling,
                 std::vector<std::unique_ptr<ScanCell>>& rising,
                 std::vector<std::unique_ptr<ScanCell>>& sorted) {
            sorted.reserve(falling.size() + rising.size());
            // Sort to reduce wire length
            OptimizeScanWirelength(falling, config_, logger_);
            OptimizeScanWirelength(rising, config_, logger_);
            // Falling edge first
            std::move(falling.begin(),
                      falling.end(),
                      std::back_inserter(sorted));
            std::move(rising.begin(), rising.end(), std::back_inserter(sorted));
          });
    }
  }
}

}  // namespace dft
