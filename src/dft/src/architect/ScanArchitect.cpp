// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2023-2025, The OpenROAD Authors

#include "ScanArchitect.hh"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <limits>
#include <map>
#include <memory>
#include <optional>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "ClockDomain.hh"
#include "ClockDomainHash.hh"
#include "ScanArchitectHeuristic.hh"
#include "utl/Logger.h"

namespace dft {

namespace {

bool CompareScanCells(const std::unique_ptr<ScanCell>& lhs,
                      const std::unique_ptr<ScanCell>& rhs)
{
  // If they have the same number of bits, then we compare the names of the
  // cells so they are ordered by name
  if (lhs->getBits() == rhs->getBits()) {
    const ClockDomain& lhs_clock_domain = lhs->getClockDomain();
    const ClockDomain& rhs_clock_domain = rhs->getClockDomain();

    if (lhs_clock_domain.getClockName() == rhs_clock_domain.getClockName()) {
      return lhs_clock_domain.getClockEdge() < rhs_clock_domain.getClockEdge();
    }

    return lhs_clock_domain.getClockName() < rhs_clock_domain.getClockName();
  }
  // Bigger elements last
  return lhs->getBits() < rhs->getBits();
}

void SortScanCells(std::vector<std::unique_ptr<ScanCell>>& scan_cells)
{
  std::stable_sort(scan_cells.begin(), scan_cells.end(), CompareScanCells);
}

}  // namespace

ScanCellsBucket::ScanCellsBucket(utl::Logger* logger) : logger_(logger)
{
}

namespace {
void showBuckets(
    utl::Logger* logger,
    const std::unordered_map<size_t, std::vector<std::unique_ptr<ScanCell>>>&
        buckets)
{
  for (auto& bucket : buckets) {
    debugPrint(logger, utl::DFT, "buckets", 2, "Bucket {}", bucket.first);
    for (auto& scan_cell : bucket.second) {
      debugPrint(logger, utl::DFT, "buckets", 2, "\t{}", scan_cell->getName());
    }
  }
}
};  // namespace

void ScanCellsBucket::init(const ScanArchitectConfig& config,
                           std::vector<std::unique_ptr<ScanCell>>& scan_cells)
{
  auto hash_fn = GetClockDomainHashFn(config, logger_);
  for (std::unique_ptr<ScanCell>& scan_cell : scan_cells) {
    buckets_[hash_fn(scan_cell->getClockDomain())].push_back(
        std::move(scan_cell));
  }

  // Sort the buckets
  for (auto& [hash_domain, scan_cells] : buckets_) {
    SortScanCells(scan_cells);
  }

  showBuckets(logger_, buckets_);
}

std::unordered_map<size_t, uint64_t>
ScanCellsBucket::getTotalBitsPerHashDomain() const
{
  std::unordered_map<size_t, uint64_t> total_bits;
  for (const auto& [hash_domain, scan_cells] : buckets_) {
    for (const std::unique_ptr<ScanCell>& scan_cell : scan_cells) {
      total_bits[hash_domain] += scan_cell->getBits();
    }
  }
  return total_bits;
}

std::unique_ptr<ScanCell> ScanCellsBucket::pop(size_t hash_domain)
{
  auto& bucket = buckets_.find(hash_domain)->second;
  std::unique_ptr<ScanCell> scan_cell = std::move(bucket.back());
  bucket.pop_back();
  return scan_cell;
}

uint64_t ScanCellsBucket::numberOfCells(size_t hash_domain) const
{
  return buckets_.find(hash_domain)->second.size();
}

std::unique_ptr<ScanArchitect> ScanArchitect::ConstructScanScanArchitect(
    const ScanArchitectConfig& config,
    std::unique_ptr<ScanCellsBucket> scan_cells_bucket,
    utl::Logger* logger)
{
  return std::make_unique<ScanArchitectHeuristic>(
      config, std::move(scan_cells_bucket), logger);
}

ScanArchitect::ScanArchitect(const ScanArchitectConfig& config,
                             std::unique_ptr<ScanCellsBucket> scan_cells_bucket)
    : config_(config), scan_cells_bucket_(std::move(scan_cells_bucket))
{
}

void ScanArchitect::init()
{
  createScanChains();
}

void ScanArchitect::inferChainCount()
{
  std::unordered_map<size_t, uint64_t> hash_domains_total_bits
      = scan_cells_bucket_->getTotalBitsPerHashDomain();

  utl::Logger* logger = scan_cells_bucket_->getLogger();
  const auto ceil_div = [](uint64_t num, uint64_t denom) -> uint64_t {
    return denom == 0 ? 0 : (num + denom - 1) / denom;
  };

  struct DomainInfo
  {
    size_t hash_domain = 0;
    uint64_t bits = 0;
    uint64_t cells = 0;
    uint64_t chains = 0;
  };

  std::vector<DomainInfo> domains;
  domains.reserve(hash_domains_total_bits.size());

  uint64_t total_bits = 0;
  uint64_t total_cells = 0;
  for (const auto& [hash_domain, bits] : hash_domains_total_bits) {
    if (bits == 0) {
      continue;
    }
    DomainInfo d;
    d.hash_domain = hash_domain;
    d.bits = bits;
    d.cells = scan_cells_bucket_->numberOfCells(hash_domain);
    domains.push_back(d);
    total_bits += bits;
    total_cells += d.cells;
  }

  hash_domain_to_limits_.clear();
  if (domains.empty()) {
    return;
  }

  const bool hard_max_length = config_.getMaxLength().has_value();
  const uint64_t hard_max_len
      = hard_max_length ? config_.getMaxLength().value() : 0;
  const uint64_t max_chain_count_limit
      = (config_.getMaxChains().has_value() ? config_.getMaxChains().value()
                                            : 0);
  const double allowed_ratio
      = 1.0 + (std::max(0.0, config_.getMaxImbalancePercent()) / 100.0);
  constexpr double kEps = 1e-12;

  struct RatioStats
  {
    double ratio = std::numeric_limits<double>::infinity();
    uint64_t min_len = 0;
    uint64_t max_len = 0;
  };

  const auto ratio_stats = [&](const std::vector<DomainInfo>& ds) -> RatioStats {
    uint64_t global_min = std::numeric_limits<uint64_t>::max();
    uint64_t global_max = 0;
    for (const DomainInfo& d : ds) {
      if (d.bits == 0) {
        continue;
      }
      if (d.chains == 0) {
        return {};
      }
      const uint64_t min_len = d.bits / d.chains;
      const uint64_t max_len = ceil_div(d.bits, d.chains);
      if (min_len == 0) {
        return {};
      }
      global_min = std::min(global_min, min_len);
      global_max = std::max(global_max, max_len);
    }
    if (global_min == std::numeric_limits<uint64_t>::max() || global_min == 0) {
      return {};
    }
    RatioStats s;
    s.min_len = global_min;
    s.max_len = global_max;
    s.ratio = static_cast<double>(global_max) / static_cast<double>(global_min);
    return s;
  };

  const auto choose_best_increment
      = [&](std::vector<DomainInfo>& ds) -> std::optional<std::size_t> {
    std::optional<std::size_t> best_idx;
    RatioStats best_stats;
    for (std::size_t i = 0; i < ds.size(); ++i) {
      if (ds[i].chains >= ds[i].cells) {
        continue;
      }
      ds[i].chains++;
      const RatioStats s = ratio_stats(ds);
      ds[i].chains--;

      if (!best_idx.has_value() || s.ratio + kEps < best_stats.ratio
          || (std::abs(s.ratio - best_stats.ratio) <= kEps
              && (s.max_len < best_stats.max_len
                  || (s.max_len == best_stats.max_len
                      && s.min_len > best_stats.min_len)))) {
        best_idx = i;
        best_stats = s;
      }
    }
    return best_idx;
  };

  // 1) Exact chain_count (total across design).
  if (auto chain_count = config_.getChainCount(); chain_count.has_value()) {
    const uint64_t target = chain_count.value();
    if (target == 0) {
      logger->error(utl::DFT, 73, "Scan architect chain_count cannot be 0.");
    }
    if (max_chain_count_limit != 0 && target > max_chain_count_limit) {
      logger->error(utl::DFT,
                    195,
                    "Scan architect constraints infeasible: chain_count={} "
                    "exceeds max_chains={}.",
                    target,
                    max_chain_count_limit);
    }
    if (target > total_cells) {
      logger->error(utl::DFT,
                    74,
                    "Scan architect constraints infeasible: chain_count={} "
                    "exceeds total_scan_cells={} (would create empty chains).",
                    target,
                    total_cells);
    }
    if (hard_max_length) {
      const uint64_t capacity = target * hard_max_len;
      if (capacity != 0 && total_bits > capacity) {
        const uint64_t required_max_len = ceil_div(total_bits, target);
        logger->error(
            utl::DFT,
            129,
            "Scan architect constraints infeasible: total_bits={} exceeds "
            "chain_count*max_length={}*{} (capacity {}); requires per-chain "
            "length at least {}.",
            total_bits,
            target,
            hard_max_len,
            capacity,
            required_max_len);
      }
    }

    uint64_t required_total = 0;
    for (DomainInfo& d : domains) {
      const uint64_t min_chains = hard_max_length ? ceil_div(d.bits, hard_max_len)
                                                  : 1;
      if (min_chains == 0) {
        d.chains = 1;
      } else {
        d.chains = min_chains;
      }
      if (d.chains > d.cells) {
        logger->error(
            utl::DFT,
            199,
            "Scan architect constraints infeasible: chain_count={} requires at "
            "least {} chain(s) in hash domain {} but it has only {} scan cell(s).",
            target,
            d.chains,
            d.hash_domain,
            d.cells);
      }
      required_total += d.chains;
    }

    if (target < required_total) {
      logger->error(utl::DFT,
                    196,
                    "Scan architect constraints infeasible: chain_count={} is "
                    "less than the minimum required {} chain(s).",
                    target,
                    required_total);
    }

    uint64_t total_chains = required_total;
    while (total_chains < target) {
      const std::optional<std::size_t> best = choose_best_increment(domains);
      if (!best.has_value()) {
        logger->error(
            utl::DFT,
            200,
            "Scan architect constraints infeasible: chain_count={} cannot be "
            "satisfied because there are only {} scan cell(s) total.",
            target,
            total_cells);
      }
      domains[best.value()].chains++;
      ++total_chains;
    }

    const RatioStats s = ratio_stats(domains);
    if (s.ratio > allowed_ratio + kEps) {
      logger->error(utl::DFT,
                    128,
                    "Scan architect constraints infeasible: chain_count={} "
                    "cannot satisfy max_imbalance={:.1f}%.",
                    target,
                    config_.getMaxImbalancePercent());
    }

    for (const DomainInfo& d : domains) {
      HashDomainLimits limits;
      limits.chain_count = d.chains;
      limits.max_length = hard_max_length ? hard_max_len : ceil_div(d.bits, d.chains);
      hash_domain_to_limits_.insert({d.hash_domain, limits});
    }
    return;
  }

  // 2) Infer chain count.
  //
  // If max_length is set, it is a hard per-chain cap.
  // If max_length is not set, we use a soft default of 200 to avoid producing
  // very long scan chains by default.
  const uint64_t soft_target_len = 200;
  uint64_t required_total = 0;
  for (DomainInfo& d : domains) {
    uint64_t c = 1;
    if (hard_max_length) {
      c = ceil_div(d.bits, hard_max_len);
    } else if (max_chain_count_limit != 0) {
      c = 1;
    } else {
      c = std::max<uint64_t>(1, ceil_div(d.bits, soft_target_len));
    }
    c = std::max<uint64_t>(1, c);
    if (c > d.cells) {
      if (hard_max_length) {
        logger->error(
            utl::DFT,
            201,
            "Scan architect constraints infeasible: max_length={} requires at "
            "least {} chain(s) in hash domain {} but it has only {} scan cell(s).",
            hard_max_len,
            c,
            d.hash_domain,
            d.cells);
      }
      c = d.cells;
    }
    d.chains = c;
    required_total += d.chains;
  }

  if (max_chain_count_limit != 0 && required_total > max_chain_count_limit) {
    if (hard_max_length) {
      logger->error(utl::DFT,
                    130,
                    "Scan architect constraints infeasible: max_length={} "
                    "requires at least {} chains for total_bits={}, but "
                    "max_chains={}.",
                    hard_max_len,
                    required_total,
                    total_bits,
                    max_chain_count_limit);
    } else {
      logger->error(utl::DFT,
                    198,
                    "Scan architect constraints infeasible: max_chains={} is "
                    "less than the minimum required {} chain(s) (one per hash "
                    "domain).",
                    max_chain_count_limit,
                    required_total);
    }
  }

  uint64_t total_chains = required_total;
  for (;;) {
    const RatioStats s = ratio_stats(domains);
    if (s.ratio <= allowed_ratio + kEps) {
      break;
    }
    if (max_chain_count_limit != 0 && total_chains >= max_chain_count_limit) {
      logger->error(
          utl::DFT,
          132,
          "Scan architect constraints infeasible: cannot satisfy "
          "max_imbalance={:.1f}% with max_length={} and max_chains={}.",
          config_.getMaxImbalancePercent(),
          (hard_max_length ? hard_max_len : 0),
          max_chain_count_limit);
    }
    const std::optional<std::size_t> best = choose_best_increment(domains);
    if (!best.has_value()) {
      logger->error(utl::DFT,
                    203,
                    "Scan architect constraints infeasible: cannot satisfy "
                    "max_imbalance={:.1f}% due to scan chain assignment limits.",
                    config_.getMaxImbalancePercent());
    }
    domains[best.value()].chains++;
    ++total_chains;
  }

  for (const DomainInfo& d : domains) {
    HashDomainLimits limits;
    limits.chain_count = d.chains;
    limits.max_length = hard_max_length ? hard_max_len : ceil_div(d.bits, d.chains);
    hash_domain_to_limits_.insert({d.hash_domain, limits});
  }
}

std::map<size_t, ScanArchitect::HashDomainLimits>
ScanArchitect::inferChainCountFromMaxLength(
    const std::unordered_map<size_t, uint64_t>& hash_domains_total_bit,
    uint64_t max_length,
    const std::optional<uint64_t>& max_chains)
{
  std::map<size_t, HashDomainLimits> hash_domain_to_limits;
  for (const auto& [hash_domain, bits] : hash_domains_total_bit) {
    if (bits == 0) {
      continue;
    }
    uint64_t domain_chain_count = 0;
    if (bits % max_length != 0) {
      // For unbalance case, we need +1 chains to hold the aditionals cells
      domain_chain_count = bits / max_length + 1;
    } else {
      domain_chain_count = bits / max_length;
    }

    if (max_chains.has_value()) {
      domain_chain_count = std::min(domain_chain_count, max_chains.value());
    }
    domain_chain_count = std::max<uint64_t>(1, domain_chain_count);

    uint64_t domain_max_length = 0;
    if (bits % domain_chain_count == 0) {
      domain_max_length = bits / domain_chain_count;
    } else {
      domain_max_length = bits / domain_chain_count + 1;
    }

    HashDomainLimits hash_domain_limits;
    hash_domain_limits.chain_count = domain_chain_count;
    hash_domain_limits.max_length = domain_max_length;

    hash_domain_to_limits.insert({hash_domain, hash_domain_limits});
  }
  return hash_domain_to_limits;
}

void ScanArchitect::createScanChains()
{
  inferChainCount();
  uint64_t total_chains = 0;
  for (const auto& [hash_domain, limits] : hash_domain_to_limits_) {
    (void) hash_domain;
    total_chains += limits.chain_count;
  }

  const std::vector<std::string>& user_names = config_.getChainNames();
	if (!user_names.empty() && user_names.size() != total_chains) {
	  utl::Logger* logger = scan_cells_bucket_->getLogger();
	  logger->error(utl::DFT,
	                133,
	                "Scan constraints specify {} chain name(s), but the DFT plan "
	                "requires {} chain(s).",
	                user_names.size(),
	                total_chains);
  }
  if (!user_names.empty()) {
    std::unordered_set<std::string_view> seen;
    seen.reserve(user_names.size() * 2);
    for (const std::string& name : user_names) {
	      if (!seen.insert(name).second) {
	        utl::Logger* logger = scan_cells_bucket_->getLogger();
	        logger->error(utl::DFT,
	                      134,
	                      "Scan constraints contain duplicate chain name '{}'.",
	                      name);
	      }
	    }
  }

  uint64_t chain_number = 0;
  for (const auto& [hash_domain, limits] : hash_domain_to_limits_) {
    for (uint64_t i = 0; i < limits.chain_count; ++i) {
      std::string chain_name;
      if (!user_names.empty()) {
        chain_name = user_names[chain_number];
      } else {
        chain_name = fmt::format("chain_{}", chain_number);
      }
      hash_domain_scan_chains_[hash_domain].push_back(
          std::make_unique<ScanChain>(chain_name));
      ++chain_number;
    }
  }
}

std::vector<std::unique_ptr<ScanChain>> ScanArchitect::getScanChains()
{
  std::vector<std::unique_ptr<ScanChain>> scan_chains_flat;
  std::size_t total = 0;
  for (const auto& [hash_domain, scan_chains] : hash_domain_scan_chains_) {
    (void) hash_domain;
    total += scan_chains.size();
  }
  scan_chains_flat.reserve(total);

  // Preserve deterministic chain numbering / user-specified name ordering by
  // iterating hash domains in creation order (hash_domain_to_limits_ is a
  // sorted map). This also avoids lexicographic issues like "chain_10" sorting
  // before "chain_2" which can misalign scan_in_N/scan_out_N port ordinals.
  for (const auto& [hash_domain, _] : hash_domain_to_limits_) {
    (void) _;
    auto it = hash_domain_scan_chains_.find(hash_domain);
    if (it == hash_domain_scan_chains_.end()) {
      continue;
    }
    auto& scan_chains = it->second;
    std::move(std::begin(scan_chains),
              std::end(scan_chains),
              std::back_inserter(scan_chains_flat));
    scan_chains.clear();
  }

  // Defensive: move any remaining chains from unexpected domains.
  for (auto& [hash_domain, scan_chains] : hash_domain_scan_chains_) {
    (void) hash_domain;
    if (scan_chains.empty()) {
      continue;
    }
    std::move(std::begin(scan_chains),
              std::end(scan_chains),
              std::back_inserter(scan_chains_flat));
    scan_chains.clear();
  }

  return scan_chains_flat;
}
}  // namespace dft
