// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2023-2025, The OpenROAD Authors

#include "ScanArchitectConfig.hh"

#include <algorithm>
#include <cctype>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <optional>
#include <sstream>
#include <string>
#include <string_view>
#include <unordered_map>
#include <unordered_set>

#include "Formatting.hh"
#include "utl/Logger.h"

namespace dft {

namespace {

bool ParseInt(const std::string& token, int& value)
{
  if (token.empty()) {
    return false;
  }
  char* end = nullptr;
  const long v = std::strtol(token.c_str(), &end, 10);
  if (end == nullptr || *end != '\0') {
    return false;
  }
  if (v < std::numeric_limits<int>::min()
      || v > std::numeric_limits<int>::max()) {
    return false;
  }
  value = static_cast<int>(v);
  return true;
}

[[maybe_unused]] bool ParsePointTokens(std::istringstream& iss,
                                       ScanArchitectConfig::Point& point)
{
  std::string sx;
  std::string sy;
  if (!(iss >> sx >> sy)) {
    return false;
  }
  int x = 0;
  int y = 0;
  if (!ParseInt(sx, x) || !ParseInt(sy, y)) {
    return false;
  }
  point.x = x;
  point.y = y;
  return true;
}

bool ParseEndpointTokens(std::istringstream& iss,
                         ScanArchitectConfig::ChainEndpoint& endpoint)
{
  std::string first;
  if (!(iss >> first)) {
    return false;
  }

  int x = 0;
  int y = 0;
  if (ParseInt(first, x)) {
    std::string second;
    if (!(iss >> second)) {
      return false;
    }
    if (!ParseInt(second, y)) {
      return false;
    }
    endpoint.type = ScanArchitectConfig::ChainEndpoint::Type::Point;
    endpoint.point.x = x;
    endpoint.point.y = y;
    endpoint.term.clear();
    return true;
  }

  endpoint.type = ScanArchitectConfig::ChainEndpoint::Type::Term;
  endpoint.term = std::move(first);
  endpoint.point = ScanArchitectConfig::Point{};
  return true;
}

bool GlobMatch(std::string_view pattern, std::string_view text)
{
  // Simple glob matcher supporting:
  // - '*' matches any substring (including empty)
  // - '?' matches any single character
  // - '\\' escapes the next character
  std::size_t p = 0;
  std::size_t t = 0;
  std::size_t star = std::string_view::npos;
  std::size_t match = 0;

  auto next_pat = [&](std::size_t idx,
                      char& out,
                      bool& escaped) -> std::size_t {
    escaped = false;
    if (idx >= pattern.size()) {
      out = '\0';
      return idx;
    }
    if (pattern[idx] == '\\' && idx + 1 < pattern.size()) {
      escaped = true;
      out = pattern[idx + 1];
      return idx + 2;
    }
    out = pattern[idx];
    return idx + 1;
  };

  while (t < text.size()) {
    char pc = '\0';
    bool esc = false;
    const std::size_t p_next = next_pat(p, pc, esc);

    if (!esc && p < pattern.size() && pc == '*') {
      star = p;
      match = t;
      p = p_next;
      continue;
    }

    if (p < pattern.size()
        && (( !esc && pc == '?') || (pc == text[t]))) {
      p = p_next;
      ++t;
      continue;
    }

    if (star != std::string_view::npos) {
      // Backtrack: extend the '*' match by one character.
      char star_c = '\0';
      bool star_esc = false;
      (void) next_pat(star, star_c, star_esc);
      p = star + 1;
      t = ++match;
      continue;
    }

    return false;
  }

  // Consume trailing '*' in pattern.
  while (p < pattern.size()) {
    char pc = '\0';
    bool esc = false;
    const std::size_t p_next = next_pat(p, pc, esc);
    if (!esc && pc == '*') {
      p = p_next;
      continue;
    }
    break;
  }

  return p >= pattern.size();
}

bool MatchesAny(const std::vector<std::string>& patterns, std::string_view text)
{
  for (const std::string& pat : patterns) {
    if (GlobMatch(pat, text)) {
      return true;
    }
  }
  return false;
}

}  // namespace

ScanArchitectConfig::ClockMixing ScanArchitectConfig::getClockMixing() const
{
  return clock_mixing_;
}

ScanArchitectConfig::PolarityMode ScanArchitectConfig::getPolarityMode() const
{
  return polarity_mode_;
}

ScanArchitectConfig::ScanOrderMetric ScanArchitectConfig::getScanOrderMetric()
    const
{
  return scan_order_metric_;
}

ScanArchitectConfig::ScanOrderSolver ScanArchitectConfig::getScanOrderSolver()
    const
{
  return scan_order_solver_;
}

std::string_view ScanArchitectConfig::getScanEnableNamePattern() const
{
  return scan_enable_name_pattern_;
}

std::string_view ScanArchitectConfig::getScanInNamePattern() const
{
  return scan_in_name_pattern_;
}

std::string_view ScanArchitectConfig::getScanOutNamePattern() const
{
  return scan_out_name_pattern_;
}

uint64_t ScanArchitectConfig::getScanOptRounds() const
{
  return scanopt_rounds_;
}

uint64_t ScanArchitectConfig::getScanOptSeed() const
{
  return scanopt_seed_;
}

double ScanArchitectConfig::getScanOptTimeLimitSeconds() const
{
  return scanopt_time_limit_seconds_;
}

bool ScanArchitectConfig::getScanOptTempControl() const
{
  return scanopt_temp_control_;
}

double ScanArchitectConfig::getScanOptTDiv() const
{
  return scanopt_t_div_;
}

double ScanArchitectConfig::getVerticalWeight() const
{
  return vertical_weight_;
}

double ScanArchitectConfig::getTimingWeightSetup() const
{
  return timing_weight_setup_;
}

double ScanArchitectConfig::getTimingWeightHold() const
{
  return timing_weight_hold_;
}

double ScanArchitectConfig::getTimingCriticalSlack() const
{
  return timing_critical_slack_;
}

const std::optional<uint64_t>& ScanArchitectConfig::getChainCount() const
{
  return chain_count_;
}

const std::optional<uint64_t>& ScanArchitectConfig::getMaxLength() const
{
  return max_length_;
}

const std::optional<uint64_t>& ScanArchitectConfig::getMaxChains() const
{
  return max_chains_;
}

void ScanArchitectConfig::setMaxImbalancePercent(double percent)
{
  if (percent >= 0.0) {
    max_imbalance_percent_ = percent;
  }
}

double ScanArchitectConfig::getMaxImbalancePercent() const
{
  return max_imbalance_percent_;
}

void ScanArchitectConfig::setClockMixing(ClockMixing clock_mixing)
{
  clock_mixing_ = clock_mixing;
}

void ScanArchitectConfig::setPolarityMode(PolarityMode mode)
{
  polarity_mode_ = mode;
}

void ScanArchitectConfig::setScanOrderMetric(ScanOrderMetric metric)
{
  scan_order_metric_ = metric;
}

void ScanArchitectConfig::setScanOrderSolver(ScanOrderSolver solver)
{
  scan_order_solver_ = solver;
}

void ScanArchitectConfig::setScanEnableNamePattern(std::string_view pattern)
{
  scan_enable_name_pattern_ = std::string(pattern);
}

void ScanArchitectConfig::setScanInNamePattern(std::string_view pattern)
{
  scan_in_name_pattern_ = std::string(pattern);
}

void ScanArchitectConfig::setScanOutNamePattern(std::string_view pattern)
{
  scan_out_name_pattern_ = std::string(pattern);
}

void ScanArchitectConfig::setScanOptRounds(uint64_t rounds)
{
  scanopt_rounds_ = rounds;
}

void ScanArchitectConfig::setScanOptSeed(uint64_t seed)
{
  scanopt_seed_ = seed;
}

void ScanArchitectConfig::setScanOptTimeLimitSeconds(double seconds)
{
  if (seconds >= 0.0) {
    scanopt_time_limit_seconds_ = seconds;
  }
}

void ScanArchitectConfig::setScanOptTempControl(bool enable)
{
  scanopt_temp_control_ = enable;
}

void ScanArchitectConfig::setScanOptTDiv(double t_div)
{
  if (t_div > 0.0) {
    scanopt_t_div_ = t_div;
  }
}

void ScanArchitectConfig::setVerticalWeight(double weight)
{
  if (weight > 0.0) {
    vertical_weight_ = weight;
  }
}

void ScanArchitectConfig::setTimingWeightSetup(double weight)
{
  if (weight >= 0.0) {
    timing_weight_setup_ = weight;
  }
}

void ScanArchitectConfig::setTimingWeightHold(double weight)
{
  if (weight >= 0.0) {
    timing_weight_hold_ = weight;
  }
}

void ScanArchitectConfig::setTimingCriticalSlack(double slack)
{
  if (slack >= 0.0) {
    timing_critical_slack_ = slack;
  }
}

void ScanArchitectConfig::setAutoExcludeShiftRegisters(bool enable)
{
  auto_exclude_shift_registers_ = enable;
}

bool ScanArchitectConfig::getAutoExcludeShiftRegisters() const
{
  return auto_exclude_shift_registers_;
}

void ScanArchitectConfig::setPreferQbarScanOut(bool enable)
{
  prefer_qbar_scan_out_ = enable;
}

bool ScanArchitectConfig::getPreferQbarScanOut() const
{
  return prefer_qbar_scan_out_;
}

void ScanArchitectConfig::setShiftRegisterMinLength(int min_length)
{
  // Keep this conservative; lengths < 2 do not represent a chain.
  shift_register_min_length_ = std::max(2, min_length);
}

int ScanArchitectConfig::getShiftRegisterMinLength() const
{
  return shift_register_min_length_;
}

void ScanArchitectConfig::clearAutoExcludedInstances()
{
  auto_excluded_instances_.clear();
}

void ScanArchitectConfig::setAutoExcludedInstances(
    std::unordered_set<std::string> instances)
{
  auto_excluded_instances_ = std::move(instances);
}

void ScanArchitectConfig::setChainCount(uint64_t chain_count)
{
  chain_count_ = chain_count;
}

void ScanArchitectConfig::setMaxChains(uint64_t max_chains)
{
  max_chains_ = max_chains;
}

void ScanArchitectConfig::setMaxLength(uint64_t max_length)
{
  max_length_ = max_length;
}

void ScanArchitectConfig::clearScanOrderConstraints()
{
  scan_order_groups_.clear();
  scan_order_fixed_edges_.clear();
  scan_order_before_.clear();
  chain_names_.clear();
  chain_endpoints_by_name_.clear();
  instance_to_chain_name_.clear();
  excluded_instances_.clear();
  excluded_instance_patterns_.clear();
  excluded_master_patterns_.clear();
}

void ScanArchitectConfig::setDefaultGroupPriority(int priority)
{
  default_group_priority_ = std::clamp(priority, 0, 127);
}

int ScanArchitectConfig::getDefaultGroupPriority() const
{
  return default_group_priority_;
}

const std::vector<ScanArchitectConfig::ScanOrderGroupConstraint>&
ScanArchitectConfig::getScanOrderGroups() const
{
  return scan_order_groups_;
}

const std::vector<ScanArchitectConfig::ScanOrderFixedEdgeConstraint>&
ScanArchitectConfig::getScanOrderFixedEdges() const
{
  return scan_order_fixed_edges_;
}

const std::vector<ScanArchitectConfig::ScanOrderBeforeConstraint>&
ScanArchitectConfig::getScanOrderBeforeConstraints() const
{
  return scan_order_before_;
}

std::optional<ScanArchitectConfig::ChainEndpoints> ScanArchitectConfig::getChainEndpoints(
    std::string_view chain_name) const
{
  const auto it = chain_endpoints_by_name_.find(std::string(chain_name));
  if (it == chain_endpoints_by_name_.end()) {
    return std::nullopt;
  }
  return it->second;
}

const std::vector<std::string>& ScanArchitectConfig::getChainNames() const
{
  return chain_names_;
}

std::optional<std::string_view> ScanArchitectConfig::getAssignedChainForInstance(
    std::string_view inst_name) const
{
  const auto it = instance_to_chain_name_.find(std::string(inst_name));
  if (it == instance_to_chain_name_.end()) {
    return std::nullopt;
  }
  return std::string_view(it->second);
}

bool ScanArchitectConfig::isInstanceExcluded(std::string_view inst_name) const
{
  if (auto_excluded_instances_.find(std::string(inst_name))
      != auto_excluded_instances_.end()) {
    return true;
  }
  if (excluded_instances_.find(std::string(inst_name)) != excluded_instances_.end()) {
    return true;
  }
  return MatchesAny(excluded_instance_patterns_, inst_name);
}

bool ScanArchitectConfig::isMasterExcluded(std::string_view master_name) const
{
  return MatchesAny(excluded_master_patterns_, master_name);
}

bool ScanArchitectConfig::isInstanceExcluded(std::string_view inst_name,
                                             std::string_view master_name) const
{
  return isInstanceExcluded(inst_name) || isMasterExcluded(master_name);
}

const std::unordered_set<std::string>& ScanArchitectConfig::getExcludedInstances() const
{
  return excluded_instances_;
}

const std::unordered_set<std::string>& ScanArchitectConfig::getAutoExcludedInstances() const
{
  return auto_excluded_instances_;
}

const std::vector<std::string>& ScanArchitectConfig::getExcludedInstancePatterns()
    const
{
  return excluded_instance_patterns_;
}

const std::vector<std::string>& ScanArchitectConfig::getExcludedMasterPatterns()
    const
{
  return excluded_master_patterns_;
}

bool ScanArchitectConfig::loadScanOrderConstraintsFile(
    const std::string& path,
    utl::Logger* logger)
{
  std::ifstream f(path);
  if (!f) {
    if (logger) {
      logger->error(utl::DFT, 135, "Couldn't open scan constraints file '{}'", path);
    }
    return false;
  }

  clearScanOrderConstraints();

  struct RawGroup
  {
    std::string name;
    int priority = 64;
    bool is_path = false;
    std::vector<std::string> members;
  };

  struct RawAssignment
  {
    std::string chain_name;
    std::vector<std::string> items;
  };

  std::vector<std::string> raw_excludes;
  std::vector<std::string> raw_exclude_instance_patterns;
  std::vector<std::string> raw_exclude_master_patterns;

  std::unordered_map<std::string, RawGroup> raw_groups;
  std::vector<std::string> raw_group_order;
  std::vector<RawAssignment> raw_assignments;

  const auto trim = [](std::string& s) {
    auto not_space = [](unsigned char c) { return !std::isspace(c); };
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), not_space));
    s.erase(std::find_if(s.rbegin(), s.rend(), not_space).base(), s.end());
  };

  std::string raw;
  int lineno = 0;
  while (std::getline(f, raw)) {
    ++lineno;
    std::string line = raw;
    if (const std::size_t hash = line.find('#'); hash != std::string::npos) {
      line = line.substr(0, hash);
    }
    trim(line);
    if (line.empty()) {
      continue;
    }

    std::istringstream iss(line);
    std::string keyword;
    iss >> keyword;
    std::transform(keyword.begin(),
                   keyword.end(),
                   keyword.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });

    if (keyword == "default_priority" || keyword == "default_group_priority") {
      int pr = 0;
      if (!(iss >> pr)) {
        if (logger) {
          logger->warn(utl::DFT,
                       136,
                       "Scan constraints parse error at {}:{}: expected integer "
                       "priority",
                       path,
                       lineno);
        }
        continue;
      }
      setDefaultGroupPriority(pr);
      continue;
    }

    if (keyword == "chain") {
      std::string name;
      if (!(iss >> name)) {
        if (logger) {
          logger->warn(utl::DFT,
                       137,
                       "Scan constraints parse error at {}:{}: expected 'chain "
                       "<name> [begin <x> <y>] [end <x> <y>]'",
                       path,
                       lineno);
        }
        continue;
      }

      if (std::find(chain_names_.begin(), chain_names_.end(), name)
          != chain_names_.end()) {
        if (logger) {
          logger->error(utl::DFT,
                        138,
                        "Scan constraints parse error at {}:{}: duplicate chain "
                        "name '{}' (file '{}')",
                        path,
                        lineno,
                        name,
                        path);
        }
        continue;
      }
      chain_names_.push_back(name);

      ChainEndpoints& endpoints = chain_endpoints_by_name_[name];
      std::string tok;
      while (iss >> tok) {
        std::string k = tok;
        std::transform(
            k.begin(), k.end(), k.begin(), [](unsigned char c) {
              return static_cast<char>(std::tolower(c));
            });

        if (k == "begin" || k == "begin_port" || k == "beginport") {
          ChainEndpoint ep;
          if (!ParseEndpointTokens(iss, ep)) {
            if (logger) {
              logger->warn(utl::DFT,
                           139,
                           "Scan constraints parse error at {}:{}: expected "
                           "'begin <x> <y>' or 'begin <port|inst/pin>'",
                           path,
                           lineno);
            }
            break;
          }
          endpoints.begin = std::move(ep);
          continue;
        }
        if (k == "end" || k == "end_port" || k == "endport") {
          ChainEndpoint ep;
          if (!ParseEndpointTokens(iss, ep)) {
            if (logger) {
              logger->warn(utl::DFT,
                           140,
                           "Scan constraints parse error at {}:{}: expected "
                           "'end <x> <y>' or 'end <port|inst/pin>'",
                           path,
                           lineno);
            }
            break;
          }
          endpoints.end = std::move(ep);
          continue;
        }

        if (logger) {
          logger->warn(utl::DFT,
                       141,
                       "Scan constraints: ignoring unknown chain token '{}' at "
                       "{}:{}",
                       tok,
                       path,
                       lineno);
        }
      }
      continue;
    }

    if (keyword == "chain_begin" || keyword == "chainbegin") {
      std::string name;
      if (!(iss >> name)) {
        if (logger) {
          logger->warn(utl::DFT,
                       142,
                       "Scan constraints parse error at {}:{}: expected "
                       "'chain_begin <name> <x> <y>' or 'chain_begin <name> "
                       "<port|inst/pin>'",
                       path,
                       lineno);
        }
        continue;
      }
      ChainEndpoint ep;
      if (!ParseEndpointTokens(iss, ep)) {
        if (logger) {
          logger->warn(utl::DFT,
                       143,
                       "Scan constraints parse error at {}:{}: expected "
                       "'chain_begin <name> <x> <y>' or 'chain_begin <name> "
                       "<port|inst/pin>'",
                       path,
                       lineno);
        }
        continue;
      }
      chain_endpoints_by_name_[name].begin = std::move(ep);
      continue;
    }

    if (keyword == "chain_end" || keyword == "chainend") {
      std::string name;
      if (!(iss >> name)) {
        if (logger) {
          logger->warn(utl::DFT,
                       144,
                       "Scan constraints parse error at {}:{}: expected "
                       "'chain_end <name> <x> <y>' or 'chain_end <name> "
                       "<port|inst/pin>'",
                       path,
                       lineno);
        }
        continue;
      }
      ChainEndpoint ep;
      if (!ParseEndpointTokens(iss, ep)) {
        if (logger) {
          logger->warn(utl::DFT,
                       145,
                       "Scan constraints parse error at {}:{}: expected "
                       "'chain_end <name> <x> <y>' or 'chain_end <name> "
                       "<port|inst/pin>'",
                       path,
                       lineno);
        }
        continue;
      }
      chain_endpoints_by_name_[name].end = std::move(ep);
      continue;
    }

    if (keyword == "group") {
      std::string first;
      if (!(iss >> first)) {
        if (logger) {
          logger->warn(
              utl::DFT,
              146,
              "Scan constraints parse error at {}:{}: expected 'group "
              "[<name>] [<priority>] <inst/group...>'",
              path,
              lineno);
        }
        continue;
      }

      int pr = default_group_priority_;
      std::string group_name;
      std::vector<std::string> members;

      int parsed_pr = 0;
      if (ParseInt(first, parsed_pr)) {
        pr = std::clamp(parsed_pr, 0, 127);
        group_name = fmt::format("__group_{}", lineno);
      } else {
        group_name = first;
        std::string second;
        if (!(iss >> second)) {
          if (logger) {
            logger->warn(
                utl::DFT,
                147,
                "Scan constraints parse error at {}:{}: group '{}' has no "
                "members",
                path,
                lineno,
                group_name);
          }
          continue;
        }
        if (ParseInt(second, parsed_pr)) {
          pr = std::clamp(parsed_pr, 0, 127);
        } else {
          members.push_back(second);
        }
      }

      std::string tok;
      while (iss >> tok) {
        members.push_back(tok);
      }
      if (members.empty()) {
        if (logger) {
          logger->warn(utl::DFT,
                       148,
                       "Scan constraints parse error at {}:{}: group has no "
                       "members",
                       path,
                       lineno);
        }
        continue;
      }

      RawGroup def;
      def.name = group_name;
      def.priority = pr;
      def.is_path = false;
      def.members = std::move(members);

      raw_groups[group_name] = std::move(def);
      raw_group_order.push_back(group_name);
      continue;
    }

    if (keyword == "path" || keyword == "strict_group" || keyword == "strict") {
      std::string first;
      if (!(iss >> first)) {
        if (logger) {
          logger->warn(
              utl::DFT,
              149,
              "Scan constraints parse error at {}:{}: expected 'path "
              "[<name>] [<priority>] <inst...>'",
              path,
              lineno);
        }
        continue;
      }

      int pr = default_group_priority_;
      std::string group_name;
      std::vector<std::string> members;

      int parsed_pr = 0;
      if (ParseInt(first, parsed_pr)) {
        pr = std::clamp(parsed_pr, 0, 127);
        group_name = fmt::format("__path_{}", lineno);
      } else {
        group_name = first;
        std::string second;
        if (!(iss >> second)) {
          if (logger) {
            logger->warn(
                utl::DFT,
                150,
                "Scan constraints parse error at {}:{}: path '{}' has no "
                "members",
                path,
                lineno,
                group_name);
          }
          continue;
        }
        if (ParseInt(second, parsed_pr)) {
          pr = std::clamp(parsed_pr, 0, 127);
        } else {
          members.push_back(second);
        }
      }

      std::string tok;
      while (iss >> tok) {
        members.push_back(tok);
      }
      if (members.size() < 2) {
        if (logger) {
          logger->warn(utl::DFT,
                       151,
                       "Scan constraints parse error at {}:{}: path '{}' must "
                       "have at least 2 instances",
                       path,
                       lineno,
                       group_name);
        }
        continue;
      }

      RawGroup def;
      def.name = group_name;
      def.priority = pr;
      def.is_path = true;
      def.members = std::move(members);

      raw_groups[group_name] = std::move(def);
      raw_group_order.push_back(group_name);
      continue;
    }

    if (keyword == "fixed_edge" || keyword == "fixededge" || keyword == "edge") {
      std::string from;
      std::string to;
      if (!(iss >> from >> to)) {
        if (logger) {
          logger->warn(
              utl::DFT,
              152,
              "Scan constraints parse error at {}:{}: expected 'fixed_edge <from> "
              "<to>'",
              path,
              lineno);
        }
        continue;
      }
      scan_order_fixed_edges_.push_back(
          {.from_inst = std::move(from), .to_inst = std::move(to)});
      continue;
    }

    if (keyword == "before") {
      std::string a;
      std::string b;
      if (!(iss >> a >> b)) {
        if (logger) {
          logger->warn(
              utl::DFT,
              153,
              "Scan constraints parse error at {}:{}: expected 'before <a> <b>'",
              path,
              lineno);
        }
        continue;
      }
      scan_order_before_.push_back({.before = std::move(a), .after = std::move(b)});
      continue;
    }

    if (keyword == "assign") {
      std::string chain_name;
      if (!(iss >> chain_name)) {
        if (logger) {
          logger->warn(utl::DFT,
                       154,
                       "Scan constraints parse error at {}:{}: expected 'assign "
                       "<chain> <inst/group...>'",
                       path,
                       lineno);
        }
        continue;
      }
      std::vector<std::string> items;
      std::string item;
      while (iss >> item) {
        items.push_back(item);
      }
      if (items.empty()) {
        if (logger) {
          logger->warn(utl::DFT,
                       155,
                       "Scan constraints parse error at {}:{}: assign has no "
                       "items",
                       path,
                       lineno);
        }
        continue;
      }
      raw_assignments.push_back(
          {.chain_name = std::move(chain_name), .items = std::move(items)});
      continue;
    }

    if (keyword == "exclude" || keyword == "exclude_instance"
        || keyword == "exclude_instances") {
      std::vector<std::string> items;
      std::string item;
      while (iss >> item) {
        items.push_back(item);
      }
      if (items.empty()) {
        if (logger) {
          logger->warn(utl::DFT,
                       186,
                       "Scan constraints parse error at {}:{}: exclude has no "
                       "items",
                       path,
                       lineno);
        }
        continue;
      }
      for (auto& v : items) {
        raw_excludes.push_back(std::move(v));
      }
      continue;
    }

    if (keyword == "exclude_instance_pattern" || keyword == "exclude_name_pattern"
        || keyword == "exclude_name_patterns") {
      std::string pat;
      bool any = false;
      while (iss >> pat) {
        raw_exclude_instance_patterns.push_back(std::move(pat));
        any = true;
      }
      if (!any && logger) {
        logger->warn(utl::DFT,
                     219,
                     "Scan constraints parse error at {}:{}: {} has no patterns",
                     path,
                     lineno,
                     keyword);
      }
      continue;
    }

    if (keyword == "exclude_master" || keyword == "exclude_master_pattern"
        || keyword == "exclude_master_patterns") {
      std::string pat;
      bool any = false;
      while (iss >> pat) {
        raw_exclude_master_patterns.push_back(std::move(pat));
        any = true;
      }
      if (!any && logger) {
        logger->warn(utl::DFT,
                     220,
                     "Scan constraints parse error at {}:{}: {} has no patterns",
                     path,
                     lineno,
                     keyword);
      }
      continue;
    }

    if (logger) {
      logger->warn(utl::DFT,
                   156,
                   "Scan constraints: ignoring unknown directive '{}' at {}:{}",
                   keyword,
                   path,
                   lineno);
    }
  }

  // Resolve hierarchical groups and apply them as "stay in one chain" constraints.
  std::unordered_map<std::string, std::vector<std::string>> expanded_groups;
  std::unordered_map<std::string, int> group_state;  // 0=unseen,1=visiting,2=done

  const auto expand_group = [&](const auto& self,
                                const std::string& name) -> const std::vector<std::string>& {
    const auto state_it = group_state.find(name);
    if (state_it != group_state.end() && state_it->second == 2) {
      return expanded_groups[name];
    }
    if (state_it != group_state.end() && state_it->second == 1) {
      if (logger) {
        logger->error(utl::DFT,
                      157,
                      "Scan constraints group cycle detected at '{}' (file '{}')",
                      name,
                      path);
      }
      // logger->error throws; return to satisfy compiler.
      return expanded_groups[name];
    }
    group_state[name] = 1;

    std::vector<std::string> out;
    std::unordered_set<std::string> seen;

    const auto it = raw_groups.find(name);
    if (it == raw_groups.end()) {
      group_state[name] = 2;
      expanded_groups[name] = out;
      return expanded_groups[name];
    }

    for (const std::string& tok : it->second.members) {
      const auto git = raw_groups.find(tok);
      if (git != raw_groups.end()) {
        const auto& sub = self(self, tok);
        for (const std::string& inst : sub) {
          if (seen.insert(inst).second) {
            out.push_back(inst);
          }
        }
      } else {
        if (seen.insert(tok).second) {
          out.push_back(tok);
        }
      }
    }

    group_state[name] = 2;
    expanded_groups[name] = std::move(out);
    return expanded_groups[name];
  };

  for (const std::string& name : raw_group_order) {
    (void) expand_group(expand_group, name);
  }

  scan_order_groups_.reserve(expanded_groups.size());
  for (const std::string& name : raw_group_order) {
    const auto it = raw_groups.find(name);
    if (it == raw_groups.end()) {
      continue;
    }
    ScanOrderGroupConstraint group;
    group.name = name;
    group.priority = std::clamp(it->second.priority, 0, 127);
    group.inst_names = expanded_groups[name];
    scan_order_groups_.push_back(std::move(group));

    // Convert strict path groups into fixed edges.
    if (it->second.is_path) {
      for (const std::string& tok : it->second.members) {
        if (raw_groups.find(tok) != raw_groups.end()) {
          if (logger) {
            logger->error(
                utl::DFT,
                158,
                "Scan constraints: path '{}' cannot reference group '{}' (file "
                "'{}')",
                name,
                tok,
                path);
          }
        }
      }
      for (std::size_t i = 1; i < it->second.members.size(); ++i) {
        scan_order_fixed_edges_.push_back(
            {.from_inst = it->second.members[i - 1],
             .to_inst = it->second.members[i]});
      }
    }
  }

  // Resolve assignment constraints to instance->chain map (expanding groups).
  for (const RawAssignment& asn : raw_assignments) {
    for (const std::string& item : asn.items) {
      const auto git = expanded_groups.find(item);
      if (git != expanded_groups.end()) {
        for (const std::string& inst : git->second) {
          const auto [it, inserted]
              = instance_to_chain_name_.emplace(inst, asn.chain_name);
          if (!inserted && it->second != asn.chain_name) {
            if (logger) {
              logger->error(
                  utl::DFT,
                  159,
                  "Scan constraints: instance '{}' assigned to multiple chains "
                  "('{}' and '{}') in file '{}'",
                  inst,
                  it->second,
                  asn.chain_name,
                  path);
            }
          }
        }
      } else {
        const auto [it, inserted]
            = instance_to_chain_name_.emplace(item, asn.chain_name);
        if (!inserted && it->second != asn.chain_name) {
          if (logger) {
            logger->error(
                utl::DFT,
                160,
                "Scan constraints: instance '{}' assigned to multiple chains "
                "('{}' and '{}') in file '{}'",
                item,
                it->second,
                asn.chain_name,
                path);
          }
        }
      }
    }
  }

  // Resolve instance exclusions (expanding groups).
  for (const std::string& item : raw_excludes) {
    const auto git = expanded_groups.find(item);
    if (git != expanded_groups.end()) {
      for (const std::string& inst : git->second) {
        excluded_instances_.insert(inst);
      }
    } else {
      excluded_instances_.insert(item);
    }
  }

  // Record pattern-based instance/master exclusions.
  excluded_instance_patterns_ = std::move(raw_exclude_instance_patterns);
  excluded_master_patterns_ = std::move(raw_exclude_master_patterns);

  return true;
}

void ScanArchitectConfig::report(utl::Logger* logger) const
{
  logger->report("Scan Architect Config:");
  logger->report("- Chain Count: {}", utils::FormatForReport(chain_count_));
  logger->report("- Max Length: {}", utils::FormatForReport(max_length_));
  logger->report("- Max Chains: {}", utils::FormatForReport(max_chains_));
  logger->report("- Clock Mixing: {}", ClockMixingName(clock_mixing_));
  logger->report("- Polarity Mode: {}", PolarityModeName(polarity_mode_));
  logger->report("- Scan Enable Name Pattern: {}", scan_enable_name_pattern_);
  logger->report("- Scan In Name Pattern: {}", scan_in_name_pattern_);
  logger->report("- Scan Out Name Pattern: {}", scan_out_name_pattern_);
  logger->report("- Scan Order Metric: {}",
                 ScanOrderMetricName(scan_order_metric_));
  logger->report("- Scan Order Solver: {}",
                 ScanOrderSolverName(scan_order_solver_));
  logger->report("- Vertical Weight: {:.3f}", vertical_weight_);
  logger->report("- Max Imbalance: {:.1f}%", max_imbalance_percent_);
  if (timing_weight_setup_ != 0.0 || timing_weight_hold_ != 0.0) {
    logger->report("- Timing Setup Weight: {:.3f}", timing_weight_setup_);
    logger->report("- Timing Hold Weight: {:.3f}", timing_weight_hold_);
    logger->report("- Timing Critical Slack: {:.3f}", timing_critical_slack_);
  }
  if (scan_order_solver_ == ScanOrderSolver::ScanOpt) {
    logger->report("- ScanOpt Rounds: {}", scanopt_rounds_);
    logger->report("- ScanOpt Seed: {}", scanopt_seed_);
    logger->report("- ScanOpt Temp Control: {}", scanopt_temp_control_);
    logger->report("- ScanOpt TDiv: {:.3f}", scanopt_t_div_);
    if (scanopt_time_limit_seconds_ > 0.0) {
      logger->report("- ScanOpt Time Limit: {:.3f}s", scanopt_time_limit_seconds_);
    }
  }
  if (!scan_order_groups_.empty() || !scan_order_fixed_edges_.empty()) {
    logger->report("- Scan Order Constraints:");
    logger->report("  - Default Priority: {}", default_group_priority_);
    logger->report("  - Groups: {}", scan_order_groups_.size());
    logger->report("  - Fixed Edges: {}", scan_order_fixed_edges_.size());
  }
  if (!scan_order_before_.empty()) {
    logger->report("  - Before Constraints: {}", scan_order_before_.size());
  }
  if (!chain_names_.empty()) {
    logger->report("- Scan Chains (user-defined): {}", chain_names_.size());
  }
  if (!instance_to_chain_name_.empty()) {
    logger->report("- Scan Chain Assignments: {}", instance_to_chain_name_.size());
  }
  if (auto_exclude_shift_registers_) {
    logger->report("- Auto Exclude Shift Registers: yes (min_length={})",
                   shift_register_min_length_);
  }
  if (prefer_qbar_scan_out_) {
    logger->report("- Prefer Qbar Scan Out: yes");
  }
  if (!auto_excluded_instances_.empty()) {
    logger->report("- Auto Excluded Instances: {}", auto_excluded_instances_.size());
  }
  if (!excluded_instances_.empty()) {
    logger->report("- Excluded Instances: {}", excluded_instances_.size());
  }
  if (!excluded_instance_patterns_.empty()) {
    logger->report("- Excluded Instance Patterns: {}",
                   excluded_instance_patterns_.size());
  }
  if (!excluded_master_patterns_.empty()) {
    logger->report("- Excluded Master Patterns: {}",
                   excluded_master_patterns_.size());
  }
}

std::string ScanArchitectConfig::ClockMixingName(
    ScanArchitectConfig::ClockMixing clock_mixing)
{
  switch (clock_mixing) {
    case ScanArchitectConfig::ClockMixing::NoMix:
      return "No Mix";
    case ScanArchitectConfig::ClockMixing::ClockMix:
      return "Clock Mix";
    default:
      return "Missing case in ClockMixingName";
  }
}

std::string ScanArchitectConfig::PolarityModeName(
    ScanArchitectConfig::PolarityMode mode)
{
  switch (mode) {
    case ScanArchitectConfig::PolarityMode::Mid:
      return "Mid";
    case ScanArchitectConfig::PolarityMode::Strict:
      return "Strict";
    default:
      return "Missing case in PolarityModeName";
  }
}

std::string ScanArchitectConfig::ScanOrderMetricName(
    ScanArchitectConfig::ScanOrderMetric metric)
{
  switch (metric) {
    case ScanArchitectConfig::ScanOrderMetric::Placement:
      return "Placement (pin-to-pin)";
    case ScanArchitectConfig::ScanOrderMetric::PinToNet:
      return "Pin-to-net (routing-aware)";
    default:
      return "Missing case in ScanOrderMetricName";
  }
}

std::string ScanArchitectConfig::ScanOrderSolverName(
    ScanArchitectConfig::ScanOrderSolver solver)
{
  switch (solver) {
    case ScanArchitectConfig::ScanOrderSolver::Heuristic:
      return "Heuristic";
    case ScanArchitectConfig::ScanOrderSolver::ScanOpt:
      return "ScanOpt";
    case ScanArchitectConfig::ScanOrderSolver::UclaScanOpt:
      return "UCLA ScanOpt";
    default:
      return "Missing case in ScanOrderSolverName";
  }
}

}  // namespace dft
