// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2024-2025, The OpenROAD Authors

#include "ClockDomainHash.hh"

#include <cstddef>
#include <cstdint>
#include <functional>
#include <string_view>

#include "utl/Logger.h"

namespace dft {

namespace {

size_t FNV1a(const uint8_t* data, size_t length)
{
  size_t hash = 0xcbf29ce484222325;
  for (size_t i = 0; i < length; i += 1) {
    hash ^= data[i];
    hash *= 0x00000100000001b3;
  }
  return hash;
}

size_t HashClockName(std::string_view name)
{
  return FNV1a(reinterpret_cast<const uint8_t*>(name.data()), name.size());
}

size_t HashClockNameAndEdge(const ClockDomain& clock_domain)
{
  size_t hash = 0xcbf29ce484222325;
  for (char c : clock_domain.getClockName()) {
    hash ^= static_cast<uint8_t>(c);
    hash *= 0x00000100000001b3;
  }
  const uint8_t edge = clock_domain.getClockEdge() == ClockEdge::Falling ? 'F' : 'R';
  hash ^= edge;
  hash *= 0x00000100000001b3;
  return hash;
}

}  // namespace

std::function<size_t(const ClockDomain&)> GetClockDomainHashFn(
    const ScanArchitectConfig& config,
    utl::Logger* logger)
{
  const auto polarity_mode = config.getPolarityMode();
  switch (config.getClockMixing()) {
    case ScanArchitectConfig::ClockMixing::NoMix:
      if (polarity_mode == ScanArchitectConfig::PolarityMode::Strict) {
        return [](const ClockDomain& clock_domain) {
          return HashClockNameAndEdge(clock_domain);
        };
      }
      return [](const ClockDomain& clock_domain) {
        return HashClockName(clock_domain.getClockName());
      };
    case ScanArchitectConfig::ClockMixing::ClockMix:
      if (polarity_mode == ScanArchitectConfig::PolarityMode::Strict) {
        return [](const ClockDomain& clock_domain) {
          return clock_domain.getClockEdge() == ClockEdge::Falling ? 2 : 1;
        };
      }
      return [](const ClockDomain& clock_domain) { return 1; };
    default:
      // Not implemented
      logger->error(utl::DFT, 4, "Clock mix config requested is not supported");
  }
}

}  // namespace dft
