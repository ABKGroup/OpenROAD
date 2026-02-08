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

}  // namespace

std::function<size_t(const ClockDomain&)> GetClockDomainHashFn(
    const ScanArchitectConfig& config,
    utl::Logger* logger)
{
  switch (config.getClockMixing()) {
    // For NoMix, every clock is a different domain. Edge polarity is handled
    // within each chain via the "mid" polarity ordering (falling before rising),
    // rather than forcing separate chains per edge.
    case ScanArchitectConfig::ClockMixing::NoMix:
      return [](const ClockDomain& clock_domain) {
        return HashClockName(clock_domain.getClockName());
      };
    case ScanArchitectConfig::ClockMixing::ClockMix:
      return [](const ClockDomain& clock_domain) { return 1; };
    default:
      // Not implemented
      logger->error(utl::DFT, 4, "Clock mix config requested is not supported");
  }
}

}  // namespace dft
