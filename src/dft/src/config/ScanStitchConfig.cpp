// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2024-2025, The OpenROAD Authors

#include "ScanStitchConfig.hh"

#include <string_view>

#include "Formatting.hh"
#include "utl/Logger.h"

namespace dft {

void ScanStitchConfig::setEnableNamePattern(
    std::string_view enable_name_pattern)
{
  enable_name_pattern_ = enable_name_pattern;
}
std::string_view ScanStitchConfig::getEnableNamePattern() const
{
  return enable_name_pattern_;
};

void ScanStitchConfig::setInNamePattern(std::string_view in_name_pattern)
{
  in_name_pattern_ = in_name_pattern;
};
std::string_view ScanStitchConfig::getInNamePattern() const
{
  return in_name_pattern_;
};

void ScanStitchConfig::setOutNamePattern(std::string_view out_name_pattern)
{
  out_name_pattern_ = out_name_pattern;
};
std::string_view ScanStitchConfig::getOutNamePattern() const
{
  return out_name_pattern_;
};

void ScanStitchConfig::setInsertLockup(bool enable)
{
  insert_lockup_ = enable;
}

bool ScanStitchConfig::getInsertLockup() const
{
  return insert_lockup_;
}

void ScanStitchConfig::setLockupCellRising(std::string_view cell_name)
{
  lockup_cell_rising_ = cell_name;
}

std::string_view ScanStitchConfig::getLockupCellRising() const
{
  return lockup_cell_rising_;
}

void ScanStitchConfig::setLockupCellFalling(std::string_view cell_name)
{
  lockup_cell_falling_ = cell_name;
}

std::string_view ScanStitchConfig::getLockupCellFalling() const
{
  return lockup_cell_falling_;
}

void ScanStitchConfig::setLockupInPin(std::string_view pin_name)
{
  lockup_in_pin_ = pin_name;
}

std::string_view ScanStitchConfig::getLockupInPin() const
{
  return lockup_in_pin_;
}

void ScanStitchConfig::setLockupOutPin(std::string_view pin_name)
{
  lockup_out_pin_ = pin_name;
}

std::string_view ScanStitchConfig::getLockupOutPin() const
{
  return lockup_out_pin_;
}

void ScanStitchConfig::setLockupClockPinRising(std::string_view pin_name)
{
  lockup_clock_pin_rising_ = pin_name;
}

std::string_view ScanStitchConfig::getLockupClockPinRising() const
{
  return lockup_clock_pin_rising_;
}

void ScanStitchConfig::setLockupClockPinFalling(std::string_view pin_name)
{
  lockup_clock_pin_falling_ = pin_name;
}

std::string_view ScanStitchConfig::getLockupClockPinFalling() const
{
  return lockup_clock_pin_falling_;
}

void ScanStitchConfig::setTimingBufferCell(std::string_view cell_name)
{
  timing_buffer_cell_ = cell_name;
}

std::string_view ScanStitchConfig::getTimingBufferCell() const
{
  return timing_buffer_cell_;
}

void ScanStitchConfig::setTimingBufferInPin(std::string_view pin_name)
{
  timing_buffer_in_pin_ = pin_name;
}

std::string_view ScanStitchConfig::getTimingBufferInPin() const
{
  return timing_buffer_in_pin_;
}

void ScanStitchConfig::setTimingBufferOutPin(std::string_view pin_name)
{
  timing_buffer_out_pin_ = pin_name;
}

std::string_view ScanStitchConfig::getTimingBufferOutPin() const
{
  return timing_buffer_out_pin_;
}

void ScanStitchConfig::report(utl::Logger* logger) const
{
  logger->report("Scan Stitch Config:");
  logger->report("- Scan Enable Name Pattern: '{}'", enable_name_pattern_);
  logger->report("- Scan In Name Pattern: '{}'", in_name_pattern_);
  logger->report("- Scan Out Name Pattern: '{}'", out_name_pattern_);
  logger->report("- Insert Lockup Elements: {}",
                 insert_lockup_ ? "true" : "false");
  if (insert_lockup_) {
    logger->report("- Lockup (rising) Cell: '{}'", lockup_cell_rising_);
    logger->report("- Lockup (rising) Clock Pin: '{}'", lockup_clock_pin_rising_);
    logger->report("- Lockup (falling) Cell: '{}'", lockup_cell_falling_);
    logger->report("- Lockup (falling) Clock Pin: '{}'", lockup_clock_pin_falling_);
    logger->report("- Lockup Data In Pin: '{}'", lockup_in_pin_);
    logger->report("- Lockup Data Out Pin: '{}'", lockup_out_pin_);
  }

  logger->report("- Timing Buffer Cell: '{}'", timing_buffer_cell_);
  if (!timing_buffer_cell_.empty()) {
    logger->report("- Timing Buffer In Pin: '{}'", timing_buffer_in_pin_);
    logger->report("- Timing Buffer Out Pin: '{}'", timing_buffer_out_pin_);
  }
}

}  // namespace dft
