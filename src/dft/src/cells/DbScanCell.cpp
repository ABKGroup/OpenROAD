// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2023-2025, The OpenROAD Authors

#include "DbScanCell.hh"

#include <algorithm>
#include <utility>

namespace dft {

DbScanCell::DbScanCell(const std::string& name,
                       std::unique_ptr<ClockDomain> clock_domain,
                       odb::dbInst* inst,
                       ScanLoad scan_in,
                       ScanLoad scan_enable,
                       ScanDriver scan_out,
                       uint64_t bits,
                       bool has_scan_enable,
                       utl::Logger* logger)
    : ScanCell(name, std::move(clock_domain), logger),
      inst_(inst),
      scan_in_(std::move(scan_in)),
      scan_enable_(std::move(scan_enable)),
      scan_out_(std::move(scan_out)),
      bits_(bits),
      has_scan_enable_(has_scan_enable)
{
}

uint64_t DbScanCell::getBits() const
{
  return bits_;
}

void DbScanCell::connectScanEnable(const ScanDriver& driver) const
{
  if (!has_scan_enable_) {
    return;
  }
  Connect(scan_enable_, driver, /*preserve=*/false);
}

void DbScanCell::connectScanIn(const ScanDriver& driver) const
{
  Connect(scan_in_, driver, /*preserve=*/false);
}

void DbScanCell::connectScanOut(const ScanLoad& load) const
{
  // The scan out usually will be connected to functional data paths already, we
  // need to preserve the connections.
  Connect(load, scan_out_, /*preserve=*/true);
}

ScanLoad DbScanCell::getScanEnable() const
{
  return scan_enable_;
}

ScanLoad DbScanCell::getScanIn() const
{
  return scan_in_;
}

ScanDriver DbScanCell::getScanOut() const
{
  return scan_out_;
}

odb::dbInst* DbScanCell::getDbInst() const
{
  return inst_;
}

odb::Point DbScanCell::getOrigin() const
{
  return inst_ != nullptr ? inst_->getLocation() : odb::Point(0, 0);
}

bool DbScanCell::isPlaced() const
{
  return inst_ != nullptr && inst_->isPlaced();
}

}  // namespace dft
