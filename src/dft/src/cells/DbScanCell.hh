// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2023-2025, The OpenROAD Authors

#pragma once

#include <cstdint>
#include <memory>

#include "ClockDomain.hh"
#include "ScanCell.hh"
#include "ScanPin.hh"
#include "odb/db.h"

namespace dft {

// A scan cell constructed from explicit DB pins (e.g. imported SCANDEF/ODB scan
// chains). Scan enable is optional: if absent, connectScanEnable() is a no-op.
class DbScanCell : public ScanCell
{
 public:
  DbScanCell(const std::string& name,
             std::unique_ptr<ClockDomain> clock_domain,
             odb::dbInst* inst,
             ScanLoad scan_in,
             ScanLoad scan_enable,
             ScanDriver scan_out,
             uint64_t bits,
             bool has_scan_enable,
             utl::Logger* logger);

  uint64_t getBits() const override;
  void connectScanEnable(const ScanDriver& driver) const override;
  void connectScanIn(const ScanDriver& driver) const override;
  void connectScanOut(const ScanLoad& load) const override;
  ScanLoad getScanEnable() const override;
  ScanLoad getScanIn() const override;
  ScanDriver getScanOut() const override;
  odb::dbInst* getDbInst() const override;

  odb::Point getOrigin() const override;
  bool isPlaced() const override;

 private:
  odb::dbInst* inst_{nullptr};
  ScanLoad scan_in_;
  ScanLoad scan_enable_;
  ScanDriver scan_out_;
  uint64_t bits_{0};
  bool has_scan_enable_{false};
};

}  // namespace dft
