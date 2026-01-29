// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2024-2025, The OpenROAD Authors

#pragma once

#include "ScanArchitectConfig.hh"
#include "ScanCell.hh"
#include "utl/Logger.h"

namespace dft {

// Order scan cells to reduce wirelength
void OptimizeScanWirelength(std::vector<std::unique_ptr<ScanCell>>& cells,
                            const ScanArchitectConfig& config,
                            utl::Logger* logger);

// Endpoint-aware ordering: includes optional begin/end wirelength in the
// objective (begin->first scan-in, last scan-out->end).
void OptimizeScanWirelength(
    std::vector<std::unique_ptr<ScanCell>>& cells,
    const ScanArchitectConfig& config,
    utl::Logger* logger,
    const std::optional<ScanArchitectConfig::ChainEndpoints>& endpoints);

}  // namespace dft
