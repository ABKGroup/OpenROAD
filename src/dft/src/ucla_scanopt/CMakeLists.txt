# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2026, The OpenROAD Authors

set(UCLA_ROOT ${PROJECT_SOURCE_DIR}/third_party/UCLApack-3-010411)

add_library(dft_ucla_scanopt_lib
  UclaScanOpt.cc

  # UCLApack-3-010411 (ScanOpt dependencies).
  ${UCLA_ROOT}/ABKCommon/Tausworthe.cxx
  ${UCLA_ROOT}/ABKCommon/abkassert.cxx
  ${UCLA_ROOT}/ABKCommon/abkcpunorm.cxx
  ${UCLA_ROOT}/ABKCommon/abkfunc.cxx
  ${UCLA_ROOT}/ABKCommon/abkio.cxx
  ${UCLA_ROOT}/ABKCommon/abkroot.cxx
  ${UCLA_ROOT}/ABKCommon/abkstring.cxx
  ${UCLA_ROOT}/ABKCommon/abktimer.cxx
  ${UCLA_ROOT}/ABKCommon/crc32.cxx
  ${UCLA_ROOT}/ABKCommon/infolines.cxx
  ${UCLA_ROOT}/ABKCommon/paramproc.cxx
  ${UCLA_ROOT}/ABKCommon/pathDelims.cxx
  ${UCLA_ROOT}/ABKCommon/platfDepend.cxx
  ${UCLA_ROOT}/ABKCommon/r1279.cxx
  ${UCLA_ROOT}/ABKCommon/r250.cxx
  ${UCLA_ROOT}/ABKCommon/seed.cxx
  ${UCLA_ROOT}/ABKCommon/verbosity.cxx
  ${UCLA_ROOT}/ABKCommon/MD5/md5.cxx

  ${UCLA_ROOT}/ScanChain/scanCells.cxx
  ${UCLA_ROOT}/ScanChain/scanChain.cxx

  ${UCLA_ROOT}/ScanOpt/kenKDTree.cxx
  ${UCLA_ROOT}/ScanOpt/optimizer.cxx
  ${UCLA_ROOT}/ScanOpt/optimizer1.cxx
  ${UCLA_ROOT}/ScanOpt/path.cxx
  ${UCLA_ROOT}/ScanOpt/randOpts.cxx
  ${UCLA_ROOT}/ScanOpt/scanTourDZ.cxx
)

target_include_directories(dft_ucla_scanopt_lib
  PUBLIC
    ${CMAKE_CURRENT_LIST_DIR}
    ${UCLA_ROOT}
    ${UCLA_ROOT}/compat
)

target_compile_options(dft_ucla_scanopt_lib
  PRIVATE
    $<$<CXX_COMPILER_ID:MSVC>:/w>
    $<$<NOT:$<CXX_COMPILER_ID:MSVC>>:-w>
)
