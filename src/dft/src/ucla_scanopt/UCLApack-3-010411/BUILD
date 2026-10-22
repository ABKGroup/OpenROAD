# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2025, The OpenROAD Authors

load("@rules_cc//cc:cc_library.bzl", "cc_library")
load("//test:regression.bzl", "message_srcs")

package(
    default_visibility = ["//:__subpackages__"],
    features = ["layering_check"],
)

cc_library(
    name = "UCLApack-3-010411",
    srcs = [
        "ABKCommon/Tausworthe.cxx",
        "ABKCommon/abkassert.cxx",
        "ABKCommon/abkcpunorm.cxx",
        "ABKCommon/abkfunc.cxx",
        "ABKCommon/abkio.cxx",
        "ABKCommon/abkroot.cxx",
        "ABKCommon/abkstring.cxx",
        "ABKCommon/abktimer.cxx",
        "ABKCommon/crc32.cxx",
        "ABKCommon/infolines.cxx",
        "ABKCommon/paramproc.cxx",
        "ABKCommon/pathDelims.cxx",
        "ABKCommon/platfDepend.cxx",
        "ABKCommon/r1279.cxx",
        "ABKCommon/r250.cxx",
        "ABKCommon/seed.cxx",
        "ABKCommon/verbosity.cxx",
        "ABKCommon/MD5/md5.cxx",
        "ScanChain/scanCells.cxx",
        "ScanChain/scanChain.cxx",
        "ScanOpt/kenKDTree.cxx",
        "ScanOpt/optimizer.cxx",
        "ScanOpt/optimizer1.cxx",
        "ScanOpt/path.cxx",
        "ScanOpt/randOpts.cxx",
        "ScanOpt/scanTourDZ.cxx",
    ],
    hdrs = glob([
        "ABKCommon/*.h",
        "ScanChain/*.h",
        "ScanOpt/*.h",
        "compat/*.h",
    ]),
    includes = [
        ".",
        "compat"
    ],
    deps = [
    ],
    copts = ["-Wno-misleading-indentation"],

)

message_srcs(
    visibility = ["//src/dft:__pkg__"],
)
