# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2025, The OpenROAD Authors

load("@rules_cc//cc:cc_library.bzl", "cc_library")
load("//test:regression.bzl", "message_srcs")

package(
    default_visibility = ["//:__subpackages__"],
    features = ["layering_check"],
)

cc_library(
    name = "formatting",
    hdrs = [
        "Formatting.hh",
    ],
    includes = [
        ".",
    ],
    deps = [
        "@spdlog",
    ],
)

cc_library(
    name = "utils",
    srcs = [
        "ScanPin.cpp",
        "ShiftRegister.cpp",
        "Utils.cpp",
    ],
    hdrs = [
        "Formatting.hh",
        "ScanPin.hh",
        "ShiftRegister.hh",
        "Utils.hh",
    ],
    includes = [
        ".",
    ],
    deps = [
        "//src/dbSta",
        "//src/dbSta:dbNetwork",
        "//src/dft/src/config",
        "//src/odb/src/db",
        "//src/sta:opensta_lib",
        "//src/utl",
        "@spdlog",
    ],
)

message_srcs(
    visibility = ["//src/dft:__pkg__"],
)
