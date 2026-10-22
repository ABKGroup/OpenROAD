# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2025, The OpenROAD Authors

load("@rules_cc//cc:cc_library.bzl", "cc_library")
load("//test:regression.bzl", "message_srcs")

package(
    default_visibility = ["//:__subpackages__"],
    features = ["layering_check"],
)

cc_library(
    name = "ucla_scanopt",
    srcs = [
        "UclaScanOpt.cc",
    ],
    hdrs = [
        "UclaScanOpt.hh",
    ],
    includes = [
        ".",
        "UCLApack-3-010411"
    ],
    deps = [
        "//src/dft/src/ucla_scanopt/UCLApack-3-010411",
    ],
)

message_srcs(
    visibility = ["//src/dft:__pkg__"],
)
