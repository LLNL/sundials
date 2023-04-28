# ---------------------------------------------------------------
# Programmer(s): Cody J. Balos @ LLNL
# ---------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2023, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# ---------------------------------------------------------------

message(STATUS "Added HIP NVECTOR benchmark")

sundials_add_nvector_benchmark(nvector_hip_benchmark
  SOURCES test_nvector_performance_hip.cpp
  SUNDIALS_TARGETS sundials_nvechip
  INSTALL_SUBDIR nvector/hip
  )

# NUM_CORES based on Corona system for 1 node
sundials_add_benchmark(nvector_hip_benchmark nvector_hip_benchmark
  NUM_CORES 8
  ENABLE_GPU
  BENCHMARK_ARGS "10 10 10 10 10 10"
  )
