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

message(STATUS "Added CUDA NVECTOR benchmark")

sundials_add_nvector_benchmark(nvector_cuda_benchmark
  SOURCES test_nvector_performance_cuda.cu
  SUNDIALS_TARGETS sundials_nveccuda
  INSTALL_SUBDIR nvector/cuda
  )

# NUM_NODES based on 1 lassen node
sundials_add_benchmark(nvector_cuda_benchmark nvector_cuda_benchmark
  NUM_NODES 4
  ENABLE_GPU
  BENCHMARK_ARGS "10 10 10 10 10 10"
  )
