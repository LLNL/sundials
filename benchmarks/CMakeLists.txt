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
# benchmarks level CMakeLists.txt for SUNDIALS
# ---------------------------------------------------------------

if(NOT (CMAKE_BUILD_TYPE STREQUAL "Release" OR
      CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo"))
  message(WARNING "SUNDIALS is not being built in a \"Release\" mode, "
          "benchmark performance will be affected")
endif()

sundials_option(BENCHMARK_NVECTOR BOOL "NVector benchmarks are on" ON)

#----------------------------------------
# Add specific benchmarks
#----------------------------------------

if(ENABLE_MPI)
add_subdirectory(diffusion_2D)
add_subdirectory(advection_reaction_3D)
endif()

# Add the nvector benchmarks
if(BENCHMARK_NVECTOR)
  add_subdirectory(nvector)
endif()
