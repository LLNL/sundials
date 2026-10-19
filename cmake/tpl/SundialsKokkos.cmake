# ------------------------------------------------------------------------------
# Programmer(s): David J. Gardner @ LLNL
# -----------------------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2025-2026, Lawrence Livermore National Security,
# University of Maryland Baltimore County, and the SUNDIALS contributors.
# Copyright (c) 2013-2025, Lawrence Livermore National Security
# and Southern Methodist University.
# Copyright (c) 2002-2013, Lawrence Livermore National Security.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# -----------------------------------------------------------------------------
# Module to find and setup Kokkos.
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Section 1: Include guard
# -----------------------------------------------------------------------------

include_guard(GLOBAL)

# -----------------------------------------------------------------------------
# Section 2: Check to make sure options are compatible
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Section 3: Find the TPL
# -----------------------------------------------------------------------------
find_package(Kokkos REQUIRED PATHS "${Kokkos_DIR}")

# SUNDIALS uses the common Kokkos 4/5 view and policy interfaces. Kokkos 5
# itself requires C++20, whereas Kokkos 4 is usable with C++17.
if(Kokkos_VERSION VERSION_LESS 4.0.0)
  message(FATAL_ERROR "SUNDIALS requires Kokkos 4.0 or newer (found ${Kokkos_VERSION})")
elseif(Kokkos_VERSION VERSION_LESS 5.0.0)
  set(SUNDIALS_KOKKOS_CXX_FEATURE cxx_std_17 CACHE INTERNAL
      "C++ feature required by the discovered Kokkos" FORCE)
else()
  set(SUNDIALS_KOKKOS_CXX_FEATURE cxx_std_20 CACHE INTERNAL
      "C++ feature required by the discovered Kokkos" FORCE)
endif()

# We should be able to use Kokkos_DEVICES directly but it seems to get removed
# or unset in some CMake versions
set(KOKKOS_EXAMPLES_BACKENDS
    "${Kokkos_DEVICES}"
    CACHE STRING "Kokkos backends to build examples with")
mark_as_advanced(FORCE KOKKOS_EXAMPLES_BACKENDS)
message(STATUS "Kokkos VERSION: ${Kokkos_VERSION}")

# -----------------------------------------------------------------------------
# Section 4: Test the TPL
# -----------------------------------------------------------------------------

if(SUNDIALS_ENABLE_KOKKOS_CHECKS)
  message(CHECK_START "Testing Kokkos")
  message(CHECK_PASS "success")
else()
  message(STATUS "Skipped Kokkos checks.")
endif()
