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
# Module to find and setup Kokkos Kernels.
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
find_package(KokkosKernels REQUIRED PATHS "${KokkosKernels_DIR}")

message(STATUS "Kokkos Kernels VERSION: ${KokkosKernels_VERSION}")

# Kokkos Kernels must use the same Kokkos major.minor release.  Keep this
# check in the build tree as well as the installed package configuration.
string(REGEX MATCH "^[0-9]+\\.[0-9]+" SUNDIALS_KOKKOS_VERSION_SERIES
       "${Kokkos_VERSION}")
string(REGEX MATCH "^[0-9]+\\.[0-9]+" SUNDIALS_KOKKOS_KERNELS_VERSION_SERIES
       "${KokkosKernels_VERSION}")
if(NOT SUNDIALS_KOKKOS_VERSION_SERIES STREQUAL
       SUNDIALS_KOKKOS_KERNELS_VERSION_SERIES)
  message(FATAL_ERROR
          "SUNDIALS requires Kokkos and KokkosKernels to have matching "
          "major.minor versions (Kokkos ${Kokkos_VERSION}, "
          "KokkosKernels ${KokkosKernels_VERSION})")
endif()

# -----------------------------------------------------------------------------
# Section 4: Test the TPL
# -----------------------------------------------------------------------------

if(SUNDIALS_ENABLE_KOKKOS_KERNELS_CHECKS)
  message(CHECK_START "Testing Kokkos Kernels")
  message(CHECK_PASS "success")
else()
  message(STATUS "Skipped Kokkos Kernels checks.")
endif()
