# -----------------------------------------------------------------------------
# Programmer(s): Cody J. Balos @ LLNL
# -----------------------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2023, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# -----------------------------------------------------------------------------
# Module to find and setup GINKGO correctly.
# Created from the SundialsTPL.cmake template.
# All SUNDIALS modules that find and setup a TPL must:
#
# 1. Check to make sure the SUNDIALS configuration and the TPL is compatible.
# 2. Find the TPL.
# 3. Check if the TPL works with SUNDIALS, UNLESS the override option
# TPL_WORKS is TRUE - in this case the tests should not be performed and it
# should be assumed that the TPL works with SUNDIALS.
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Section 1: Include guard
# -----------------------------------------------------------------------------

if(NOT DEFINED SUNDIALS_GINKGO_INCLUDED)
  set(SUNDIALS_GINKGO_INCLUDED)
else()
  return()
endif()

# -----------------------------------------------------------------------------
# Section 2: Check to make sure options are compatible
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# Section 3: Find the TPL
# -----------------------------------------------------------------------------
find_package(Ginkgo REQUIRED
             HINTS "${Ginkgo_DIR}"
             NO_DEFAULT_PATH)

message(STATUS "GINKGO VERSION:     ${GINKGO_PROJECT_VERSION}")
message(STATUS "GINKGO BUILD TYPE:  ${GINKGO_BUILD_TYPE}")
message(STATUS "GINKGO LIBRARIES:   ${GINKGO_INTERFACE_LINK_LIBRARIES}")
message(STATUS "GINKGO LINK FLAGS:  ${GINKGO_INTERFACE_LINK_FLAGS}")
message(STATUS "GINKGO CXX FLAGS:   ${GINKGO_INTERFACE_CXX_FLAGS}")

# -----------------------------------------------------------------------------
# Section 4: Test the TPL
# -----------------------------------------------------------------------------
if(Ginkgo_FOUND AND (NOT GINKGO_WORKS))
  if(SUNDIALS_PRECISION MATCHES "extended|EXTENDED")
    print_error("SUNDIALS GINKGO interface is not compatible with extended precision")
  endif()

  if(SUNDIALS_GINKGO_BACKENDS MATCHES "CUDA" AND NOT ENABLE_CUDA)
    print_error("SUNDIALS_GINKGO_BACKENDS includes CUDA but CUDA is not enabled. Set ENABLE_CUDA=ON or change the backend.")
  endif()

  if(SUNDIALS_GINKGO_BACKENDS MATCHES "HIP" AND NOT ENABLE_HIP)
    print_error("SUNDIALS_GINKGO_BACKENDS includes HIP but HIP is not enabled. Set ENABLE_HIP=ON or change the backend.")
  endif()

  if(SUNDIALS_GINKGO_BACKENDS MATCHES "DPCPP" AND NOT ENABLE_SYCL)
    print_error("SUNDIALS_GINKGO_BACKENDS includes DPC++ but SYCL/DPC++ is not enabled. Set ENABLE_SYCL=ON or change the backend.")
  endif()

  if(SUNDIALS_GINKGO_BACKENDS MATCHES "OMP" AND NOT ENABLE_OPENMP)
    print_error("SUNDIALS_GINKGO_BACKENDS includes OMP but OpenMP is not enabled. Set ENABLE_OPENMP=ON or change the backend.")
  endif()

  message(STATUS "Checking if GINKGO works... OK")
  set(GINKGO_WORKS TRUE CACHE BOOL "GINKGO works with SUNDIALS as configured" FORCE)
elseif(Ginkgo_FOUND AND GINKGO_WORKS)
  message(STATUS "Skipped GINKGO tests, assuming GINKGO works with SUNDIALS.")
endif()
