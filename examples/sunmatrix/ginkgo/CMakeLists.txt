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

# Example lists are tuples "name\;args\;type" where the type is
# 'develop' for examples excluded from 'make test' in releases
set(examples
  "test_sunmatrix_ginkgo.cpp\;100 100 0\;"
  "test_sunmatrix_ginkgo.cpp\;100 10 0\;"
  "test_sunmatrix_ginkgo.cpp\;10 100 0\;"
  "test_sunmatrix_ginkgo.cpp\;100 100 1\;"
  "test_sunmatrix_ginkgo.cpp\;100 10 1\;"
  "test_sunmatrix_ginkgo.cpp\;10 100 1\;"
)

# Add source directory to include directories
include_directories(..)

sundials_add_examples_ginkgo(examples
  TARGETS test_sunmatrix_obj sundials_sunmatrixdense
  BACKENDS REF OMP CUDA HIP DPCPP
  UNIT_TEST)

# Install the targets
if(EXAMPLES_INSTALL)

  if(SUNDIALS_GINKGO_BACKENDS MATCHES "CUDA")
    list(APPEND vectors nveccuda)
  endif()
  if(SUNDIALS_GINKGO_BACKENDS MATCHES "HIP")
    list(APPEND vectors nvechip)
  endif()
  if(SUNDIALS_GINKGO_BACKENDS MATCHES "DPCPP")
    list(APPEND vectors nvecsycl)
  endif()
  if(SUNDIALS_GINKGO_BACKENDS MATCHES "OMP")
    list(APPEND vectors nvecopenmp)
  endif()
  if(SUNDIALS_GINKGO_BACKENDS MATCHES "REF")
    list(APPEND vectors nvecserial)
  endif()

  sundials_install_examples_ginkgo(sunmatrix
    CPU_GPU_EXAMPLES_VAR
      examples
    SUNDIALS_COMPONENTS
      sunmatrixdense
      ${vectors}
    SUNDIALS_TARGETS
      sunmatrixdense
    DEPENDENCIES
      ${SUNDIALS_SOURCE_DIR}/examples/sunmatrix/test_sunmatrix.c
      ${SUNDIALS_SOURCE_DIR}/examples/sunmatrix/test_sunmatrix.h
    DESTINATION
      sunmatrix/ginkgo
    )

endif()
