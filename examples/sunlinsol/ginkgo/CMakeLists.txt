# ------------------------------------------------------------------------------
# Programmer(s): Cody J. Balos @ LLNL
# ------------------------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2023, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# ------------------------------------------------------------------------------

# Example lists are tuples "name\;args\;type" where the type is
# 'develop' for examples excluded from 'make test' in releases
set(examples
  "test_sunlinsol_ginkgo.cpp\;bicg csr 100 0 100 0\;"
#  "test_sunlinsol_ginkgo.cpp\;bicg dense 100 0 100 0\;"
  "test_sunlinsol_ginkgo.cpp\;bicgstab csr 100 0 100 0\;"
#  "test_sunlinsol_ginkgo.cpp\;bicgstab dense 100 0 100 0\;"
  "test_sunlinsol_ginkgo.cpp\;cgs csr 100 0 100 0\;"
#  "test_sunlinsol_ginkgo.cpp\;cgs dense 100 0 100 0\;"
  "test_sunlinsol_ginkgo.cpp\;gmres csr 100 0 100 0\;"
#  "test_sunlinsol_ginkgo.cpp\;gmres dense 100 0 100 0\;"
  "test_sunlinsol_ginkgo.cpp\;idr csr 100 0 100 0\;"
#  "test_sunlinsol_ginkgo.cpp\;idr dense 100 0 100 0\;"
)

# Add source directory to include directories
include_directories(..)

sundials_add_examples_ginkgo(examples
  TARGETS test_sunlinsol_obj
  BACKENDS REF OMP CUDA HIP
  UNIT_TEST)

# Install the targets
if(EXAMPLES_INSTALL)

  if(SUNDIALS_GINKGO_BACKENDS MATCHES "CUDA")
    list(APPEND vectors nveccuda)
  endif()
  if(SUNDIALS_GINKGO_BACKENDS MATCHES "HIP")
    list(APPEND vectors nvechip)
  endif()
  if(SUNDIALS_GINKGO_BACKENDS MATCHES "OMP")
    list(APPEND vectors nvecopenmp)
  endif()
  if(SUNDIALS_GINKGO_BACKENDS MATCHES "REF")
    list(APPEND vectors nvecserial)
  endif()

  sundials_install_examples_ginkgo(sunlinsol
    CPU_GPU_EXAMPLES_VAR
      examples
    SUNDIALS_COMPONENTS
      ${vectors}
    DEPENDENCIES
      ${SUNDIALS_SOURCE_DIR}/examples/sunlinsol/test_sunlinsol.c
      ${SUNDIALS_SOURCE_DIR}/examples/sunlinsol/test_sunlinsol.h
    DESTINATION
      sunlinsol/ginkgo
    )

endif()
