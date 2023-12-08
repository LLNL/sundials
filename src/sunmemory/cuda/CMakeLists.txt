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

# Create a library out of the generic sundials modules
sundials_add_library(sundials_sunmemcuda
  SOURCES
    sundials_cuda_memory.cu
  HEADERS
    ${SUNDIALS_SOURCE_DIR}/include/sunmemory/sunmemory_cuda.h
  LINK_LIBRARIES
    PUBLIC sundials_core
  INCLUDE_SUBDIR
    sunmemory
  OBJECT_LIB_ONLY
)
