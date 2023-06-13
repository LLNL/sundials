# ---------------------------------------------------------------
# Programmer(s): Daniel R. Reynolds @ SMU
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
sundials_add_library(sundials_sunheuristicsdefault
  SOURCES
    sunheuristics_default.c
  HEADERS
    ${SUNDIALS_SOURCE_DIR}/include/sunheuristics/sunheuristics_default.h
  INCLUDE_SUBDIR
    sunheuristics
  OBJECT_LIB_ONLY
)
