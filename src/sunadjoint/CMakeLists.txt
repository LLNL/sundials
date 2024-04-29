# ---------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2024, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# ---------------------------------------------------------------

# Create a library out of the generic sundials modules
sundials_add_library(sundials_adjoint
  SOURCES
    sunadjoint_solver.c
  HEADERS
    ${SUNDIALS_SOURCE_DIR}/include/sunadjoint/sunadjoint_solver.h
  LINK_LIBRARIES
    PUBLIC sundials_core
  INCLUDE_SUBDIR
    sunadjoint
)
