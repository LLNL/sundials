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

sundials_add_f2003_library(
  sundials_fsunadjoint_mod
  SOURCES fsunadjointstepper_mod.f90
          fsunadjointstepper_mod.c
          fsunadjointcheckpointscheme_mod.f90
          fsunadjointcheckpointscheme_mod.c
          fsunadjointcheckpointscheme_basic_mod.f90
          fsunadjointcheckpointscheme_basic_mod.c
  LINK_LIBRARIES PUBLIC sundials_fcore_mod
  OUTPUT_NAME sundials_fsunadjoint_mod OBJECT_LIB_ONLY)

message(STATUS "Added SUNAdjoint interface")
