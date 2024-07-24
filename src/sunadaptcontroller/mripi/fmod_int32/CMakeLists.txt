# ---------------------------------------------------------------
# Programmer(s): Daniel R. Reynolds @ SMU
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

sundials_add_f2003_library(sundials_fsunadaptcontrollermripi_mod
  SOURCES
    fsunadaptcontroller_mripi_mod.f90 fsunadaptcontroller_mripi_mod.c
  LINK_LIBRARIES
    PUBLIC sundials_fcore_mod
  OBJECT_LIBRARIES
  OUTPUT_NAME
    sundials_fsunadaptcontrollermripi_mod
  OBJECT_LIB_ONLY
)

message(STATUS "Added SUNAdaptController_MRIPI F2003 interface")
