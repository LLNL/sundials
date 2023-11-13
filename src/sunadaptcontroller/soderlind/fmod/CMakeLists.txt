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

sundials_add_f2003_library(sundials_fsunadaptcontrollersoderlind_mod
  SOURCES
    fsunadaptcontroller_soderlind_mod.f90 fsunadaptcontroller_soderlind_mod.c
  OBJECT_LIBRARIES
    sundials_fcore_mod_obj
  OUTPUT_NAME
    sundials_fsunadaptcontrollersoderlind_mod
  OBJECT_LIB_ONLY
)

message(STATUS "Added SUNAdaptController_Soderlind F2003 interface")
