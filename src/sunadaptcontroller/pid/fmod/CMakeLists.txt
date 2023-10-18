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

sundials_add_f2003_library(sundials_fsunadaptcontrollerpid_mod
  SOURCES
    fsunadaptcontroller_pid_mod.f90 fsunadaptcontroller_pid_mod.c
  OBJECT_LIBRARIES
    sundials_fgeneric_mod_obj
  OUTPUT_NAME
    sundials_fsunadaptcontrollerpid_mod
  OBJECT_LIB_ONLY
)

message(STATUS "Added SUNAdaptController_PID F2003 interface")
