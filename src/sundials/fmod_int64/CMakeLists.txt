# ---------------------------------------------------------------
# Programmer(s): Cody J. Balos @ LLNL
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
# CMakeLists.txt file for the F2003 SUNDIALS object library
# ---------------------------------------------------------------

set(sundials_SOURCES fsundials_core_mod.c fsundials_core_mod.f90)

sundials_add_f2003_library(sundials_fcore_mod SOURCES ${sundials_SOURCES})
