# ------------------------------------------------------------------------------
# Programmer(s): Daniel R. Reynolds @ SMU
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
# controller level CMakeLists.txt for SUNDIALS
# ------------------------------------------------------------------------------

# required native matrices
#add_subdirectory(exp_gus)
add_subdirectory(fixed)
add_subdirectory(i_control)
#add_subdirectory(imex_gus)
#add_subdirectory(imp_gus)
add_subdirectory(pi_control)
add_subdirectory(pid_control)
