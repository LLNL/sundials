# ------------------------------------------------------------------------------
# Programmer(s): Daniel R. Reynolds @ SMU
# ------------------------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2025, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# ------------------------------------------------------------------------------
# adaptivity controller level CMakeLists.txt for SUNDIALS
# ------------------------------------------------------------------------------

# required native adaptivity controllers
add_subdirectory(imexgus)
add_subdirectory(soderlind)
add_subdirectory(mrihtol)
