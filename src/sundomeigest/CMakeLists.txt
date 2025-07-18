# ------------------------------------------------------------------------------
# Programmer(s): Mustafa Aggul @ SMU
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
# dominant eigenvalue estimator level CMakeLists.txt for SUNDIALS
# ------------------------------------------------------------------------------

# required native dominant eigenvalue estimators
add_subdirectory(Power)

if(ENABLE_LAPACK)
  add_subdirectory(Arnoldi)
endif()
