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
# unit_tests/arkode level CMakeLists.txt for SUNDIALS
# ---------------------------------------------------------------

# C unit tests
add_subdirectory(C_serial)

# C++ unit tests
if(CXX_FOUND)
  add_subdirectory(CXX_serial)
  if(ENABLE_MPI AND MPI_CXX_FOUND)
    add_subdirectory(CXX_parallel)
  endif()
endif()

# F2003 unit tests
if(BUILD_FORTRAN_MODULE_INTERFACE)
  add_subdirectory(F2003_serial)
endif()

