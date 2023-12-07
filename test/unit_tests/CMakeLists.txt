# ---------------------------------------------------------------
# Programmer(s): David J. Gardner @ LLNL
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
# unit_tests level CMakeLists.txt for SUNDIALS
# ---------------------------------------------------------------

message("Adding units tests")

add_subdirectory(sundials)

if(BUILD_ARKODE)
  add_subdirectory(arkode)
endif()

if(BUILD_CVODE)
  add_subdirectory(cvode)
endif()

if(BUILD_CVODES)
  add_subdirectory(cvodes)
endif()

if(BUILD_IDA)
  add_subdirectory(ida)
endif()

if(BUILD_IDAS)
  add_subdirectory(idas)
endif()

if(BUILD_KINSOL)
  add_subdirectory(kinsol)
endif()

if(CXX_FOUND)
  add_subdirectory(sunmemory)
  if(SUNDIALS_BUILD_WITH_PROFILING)
    add_subdirectory(profiling)
  endif()
endif()
