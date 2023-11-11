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

# Create a library out of the generic sundials modules
sundials_add_library(sundials_sunadaptcontrollersoderlind
  SOURCES
    sunadaptcontroller_soderlind.c
  HEADERS
    ${SUNDIALS_SOURCE_DIR}/include/sunadaptcontroller/sunadaptcontroller_soderlind.h
  INCLUDE_SUBDIR
    sunadaptcontroller
  OBJECT_LIB_ONLY
)

# Add F2003 module if the interface is enabled
if(BUILD_FORTRAN_MODULE_INTERFACE)
  add_subdirectory(fmod)
endif()
