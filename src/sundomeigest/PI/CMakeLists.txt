# ---------------------------------------------------------------
# Programmer(s): Mustafa Aggul @ SMU
# ---------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2025, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# ---------------------------------------------------------------
# CMakeLists.txt file for the PI SUNDomEigEstimator library
# ---------------------------------------------------------------

install(CODE "MESSAGE(\"\nInstall SUNDOMEIGEST_PI\n\")")

# Add the sundomeigest_pi library
sundials_add_library(
  sundials_sundomeigestpi
  SOURCES sundomeigest_pi.c
  HEADERS ${SUNDIALS_SOURCE_DIR}/include/sundomeigest/sundomeigest_pi.h
  INCLUDE_SUBDIR sundomeigest
  LINK_LIBRARIES PUBLIC sundials_core
  OBJECT_LIBRARIES
  OUTPUT_NAME sundials_sundomeigestpi
  VERSION ${sundomeigestlib_VERSION}
  SOVERSION ${sundomeigestlib_SOVERSION})

message(STATUS "Added SUNDOMEIGEST_PI module")
