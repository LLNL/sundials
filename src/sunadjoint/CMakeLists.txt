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

# Create a library out of the generic sundials modules
sundials_add_library(
  sundials_adjoint
  SOURCES sunadjoint_stepper.c sunadjoint_checkpointscheme.c
          sunadjoint_checkpointscheme_basic.c
  HEADERS
    ${SUNDIALS_SOURCE_DIR}/include/sunadjoint/sunadjoint_stepper.h
    ${SUNDIALS_SOURCE_DIR}/include/sunadjoint/sunadjoint_checkpointscheme.h
    ${SUNDIALS_SOURCE_DIR}/include/sunadjoint/sunadjoint_checkpointscheme_basic.h
  LINK_LIBRARIES PUBLIC sundials_core
  OBJECT_LIBRARIES sundials_sunmemsys_obj
  INCLUDE_SUBDIR sunadjoint)

# Add F2003 module if the interface is enabled
if(BUILD_FORTRAN_MODULE_INTERFACE)
  add_subdirectory("fmod_int${SUNDIALS_INDEX_SIZE}")
endif()
