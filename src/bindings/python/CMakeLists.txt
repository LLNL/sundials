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

# Add the source file
set(sundials_SOURCES
    pysundials.cpp
    sundials/pysundials_adaptcontroller.cpp
    sundials/pysundials_adjointcheckpointscheme.cpp
    sundials/pysundials_adjointstepper.cpp
    sundials/pysundials_context.cpp
    sundials/pysundials_core.cpp
    sundials/pysundials_linearsolver.cpp
    sundials/pysundials_logger.cpp
    sundials/pysundials_matrix.cpp
    sundials/pysundials_memory.cpp
    sundials/pysundials_nonlinearsolver.cpp
    sundials/pysundials_nvector.cpp
    sundials/pysundials_profiler.cpp
    sundials/pysundials_stepper.cpp
    arkode/arkode.cpp
    arkode/arkode_arkstep.cpp
    arkode/arkode_erkstep.cpp
    nvector/nvector_serial.cpp
    sunlinsol/sunlinsol_spgmr.cpp
    sunlinsol/sunlinsol_dense.cpp
    sunlinsol/sunlinsol_band.cpp
    sunlinsol/sunlinsol_spbcgs.cpp
    sunlinsol/sunlinsol_spfgmr.cpp
    sunlinsol/sunlinsol_sptfqmr.cpp
    sunlinsol/sunlinsol_pcg.cpp
)

# Create the Python sundials library
nanobind_add_module(pysundials ${sundials_SOURCES})

# Include private header locations
target_include_directories(pysundials PRIVATE ${SUNDIALS_SOURCE_DIR}/src
                                              ${SUNDIALS_SOURCE_DIR}/src/arkode)

# Link against sundials libraries
target_link_libraries(pysundials PRIVATE sundials_arkode sundials_nvecserial
                                         sundials_core)

# Copy tests to the build directory
add_custom_target(
  copy-python-tests ALL
  COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_LIST_DIR}/tests
          ${CMAKE_CURRENT_BINARY_DIR})

add_dependencies(pysundials copy-python-tests)

add_test(NAME python-tests
  COMMAND ${Python_EXECUTABLE} -m pytest
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)
