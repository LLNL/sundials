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

# Setup nanobind
find_package(
  Python 3.8
  COMPONENTS Interpreter Development.Module
  REQUIRED)
execute_process(
  COMMAND "${Python_EXECUTABLE}" -m nanobind --cmake_dir
  OUTPUT_STRIP_TRAILING_WHITESPACE
  OUTPUT_VARIABLE nanobind_ROOT)
find_package(nanobind CONFIG REQUIRED)

# Add the source file
set(sundials_SOURCES
    pysundials.cpp
    sundials/sundials_core.cpp
    sundials/sundials_adaptcontroller.cpp
    # sundials/sundials_adjointcheckpointscheme.cpp
    sundials/sundials_adjointstepper.cpp
    sundials/sundials_matrix.cpp
    # sundials/sundials_memory.cpp
    sundials/sundials_nonlinearsolver.cpp
    sundials/sundials_linearsolver.cpp
    sundials/sundials_nvector.cpp
    sundials/sundials_stepper.cpp
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
