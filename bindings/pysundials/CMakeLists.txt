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

# Warn if the user invokes CMake directly
if(NOT SKBUILD)
  message(
    WARNING
      "\
  This CMake file is meant to be executed using 'scikit-build-core'.
  Running it directly will almost certainly not produce the desired
  result. If you are a user trying to install this package, use the
  command below, which will install all necessary build dependencies,
  compile the package in an isolated environment, and then install it.
  =====================================================================
   $ pip install .
  =====================================================================
  If you are a software developer, and this is your own package, then
  it is usually much more efficient to install the build dependencies
  in your environment once and use the following command that avoids
  a costly creation of a new virtual environment at every compilation:
  =====================================================================
   $ pip install nanobind scikit-build-core[pyproject]
   $ pip install --no-build-isolation -ve .
  =====================================================================
  You may optionally add -Ceditable.rebuild=true to auto-rebuild when
  the package is imported. Otherwise, you need to rerun the above
  after editing C++ files.")
endif()

# nanobind needs the Python Interpreter and Development components
find_package(
  Python 3.8
  COMPONENTS Interpreter Development.Module
  OPTIONAL_COMPONENTS Development.SABIModule
  REQUIRED)

# Determine location of nanobind cmake config file
execute_process(
  COMMAND "${Python_EXECUTABLE}" -m nanobind --cmake_dir
  OUTPUT_STRIP_TRAILING_WHITESPACE
  OUTPUT_VARIABLE nanobind_ROOT)

# nanobind must already be installed as a Python module (e.g., with pip)
find_package(nanobind CONFIG REQUIRED)

# Add the source files for the bindings
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
    nvector/nvector_serial.cpp
    sunlinsol/sunlinsol_spgmr.cpp
    sunlinsol/sunlinsol_dense.cpp
    sunlinsol/sunlinsol_band.cpp
    sunlinsol/sunlinsol_spbcgs.cpp
    sunlinsol/sunlinsol_spfgmr.cpp
    sunlinsol/sunlinsol_sptfqmr.cpp
    sunlinsol/sunlinsol_pcg.cpp
    arkode/pysundials_arkode.cpp
    arkode/pysundials_arkode_arkstep.cpp
    arkode/pysundials_arkode_erkstep.cpp
    arkode/pysundials_arkode_sprkstep.cpp
    arkode/pysundials_arkode_lsrkstep.cpp
    arkode/pysundials_arkode_mristep.cpp
    arkode/pysundials_arkode_forcingstep.cpp
    arkode/pysundials_arkode_splittingstep.cpp
    cvodes/pysundials_cvodes.cpp
    idas/pysundials_idas.cpp)

# Create the Python sundials library
nanobind_add_module(pysundials STABLE_ABI NB_STATIC ${sundials_SOURCES})

# Include private header locations
target_include_directories(pysundials PRIVATE helpers
                                              ${SUNDIALS_SOURCE_DIR}/src
                                              ${SUNDIALS_SOURCE_DIR}/src/arkode)

# Link against sundials libraries
target_link_libraries(pysundials PRIVATE sundials_arkode sundials_cvodes sundials_idas
                                         sundials_nvecserial sundials_core)

# Install directive for scikit-build-core
install(TARGETS pysundials LIBRARY DESTINATION .)

# # Setup tests find_program(PYTEST_EXECUTABLE pytest)

# # Copy tests to the build directory add_custom_target( copy-python-tests ALL
# COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_LIST_DIR}/tests
# ${CMAKE_CURRENT_BINARY_DIR})

# add_dependencies(pysundials copy-python-tests)

# add_test(NAME python-tests COMMAND ${Python_EXECUTABLE} -m pytest
# WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} )
