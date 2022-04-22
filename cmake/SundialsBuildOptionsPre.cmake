# ---------------------------------------------------------------
# Programmer(s): Cody J. Balos @ LLNL
# ---------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2022, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# ---------------------------------------------------------------
# SUNDIALS build options that are interepreted prior to any
# other CMake configuration.
# ---------------------------------------------------------------

# ---------------------------------------------------------------
# xSDK specific options and defaults
# ---------------------------------------------------------------

# always show the option to turn on xSDK defaults
sundials_option(USE_XSDK_DEFAULTS BOOL "Enable default xSDK settings" OFF)

if(USE_XSDK_DEFAULTS)
  message(STATUS "Enabling xSDK defaults:")

  # set the CMake build type, SUNDIALS does not set a build type by default
  if(NOT CMAKE_BUILD_TYPE)
    message(STATUS "  Setting build type to Debug")
    set(DOCSTR "Choose the type of build: None Debug Release RelWithDebInfo MinSizeRel")
    set(CMAKE_BUILD_TYPE "Debug" CACHE STRING "${DOCSTR}" FORCE)
  endif()
endif()

# ---------------------------------------------------------------
# Option to specify precision (realtype)
# ---------------------------------------------------------------

set(DOCSTR "single, double, or extended")
sundials_option(SUNDIALS_PRECISION STRING "${DOCSTR}" "DOUBLE")
string(TOUPPER ${SUNDIALS_PRECISION} _upper_SUNDIALS_PRECISION)
force_variable(SUNDIALS_PRECISION STRING "${DOCSTR}" ${_upper_SUNDIALS_PRECISION})

# ---------------------------------------------------------------
# Option to specify index type
# ---------------------------------------------------------------

# set the index size, SUNDIALS_INDEX_SIZE defaults to 64
set(DOCSTR "Signed 64-bit (64) or signed 32-bit (32) integer")
if(USE_XSDK_DEFAULTS)
  sundials_option(SUNDIALS_INDEX_SIZE STRING "${DOCSTR}" "32")
else()
  sundials_option(SUNDIALS_INDEX_SIZE STRING "${DOCSTR}" "64")
endif()

set(DOCSTR "Integer type to use for indices in SUNDIALS")
sundials_option(SUNDIALS_INDEX_TYPE STRING "${DOCSTR}" "" ADVANCED)

# ---------------------------------------------------------------
# Option to enable monitoring
# ---------------------------------------------------------------

set(DOCSTR "Build with simulation monitoring capabilities enabled")
sundials_option(SUNDIALS_BUILD_WITH_MONITORING BOOL "${DOCSTR}" OFF)

# ---------------------------------------------------------------
# Option to enable profiling
# ---------------------------------------------------------------

set(DOCSTR "Build with simulation profiling capabilities enabled")
sundials_option(SUNDIALS_BUILD_WITH_PROFILING BOOL "${DOCSTR}" OFF)

# ---------------------------------------------------------------
# Option to enable logging
# ---------------------------------------------------------------

set(DOCSTR "Build with logging capabilities enabled (0 = no logging, 1 = errors, 2 = +warnings, 3 = +info, 4 = +debug, 5 = +extras")
sundials_option(SUNDIALS_LOGGING_LEVEL STRING "${DOCSTR}" 0
                OPTIONS "0;1;2;3;4;5")

if(SUNDIALS_LOGGING_LEVEL GREATER_EQUAL 1)
  message(STATUS "SUNDIALS logging level set to ${SUNDIALS_LOGGING_LEVEL}")
  message(WARNING "SUNDIALS built with logging turned on, performance may be affected.")
endif()

set(DOCSTR "Build SUNDIALS logging with MPI support")
sundials_option(SUNDIALS_LOGGING_ENABLE_MPI BOOL "${DOCSTR}" "${ENABLE_MPI}"
                DEPENDS_ON ENABLE_MPI)

# ---------------------------------------------------------------
# Option to use the generic math libraries (UNIX only)
# ---------------------------------------------------------------

if(UNIX)
  sundials_option(USE_GENERIC_MATH BOOL "Use generic (std-c) math libraries" ON)
  # all executables will be linked against -lm
  set(EXTRA_LINK_LIBS -lm)
endif()

# ---------------------------------------------------------------
# Options to enable static and/or shared libraries
# ---------------------------------------------------------------

sundials_option(BUILD_STATIC_LIBS BOOL "Build static libraries" ON)
sundials_option(BUILD_SHARED_LIBS BOOL "Build shared libraries" ON)

# Make sure we build at least one type of libraries
if(NOT BUILD_STATIC_LIBS AND NOT BUILD_SHARED_LIBS)
  print_error("Both static and shared library generation were disabled.")
endif()

# ---------------------------------------------------------------
# Options to enable SUNDIALS packages and modules
# ---------------------------------------------------------------

# For each SUNDIALS package available (i.e. for which we have the
# sources), give the user the option of enabling/disabling it.

if(IS_DIRECTORY "${SUNDIALS_SOURCE_DIR}/src/arkode")
  sundials_option(BUILD_ARKODE BOOL "Build the ARKODE library" ON)
  list(APPEND SUNDIALS_BUILD_LIST "BUILD_ARKODE")
else()
  set(BUILD_ARKODE OFF)
endif()

if(IS_DIRECTORY "${SUNDIALS_SOURCE_DIR}/src/cvode")
  sundials_option(BUILD_CVODE BOOL "Build the CVODE library" ON)
  list(APPEND SUNDIALS_BUILD_LIST "BUILD_CVODE")
else()
  set(BUILD_CVODE OFF)
endif()

if(IS_DIRECTORY "${SUNDIALS_SOURCE_DIR}/src/cvodes")
  sundials_option(BUILD_CVODES BOOL "Build the CVODES library" ON)
  list(APPEND SUNDIALS_BUILD_LIST "BUILD_CVODES")
else()
  set(BUILD_CVODES OFF)
endif()

if(IS_DIRECTORY "${SUNDIALS_SOURCE_DIR}/src/ida")
  sundials_option(BUILD_IDA BOOL "Build the IDA library" ON)
  list(APPEND SUNDIALS_BUILD_LIST "BUILD_IDA")
else()
  set(BUILD_IDA OFF)
endif()

if(IS_DIRECTORY "${SUNDIALS_SOURCE_DIR}/src/idas")
  sundials_option(BUILD_IDAS BOOL "Build the IDAS library" ON)
  list(APPEND SUNDIALS_BUILD_LIST "BUILD_IDAS")
else()
  set(BUILD_IDAS OFF)
endif()

if(IS_DIRECTORY "${SUNDIALS_SOURCE_DIR}/src/kinsol")
  sundials_option(BUILD_KINSOL BOOL "Build the KINSOL library" ON)
  list(APPEND SUNDIALS_BUILD_LIST "BUILD_KINSOL")
else()
  set(BUILD_KINSOL OFF)
endif()

# ---------------------------------------------------------------
# Options to enable Fortran interfaces.
# ---------------------------------------------------------------

# Fortran 2003 interface is disabled by default
set(DOCSTR "Enable Fortran 2003 modules")
sundials_option(BUILD_FORTRAN_MODULE_INTERFACE BOOL "${DOCSTR}" OFF)

if(BUILD_FORTRAN_MODULE_INTERFACE)
  # F2003 interface only supports double precision
  if(NOT (SUNDIALS_PRECISION MATCHES "DOUBLE"))
    print_error("F2003 interface is not compatible with ${SUNDIALS_PRECISION} precision")
  endif()

  # F2003 interface only supports 64-bit indices
  if(NOT (SUNDIALS_INDEX_SIZE MATCHES "64"))
    print_error("F2003 interface is not compatible with ${SUNDIALS_INDEX_SIZE}-bit indicies")
  endif()

  # Allow a user to set where the Fortran modules will be installed
  set(DOCSTR "Directory where Fortran module files are installed")
  sundials_option(Fortran_INSTALL_MODDIR STRING "${DOCSTR}" "fortran")
endif()

# ---------------------------------------------------------------
# Options for benchmark suite
# ---------------------------------------------------------------

sundials_option(BUILD_BENCHMARKS BOOL "Build the SUNDIALS benchmark suite" OFF)

# ---------------------------------------------------------------
# Options for CMake config installation
# ---------------------------------------------------------------

set(DOCSTR "Path to SUNDIALS cmake files")
sundials_option(SUNDIALS_INSTALL_CMAKEDIR STRING "${DOCSTR}"
                "${CMAKE_INSTALL_LIBDIR}/cmake/sundials")

# ---------------------------------------------------------------
# Options to enable compiler warnings, address sanitizer
# ---------------------------------------------------------------

sundials_option(ENABLE_ALL_WARNINGS BOOL
  "Enable all compiler warnings" OFF ADVANCED)

sundials_option(ENABLE_WARNINGS_AS_ERRORS BOOL
  "Enable compiler warnings as errors" OFF ADVANCED)

sundials_option(ENABLE_ADDRESS_SANITIZER BOOL
  "Enable address sanitizer" OFF ADVANCED)

# ---------------------------------------------------------------
# Options to enable SUNDIALS debugging
# ---------------------------------------------------------------

# List of debugging options (used to add preprocessor directives)
set(_SUNDIALS_DEBUG_OPTIONS
  SUNDIALS_DEBUG
  SUNDIALS_DEBUG_ASSERT
  SUNDIALS_DEBUG_CUDA_LASTERROR
  SUNDIALS_DEBUG_HIP_LASTERROR
  SUNDIALS_DEBUG_PRINTVEC)

sundials_option(SUNDIALS_DEBUG BOOL
  "Enable additional debugging output and options" OFF
  ADVANCED)

if(SUNDIALS_DEBUG AND SUNDIALS_LOGGING_LEVEL LESS 4)
  set(DOCSTR "SUNDIALS_DEBUG=ON forced the logging level to 4")
  message(STATUS "${DOCSTR}")
  set(SUNDIALS_LOGGING_LEVEL "4" CACHE STRING "${DOCSTR}" FORCE)
endif()

sundials_option(SUNDIALS_DEBUG_ASSERT BOOL
  "Enable assert when debugging" OFF
  DEPENDS_ON SUNDIALS_DEBUG
  ADVANCED)

sundials_option(SUNDIALS_DEBUG_CUDA_LASTERROR BOOL
  "Enable CUDA last error checks when debugging" OFF
  DEPENDS_ON SUNDIALS_DEBUG ENABLE_CUDA
  ADVANCED)

sundials_option(SUNDIALS_DEBUG_HIP_LASTERROR BOOL
  "Enable HIP last error checks when debugging" OFF
  DEPENDS_ON SUNDIALS_DEBUG ENABLE_HIP
  ADVANCED)

sundials_option(SUNDIALS_DEBUG_PRINTVEC BOOL
  "Enable vector printing when debugging" OFF
  DEPENDS_ON SUNDIALS_DEBUG
  ADVANCED)

if(SUNDIALS_DEBUG_PRINTVEC AND SUNDIALS_LOGGING_LEVEL LESS 5)
  set(DOCSTR "SUNDIALS_DEBUG_PRINTVEC=ON forced the logging level to 5")
  message(STATUS "${DOCSTR}")
  set(SUNDIALS_LOGGING_LEVEL "5" CACHE STRING "${DOCSTR}" FORCE)
endif()

# ---------------------------------------------------------------
# Options for SUNDIALS testing
# ---------------------------------------------------------------

sundials_option(SUNDIALS_TEST_FLOAT_PRECISION STRING
  "Precision for floating point comparisons (number of digits)" "-1" ADVANCED)

sundials_option(SUNDIALS_TEST_INTEGER_PRECISION STRING
  "Precision for integer comparisons (percent difference)" "-1" ADVANCED)

sundials_option(SUNDIALS_TEST_OUTPUT_DIR PATH
  "Location to write testing output files" "" ADVANCED)

sundials_option(SUNDIALS_TEST_ANSWER_DIR PATH
  "Location of testing answer files" "" ADVANCED)
