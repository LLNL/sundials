# ---------------------------------------------------------------
# $Revision: 4958 $
# $Date: 2016-09-23 14:02:13 -0700 (Fri, 23 Sep 2016) $
# ---------------------------------------------------------------
# Programmer:  Radu Serban @ LLNL
# ---------------------------------------------------------------
# LLNS Copyright Start
# Copyright (c) 2014, Lawrence Livermore National Security
# This work was performed under the auspices of the U.S. Department 
# of Energy by Lawrence Livermore National Laboratory in part under 
# Contract W-7405-Eng-48 and in part under Contract DE-AC52-07NA27344.
# Produced at the Lawrence Livermore National Laboratory.
# All rights reserved.
# For details, see the LICENSE file.
# LLNS Copyright End
# ---------------------------------------------------------------
# Top level CMakeLists.txt for SUNDIALS (for cmake build system)


# -------------------------------------------------------------
# Initial commands
# -------------------------------------------------------------

# Require a fairly recent cmake version

CMAKE_MINIMUM_REQUIRED(VERSION 2.8.1)

# Set CMake policy to allow examples to build
if(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

# Project SUNDIALS (initially only C supported)

PROJECT(sundials C)

# Set some variables with info on the SUNDIALS project

SET(PACKAGE_BUGREPORT "woodward6@llnl.gov")
SET(PACKAGE_NAME "SUNDIALS")
SET(PACKAGE_STRING "SUNDIALS 2.7.0")
SET(PACKAGE_TARNAME "sundials")
SET(PACKAGE_VERSION "2.7.0")

# 
SET_PROPERTY(GLOBAL PROPERTY USE_FOLDERS ON)

# Prohibit in-source build

IF("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_BINARY_DIR}")
  MESSAGE(FATAL_ERROR "In-source build prohibited.")
ENDIF("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_BINARY_DIR}")

# Hide/show some cache variables

MARK_AS_ADVANCED(EXECUTABLE_OUTPUT_PATH LIBRARY_OUTPUT_PATH)
MARK_AS_ADVANCED(CLEAR
  CMAKE_C_COMPILER
  CMAKE_C_FLAGS
  CMAKE_C_FLAGS_DEBUG
  CMAKE_C_FLAGS_MINSIZEREL
  CMAKE_C_FLAGS_RELEASE
  CMAKE_C_FLAGS_RELWITHDEB)

# Specify the VERSION and SOVERSION for shared libraries

SET(arkodelib_VERSION "1.1.0")
SET(arkodelib_SOVERSION "1")

SET(cvodelib_VERSION "2.9.0")
SET(cvodelib_SOVERSION "2")

SET(cvodeslib_VERSION "2.9.0")
SET(cvodeslib_SOVERSION "2")

SET(idalib_VERSION "2.9.0")
SET(idalib_SOVERSION "2")

SET(idaslib_VERSION "1.3.0")
SET(idaslib_SOVERSION "1")

SET(kinsollib_VERSION "2.9.0")
SET(kinsollib_SOVERSION "2")

SET(cpodeslib_VERSION "0.0.0")
SET(cpodeslib_SOVERSION "0")

SET(nveclib_VERSION "2.7.0")
SET(nveclib_SOVERSION "2")

# Specify the location of additional CMAKE modules

SET(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/config)

# -------------------------------------------------------------
# MACRO definitions
# -------------------------------------------------------------

# Macros to hide/show cached variables.
# These two macros can be used to "hide" or "show" in the
# list of cached variables various variables and/or options 
# that depend on other options.
# Note that once a variable is modified, it will preserve its
# value (hidding it merely makes it internal)

MACRO(HIDE_VARIABLE var)
  IF(DEFINED ${var})
    SET(${var} "${${var}}" CACHE INTERNAL "")
  ENDIF(DEFINED ${var})
ENDMACRO(HIDE_VARIABLE)

MACRO(SHOW_VARIABLE var type doc default)
  IF(DEFINED ${var})
    SET(${var} "${${var}}" CACHE "${type}" "${doc}" FORCE)
  ELSE(DEFINED ${var})
    SET(${var} "${default}" CACHE "${type}" "${doc}")
  ENDIF(DEFINED ${var})
ENDMACRO(SHOW_VARIABLE)

# Macros to append a common suffix or prefix to the elements of a list

MACRO(ADD_SUFFIX rootlist suffix)
  SET(outlist )
  FOREACH(root ${${rootlist}})
    LIST(APPEND outlist ${root}${suffix})
  ENDFOREACH(root)
  SET(${rootlist} ${outlist})
ENDMACRO(ADD_SUFFIX)

MACRO(ADD_PREFIX prefix rootlist)
  SET(outlist )
  FOREACH(root ${${rootlist}})
    LIST(APPEND outlist ${prefix}${root})
  ENDFOREACH(root)
  SET(${rootlist} ${outlist})
ENDMACRO(ADD_PREFIX)

# Macro to print warning that some features will be disabled
# due to some failure.

MACRO(PRINT_WARNING message action)
  MESSAGE("WARNING: ${message}.\n   ${action}.")
ENDMACRO(PRINT_WARNING)

# Returns an unquoted string.  Note that CMake will readily turn such
# strings back into lists, due to the duality of lists and
# semicolon-separated strings.  So be careful how you use it.

MACRO(LIST2STRING alist astring)
  FOREACH(elem ${${alist}})
   SET(${astring} "${${astring}} ${elem}")
  ENDFOREACH(elem)
ENDMACRO(LIST2STRING)

# -------------------------------------------------------------
# Which modules to build?
# -------------------------------------------------------------

# For each SUNDIALS solver available (i.e. for which we have the
# sources), give the user the option of enabling/disabling it.

IF(IS_DIRECTORY "${sundials_SOURCE_DIR}/src/arkode")
  OPTION(BUILD_ARKODE  "Build the ARKODE library"  ON)
ELSE(IS_DIRECTORY "${sundials_SOURCE_DIR}/src/arkode")
  SET(BUILD_ARKODE OFF)
ENDIF(IS_DIRECTORY "${sundials_SOURCE_DIR}/src/arkode")

IF(IS_DIRECTORY "${sundials_SOURCE_DIR}/src/cvode")
  OPTION(BUILD_CVODE  "Build the CVODE library"  ON)
ELSE(IS_DIRECTORY "${sundials_SOURCE_DIR}/src/cvode")
  SET(BUILD_CVODE OFF)
ENDIF(IS_DIRECTORY "${sundials_SOURCE_DIR}/src/cvode")

IF(IS_DIRECTORY "${sundials_SOURCE_DIR}/src/cvodes")
  OPTION(BUILD_CVODES  "Build the CVODES library"  ON)
ELSE(IS_DIRECTORY "${sundials_SOURCE_DIR}/src/cvodes")
  SET(BUILD_CVODES OFF)
ENDIF(IS_DIRECTORY "${sundials_SOURCE_DIR}/src/cvodes")

IF(IS_DIRECTORY "${sundials_SOURCE_DIR}/src/ida")
  OPTION(BUILD_IDA  "Build the IDA library"  ON)
ELSE(IS_DIRECTORY "${sundials_SOURCE_DIR}/src/ida")
  SET(BUILD_IDA OFF)
ENDIF(IS_DIRECTORY "${sundials_SOURCE_DIR}/src/ida")

IF(IS_DIRECTORY "${sundials_SOURCE_DIR}/src/idas")
  OPTION(BUILD_IDAS  "Build the IDAS library"  ON)
ELSE(IS_DIRECTORY "${sundials_SOURCE_DIR}/src/idas")
  SET(BUILD_IDAS OFF)
ENDIF(IS_DIRECTORY "${sundials_SOURCE_DIR}/src/idas")

IF(IS_DIRECTORY "${sundials_SOURCE_DIR}/src/kinsol")
  OPTION(BUILD_KINSOL  "Build the KINSOL library"  ON)
ELSE(IS_DIRECTORY "${sundials_SOURCE_DIR}/src/kinsol")
  SET(BUILD_KINSOL OFF)
ENDIF(IS_DIRECTORY "${sundials_SOURCE_DIR}/src/kinsol")

# CPODES is always OFF for now.  (commented out for Release); ToDo: better way to do this?
#IF(IS_DIRECTORY "${sundials_SOURCE_DIR}/src/cpodes")
#  OPTION(BUILD_CPODES  "Build the CPODES library"  OFF)
#ELSE(IS_DIRECTORY "${sundials_SOURCE_DIR}/src/cpodes")
#  SET(BUILD_CPODES OFF)
#ENDIF(IS_DIRECTORY "${sundials_SOURCE_DIR}/src/cpodes")

# -------------------------------------------------------------
# Other configuration options
# -------------------------------------------------------------

# Option that allows users to build static and/or shared libraries
# ----------------------------------------------------------------

OPTION(BUILD_STATIC_LIBS "Build static libraries" ON)
OPTION(BUILD_SHARED_LIBS "Build shared libraries" ON)

# Prepare substitution variable SUNDIALS_EXPORT for sundials_config.h
# When building shared SUNDIALS libraries under Windows, use
#      #define SUNDIALS_EXPORT __declspec(dllexport)
# When linking to shared SUNDIALS libraries under Windows, use
#      #define SUNDIALS_EXPORT __declspec(dllimport)
# In all other cases (other platforms or static libraries
# under Windows), the SUNDIALS_EXPORT macro is empty

IF(BUILD_SHARED_LIBS AND WIN32)
  SET(SUNDIALS_EXPORT
    "#ifdef BUILD_SUNDIALS_LIBRARY
#define SUNDIALS_EXPORT __declspec(dllexport)
#else
#define SUNDIALS_EXPORT __declspec(dllimport)
#endif")
ELSE(BUILD_SHARED_LIBS AND WIN32)
  SET(SUNDIALS_EXPORT "#define SUNDIALS_EXPORT")
ENDIF(BUILD_SHARED_LIBS AND WIN32)

# Make sure we build at least one type of libraries
IF(NOT BUILD_STATIC_LIBS AND NOT BUILD_SHARED_LIBS)
  PRINT_WARNING("Both static and shared library generation were disabled"
    "Building static libraries was re-enabled")
  SET(BUILD_STATIC_LIBS ON CACHE BOOL "Build static libraries" FORCE)
ENDIF(NOT BUILD_STATIC_LIBS AND NOT BUILD_SHARED_LIBS)

# Option to specify precision
# ---------------------------

SET(SUNDIALS_PRECISION "double" CACHE STRING "double, single or extended")

# prepare substitution variable PRECISION_LEVEL for sundials_config.h
STRING(TOUPPER ${SUNDIALS_PRECISION} SUNDIALS_PRECISION)
SET(PRECISION_LEVEL "#define SUNDIALS_${SUNDIALS_PRECISION}_PRECISION 1")

# Option to use the generic math libraries (UNIX only)
# ----------------------------------------------------

IF(UNIX)
  OPTION(USE_GENERIC_MATH "Use generic (std-c) math libraries" ON)
  IF(USE_GENERIC_MATH)
    # executables will be linked against -lm 
    SET(EXTRA_LINK_LIBS -lm)
    # prepare substitution variable for sundials_config.h
    SET(SUNDIALS_USE_GENERIC_MATH TRUE)
  ENDIF(USE_GENERIC_MATH)
ENDIF(UNIX)

## clock-monotonic, see if we need to link with rt
include(CheckSymbolExists)
set(CMAKE_REQUIRED_LIBRARIES_SAVE ${CMAKE_REQUIRED_LIBRARIES})
set(CMAKE_REQUIRED_LIBRARIES rt)
CHECK_SYMBOL_EXISTS(_POSIX_TIMERS "unistd.h;time.h" SUNDIALS_POSIX_TIMERS)
set(CMAKE_REQUIRED_LIBRARIES ${CMAKE_REQUIRED_LIBRARIES_SAVE})
if(SUNDIALS_POSIX_TIMERS)
  find_library(SUNDIALS_RT_LIBRARY NAMES rt)
  mark_as_advanced(SUNDIALS_RT_LIBRARY)
  if(SUNDIALS_RT_LIBRARY)
    # sundials_config.h symbol
    SET(SUNDIALS_HAVE_POSIX_TIMERS TRUE)
    set(EXTRA_LINK_LIBS ${EXTRA_LINK_LIBS} ${SUNDIALS_RT_LIBRARY})
  endif()
endif()

# -------------------------------------------------------------
# Enable Fortran support?
# -------------------------------------------------------------

# FCMIX support is an option only if at least one solver that 
# provides such an interface is built.

IF(BUILD_ARKODE OR BUILD_CVODE OR BUILD_IDA OR BUILD_KINSOL)
  SHOW_VARIABLE(FCMIX_ENABLE BOOL "Enable Fortran-C support" OFF)
ELSE(BUILD_ARKODE OR BUILD_CVODE OR BUILD_IDA OR BUILD_KINSOL)
  HIDE_VARIABLE(FCMIX_ENABLE)
ENDIF(BUILD_ARKODE OR BUILD_CVODE OR BUILD_IDA OR BUILD_KINSOL)

# -------------------------------------------------------------
# Enable Fortran90 support?
# -------------------------------------------------------------

# F90 support is an option only if ARKode is built and 
# examples are enabled

IF(BUILD_ARKODE AND EXAMPLES_ENABLE AND FCMIX_ENABLE)
  SHOW_VARIABLE(F90_ENABLE BOOL "Enable Fortran90 ARKode examples" OFF)
ELSE(BUILD_ARKODE AND EXAMPLES_ENABLE AND FCMIX_ENABLE)
  # # set back to OFF (in case was ON)
  # SET(F90_ENABLE OFF)
  # HIDE_VARIABLE(F90_ENABLE AND FCMIX_ENABLE)
  HIDE_VARIABLE(F90_ENABLE)
ENDIF(BUILD_ARKODE AND EXAMPLES_ENABLE AND FCMIX_ENABLE)

# -------------------------------------------------------------
# Enable C++ support?
# -------------------------------------------------------------

# C++ support is an option only if ARKode is built and 
# examples are enabled

IF(BUILD_ARKODE AND EXAMPLES_ENABLE)
  SHOW_VARIABLE(CXX_ENABLE BOOL "Enable C++ ARKode examples" OFF)
ELSE(BUILD_ARKODE AND EXAMPLES_ENABLE)
  HIDE_VARIABLE(CXX_ENABLE)
ENDIF(BUILD_ARKODE AND EXAMPLES_ENABLE)

# -------------------------------------------------------------
# Enable BLAS/LAPACK support?
# -------------------------------------------------------------

OPTION(LAPACK_ENABLE "Enable Lapack support" OFF)

IF(NOT LAPACK_ENABLE)
  HIDE_VARIABLE(SUNDIALS_F77_FUNC_CASE)
  HIDE_VARIABLE(SUNDIALS_F77_FUNC_UNDERSCORES)
  HIDE_VARIABLE(LAPACK_LIBRARIES)
ENDIF(NOT LAPACK_ENABLE)

# -------------------------------------------------------------
# Enable SUPERLUMT support?
# -------------------------------------------------------------

OPTION(SUPERLUMT_ENABLE "Enable SUPERLUMT support" OFF)

IF(NOT SUPERLUMT_ENABLE)
  HIDE_VARIABLE(SUPERLUMT_THREAD_TYPE)
  HIDE_VARIABLE(SUPERLUMT_LIBRARY_DIR)
  HIDE_VARIABLE(SUPERLUMT_INCLUDE_DIR)
  SET (SUPERLUMT_DISABLED TRUE CACHE INTERNAL "GUI - return when first set")
ENDIF(NOT SUPERLUMT_ENABLE)

# -------------------------------------------------------------
# Enable KLU support?
# -------------------------------------------------------------

OPTION(KLU_ENABLE "Enable KLU support" OFF)

IF(NOT KLU_ENABLE)
  HIDE_VARIABLE(KLU_LIBRARY_DIR)
  HIDE_VARIABLE(KLU_INCLUDE_DIR)
  SET (KLU_DISABLED TRUE CACHE INTERNAL "GUI - return when first set")
ENDIF(NOT KLU_ENABLE)

# -------------------------------------------------------------
# Enable MPI support?
# -------------------------------------------------------------

OPTION(MPI_ENABLE "Enable MPI support" OFF)

IF(NOT MPI_ENABLE)
  HIDE_VARIABLE(MPI_INCLUDE_PATH)
  HIDE_VARIABLE(MPI_LIBRARIES)
  HIDE_VARIABLE(MPI_EXTRA_LIBRARIES)
  HIDE_VARIABLE(MPI_MPICC)
  HIDE_VARIABLE(MPI_MPICXX)
  HIDE_VARIABLE(MPI_MPIF77)
  HIDE_VARIABLE(MPI_MPIF90)
ENDIF(NOT MPI_ENABLE)

# -------------------------------------------------------------
# Enable HYPRE Vector support?
# -------------------------------------------------------------

OPTION(HYPRE_ENABLE "Enable Hypre support" OFF)
IF(NOT MPI_ENABLE)
  HIDE_VARIABLE(HYPRE_INCLUDE_DIR)
  HIDE_VARIABLE(HYPRE_LIBRARY_DIR)
  SET (HYPRE_DISABLED TRUE CACHE INTERNAL "GUI - return when first set")
ENDIF(NOT MPI_ENABLE)

# -------------------------------------------------------------
# Enable OpenMP support?
# -------------------------------------------------------------

OPTION(OPENMP_ENABLE "Enable OpenMP support" OFF)

IF(OPENMP_ENABLE)
  FIND_PACKAGE(OpenMP)
  IF(NOT OPENMP_FOUND)
    message(STATUS "Disabling OpenMP support, could not determine compiler flags")
  ENDIF(NOT OPENMP_FOUND)
ENDIF(OPENMP_ENABLE)

# -------------------------------------------------------------
# Enable Pthread support?
# -------------------------------------------------------------

OPTION(PTHREAD_ENABLE "Enable Pthreads support" OFF)

IF(PTHREAD_ENABLE)
  FIND_PACKAGE(Threads)
  IF(CMAKE_USE_PTHREADS_INIT)
    message(STATUS "Using Pthreads")
    SET(PTHREADS_FOUND TRUE)
    # SGS    
  ELSE()
    message(STATUS "Disabling Pthreads support, could not determine compiler flags")
  endif()
ENDIF(PTHREAD_ENABLE)

# -------------------------------------------------------------
# Enable PETSc support?
# -------------------------------------------------------------

OPTION(PETSC_ENABLE "Enable PETSC support" OFF)

IF(NOT PETSC_ENABLE)
  HIDE_VARIABLE(PETSC_LIBRARY_DIR)
  HIDE_VARIABLE(PETSC_INCLUDE_DIR)
  SET (PETSC_DISABLED TRUE CACHE INTERNAL "GUI - return when first set")
ENDIF(NOT PETSC_ENABLE)

# -------------------------------------------------------------
# Enable examples?
# -------------------------------------------------------------

OPTION(EXAMPLES_ENABLE "Build the SUNDIALS examples" ON)

IF(EXAMPLES_ENABLE)

  # If examples are enabled, set different options

  # The examples will be linked with the library corresponding to the build type.
  # Whenever building shared libraries, use them to link the examples.
  IF(BUILD_SHARED_LIBS)
    SET(LINK_LIBRARY_TYPE "shared")
  ELSE(BUILD_SHARED_LIBS)
    SET(LINK_LIBRARY_TYPE "static")
  ENDIF(BUILD_SHARED_LIBS)

  # Check if example files are to be exported
  SHOW_VARIABLE(EXAMPLES_INSTALL BOOL "Install example files" ON)

  # If examples are to be exported, check where we should install them.
  IF(EXAMPLES_INSTALL)

    SHOW_VARIABLE(EXAMPLES_INSTALL_PATH PATH
      "Output directory for installing example files" "${CMAKE_INSTALL_PREFIX}/examples")
    
    IF(NOT EXAMPLES_INSTALL_PATH)
      PRINT_WARNING("The example installation path is empty"
        "Example installation path was reset to its default value")
      SET(EXAMPLES_INSTALL_PATH "${CMAKE_INSTALL_PREFIX}/examples" CACHE STRING
        "Output directory for installing example files" FORCE)
    ENDIF(NOT EXAMPLES_INSTALL_PATH)

  ELSE(EXAMPLES_INSTALL)

    HIDE_VARIABLE(EXAMPLES_INSTALL_PATH)

  ENDIF(EXAMPLES_INSTALL)

ELSE(EXAMPLES_ENABLE)

  # If examples are disabled, hide all options related to
  # building and installing the SUNDIALS examples

  HIDE_VARIABLE(EXAMPLES_INSTALL)
  HIDE_VARIABLE(EXAMPLES_INSTALL_PATH)
  
ENDIF(EXAMPLES_ENABLE)

# -------------------------------------------------------------
# Add any other necessary compiler flags & definitions
# -------------------------------------------------------------

# Under Windows, add compiler directive to inhibit warnings
# about use of unsecure functions

IF(WIN32)
  ADD_DEFINITIONS(-D_CRT_SECURE_NO_WARNINGS)
ENDIF(WIN32)

IF(APPLE)
  SET(CMAKE_SHARED_LIBRARY_CREATE_C_FLAGS "${CMAKE_SHARED_LIBRARY_CREATE_C_FLAGS} -undefined dynamic_lookup")
ENDIF(APPLE)

# -------------------------------------------------------------
# A Fortran compiler is needed if:
# (a) FCMIX is enabled
# (b) LAPACK is enabled (for the name-mangling scheme)
# -------------------------------------------------------------

IF(FCMIX_ENABLE OR LAPACK_ENABLE)
  INCLUDE(SundialsFortran)
  IF(NOT F77_FOUND AND FCMIX_ENABLE)
    PRINT_WARNING("Fortran compiler not functional"
      "FCMIX support will not be provided")
  ENDIF(NOT F77_FOUND AND FCMIX_ENABLE)
ENDIF(FCMIX_ENABLE OR LAPACK_ENABLE)

# -------------------------------------------------------------
# A Fortran90 compiler is needed if:
# (a) F90 is enabled, and
# (b) ARKODE examples are enabled
# -------------------------------------------------------------

IF(F90_ENABLE AND BUILD_ARKODE)
  INCLUDE(SundialsFortran90)
  IF(NOT F90_FOUND AND F90_ENABLE)
    PRINT_WARNING("Fortran90 compiler not functional"
      "F90 support will not be provided")
  ENDIF(NOT F90_FOUND AND F90_ENABLE)
ENDIF(F90_ENABLE AND BUILD_ARKODE)

# -------------------------------------------------------------
# A C++ compiler is needed if:
# (a) C++ is enabled, and
# (b) ARKODE examples are enabled
# -------------------------------------------------------------

IF(CXX_ENABLE AND BUILD_ARKODE)
  INCLUDE(SundialsCXX)
  IF(NOT CXX_FOUND AND CXX_ENABLE)
    PRINT_WARNING("C++ compiler not functional"
      "C++ support will not be provided")
  ENDIF(NOT CXX_FOUND AND CXX_ENABLE)
ENDIF(CXX_ENABLE AND BUILD_ARKODE)

# -------------------------------------------------------------
# Check if we need an alternate way of specifying the Fortran
# name-mangling scheme if we were unable to infer it using a
# compiler. 
# Ask the user to specify the case and number of appended underscores
# corresponding to the Fortran name-mangling scheme of symbol names 
# that do not themselves contain underscores (recall that this is all
# we really need for the interfaces to LAPACK).
# Note: the default scheme is lower case - one underscore
# -------------------------------------------------------------

IF(LAPACK_ENABLE AND NOT F77SCHEME_FOUND)
  # Specify the case for the Fortran name-mangling scheme
  SHOW_VARIABLE(SUNDIALS_F77_FUNC_CASE STRING
    "case of Fortran function names (lower/upper)"
    "lower")
  # Specify the number of appended underscores for the Fortran name-mangling scheme
  SHOW_VARIABLE(SUNDIALS_F77_FUNC_UNDERSCORES STRING 
    "number of underscores appended to Fortran function names"
    "one")  
  # Based on the given case and number of underscores,
  # set the C preprocessor macro definition
  IF(${SUNDIALS_F77_FUNC_CASE} MATCHES "lower")
    IF(${SUNDIALS_F77_FUNC_UNDERSCORES} MATCHES "none")
      SET(CMAKE_Fortran_SCHEME_NO_UNDERSCORES "mysub")
    ENDIF(${SUNDIALS_F77_FUNC_UNDERSCORES} MATCHES "none")
    IF(${SUNDIALS_F77_FUNC_UNDERSCORES} MATCHES "one")
      SET(CMAKE_Fortran_SCHEME_NO_UNDERSCORES "mysub_")
    ENDIF(${SUNDIALS_F77_FUNC_UNDERSCORES} MATCHES "one")
    IF(${SUNDIALS_F77_FUNC_UNDERSCORES} MATCHES "two")
      SET(CMAKE_Fortran_SCHEME_NO_UNDERSCORES "mysub__")
    ENDIF(${SUNDIALS_F77_FUNC_UNDERSCORES} MATCHES "two")
  ELSE(${SUNDIALS_F77_FUNC_CASE} MATCHES "lower")
    IF(${SUNDIALS_F77_FUNC_UNDERSCORES} MATCHES "none")
      SET(CMAKE_Fortran_SCHEME_NO_UNDERSCORES "MYSUB")
    ENDIF(${SUNDIALS_F77_FUNC_UNDERSCORES} MATCHES "none")
    IF(${SUNDIALS_F77_FUNC_UNDERSCORES} MATCHES "one")
      SET(CMAKE_Fortran_SCHEME_NO_UNDERSCORES "MYSUB_")
    ENDIF(${SUNDIALS_F77_FUNC_UNDERSCORES} MATCHES "one")
    IF(${SUNDIALS_F77_FUNC_UNDERSCORES} MATCHES "two")
      SET(CMAKE_Fortran_SCHEME_NO_UNDERSCORES "MYSUB__")
    ENDIF(${SUNDIALS_F77_FUNC_UNDERSCORES} MATCHES "two")
  ENDIF(${SUNDIALS_F77_FUNC_CASE} MATCHES "lower")
  # Since the SUNDIALS codes never use symbol names containing
  # underscores, set a default scheme (probably wrong) for symbols
  # with underscores.
  SET(CMAKE_Fortran_SCHEME_WITH_UNDERSCORES "my_sub_")
  # We now "have" a scheme.
  SET(F77SCHEME_FOUND TRUE)
ENDIF(LAPACK_ENABLE AND NOT F77SCHEME_FOUND)

# -------------------------------------------------------------
# If we have a name-mangling scheme (either automatically
# inferred or provided by the user), set the SUNDIALS 
# compiler preprocessor macro definitions.
# -------------------------------------------------------------

SET(F77_MANGLE_MACRO1 "")
SET(F77_MANGLE_MACRO2 "")

IF(F77SCHEME_FOUND)
  # Symbols WITHOUT underscores
  IF(${CMAKE_Fortran_SCHEME_NO_UNDERSCORES} MATCHES "mysub")
    SET(F77_MANGLE_MACRO1 "#define SUNDIALS_F77_FUNC(name,NAME) name")
  ENDIF(${CMAKE_Fortran_SCHEME_NO_UNDERSCORES} MATCHES "mysub")
  IF(${CMAKE_Fortran_SCHEME_NO_UNDERSCORES} MATCHES "mysub_")
    SET(F77_MANGLE_MACRO1 "#define SUNDIALS_F77_FUNC(name,NAME) name ## _")
  ENDIF(${CMAKE_Fortran_SCHEME_NO_UNDERSCORES} MATCHES "mysub_")
  IF(${CMAKE_Fortran_SCHEME_NO_UNDERSCORES} MATCHES "mysub__")
    SET(F77_MANGLE_MACRO1 "#define SUNDIALS_F77_FUNC(name,NAME) name ## __")
  ENDIF(${CMAKE_Fortran_SCHEME_NO_UNDERSCORES} MATCHES "mysub__")
  IF(${CMAKE_Fortran_SCHEME_NO_UNDERSCORES} MATCHES "MYSUB")
    SET(F77_MANGLE_MACRO1 "#define SUNDIALS_F77_FUNC(name,NAME) NAME")
  ENDIF(${CMAKE_Fortran_SCHEME_NO_UNDERSCORES} MATCHES "MYSUB")
  IF(${CMAKE_Fortran_SCHEME_NO_UNDERSCORES} MATCHES "MYSUB_")
    SET(F77_MANGLE_MACRO1 "#define SUNDIALS_F77_FUNC(name,NAME) NAME ## _")
  ENDIF(${CMAKE_Fortran_SCHEME_NO_UNDERSCORES} MATCHES "MYSUB_")
  IF(${CMAKE_Fortran_SCHEME_NO_UNDERSCORES} MATCHES "MYSUB__")
    SET(F77_MANGLE_MACRO1 "#define SUNDIALS_F77_FUNC(name,NAME) NAME ## __")
  ENDIF(${CMAKE_Fortran_SCHEME_NO_UNDERSCORES} MATCHES "MYSUB__")
  # Symbols with underscores 
  IF(${CMAKE_Fortran_SCHEME_NO_UNDERSCORES} MATCHES "my_sub")
    SET(F77_MANGLE_MACRO2 "#define SUNDIALS_F77_FUNC_(name,NAME) name")
  ENDIF(${CMAKE_Fortran_SCHEME_NO_UNDERSCORES} MATCHES "my_sub")
  IF(${CMAKE_Fortran_SCHEME_NO_UNDERSCORES} MATCHES "my_sub_")
    SET(F77_MANGLE_MACRO2 "#define SUNDIALS_F77_FUNC_(name,NAME) name ## _")
  ENDIF(${CMAKE_Fortran_SCHEME_NO_UNDERSCORES} MATCHES "my_sub_")
  IF(${CMAKE_Fortran_SCHEME_NO_UNDERSCORES} MATCHES "my_sub__")
    SET(F77_MANGLE_MACRO2 "#define SUNDIALS_F77_FUNC_(name,NAME) name ## __")
  ENDIF(${CMAKE_Fortran_SCHEME_NO_UNDERSCORES} MATCHES "my_sub__")
  IF(${CMAKE_Fortran_SCHEME_NO_UNDERSCORES} MATCHES "MY_SUB")
    SET(F77_MANGLE_MACRO2 "#define SUNDIALS_F77_FUNC_(name,NAME) NAME")
  ENDIF(${CMAKE_Fortran_SCHEME_NO_UNDERSCORES} MATCHES "MY_SUB")
  IF(${CMAKE_Fortran_SCHEME_NO_UNDERSCORES} MATCHES "MY_SUB_")
    SET(F77_MANGLE_MACRO2 "#define SUNDIALS_F77_FUNC_(name,NAME) NAME ## _")
  ENDIF(${CMAKE_Fortran_SCHEME_NO_UNDERSCORES} MATCHES "MY_SUB_")
  IF(${CMAKE_Fortran_SCHEME_NO_UNDERSCORES} MATCHES "MY_SUB__")
    SET(F77_MANGLE_MACRO2 "#define SUNDIALS_F77_FUNC_(name,NAME) NAME ## __")
  ENDIF(${CMAKE_Fortran_SCHEME_NO_UNDERSCORES} MATCHES "MY_SUB__")
ENDIF(F77SCHEME_FOUND)

# -------------------------------------------------------------
# Find (and test) the Lapack libraries
# -------------------------------------------------------------

# If LAPACK is needed, first try to find the appropriate
# libraries and linker flags needed to link against them.

IF(LAPACK_ENABLE)

  # find BLAS and LAPACK Libraries
  INCLUDE(SundialsLapack)

  SHOW_VARIABLE(LAPACK_LIBRARIES STRING "Blas and Lapack libraries" "${LAPACK_LIBRARIES}")

  IF(LAPACK_LIBRARIES AND NOT LAPACK_FOUND)
    PRINT_WARNING("LAPACK not functional"
      "Blas/Lapack support will not be provided")
  ELSE(LAPACK_LIBRARIES AND NOT LAPACK_FOUND)
    #set sundials_config.h symbol via sundials_config.in
    SET(SUNDIALS_BLAS_LAPACK TRUE)
  ENDIF(LAPACK_LIBRARIES AND NOT LAPACK_FOUND)
  

ELSE(LAPACK_ENABLE)

  HIDE_VARIABLE(LAPACK_LIBRARIES)

ENDIF(LAPACK_ENABLE)

# -------------------------------------------------------------
# Find (and test) the SUPERLUMT libraries
# -------------------------------------------------------------

# If SUPERLUMT is needed, first try to find the appropriate
# libraries to link against them.

IF(SUPERLUMT_ENABLE)

  # get thread type for SUPERLUMT (OpenMP or Pthreads(default))
  SET(SUPERLUMT_THREAD_TYPE "Pthread" CACHE STRING "OpenMP or Pthread")
  SHOW_VARIABLE(SUPERLUMT_THREAD_TYPE STRING "SUPERLUMT threading type: OpenMP or Pthread" "${SUPERLUMT_THREAD_TYPE}")
  SHOW_VARIABLE(SUPERLUMT_INCLUDE_DIR PATH "SUPERLUMT include directory" "${SUPERLUMT_INCLUDE_DIR}")
  SHOW_VARIABLE(SUPERLUMT_LIBRARY_DIR PATH "SUPERLUMT library directory" "${SUPERLUMT_LIBRARY_DIR}")

  INCLUDE(SundialsSuperLUMT)

  IF(SUPERLUMT_FOUND)
    # sundials_config.h symbols
    SET(SUNDIALS_SUPERLUMT TRUE)
    SET(SUNDIALS_SUPERLUMT_THREAD_TYPE ${SUPERLUMT_THREAD_TYPE})
  ENDIF(SUPERLUMT_FOUND)

  IF(SUPERLUMT_LIBRARIES AND NOT SUPERLUMT_FOUND)
    PRINT_WARNING("SUPERLUMT not functional - support will not be provided"
      "Double check spelling specified libraries (search is case sensitive)")
  ENDIF(SUPERLUMT_LIBRARIES AND NOT SUPERLUMT_FOUND)

ELSE(SUPERLUMT_ENABLE)

  HIDE_VARIABLE(SUPERLUMT_THREAD_TYPE)
  HIDE_VARIABLE(SUPERLUMT_INCLUDE_DIR)
  HIDE_VARIABLE(SUPERLUMT_LIBRARY_DIR)

ENDIF(SUPERLUMT_ENABLE)

# -------------------------------------------------------------
# Find (and test) the KLU libraries
# -------------------------------------------------------------

# If KLU is requested, first try to find the appropriate libraries to
# link against them.

IF(KLU_ENABLE)

  SHOW_VARIABLE(KLU_INCLUDE_DIR PATH "KLU include directory"
    "${KLU_INCLUDE_DIR}") 
  SHOW_VARIABLE(KLU_LIBRARY_DIR PATH
    "Klu library directory" "${KLU_LIBRARY_DIR}")

  INCLUDE(SundialsKLU)

  IF(KLU_FOUND) 
    # sundials_config.h symbol
    SET(SUNDIALS_KLU TRUE)
    INCLUDE_DIRECTORIES(${KLU_INCLUDE_DIR})
  ENDIF(KLU_FOUND)

  IF(KLU_LIBRARIES AND NOT KLU_FOUND) 
    PRINT_WARNING("KLU not functional - support will not be provided" 
      "Double check spelling of include path and specified libraries (search is case sensitive)")
  ENDIF(KLU_LIBRARIES AND NOT KLU_FOUND)

ENDIF(KLU_ENABLE)

# -------------------------------------------------------------
# Decide how to compile MPI codes.
# -------------------------------------------------------------

IF(MPI_ENABLE)
  # get run command for mpi (openmpi or OpenMP or Pthreads(default))
  SET(MPI_RUN_COMMAND "mpirun" CACHE STRING "mpirun or srun")
  SHOW_VARIABLE(MPI_RUN_COMMAND STRING "MPI run command: mpirun or srun" "${MPI_RUN_COMMAND}")
  INCLUDE(SundialsMPIC)
  IF(NOT MPIC_FOUND)
    PRINT_WARNING("MPI not functional"
      "Parallel support will not be provided")
  ENDIF(NOT MPIC_FOUND)
  IF(MPIC_MPI2)
    SET(F77_MPI_COMM_F2C "#define SUNDIALS_MPI_COMM_F2C 1")
  ELSE(MPIC_MPI2)
    SET(F77_MPI_COMM_F2C "#define SUNDIALS_MPI_COMM_F2C 0")
  ENDIF(MPIC_MPI2)
  IF(MPIC_FOUND AND FCMIX_ENABLE)
    INCLUDE(SundialsMPIF)
  ENDIF(MPIC_FOUND AND FCMIX_ENABLE)
  IF(MPIC_FOUND AND CXX_ENABLE)
    INCLUDE(SundialsMPICXX)
  ENDIF(MPIC_FOUND AND CXX_ENABLE)
  IF(MPIC_FOUND AND F90_ENABLE)
    INCLUDE(SundialsMPIF90)
  ENDIF(MPIC_FOUND AND F90_ENABLE)
ENDIF(MPI_ENABLE)

#-------------------------------------------------------------
# Find the Hypre libraries
#-------------------------------------------------------------
IF(HYPRE_ENABLE AND NOT MPI_ENABLE)
    PRINT_WARNING("MPI not enabled - HYPRE support will not be provided" 
      "Set MPI_ENABLE to ON to use parhyp")
ENDIF(HYPRE_ENABLE AND NOT MPI_ENABLE)

IF(HYPRE_ENABLE AND MPI_ENABLE)
  SHOW_VARIABLE(HYPRE_INCLUDE_DIR PATH "HYPRE include directory"
    "${HYPRE_INCLUDE_DIR}") 
  SHOW_VARIABLE(HYPRE_LIBRARY_DIR PATH
    "HYPRE library directory" "${HYPRE_LIBRARY_DIR}")

  INCLUDE(SundialsHypre)

  IF(HYPRE_FOUND) 
    # sundials_config.h symbol
    SET(SUNDIALS_HYPRE TRUE)
    INCLUDE_DIRECTORIES(${HYPRE_INCLUDE_DIR})
  ENDIF(HYPRE_FOUND)

  IF(HYPRE_LIBRARIES AND NOT HYPRE_FOUND) 
    PRINT_WARNING("HYPRE not functional - support will not be provided" 
      "Found hypre library, test code does not work")
  ENDIF(HYPRE_LIBRARIES AND NOT HYPRE_FOUND)
  
ENDIF(HYPRE_ENABLE AND MPI_ENABLE)

# -------------------------------------------------------------
# Find (and test) the PETSc libraries
# -------------------------------------------------------------
IF(PETSC_ENABLE AND NOT MPI_ENABLE)
    PRINT_WARNING("MPI not enabled - PETSc support will not be provided" 
      "Set MPI_ENABLE to ON to use PETSc.")
ENDIF(PETSC_ENABLE AND NOT MPI_ENABLE)

IF(PETSC_ENABLE AND MPI_ENABLE)
  SHOW_VARIABLE(PETSC_INCLUDE_DIR PATH "PETSc include directory"
    "${PETSC_INCLUDE_DIR}") 
  SHOW_VARIABLE(PETSC_LIBRARY_DIR PATH
    "PETSc library directory" "${PETSC_LIBRARY_DIR}")

  INCLUDE(SundialsPETSc)

  IF(PETSC_FOUND) 
    # sundials_config.h symbol
    SET(SUNDIALS_PETSC TRUE)
    INCLUDE_DIRECTORIES(${PETSC_INCLUDE_DIR})
  ENDIF(PETSC_FOUND)

  IF(PETSC_LIBRARIES AND NOT PETSC_FOUND) 
    PRINT_WARNING("PETSC not functional - support will not be provided" 
      "Double check spelling specified libraries (search is case sensitive)")
  ENDIF(PETSC_LIBRARIES AND NOT PETSC_FOUND)

ENDIF(PETSC_ENABLE AND MPI_ENABLE)


# -------------------------------------------------------------
# If using MPI with C++, disable C++ extensions (for known wrappers)
# -------------------------------------------------------------

# IF(MPICXX_FOUND)
#   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DMPICH_SKIP_MPICXX -DOMPI_SKIP_MPICXX -DLAM_BUILDING")
# ENDIF(MPICXX_FOUND)

# -------------------------------------------------------------
# Configure the header file sundials_config.h
# -------------------------------------------------------------

# All required substitution variables should be available at this point.
# Generate the header file and place it in the binary dir.
CONFIGURE_FILE(
  ${PROJECT_SOURCE_DIR}/include/sundials/sundials_config.in
  ${PROJECT_BINARY_DIR}/include/sundials/sundials_config.h
  )

# Add the include directory in the source tree and the one in
# the binary tree (for the header file sundials_config.h)
IF(SUPERLUMT_FOUND)
  INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/include ${PROJECT_BINARY_DIR}/include ${SUPERLUMT_INCLUDE_DIR})
ELSE(SUPERLUMT_FOUND)
  INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/include ${PROJECT_BINARY_DIR}/include)
ENDIF(SUPERLUMT_FOUND)

# -------------------------------------------------------------
# Add selected modules to the build system
# -------------------------------------------------------------

# Shared components

ADD_SUBDIRECTORY(src/sundials)
ADD_SUBDIRECTORY(src/nvec_ser)
IF(MPIC_FOUND)
  ADD_SUBDIRECTORY(src/nvec_par)
ENDIF(MPIC_FOUND)

IF(SUNDIALS_HYPRE)
  ADD_SUBDIRECTORY(src/nvec_parhyp)
ENDIF(SUNDIALS_HYPRE)

IF(OPENMP_FOUND)
  ADD_SUBDIRECTORY(src/nvec_openmp)	
ENDIF(OPENMP_FOUND)

IF(PTHREADS_FOUND)
  ADD_SUBDIRECTORY(src/nvec_pthreads)   
ENDIF(PTHREADS_FOUND)

IF(PETSC_FOUND)
  ADD_SUBDIRECTORY(src/nvec_petsc)   
ENDIF(PETSC_FOUND)

# ARKODE library

IF(BUILD_ARKODE)	
  ADD_SUBDIRECTORY(src/arkode)
  IF(FCMIX_ENABLE AND F77_FOUND)
    ADD_SUBDIRECTORY(src/arkode/fcmix)
  ENDIF(FCMIX_ENABLE AND F77_FOUND)
ENDIF(BUILD_ARKODE)

# CVODE library

IF(BUILD_CVODE)	
  ADD_SUBDIRECTORY(src/cvode)
  IF(FCMIX_ENABLE AND F77_FOUND)
    ADD_SUBDIRECTORY(src/cvode/fcmix)
  ENDIF(FCMIX_ENABLE AND F77_FOUND)
ENDIF(BUILD_CVODE)

# CVODES library

IF(BUILD_CVODES)	
  ADD_SUBDIRECTORY(src/cvodes)
ENDIF(BUILD_CVODES)

# IDA library

IF(BUILD_IDA)	
  ADD_SUBDIRECTORY(src/ida)
  IF(FCMIX_ENABLE AND F77_FOUND)
    ADD_SUBDIRECTORY(src/ida/fcmix)
  ENDIF(FCMIX_ENABLE AND F77_FOUND)
ENDIF(BUILD_IDA)

# IDAS library

IF(BUILD_IDAS)	
  ADD_SUBDIRECTORY(src/idas)
ENDIF(BUILD_IDAS)

# KINSOL library

IF(BUILD_KINSOL)	
  ADD_SUBDIRECTORY(src/kinsol)
  IF(FCMIX_ENABLE AND F77_FOUND)
    ADD_SUBDIRECTORY(src/kinsol/fcmix)
  ENDIF(FCMIX_ENABLE AND F77_FOUND)
ENDIF(BUILD_KINSOL)

# CPODES library

IF(BUILD_CPODES)	
  ADD_SUBDIRECTORY(src/cpodes)
ENDIF(BUILD_CPODES)

# -------------------------------------------------------------
# Include the subdirectories corresponding to various examples
# -------------------------------------------------------------

# If building and installing the examples is enabled, include
# the subdirectories for those examples that will be built.
# Also, if we will generate exported example Makefiles, set 
# variables needed in generating them from templates.

# For now, TestRunner is not being distributed.
# So:
#  - Don't show TESTRUNNER variable
#  - Don't enable testing if TestRunner if not found.
#  - There will be no 'make test' target

INCLUDE(SundialsAddTest)
HIDE_VARIABLE(TESTRUNNER)
IF(EXAMPLES_ENABLE)
  IF(TESTRUNNER)
    ENABLE_TESTING()
  ENDIF(TESTRUNNER)

  IF(EXAMPLES_INSTALL)
    SET(SHELL "sh")
    SET(prefix "${CMAKE_INSTALL_PREFIX}")
    SET(exec_prefix "${CMAKE_INSTALL_PREFIX}")
    SET(includedir "${prefix}/include")
    SET(libdir "${exec_prefix}/lib")
    SET(CPP "${CMAKE_C_COMPILER}")
    SET(CC "${CMAKE_C_COMPILER}")
    SET(CPPFLAGS "${CMAKE_C_FLAGS_RELEASE}")
    SET(CFLAGS "${CMAKE_C_FLAGS_RELEASE}")
    SET(LDFLAGS "${CMAKE_EXE_LINKER_FLAGS_RELEASE}")
    LIST2STRING(EXTRA_LINK_LIBS LIBS)
    IF(CXX_FOUND)
      SET(CXX "${CMAKE_CXX_COMPILER}")
      SET(CXX_LNKR "${CMAKE_CXX_COMPILER}")
      SET(CXXFLAGS "${CMAKE_CXX_FLAGS_RELEASE}")
      SET(CXX_LDFLAGS "${CMAKE_CXX_FLAGS_RELEASE}")
      LIST2STRING(EXTRA_LINK_LIBS CXX_LIBS)
    ENDIF(CXX_FOUND)
    IF(F77_FOUND)
      SET(F77 "${CMAKE_Fortran_COMPILER}")
      SET(F77_LNKR "${CMAKE_Fortran_COMPILER}")
      SET(FFLAGS "${CMAKE_Fortran_FLAGS_RELEASE}")
      SET(F77_LDFLAGS "${CMAKE_Fortran_FLAGS_RELEASE}")
      LIST2STRING(EXTRA_LINK_LIBS F77_LIBS)
    ENDIF(F77_FOUND)
    IF(F90_FOUND)
      SET(F90 "${CMAKE_Fortran_COMPILER}")
      SET(F90_LNKR "${CMAKE_Fortran_COMPILER}")
      SET(F90FLAGS "${CMAKE_Fortran_FLAGS_RELEASE}")
      SET(F90_LDFLAGS "${CMAKE_Fortran_FLAGS_RELEASE}")
      LIST2STRING(EXTRA_LINK_LIBS F90_LIBS)
    ENDIF(F90_FOUND)
    IF(SUPERLUMT_FOUND)
      LIST2STRING(SUPERLUMT_LIBRARIES SUPERLUMT_LIBS)
      SET(SUPERLUMT_LIBS "${SUPERLUMT_LINKER_FLAGS} ${SUPERLUMT_LIBS}")
    ENDIF(SUPERLUMT_FOUND)
    IF(KLU_FOUND)
      LIST2STRING(KLU_LIBRARIES KLU_LIBS)
      SET(KLU_LIBS "${KLU_LINKER_FLAGS} ${KLU_LIBS}")
    ENDIF(KLU_FOUND)
    IF(LAPACK_FOUND)
      LIST2STRING(LAPACK_LIBRARIES BLAS_LAPACK_LIBS)
    ENDIF(LAPACK_FOUND)
    IF(MPIC_FOUND)
      IF(MPI_MPICC)
        SET(MPICC "${MPI_MPICC}")
        SET(MPI_INC_DIR ".")
        SET(MPI_LIB_DIR ".")
        SET(MPI_LIBS "")
        SET(MPI_FLAGS "")
      ELSE(MPI_MPICC)
        SET(MPICC "${CMAKE_C_COMPILER}")
        SET(MPI_INC_DIR "${MPI_INCLUDE_PATH}")
        SET(MPI_LIB_DIR ".")
        LIST2STRING(MPI_LIBRARIES MPI_LIBS)
      ENDIF(MPI_MPICC)
      SET(HYPRE_INC_DIR "${HYPRE_INCLUDE_DIR}")
      SET(HYPRE_LIB_DIR "${HYPRE_LIBRARY_DIR}")
      SET(HYPRE_LIBS "${HYPRE_LIBRARIES}")
    ENDIF(MPIC_FOUND)
    IF(MPICXX_FOUND)
      IF(MPI_MPICXX)
        SET(MPICXX "${MPI_MPICXX}")
      ELSE(MPI_MPICXX)
        SET(MPICXX "${CMAKE_CXX_COMPILER}")
        LIST2STRING(MPI_LIBRARIES MPI_LIBS)
      ENDIF(MPI_MPICXX)
    ENDIF(MPICXX_FOUND)
    IF(MPIF_FOUND)
      IF(MPI_MPIF77)
        SET(MPIF77 "${MPI_MPIF77}")
        SET(MPIF77_LNKR "${MPI_MPIF77}")
      ELSE(MPI_MPIF77)
        SET(MPIF77 "${CMAKE_Fortran_COMPILER}")
        SET(MPIF77_LNKR "${CMAKE_Fortran_COMPILER}")
        SET(MPI_INC_DIR "${MPI_INCLUDE_PATH}")
        SET(MPI_LIB_DIR ".")
        LIST2STRING(MPI_LIBRARIES MPI_LIBS)
      ENDIF(MPI_MPIF77)
    ENDIF(MPIF_FOUND)
    IF(MPIF90_FOUND)
      IF(MPI_MPIF90)
        SET(MPIF90 "${MPI_MPIF90}")
        SET(MPIF90_LNKR "${MPI_MPIF90}")
      ELSE(MPI_MPIF90)
        SET(MPIF90 "${CMAKE_Fortran_COMPILER}")
        SET(MPIF90_LNKR "${CMAKE_Fortran_COMPILER}")
        LIST2STRING(MPI_LIBRARIES MPI_LIBS)
      ENDIF(MPI_MPIF90)
    ENDIF(MPIF90_FOUND)
  ENDIF(EXAMPLES_INSTALL)

  IF(BUILD_ARKODE)
    ADD_SUBDIRECTORY(examples/arkode/C_serial)
    IF(OPENMP_FOUND)
      ADD_SUBDIRECTORY(examples/arkode/C_openmp)
    ENDIF(OPENMP_FOUND)
    IF(CXX_ENABLE AND CXX_FOUND)
      ADD_SUBDIRECTORY(examples/arkode/CXX_serial)
    ENDIF(CXX_ENABLE AND CXX_FOUND)
    IF(FCMIX_ENABLE AND F77_FOUND)
      ADD_SUBDIRECTORY(examples/arkode/F77_serial)
    ENDIF(FCMIX_ENABLE AND F77_FOUND)
    IF(F90_ENABLE AND F90_FOUND)
      ADD_SUBDIRECTORY(examples/arkode/F90_serial)
    ENDIF(F90_ENABLE AND F90_FOUND)
    IF(MPIC_FOUND)
      ADD_SUBDIRECTORY(examples/arkode/C_parallel)
    ENDIF(MPIC_FOUND)
    IF(MPICXX_FOUND)
      ADD_SUBDIRECTORY(examples/arkode/CXX_parallel)
    ENDIF(MPICXX_FOUND)
    IF(MPIF_FOUND)
      ADD_SUBDIRECTORY(examples/arkode/F77_parallel)
    ENDIF(MPIF_FOUND)
    IF(MPIF90_FOUND)
      ADD_SUBDIRECTORY(examples/arkode/F90_parallel)
    ENDIF(MPIF90_FOUND)
    IF(HYPRE_ENABLE AND HYPRE_FOUND)
      ADD_SUBDIRECTORY(examples/arkode/C_parhyp)
    ENDIF(HYPRE_ENABLE AND HYPRE_FOUND)
    
  ENDIF(BUILD_ARKODE)
  
  IF(BUILD_CVODE)
    ADD_SUBDIRECTORY(examples/cvode/serial)
    IF(FCMIX_ENABLE AND F77_FOUND)
      ADD_SUBDIRECTORY(examples/cvode/fcmix_serial)
    ENDIF(FCMIX_ENABLE AND F77_FOUND)
    IF(MPIC_FOUND)
      ADD_SUBDIRECTORY(examples/cvode/parallel)
    ENDIF(MPIC_FOUND)
    IF(MPIF_FOUND)
      ADD_SUBDIRECTORY(examples/cvode/fcmix_parallel)
    ENDIF(MPIF_FOUND)
    IF(OPENMP_FOUND)
      ADD_SUBDIRECTORY(examples/cvode/C_openmp)
    ENDIF(OPENMP_FOUND)
    IF(HYPRE_ENABLE AND HYPRE_FOUND)
      ADD_SUBDIRECTORY(examples/cvode/parhyp)
    ENDIF(HYPRE_ENABLE AND HYPRE_FOUND)
  ENDIF(BUILD_CVODE)
  
  IF(BUILD_CVODES)	
    ADD_SUBDIRECTORY(examples/cvodes/serial)
    IF(MPIC_FOUND)
      ADD_SUBDIRECTORY(examples/cvodes/parallel)
    ENDIF(MPIC_FOUND)
    IF(OPENMP_FOUND)
      ADD_SUBDIRECTORY(examples/cvodes/C_openmp)
    ENDIF(OPENMP_FOUND)
  ENDIF(BUILD_CVODES)
  
  IF(BUILD_IDA)
    ADD_SUBDIRECTORY(examples/ida/serial)
    IF(FCMIX_ENABLE AND F77_FOUND)
      ADD_SUBDIRECTORY(examples/ida/fcmix_serial)
      IF(OPENMP_FOUND)
        ADD_SUBDIRECTORY(examples/ida/fcmix_openmp)
      ENDIF(OPENMP_FOUND)
      IF(PTHREADS_FOUND)
	ADD_SUBDIRECTORY(examples/ida/fcmix_pthreads)
      ENDIF(PTHREADS_FOUND)
    ENDIF(FCMIX_ENABLE AND F77_FOUND)
    IF(OPENMP_FOUND)
      ADD_SUBDIRECTORY(examples/ida/C_openmp)
    ENDIF(OPENMP_FOUND)
    IF(MPIC_FOUND)
      ADD_SUBDIRECTORY(examples/ida/parallel)
    ENDIF(MPIC_FOUND)
    IF(PETSC_FOUND)
      ADD_SUBDIRECTORY(examples/ida/petsc)
    ENDIF(PETSC_FOUND)
    IF(MPIF_FOUND)
      ADD_SUBDIRECTORY(examples/ida/fcmix_parallel)
    ENDIF(MPIF_FOUND)
  ENDIF(BUILD_IDA)
  
  IF(BUILD_IDAS)	
    ADD_SUBDIRECTORY(examples/idas/serial)
    IF(OPENMP_FOUND)
      ADD_SUBDIRECTORY(examples/idas/C_openmp)
    ENDIF(OPENMP_FOUND)
    IF(MPIC_FOUND)
      ADD_SUBDIRECTORY(examples/idas/parallel)
    ENDIF(MPIC_FOUND)
  ENDIF(BUILD_IDAS)

  IF(BUILD_KINSOL)	
    ADD_SUBDIRECTORY(examples/kinsol/serial)
    IF(OPENMP_FOUND)
      ADD_SUBDIRECTORY(examples/kinsol/C_openmp) # the only example here need special handling from testrunner (not yet implemented)
    ENDIF(OPENMP_FOUND)
    IF(FCMIX_ENABLE AND F77_FOUND)
      ADD_SUBDIRECTORY(examples/kinsol/fcmix_serial)
    ENDIF(FCMIX_ENABLE AND F77_FOUND)
    IF(MPIC_FOUND)
      ADD_SUBDIRECTORY(examples/kinsol/parallel)
    ENDIF(MPIC_FOUND)
    IF(MPIF_FOUND)
      ADD_SUBDIRECTORY(examples/kinsol/fcmix_parallel)
    ENDIF(MPIF_FOUND)
  ENDIF(BUILD_KINSOL)

  IF(BUILD_CPODES)	
    ADD_SUBDIRECTORY(examples/cpodes/serial)
    IF(MPIC_FOUND)
      ADD_SUBDIRECTORY(examples/cpodes/parallel)
    ENDIF(MPIC_FOUND)
  ENDIF(BUILD_CPODES)

  # Always add the nvector serial examples
  ADD_SUBDIRECTORY(examples/nvector/serial)

  IF(MPIC_FOUND)
      ADD_SUBDIRECTORY(examples/nvector/parallel)
  ENDIF(MPIC_FOUND)
  
  IF(SUNDIALS_HYPRE)
      ADD_SUBDIRECTORY(examples/nvector/parhyp)
  ENDIF(SUNDIALS_HYPRE)

  IF(PTHREADS_FOUND)
      ADD_SUBDIRECTORY(examples/nvector/pthreads)
  ENDIF(PTHREADS_FOUND)

  IF(OPENMP_FOUND)
      ADD_SUBDIRECTORY(examples/nvector/C_openmp)
  ENDIF(OPENMP_FOUND)
  
  IF(PETSC_FOUND)
      ADD_SUBDIRECTORY(examples/nvector/petsc)
  ENDIF(PETSC_FOUND)
  
ENDIF(EXAMPLES_ENABLE)

#----------------------------------
# Install configuration header file
#----------------------------------

# install configured header file
INSTALL(
  FILES ${PROJECT_BINARY_DIR}/include/sundials/sundials_config.h
  DESTINATION include/sundials
  )

