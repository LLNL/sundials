# ---------------------------------------------------------------
# Programmer:  Cody Balos @ LLNL
# ---------------------------------------------------------------
# Based on the FindPETSC module by Jed Brown.
# ---------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2019, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# ---------------------------------------------------------------
# Copyright Jed Brown
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions
# are met:
#
# * Redistributions of source code must retain the above copyright
#   notice, this list of conditions and the following disclaimer.
#
# * Redistributions in binary form must reproduce the above copyright
#   notice, this list of conditions and the following disclaimer in
#   the documentation and/or other materials provided with the
#   distribution.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
# FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
# COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
# INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
# BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS
# OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR
# TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
# USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# ---------------------------------------------------------------
# Try to find PETSC.
# Once done this will define:
#
#  PETSC::ALL         - a CMake target for all of PETSc
#  PETSC::SYS         - a CMake target for the main PETSc library
#  PETSC::VEC         - a CMake target for the PETSc vector library
#  PETSC::MAT         - a CMake target for the PETSc matrix library
#  PETSC::DM          - a CMake target for the PETSc DM library
#  PETSC::KSP         - a CMake target for the PETSc KSP library
#  PETSC::SNES        - a CMake target for the PETSc SNES library
#  PETSC::TS          - a CMake target for the PETSc TS library
#  PETSC_INDEX_SIZE   - (internal) the size of indices in PETSC
#  PETSC_PRECISION    - (internal) the real type precision in PETSC
#  PETSC_FOUND        - (internal) system has PETSC
#  PETSC_INCLUDES     - (advanced) the PETSC include directories
#  PETSC_LIBRARIES    - (advanced) Link these to use PETSC
#  PETSC_COMPILER     - (advanced) Compiler used by PETSC, helpful to find a compatible MPI
#  PETSC_DEFINITIONS  - (advanced) Compiler switches for using PETSC
#  PETSC_MPIEXEC      - (advanced) Executable for running MPI programs
#  PETSC_VERSION      - (internal) Version string (MAJOR.MINOR.SUBMINOR)
#
# Usage:
#  find_package(PETSC COMPONENTS CXX)  - required if build --with-clanguage=C++ --with-c-support=0
#  find_package(PETSC COMPONENTS C)    - standard behavior of checking build using a C compiler
#  find_package(PETSC)                 - same as above
#
# Setting these changes the behavior of the search
#  PETSC_DIR             - directory in which PETSC resides
#  PETSC_ARCH            - build architecture
#  PETSC_CURRENT         - (advanced) redo the executable tests
#  PETSC_EXECUTABLE_RUNS - (advanced) set to ON to ignore the output of the executable tests (not recommended)
#
# Redistribution and use is allowed according to the terms of the BSD license.
# ---------------------------------------------------------------

# cmake_policy(VERSION 3.3)

#-------- helper macros and functions

function (petsc_get_version)
  if (EXISTS "${PETSC_DIR}/include/petscversion.h")
    file (STRINGS "${PETSC_DIR}/include/petscversion.h" vstrings REGEX "#define PETSC_VERSION_(RELEASE|MAJOR|MINOR|SUBMINOR|PATCH) ")
    foreach (line ${vstrings})
      string (REGEX REPLACE " +" ";" fields ${line}) # break line into three fields (the first is always "#define")
      list (GET fields 1 var)
      list (GET fields 2 val)
      set (${var} ${val} PARENT_SCOPE)
      set (${var} ${val})         # Also in local scope so we have access below
    endforeach ()
    if (PETSC_VERSION_RELEASE)
      if ($(PETSC_VERSION_PATCH) GREATER 0)
        set (PETSC_VERSION "${PETSC_VERSION_MAJOR}.${PETSC_VERSION_MINOR}.${PETSC_VERSION_SUBMINOR}p${PETSC_VERSION_PATCH}" CACHE INTERNAL "PETSC version")
      else ()
        set (PETSC_VERSION "${PETSC_VERSION_MAJOR}.${PETSC_VERSION_MINOR}.${PETSC_VERSION_SUBMINOR}" CACHE INTERNAL "PETSC version")
      endif ()
    else ()
      # make dev version compare higher than any patch level of a released version
      set (PETSC_VERSION "${PETSC_VERSION_MAJOR}.${PETSC_VERSION_MINOR}.${PETSC_VERSION_SUBMINOR}.99" CACHE INTERNAL "PETSC version")
    endif ()
  else ()
    message (SEND_ERROR "PETSC_DIR can not be used, ${PETSC_DIR}/include/petscversion.h does not exist")
  endif ()
endfunction ()

macro (PETSC_GET_VARIABLE name var)
  set (${var} "NOTFOUND" CACHE INTERNAL "Cleared" FORCE)
  execute_process (COMMAND ${MAKE_EXECUTABLE} --no-print-directory -f ${petsc_config_makefile} show VARIABLE=${name}
    OUTPUT_VARIABLE ${var}
    RESULT_VARIABLE petsc_return)
endmacro (PETSC_GET_VARIABLE)

macro (PETSC_TEST_RUNS includes libraries runs)
  if (PETSC_VERSION VERSION_GREATER 3.1)
    set (_PETSC_TSDestroy "TSDestroy(&ts)")
  else ()
    set (_PETSC_TSDestroy "TSDestroy(ts)")
  endif ()

  set(_PETSC_TEST_SOURCE "
static const char help[] = \"PETSC test program.\";
#include <petscts.h>
int main(int argc,char *argv[]) {
  PetscErrorCode ierr;
  TS ts;

  ierr = PetscInitialize(&argc,&argv,0,help);CHKERRQ(ierr);
  ierr = TSCreate(PETSC_COMM_WORLD,&ts);CHKERRQ(ierr);
  ierr = TSSetFromOptions(ts);CHKERRQ(ierr);
  ierr = ${_PETSC_TSDestroy};CHKERRQ(ierr);
  ierr = PetscFinalize();CHKERRQ(ierr);
  return 0;
}
")

  multipass_source_runs ("${includes}" "${libraries}" "${_PETSC_TEST_SOURCE}" ${runs} "${PETSC_LANGUAGE_BINDINGS}")
  if (${${runs}})
    string (CONCAT docstr 
            "Can the system successfully run a PETSC executable? "
            "This variable can be manually set to \"YES\" to force CMake to accept a given PETSC "
            "configuration, but this will almost always result in a broken build. "
            "If you change PETSC_DIR, PETSC_ARCH, or PETSC_CURRENT you would have to reset this variable.")
    set (PETSC_EXECUTABLE_RUNS "YES" CACHE BOOL ${docstr} FORCE)
  endif (${${runs}})
endmacro (PETSC_TEST_RUNS)

macro (PETSC_FIND_LIBRARY suffix name)
  # Clear any stale value, if we got here, we need to find it again
  set (PETSC_LIBRARY_${suffix} "NOTFOUND" CACHE INTERNAL "Cleared" FORCE)
  if (WIN32)
    set (libname lib${name}) # windows expects "libfoo", linux expects "foo"
  else (WIN32)
    set (libname ${name})
  endif (WIN32)
  
  find_library (PETSC_LIBRARY_${suffix} NAMES ${libname} HINTS ${petsc_lib_dir} NO_DEFAULT_PATH)
  set (PETSC_LIBRARIES_${suffix} "${PETSC_LIBRARY_${suffix}}")
  mark_as_advanced (PETSC_LIBRARY_${suffix}) 
endmacro (PETSC_FIND_LIBRARY suffix name)

macro (PETSC_JOIN libs deps)
  list (APPEND PETSC_LIBRARIES_${libs} ${PETSC_LIBRARIES_${deps}})
endmacro (PETSC_JOIN libs deps)
 

#----------------- FindPETSC -----------------------

set(PETSC_VALID_COMPONENTS C CXX)

if(NOT PETSC_FIND_COMPONENTS)

  get_property (_enabled_langs GLOBAL PROPERTY ENABLED_LANGUAGES)
  list(FIND _enabled_langs "C" _c_index)
  if (${_c_index} GREATER -1)
    set(PETSC_LANGUAGE_BINDINGS "C")
  else ()
    set(PETSC_LANGUAGE_BINDINGS "CXX")
  endif ()

else()

  # Right now, this is designed for compatability with the --with-clanguage option, so
  # only allow one item in the components list.
  list(LENGTH ${PETSC_FIND_COMPONENTS} components_length)
  if(${components_length} GREATER 1)
    message(FATAL_ERROR "Only one component for PETSC is allowed to be specified")
  endif()
  # This is a stub for allowing multiple components should that time ever come. Perhaps
  # to also test Fortran bindings?
  foreach(component ${PETSC_FIND_COMPONENTS})
    list(FIND PETSC_VALID_COMPONENTS ${component} component_location)
    if(${component_location} EQUAL -1)
      message(FATAL_ERROR "\"${component}\" is not a valid PETSC component.")
    else()
      list(APPEND PETSC_LANGUAGE_BINDINGS ${component})
    endif()
  endforeach()

endif()

# if PETSC_DIR was not explicity set, try and find PETSc
if (NOT PETSC_DIR)

  # Debian uses versioned paths e.g /usr/lib/petscdir/3.5/
  file (GLOB DEB_PATHS "/usr/lib/petscdir/*")

  find_path (PETSC_DIR include/petsc.h
    HINTS ENV PETSC_DIR
    PATHS
    /usr/lib/petsc
    # Debian paths
    ${DEB_PATHS}
    # Arch Linux path
    /opt/petsc/linux-c-opt
    # MacPorts path
    /opt/local/lib/petsc
    $ENV{HOME}/petsc
    DOC "PETSC Directory")

    if (PETSC_DIR AND NOT PETSC_ARCH)

      set (_petsc_arches
        $ENV{PETSC_ARCH}                   # If set, use environment variable first
        linux-gnu-c-debug linux-gnu-c-opt  # Debian defaults
        x86_64-unknown-linux-gnu i386-unknown-linux-gnu)
      set (petscconf "NOTFOUND" CACHE FILEPATH "Cleared" FORCE)
      foreach (arch ${_petsc_arches})
        if (NOT PETSC_ARCH)
          find_path (petscconf petscconf.h
            HINTS ${PETSC_DIR}
            PATH_SUFFIXES ${arch}/include bmake/${arch}
            NO_DEFAULT_PATH)
          if (petscconf)
            set (PETSC_ARCH "${arch}" CACHE STRING "PETSC build architecture")
          endif (petscconf)
        endif (NOT PETSC_ARCH)
      endforeach (arch)
      set (petscconf "NOTFOUND" CACHE INTERNAL "Scratch variable" FORCE)

    endif (PETSC_DIR AND NOT PETSC_ARCH)

endif()

find_program (MAKE_EXECUTABLE NAMES make gmake)

set (petsc_slaves LIBRARIES_SYS LIBRARIES_VEC LIBRARIES_MAT LIBRARIES_DM LIBRARIES_KSP
                  LIBRARIES_SNES LIBRARIES_TS INCLUDE_DIR INCLUDE_CONF)

include (FindPackageMultipass)
find_package_multipass (PETSC petsc_config_current
                        STATES DIR ARCH DEPENDENTS INCLUDES LIBRARIES COMPILER MPIEXEC
                        ${petsc_slaves})

# Determine whether the PETSC layout is old-style (through 2.3.3) or
# new-style (>= 3.0.0)
if (EXISTS "${PETSC_DIR}/${PETSC_ARCH}/lib/petsc/conf/petscvariables") # > 3.5
  set (petsc_conf_rules "${PETSC_DIR}/lib/petsc/conf/rules")
  set (petsc_conf_variables "${PETSC_DIR}/lib/petsc/conf/variables")
elseif (EXISTS "${PETSC_DIR}/${PETSC_ARCH}/include/petscconf.h")   # > 2.3.3
  set (petsc_conf_rules "${PETSC_DIR}/conf/rules")
  set (petsc_conf_variables "${PETSC_DIR}/conf/variables")
elseif (EXISTS "${PETSC_DIR}/bmake/${PETSC_ARCH}/petscconf.h") # <= 2.3.3
  set (petsc_conf_rules "${PETSC_DIR}/bmake/common/rules")
  set (petsc_conf_variables "${PETSC_DIR}/bmake/common/variables")
elseif (PETSC_DIR)
  message (SEND_ERROR "The pair PETSC_DIR=${PETSC_DIR} PETSC_ARCH=${PETSC_ARCH} do not specify a valid PETSC installation")
endif ()

if (petsc_conf_rules AND petsc_conf_variables AND NOT petsc_config_current)
  petsc_get_version()

  # Put variables into environment since they are needed to get
  # configuration (petscvariables) in the PETSC makefile
  set (ENV{PETSC_DIR} "${PETSC_DIR}")
  set (ENV{PETSC_ARCH} "${PETSC_ARCH}")

  # A temporary makefile to probe the PETSC configuration
  set (petsc_config_makefile "${PROJECT_BINARY_DIR}/Makefile.petsc")
  file (WRITE "${petsc_config_makefile}"
"## This file was autogenerated by FindPETSC.cmake
# PETSC_DIR  = ${PETSC_DIR}
# PETSC_ARCH = ${PETSC_ARCH}
include ${petsc_conf_rules}
include ${petsc_conf_variables}
show :
\t-@echo -n \${\${VARIABLE}}
")

  petsc_get_variable (PETSC_LIB_DIR            petsc_lib_dir)
  petsc_get_variable (PETSC_EXTERNAL_LIB_BASIC petsc_libs_external)
  petsc_get_variable (PETSC_CCPPFLAGS          petsc_cpp_line)
  petsc_get_variable (PETSC_INCLUDE            petsc_include)
  petsc_get_variable (PCC                      petsc_cc)
  petsc_get_variable (PCC_FLAGS                petsc_cc_flags)
  petsc_get_variable (MPIEXEC                  petsc_mpiexec)
  petsc_get_variable (PETSC_INDEX_SIZE         petsc_index_size)
  petsc_get_variable (PETSC_PRECISION          petsc_precision)
  # We are done with the temporary Makefile, calling PETSC_GET_VARIABLE after this point is invalid!
  file (REMOVE ${petsc_config_makefile})

  include (ResolveCompilerPaths)
  # Extract include paths and libraries from compile command line
  resolve_includes (petsc_includes_all "${petsc_cpp_line}")

  #on windows we need to make sure we're linking against the right
  #runtime library
  if (WIN32)
    if (petsc_cc_flags MATCHES "-MT")

      set(using_md False)
      foreach(flag_var
          CMAKE_C_FLAGS CMAKE_C_FLAGS_DEBUG CMAKE_C_FLAGS_RELEASE
          CMAKE_C_FLAGS_MINSIZEREL CMAKE_C_FLAGS_RELWITHDEBINFO
          CMAKE_CXX_FLAGS CMAKE_CXX_FLAGS_DEBUG CMAKE_CXX_FLAGS_RELEASE
          CMAKE_CXX_FLAGS_MINSIZEREL CMAKE_CXX_FLAGS_RELWITHDEBINFO)
        if(${flag_var} MATCHES "/MD")
          set(using_md True)
        endif(${flag_var} MATCHES "/MD")
      endforeach(flag_var)
      if(${using_md} MATCHES "True")
        string(CONCAT msg "PETSC was built with /MT, but /MD is currently set.\n"
                          "See http://www.cmake.org/Wiki/CMake_FAQ#How_can_I_build_my_MSVC_application_with_a_static_runtime.3F")
        message(WARNING ${msg})
      endif(${using_md} MATCHES "True")

    endif (petsc_cc_flags MATCHES "-MT")
  endif (WIN32)

  include (CorrectWindowsPaths)
  convert_cygwin_path(petsc_lib_dir)
  message (STATUS "petsc_lib_dir ${petsc_lib_dir}")

  # Look for petscvec first, if it doesn't exist, we must be using single-library
  petsc_find_library (VEC petscvec)
  if (PETSC_LIBRARY_VEC)

    petsc_find_library (SYS  "petscsys;petsc") # libpetscsys is called libpetsc prior to 3.1 (when single-library was introduced)
    petsc_find_library (MAT  petscmat)
    petsc_find_library (DM   petscdm)
    petsc_find_library (KSP  petscksp)
    petsc_find_library (SNES petscsnes)
    petsc_find_library (TS   petscts)
    petsc_join (VEC  SYS)
    petsc_join (MAT  VEC)
    petsc_join (DM   MAT)
    petsc_join (KSP  DM)
    petsc_join (SNES KSP)
    petsc_join (TS   SNES)
    petsc_join (ALL  TS)

  else ()

    set (PETSC_LIBRARY_VEC "NOTFOUND" CACHE INTERNAL "Cleared" FORCE) # There is no libpetscvec
    petsc_find_library (SINGLE petsc)
    # Debian 9/Ubuntu 16.04 uses _real and _complex extensions when using libraries in /usr/lib/petsc.
    if (NOT PETSC_LIBRARY_SINGLE)
      petsc_find_library (SINGLE petsc_real)
    endif()
    if (NOT PETSC_LIBRARY_SINGLE)
      petsc_find_library (SINGLE petsc_complex)
    endif()
    foreach (pkg SYS VEC MAT DM KSP SNES TS ALL)
      set (PETSC_LIBRARIES_${pkg} "${PETSC_LIBRARY_SINGLE}")
    endforeach ()

  endif ()

  if (PETSC_LIBRARY_TS)
    message (STATUS "Recognized PETSC install with separate libraries for each package")
  else ()
    message (STATUS "Recognized PETSC install with single library for all packages")
  endif ()

  find_path (PETSC_INCLUDE_DIR petscts.h HINTS "${PETSC_DIR}" PATH_SUFFIXES include NO_DEFAULT_PATH)
  find_path (PETSC_INCLUDE_CONF petscconf.h HINTS "${PETSC_DIR}" PATH_SUFFIXES "${PETSC_ARCH}/include" "bmake/${PETSC_ARCH}" NO_DEFAULT_PATH)
  mark_as_advanced (PETSC_INCLUDE_DIR PETSC_INCLUDE_CONF)
  set (petsc_includes_minimal ${PETSC_INCLUDE_CONF} ${PETSC_INCLUDE_DIR})

  include(Check${PETSC_LANGUAGE_BINDINGS}SourceRuns)
  petsc_test_runs ("${petsc_includes_minimal}" "${PETSC_LIBRARIES_TS}" petsc_works_minimal)
  if (petsc_works_minimal)

    message (STATUS "Minimal PETSC includes and libraries work.  This probably means we are building with shared libs.")
    set (petsc_includes_needed "${petsc_includes_minimal}")

  else (petsc_works_minimal)     # Minimal includes fail, see if just adding full includes fixes it

    petsc_test_runs ("${petsc_includes_all}" "${PETSC_LIBRARIES_TS}" petsc_works_allincludes)
    if (petsc_works_allincludes) # It does, we just need all the includes (

      string (CONCAT msg "PETSC requires extra include paths, but links correctly with only interface libraries.\n"
                         "This is an unexpected configuration (but it seems to work fine).")
      message (STATUS ${msg})
      set (petsc_includes_needed ${petsc_includes_all})

    else (petsc_works_allincludes) # We are going to need to link the external libs explicitly

      resolve_libraries (petsc_libraries_external "${petsc_libs_external}")
      foreach (pkg SYS VEC MAT DM KSP SNES TS ALL)
        list (APPEND PETSC_LIBRARIES_${pkg}  ${petsc_libraries_external})
      endforeach (pkg)

      petsc_test_runs ("${petsc_includes_minimal}" "${PETSC_LIBRARIES_TS}" petsc_works_alllibraries)
      if (petsc_works_alllibraries)

        string (CONCAT msg "PETSC only need minimal includes, but requires explicit linking to all dependencies.\n"
                           "This is expected when PETSC is built with static libraries.")
        message(STATUS ${msg})
        set (petsc_includes_needed ${petsc_includes_minimal})

      else (petsc_works_alllibraries)

        # It looks like we really need everything, should have listened to Matt
        set (petsc_includes_needed ${petsc_includes_all})
        petsc_test_runs ("${petsc_includes_all}" "${PETSC_LIBRARIES_TS}" petsc_works_all)
        if (petsc_works_all) # We fail anyways
          string (CONCAT msg "PETSC requires extra include paths and explicit linking to all dependencies.\n"
                             "This probably means you have static libraries and something unexpected in PETSC headers.")
          message (STATUS ${msg})
        else (petsc_works_all) # We fail anyways
          message (STATUS "PETSC could not be used, maybe the install is broken.")
        endif (petsc_works_all)

      endif (petsc_works_alllibraries)

    endif (petsc_works_allincludes)

  endif (petsc_works_minimal)

  # We do an out-of-source build so __FILE__ will be an absolute path, hence __INSDIR__ is superfluous
  if (${PETSC_VERSION} VERSION_LESS 3.1)
    set (PETSC_DEFINITIONS "-D__SDIR__=\"\"" CACHE STRING "PETSC definitions" FORCE)
  else ()
    set (PETSC_DEFINITIONS "-D__INSDIR__=" CACHE STRING "PETSC definitions" FORCE)
  endif ()
  
  set (PETSC_INDEX_SIZE ${petsc_index_size} CACHE INTERNAL "PETSC index size" FORCE)
  set (PETSC_PRECISION ${petsc_precision} CACHE INTERNAL "PETSC real type precision" FORCE)
  
  # Sometimes this can be used to assist FindMPI.cmake
  set (PETSC_MPIEXEC ${petsc_mpiexec} CACHE FILEPATH "Executable for running PETSC MPI programs" FORCE)
  set (PETSC_INCLUDES ${petsc_includes_needed} CACHE STRING "PETSC include path" FORCE)
  set (PETSC_LIBRARIES ${PETSC_LIBRARIES_ALL} CACHE STRING "PETSC libraries" FORCE)
  set (PETSC_COMPILER ${petsc_cc} CACHE FILEPATH "PETSC compiler" FORCE)
  
  # Note that we have forced values for all these choices.  If you
  # change these, you are telling the system to trust you that they
  # work.  It is likely that you will end up with a broken build.
  mark_as_advanced (PETSC_CURRENT PETSC_INCLUDES PETSC_LIBRARIES PETSC_COMPILER PETSC_DEFINITIONS
    PETSC_MPIEXEC PETSC_EXECUTABLE_RUNS)
endif ()

# Create targets
if (PETSC_LIBRARY_SINGLE)
  foreach (suffix SYS VEC MAT DM KSP SNES TS ALL)
    if (NOT TARGET PETSC::${suffix})
      add_library (PETSC::${suffix} UNKNOWN IMPORTED)
      set_target_properties (PETSC::${suffix} PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES ${PETSC_INCLUDES}
        INTERFACE_LINK_LIBRARIES ${PETSC_LIBRARIES}
        INTERFACE_COMPILE_OPTIONS ${PETSC_DEFINTIONS}
        IMPORTED_LOCATION ${PETSC_LIBRARY_SINGLE})
    endif ()
  endforeach ()
else ()
  foreach (suffix SYS VEC MAT DM KSP SNES TS ALL)
    if (PETSC_LIBRARY_${suffix} AND (NOT TARGET PETSC::${suffix}))
      add_library (PETSC::${suffix} UNKNOWN IMPORTED)
      set_target_properties (PETSC::${suffix} PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES ${PETSC_INCLUDES}
        INTERFACE_LINK_LIBRARIES ${PETSC_LIBRARIES_${suffix}}
        INTERFACE_COMPILE_OPTIONS ${PETSC_DEFINTIONS}
        IMPORTED_LOCATION ${PETSC_LIBRARY_${suffix}})
    endif ()
  endforeach ()
endif ()

include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (PETSC
  REQUIRED_VARS PETSC_INCLUDES PETSC_LIBRARIES PETSC_EXECUTABLE_RUNS
  VERSION_VAR PETSC_VERSION
  FAIL_MESSAGE "PETSC could not be found. Be sure to set PETSC_DIR and PETSC_ARCH.")

