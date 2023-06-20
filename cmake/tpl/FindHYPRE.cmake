# ---------------------------------------------------------------
# Programmer(s): Eddy Banks, Slaven Peles, Cody J. Balos, and
#                Jean Sexton @ LLNL
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
# HYPRE find module that creates an imported target for HYPRE.
# The target is SUNDIALS::HYPRE.
#
# The variable HYPRE_LIBRARY_DIR can be used to control
# where the module looks for the library.
#
# The variable HYPRE_INCLUDE_DIR can be used to set the
# include path for the library.
#
# This module also defines variables, but it is best to use
# the defined target to ensure includes and compile/link
# options are correctly passed to consumers.
#
#   HYPRE_FOUND       - system has HYPRE library
#   HYPRE_DIR
#   HYPRE_LIBRARY     - the HYPRE library
#   HYPRE_INCLUDE_DIR - the HYPRE include path
#   HYPRE_LIBRARIES   - all of the libraries needed for HYPRE
# ---------------------------------------------------------------

# --- Find hypre include dir ---
find_path(temp_HYPRE_INCLUDE_DIR
          NAMES HYPRE.h hypre.h
          HINTS "${HYPRE_DIR}" "${HYPRE_DIR}/include" "${HYPRE_INCLUDE_DIR}")
if (temp_HYPRE_INCLUDE_DIR)
    set(HYPRE_INCLUDE_DIR "${temp_HYPRE_INCLUDE_DIR}" CACHE PATH "" FORCE)
endif()
unset(temp_HYPRE_INCLUDE_DIR CACHE)

# --- Find hypre library ---
if (HYPRE_LIBRARY)
    # We have (or were given) HYPRE_LIBRARY - get path to use for any related libs
    get_filename_component(HYPRE_LIBRARY_DIR ${HYPRE_LIBRARY} PATH)

    # force CACHE update to show user DIR that will be used
    set(HYPRE_LIBRARY_DIR ${HYPRE_LIBRARY_DIR} CACHE PATH "" FORCE)
else ()
    # find library with user provided directory path
    set(HYPRE_LIBRARY_NAMES hypre HYPRE)
    find_library(HYPRE_LIBRARY
      NAMES ${HYPRE_LIBRARY_NAMES}
      HINTS "${HYPRE_DIR}" "${HYPRE_DIR}/lib" "${HYPRE_DIR}/lib64" "${HYPRE_LIBRARY_DIR}"
      NO_DEFAULT_PATH
      )
endif ()
mark_as_advanced(HYPRE_LIBRARY)

# --- Append found library to HYPRE_LIBRARIES ---
list(FIND HYPRE_LIBRARIES ${HYPRE_LIBRARY} _idx)
if (_idx EQUAL -1)
  # Automatically overwrite "HYPRE_LIBRARY-NOTFOUND" entry if present
  if (HYPRE_LIBRARIES MATCHES "HYPRE_LIBRARY-NOTFOUND")
    set(HYPRE_LIBRARIES "${HYPRE_LIBRARY}" CACHE STRING "" FORCE)
  # ...otherwise append HYPRE_LIBRARY entry.
  else ()
    set(HYPRE_LIBRARIES "${HYPRE_LIBRARY};${HYPRE_LIBRARIES}" CACHE STRING "" FORCE)
  endif ()
endif ()

# --- Set a more informative error message in case the library was not found ---
set(HYPRE_CONFIG_NOT_FOUND_MESSAGE "\
************************************************************************\n\
ERROR: Could not find hypre library configuration file (HYPRE_config.h).\n\
       Please specify HYPRE_DIR and ensure that it contains\n\
       \"include\" and \"lib\" or \"lib64\" subdirectories.\n\
       (e.g. \".../hypre/src/hypre\")\n\
************************************************************************")
set(HYPRE_NOT_FOUND_MESSAGE "\
************************************************************************\n\
ERROR: Could not find hypre. Please check the variables:\n\
       HYPRE_INCLUDE_DIR and HYPRE_LIBRARY_DIR\n\
************************************************************************")

# --- Find the hypre library configuration file (HYPRE_config.h) ---
find_file(HYPRE_CONFIGH_PATH HYPRE_config.h
          HINTS "${HYPRE_DIR}"
          PATH_SUFFIXES include
          NO_DEFAULT_PATH)
mark_as_advanced(FORCE HYPRE_CONFIGH_PATH)
if (HYPRE_CONFIGH_PATH)
  message(STATUS "hypre library configuration file found. Parsing for version and backends...")
else ()
  message(ERROR "${HYPRE_CONFIG_NOT_FOUND_MESSAGE}")
endif ()

# --- Parse config for hypre version ---
file(READ "${HYPRE_CONFIGH_PATH}" _hypre_config_file_text)
string(REGEX MATCH "[0-9]+\.[0-9]+\.[0-9]+" _hypre_release_version "${_hypre_config_file_text}")
message(STATUS "hypre Version: ${_hypre_release_version}")

# --- Parse config for hypre backends ---
foreach(_backend CUDA HIP)
  file(STRINGS "${HYPRE_CONFIGH_PATH}" _hypre_has_backend REGEX "^#define HYPRE_USING_${_backend}")
  if(_hypre_has_backend)
    set(HYPRE_BACKENDS "${_backend};${HYPRE_BACKENDS}")
    message(STATUS "hypre built with ${_backend} backend? - YES")
  else()
    message(STATUS "hypre built with ${_backend} backend? - NO")
  endif()
endforeach()

# --- Parse config for CUDA Unified Memory ---
file(STRINGS "${HYPRE_CONFIGH_PATH}" _hypre_using_unified_memory REGEX "^#define HYPRE_USING_UNIFIED_MEMORY")
if(_hypre_using_unified_memory)
  set(SUNDIALS_HYPRE_USING_UNIFIED_MEMORY TRUE)
  message(STATUS "hypre using CUDA Unified Memory? - YES")
else()
  message(STATUS "hypre using CUDA Unified Memory? - NO")
endif()

# --- Add libraries for backend support ---
if(SUNDIALS_HYPRE_BACKENDS MATCHES "CUDA")
  list(APPEND HYPRE_LIBRARIES CUDA::cudart CUDA::cublas)
#   find_package(CUDA REQUIRED)
#   include_directories(${CUDA_INCLUDE_DIRS})
endif()

# --- Set package variables including HYPRE_FOUND ---
find_package_handle_standard_args(HYPRE
  REQUIRED_VARS
    HYPRE_LIBRARY
    HYPRE_LIBRARIES
    HYPRE_INCLUDE_DIR
  FAIL_MESSAGE
    "${HYPRE_NOT_FOUND_MESSAGE}"
  )

# --- Create target for HYPRE ---
if(HYPRE_FOUND)

  if(NOT TARGET SUNDIALS::HYPRE)
    add_library(SUNDIALS::HYPRE UNKNOWN IMPORTED)
  endif()

  set_target_properties(SUNDIALS::HYPRE PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${HYPRE_INCLUDE_DIR}"
    INTERFACE_LINK_LIBRARIES "${HYPRE_LIBRARIES}"
    IMPORTED_LOCATION "${HYPRE_LIBRARY}")

endif()
