# ---------------------------------------------------------------
# Programmer:  Cody Balos @ LLNL
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
# SuperLUDIST find module that creates an imported target for
# SuperLU_DIST. The target is SuperLU_DIST::SuperLU_DIST.
#
# The variable SUPERLUDIST_LIBRARY_DIR can be used to control
# where the module looks for the library.
#
# The variable SUPERLUDIST_INCLUDE_DIR can be used to set the
# include path for the library.
#
# Additional libraries can be passed in SUPERLUDIST_LIBRARIES.
#
# This module also defines variables, but it is best to use
# the defined target to ensure includes and compile/link
# options are correctly passed to consumers.
#
#   SUPERLUDIST_FOUND      - system has SuperLU_DIST library
#   SUPERLUDIST_LIBRARY    - the SuperLU_DIST library
#   SUPERLUDIST_INDEX_SIZE - the bit width of indices in SUPERLUDIST
# ---------------------------------------------------------------

# Check if SUPERLUDIST_LIBRARIES contains the superlu_dist
# library as well as TPLs. If so, extract it into the
# SUPERLUDIST_LIBRARY variable.
if (SUPERLUDIST_LIBRARIES MATCHES "superlu_dist")
  foreach(lib ${SUPERLUDIST_LIBRARIES})
    if(lib MATCHES "superlu_dist")
      set(SUPERLUDIST_LIBRARY ${lib})
    endif()
  endforeach()
endif()

# find library with user provided directory path
if(NOT SUPERLUDIST_LIBRARY)
  find_library(SUPERLUDIST_LIBRARY superlu_dist PATHS ${SUPERLUDIST_LIBRARY_DIR} NO_DEFAULT_PATH)
  # if user didnt provide path, search anywhere
  if (NOT (SUPERLUDIST_LIBRARY_DIR OR SUPERLUDIST_LIBRARY))
    find_library(SUPERLUDIST_LIBRARY superlu_dist)
  endif()
  mark_as_advanced(SUPERLUDIST_LIBRARY)
endif()

# set the library dir option if it wasnt preset
if(SUPERLUDIST_LIBRARY AND (NOT SUPERLUDIST_LIBRARY_DIR))
  get_filename_component(SUPERLUDIST_LIBRARY_DIR ${SUPERLUDIST_LIBRARY} DIRECTORY)
  set(SUPERLUDIST_LIBRARY_DIR ${SUPERLUDIST_LIBRARY_DIR} CACHE STRING "" FORCE)
endif()

# set the include dir option if it wasnt preset
if(SUPERLUDIST_LIBRARY AND (NOT SUPERLUDIST_INCLUDE_DIR))
  get_filename_component(SUPERLUDIST_INCLUDE_DIR ${SUPERLUDIST_LIBRARY_DIR} DIRECTORY)
  set(SUPERLUDIST_INCLUDE_DIR "${SUPERLUDIST_INCLUDE_DIR}/include" CACHE STRING "" FORCE)
endif()

if(SUPERLUDIST_LIBRARY AND SUPERLUDIST_INCLUDE_DIR)
  find_file(SUPERLUDIST_CONFIG_PATH superlu_dist_config.h PATHS ${SUPERLUDIST_INCLUDE_DIR})
  file(STRINGS ${SUPERLUDIST_CONFIG_PATH} _strings_with_index_size REGEX "XSDK_INDEX_SIZE")
  list(GET _strings_with_index_size 0 _index_size_string)
  string(REGEX MATCHALL "[0-9][0-9]" SUPERLUDIST_INDEX_SIZE "${_index_size_string}")
  mark_as_advanced(FORCE SUPERLUDIST_CONFIG_PATH)
endif()

# Create target for SuperLU_DIST
if((NOT TARGET SuperLU_DIST::SuperLU_DIST) AND SUPERLUDIST_LIBRARY)
  add_library(SuperLU_DIST::SuperLU_DIST UNKNOWN IMPORTED)
  set_property(TARGET SuperLU_DIST::SuperLU_DIST PROPERTY INTERFACE_INCLUDE_DIRECTORIES ${SUPERLUDIST_INCLUDE_DIR})
  set_property(TARGET SuperLU_DIST::SuperLU_DIST PROPERTY INTERFACE_LINK_LIBRARIES ${SUPERLUDIST_LIBRARIES})
  set_property(TARGET SuperLU_DIST::SuperLU_DIST PROPERTY IMPORTED_LOCATION ${SUPERLUDIST_LIBRARY})
endif()

find_package_handle_standard_args(SUPERLUDIST DEFAULT_MSG SUPERLUDIST_LIBRARY SUPERLUDIST_INCLUDE_DIR SUPERLUDIST_INDEX_SIZE)
