# ---------------------------------------------------------------
# Programmer(s): Eddy Banks and David J. Gardner @ LLNL
# ---------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2021, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# ---------------------------------------------------------------
# SuperLU find module that creates an imported target for
# SuperLU. The target is SUNDIALS::SUPERLU.
#
# The variable SUPERLU_LIBRARY_DIR can be used to control
# where the module looks for the library.
#
# The variable SUPERLU_INCLUDE_DIR can be used to set the
# include path for the library.
#
# Additional libraries can be passed in SUPERLU_LIBRARIES.
#
# This module also defines variables, but it is best to use
# the defined target to ensure includes and compile/link
# options are correctly passed to consumers.
#
#   SUPERLU_FOUND       - system has SuperLU library
#   SUPERLU_LIBRARY     - the SuperLU library
# ---------------------------------------------------------------

# Set SuperLU library name with thread type postfix
set(SUPERLU_LIBRARY_NAME superlu)

if(MSVC)
  set(CMAKE_FIND_LIBRARY_PREFIXES lib ${CMAKE_FIND_LIBRARY_PREFIXES})
endif()

# Check if SUPERLU_LIBRARIES contains the superlu
# library as well as TPLs. If so, extract it into the
# SUPERLU_LIBRARY variable.
if(SUPERLU_LIBRARIES MATCHES "${SUPERLU_LIBRARY_NAME}")
  foreach(lib ${SUPERLU_LIBRARIES})
    if(lib MATCHES "${SUPERLU_LIBRARY_NAME}")
      set(SUPERLU_LIBRARY ${lib})
    endif()
  endforeach()
endif()

# find library
if(NOT SUPERLU_LIBRARY)
  # search user provided directory path
  find_library(SUPERLU_LIBRARY ${SUPERLU_LIBRARY_NAME}
    PATHS ${SUPERLU_LIBRARY_DIR} NO_DEFAULT_PATH)
  # if user didn't provide a path, search anywhere
  if(NOT (SUPERLU_LIBRARY_DIR OR SUPERLU_LIBRARY))
    find_library(SUPERLU_LIBRARY ${SUPERLU_LIBRARY_NAME})
  endif()
  mark_as_advanced(SUPERLU_LIBRARY)
endif()

# set the libraries, stripping out 'NOTFOUND' from previous attempts
string(REPLACE "SUPERLU_LIBRARY-NOTFOUND" "" SUPERLU_LIBRARIES "${SUPERLU_LIBRARIES}")
set(SUPERLU_LIBRARIES "${SUPERLU_LIBRARY};${SUPERLU_LIBRARIES}" CACHE STRING "" FORCE)

# set the library dir option if it wasn't preset
if(SUPERLU_LIBRARY AND (NOT SUPERLU_LIBRARY_DIR))
  get_filename_component(SUPERLU_LIBRARY_DIR ${SUPERLU_LIBRARY} DIRECTORY)
  set(SUPERLU_LIBRARY_DIR ${SUPERLU_LIBRARY_DIR} CACHE PATH "" FORCE)
endif()

# set the include dir option if it wasn't preset
if(SUPERLU_LIBRARY AND (NOT SUPERLU_INCLUDE_DIR))
  get_filename_component(SUPERLU_INCLUDE_DIR ${SUPERLU_LIBRARY_DIR} DIRECTORY)
  set(SUPERLU_INCLUDE_DIR "${SUPERLU_INCLUDE_DIR}/include" CACHE PATH "" FORCE)
endif()

# set a more informative error message in case the library was not found
set(SUPERLU_NOT_FOUND_MESSAGE "\
************************************************************************\n\
ERROR: Could not find SuperLU. Please check the variables:\n\
       SUPERLU_INCLUDE_DIR and SUPERLU_LIBRARY_DIR\n\
************************************************************************")

# set package variables including SUPERLU_FOUND
find_package_handle_standard_args(SUPERLU
  REQUIRED_VARS
    SUPERLU_LIBRARY
    SUPERLU_LIBRARIES
    SUPERLU_INCLUDE_DIR
  FAIL_MESSAGE
    "${SUPERLU_NOT_FOUND_MESSAGE}"
  )

# Create target for SuperLU
if(SUPERLU_FOUND)

  if(NOT TARGET SUNDIALS::SUPERLU)
    add_library(SUNDIALS::SUPERLU UNKNOWN IMPORTED)
  endif()

  set_target_properties(SUNDIALS::SUPERLU PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${SUPERLU_INCLUDE_DIR}"
    INTERFACE_LINK_LIBRARIES "${SUPERLU_LIBRARIES}"
    IMPORTED_LOCATION "${SUPERLU_LIBRARY}")

endif()
