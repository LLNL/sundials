# ------------------------------------------------------------------------------
# Programmer(s): David J. Gardner @ LLNL
# ------------------------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2025, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# ------------------------------------------------------------------------------
# The macro:
#
#   sundials_add_example(EXAMPLE_SRC)
#
# ------------------------------------------------------------------------------

function(sundials_add_example SOURCE)

  set(options)
  set(singleValueArgs LANGUAGE BACKEND)
  set(multiValueArgs INCLUDE_DIRECTORIES LINK_LIBRARIES ONLY EXCLUDE)

  # parse inputs and create variables arg_<keyword>
  cmake_parse_arguments(arg "${options}" "${oneValueArgs}" "${multiValueArgs}"
                        ${ARGN})

  # examples can also be included or excluded using an "if" statement around the
  # add_subdirectory call for the example -- the ONLY and EXCLUDE options can be
  # helpful when transitioning to the one directory for each example structure
  foreach(_index_size "32" "64")
    if(("${_index_size}" IN_LIST arg_ONLY) AND (NOT SUNDIALS_INDEX_SIZE MATCHES "${_index_size}"))
      return()
    endif()
    if(("${_index_size}" IN_LIST arg_EXCLUDE) AND (SUNDIALS_INDEX_SIZE MATCHES "${_index_size}"))
      return()
    endif()
  endforeach()

  foreach(_precision "single" "double" "extended")
    if(("${_precision}" IN_LIST arg_ONLY) AND (NOT SUNDIALS_PRECISION MATCHES "${_precision}"))
      return()
    endif()
    if(("${_precision}" IN_LIST arg_EXCLUDE) AND (SUNDIALS_PRECISION MATCHES "${_precision}"))
      return()
    endif()
  endforeach()

  foreach(_scalar_type "real" "complex")
    if(("${_scalar_type}" IN_LIST arg_ONLY) AND (NOT SUNDIALS_SCALAR_TYPE MATCHES "${_scalar_type}"))
      return()
    endif()
    if(("${_scalar_type}" IN_LIST arg_EXCLUDE) AND (SUNDIALS_SCALAR_TYPE MATCHES "${_scalar_type}"))
      return()
    endif()
  endforeach()

  # extract the file name without extension and create the target name
  get_filename_component(_name ${SOURCE} NAME_WE)

  if(arg_BACKEND)
    set(_target "${_name}.${arg_BACKEND}")
  else()
    set(_target "${_name}")
  endif()

  if(arg_LANGUAGE)
    set_source_files_properties(${SOURCE} PROPERTIES LANGUAGE ${arg_LANGUAGE})
  endif()

  # create target
  add_executable(${_target} ${SOURCE})

  # folder for IDEs
  set_target_properties(${_target} PROPERTIES FOLDER "Examples")

  # which backend to use
  if(arg_BACKEND)
    target_compile_definitions(${_target} PRIVATE USE_${arg_BACKEND})
  endif()

  # directories to include
  target_include_directories(${_target} PRIVATE ${arg_INCLUDE_DIRECTORIES})

  # libraries to link against
  target_link_libraries(${_target} PRIVATE ${arg_LINK_LIBRARIES})

endfunction()
