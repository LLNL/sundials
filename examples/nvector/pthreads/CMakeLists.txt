# ---------------------------------------------------------------
# Programmer(s): Steven Smith @ LLNL
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
# CMakeLists.txt file for PThread nvector examples
# ---------------------------------------------------------------

# Example lists are tuples "name\;args\;type" where the type is
# 'develop' for examples excluded from 'make test' in releases

# Examples using SUNDIALS pthread nvector
set(nvector_pthreads_examples
  "test_nvector_pthreads.c\; \"1000 1 0\" \; \"\" "
  "test_nvector_pthreads.c\; \"1000 2 0\" \; \"\" "
  "test_nvector_pthreads.c\; \"1000 4 0\" \; \"\" "
  "test_nvector_pthreads.c\; \"10000 1 0\" \; \"\" "
  "test_nvector_pthreads.c\; \"10000 2 0\" \; \"\" "
  "test_nvector_pthreads.c\; \"10000 4 0\" \; \"\" "
  )

# Dependencies for nvector examples
set(nvector_examples_dependencies
  test_nvector.c
  test_nvector.h
  )
# Adding extra files 
set(nvector_examples_files
  ../test_nvector.c
  ../test_nvector.h
)

# If building F2003 tests
if (BUILD_FORTRAN_MODULE_INTERFACE)
  set(nvector_pthreads_fortran_examples
    "test_fnvector_pthreads_mod\;\;")
endif()

# Add source directory to include directories
include_directories(. ..)

# Specify libraries to link against
set(NVECS_LIB sundials_nvecpthreads)
if (BUILD_FORTRAN_MODULE_INTERFACE)
  list(APPEND NVECS_LIB sundials_fnvecpthreads_mod)
endif()

# Set-up linker flags and link libraries
set(SUNDIALS_LIBS ${NVECS_LIB} ${EXE_EXTRA_LINK_LIBS})

# Add the build and install targets for each example
foreach(example_tuple ${nvector_pthreads_examples})
  
  # parse the example tuple using sundials_parse_tuple macro
  sundials_parse_tuple( nvector_pthreads_examples
    FILENAME_VAR example
    PARAMS_VAR example_args
    TYPE_VAR example_type
  ) 

  # check if this example has already been added, only need to add
  # example source files once for testing with different inputs
  if(NOT TARGET ${example})
    # example source files
    add_executable(${example} ${example})

    # link vector test utilties
    target_link_libraries(${example} PRIVATE test_nvector_obj)

    # folder to organize targets in an IDE
    set_target_properties(${example} PROPERTIES FOLDER "Examples")

    # libraries to link against
    target_link_libraries(${example} PRIVATE ${SUNDIALS_LIBS} ${CMAKE_THREAD_LIBS_INIT})
  endif()

  # check if example args are provided and set the test name
  if("${example_args}" STREQUAL "")
    set(test_name ${example})
  else()
    string(REGEX REPLACE " " "_" test_name ${example}_${example_args})
  endif()

  # add example to regression tests
  sundials_add_test(${test_name} ${example}
    TEST_ARGS ${example_args}
    EXAMPLE_TYPE ${example_type}
    NODIFF)

endforeach(example_tuple ${nvector_pthreads_examples})

# Add the build and install targets for each example
foreach(example_tuple ${nvector_pthreads_fortran_examples})

  # parse the example tuple using sundials_parse_tuple macro
  sundials_parse_tuple( nvector_openmp_fortran_examples
    FILENAME_VAR example
    PARAMS_VAR example_args
    TYPE_VAR example_type
  ) 

  # check if this example has already been added, only need to add
  # example source files once for testing with different inputs
  if(NOT TARGET ${example})
    # example source files
    add_executable(${example} ${example}.f90
      ${SUNDIALS_SOURCE_DIR}/examples/utilities/test_utilities.f90
      ${SUNDIALS_SOURCE_DIR}/examples/nvector/test_nvector.f90)

    # folder to organize targets in an IDE
    set_target_properties(${example} PROPERTIES FOLDER "Examples")
    set_target_properties(${example} PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

    # libraries to link against
    target_link_libraries(${example} ${SUNDIALS_LIBS} ${CMAKE_THREAD_LIBS_INIT})
  endif()

  # check if example args are provided and set the test name
  if("${example_args}" STREQUAL "")
    set(test_name ${example})
  else()
    string(REGEX REPLACE " " "_" test_name ${example}_${example_args})
  endif()

  # add example to regression tests
  sundials_add_test(${test_name} ${example}
    TEST_ARGS ${example_args}
    EXAMPLE_TYPE ${example_type}
    NODIFF)

endforeach(example_tuple ${nvector_pthreads_fortran_examples})

# create Makfile and CMakeLists.txt for examples
if(EXAMPLES_INSTALL)

  sundials_install_examples(nvector nvector_pthreads_examples
    CMAKE_TEMPLATE
      cmakelists_C_pthreads_ex.in
    MAKE_TEMPLATE
      makefile_pthreads_C_ex.in
    SUNDIALS_TARGETS
      nvecpthreads
    DESTINATION
      nvector/pthreads
    EXAMPLES_DEPENDENCIES
      ${nvector_examples_dependencies}
    EXTRA_FILES
      ${nvector_examples_files}
    TEST_INSTALL
      pthreads
  )

endif()
