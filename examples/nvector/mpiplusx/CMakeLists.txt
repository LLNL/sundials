# ---------------------------------------------------------------
# Programmer(s): Cody J. Balos @ LLNL
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
# CMakeLists.txt file for parallel MPIPlusX nvector examples
# ---------------------------------------------------------------

# Example lists are tuples "name\;args\;nodes\;tasks\;type" where the
# type is develop for examples excluded from 'make test' in releases

# Examples using SUNDIALS MPI nvector
set(nvector_mpiplusx_examples
  "test_nvector_mpiplusx.c\; \"1000 9\" \;1\;4\;passfail"  # run parallel on 4 procs
  )

# If building F2003 tests
if (BUILD_FORTRAN_MODULE_INTERFACE)
  set(nvector_mpiplusx_fortran_examples
    "test_fnvector_mpiplusx_mod\;\;1\;1\;"
    "test_fnvector_mpiplusx_mod\;\;1\;4\;")
endif()

# Dependencies for nvector examples
set(nvector_examples_dependencies
  test_nvector.c
  test_mpinvector.c
  test_nvector.h
  )

# Adding extra files 
set(nvector_dependencies
  ../test_nvector.c
  ../test_nvector.h
  ../test_mpinvector.c
) 
  
# Add source directory to include directories
include_directories(. ..)

# Specify libraries to link against
set(NVECS_LIB sundials_nvecmpiplusx
              sundials_nvecparallel
              sundials_nvecserial)
if (BUILD_FORTRAN_MODULE_INTERFACE)
  list(APPEND NVECS_LIB sundials_fnvecmpiplusx_mod
                        sundials_fnvecserial_mod)
endif()

# Set-up linker flags and link libraries
set(SUNDIALS_LIBS ${NVECS_LIB} ${EXE_EXTRA_LINK_LIBS})

if(MPI_C_COMPILER)
  # use MPI wrapper as the compiler
  set(CMAKE_C_COMPILER ${MPI_C_COMPILER})
else()
  # add MPI_INCLUDE_PATH to include directories
  include_directories(${MPI_INCLUDE_PATH})
endif()

if (BUILD_FORTRAN_MODULE_INTERFACE)
  if(MPI_Fortran_COMPILER)
    # use MPI wrapper as the compiler
    set(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
  else()
    # add MPI_INCLUDE_PATH to include directories
    include_directories(${MPI_INCLUDE_PATH})
  endif()
endif()

# Add the build and install targets for each example
foreach(example_tuple ${nvector_mpiplusx_examples})

  # parse the example tuple using sundials_parse_tuple macro
  sundials_parse_tuple( nvector_mpiplusx_examples 
    FILENAME_VAR example
    PARAMS_VAR example_args
    NODES_VAR number_of_nodes
    NPROCS_VAR number_of_tasks
    TYPE_VAR example_type
  )

  # check if this example has already been added, only need to add
  # example source files once for testing with different inputs
  if(NOT TARGET ${example})
    # example source files
    add_executable(${example} ${example})

    # link vector test utilties
    target_link_libraries(${example} PRIVATE test_nvector_obj test_nvectormpi_obj)

    # folder to organize targets in an IDE
    set_target_properties(${example} PROPERTIES FOLDER "Examples")

    # libraries to link against
    target_link_libraries(${example} PRIVATE ${SUNDIALS_LIBS})

    if(NOT MPI_C_COMPILER)
      target_link_libraries(${example} PRIVATE ${MPI_LIBRARIES})
    endif()
  endif()

  # check if example args are provided and set the test name
  if("${example_args}" STREQUAL "")
    if("${number_of_tasks}" STREQUAL "")
      set(test_name ${example})
    else()
      string(REGEX REPLACE " " "_" test_name ${example}_${number_of_tasks})
    endif()
  else()
    if("${number_of_tasks}" STREQUAL "")
      string(REGEX REPLACE " " "_" test_name ${example}_${example_args})
    else()
      string(REGEX REPLACE " " "_" test_name ${example}_${number_of_tasks}_${example_args})
    endif()
  endif()

  # add example to regression tests
  sundials_add_test(${test_name} ${example}
    TEST_ARGS ${example_args}
    MPI_NPROCS ${number_of_tasks}
    EXAMPLE_TYPE ${example_type}
    NODIFF)


endforeach(example_tuple ${nvector_mpiplusx_examples})

# Add the build and install targets for each example
foreach(example_tuple ${nvector_mpiplusx_fortran_examples})

  # parse the example tuple using sundials_parse_tuple macro
  sundials_parse_tuple( nvector_mpiplusx_fortran_examples 
    FILENAME_VAR example
    PARAMS_VAR example_args
    NODES_VAR number_of_nodes
    NPROCS_VAR number_of_tasks
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
    target_link_libraries(${example} ${SUNDIALS_LIBS})

    if(NOT MPI_Fortran_COMPILER)
      target_link_libraries(${example} ${MPI_LIBRARIES})
    endif()
  endif()

  # check if example args are provided and set the test name
  if("${example_args}" STREQUAL "")
    if("${number_of_tasks}" STREQUAL "")
      set(test_name ${example})
    else()
      string(REGEX REPLACE " " "_" test_name ${example}_${number_of_tasks})
    endif()
  else()
    if("${number_of_tasks}" STREQUAL "")
      string(REGEX REPLACE " " "_" test_name ${example}_${example_args})
    else()
      string(REGEX REPLACE " " "_" test_name ${example}_${number_of_tasks}_${example_args})
    endif()
  endif()

  # add example to regression tests
  sundials_add_test(${test_name} ${example}
    TEST_ARGS ${example_args}
    MPI_NPROCS ${number_of_tasks}
    EXAMPLE_TYPE ${example_type}
    NODIFF)

endforeach(example_tuple ${nvector_mpiplusx_fortran_examples})

# create Makfile and CMakeLists.txt for examples
if(EXAMPLES_INSTALL)

  sundials_install_examples(nvector nvector_mpiplusx_examples
    HAS_MPI
    CMAKE_TEMPLATE
      cmakelists_C_MPI_ex.in
    MAKE_TEMPLATE
      makefile_parallel_C_ex.in
    SUNDIALS_TARGETS
      nvecparallel
      nvecserial
      nvecmpiplusx
    DESTINATION
      nvector/mpiplusx
    EXAMPLES_DEPENDENCIES
      ${nvector_examples_dependencies}
    EXTRA_FILES
      #README
      #${NVECTOR_extras}
      #${NVECTOR_headers}
      ${nvector_dependencies}
    TEST_INSTALL
      mpiplusx
  )

endif()