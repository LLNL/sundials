# ---------------------------------------------------------------
# Programmer(s): Daniel R. Reynolds @ SMU
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
# CMakeLists.txt file for sunlinsol SPFGMR examples
# ---------------------------------------------------------------

# Set tolerance for linear solver test based on Sundials precision
if(SUNDIALS_PRECISION MATCHES "SINGLE")
  set(TOL "1e-5")
elseif(SUNDIALS_PRECISION MATCHES "DOUBLE")
  set(TOL "1e-13")
else()
  set(TOL "1e-14")
endif()

# Example lists are tuples "name\;args\;type" where the type is
# 'develop' for examples excluded from 'make test' in releases

# Examples using SUNDIALS SPFGMR linear solver
set(sunlinsol_spfgmr_examples
  "test_sunlinsol_spfgmr_serial.c\;100 1 100 ${TOL} 0\;passfail"
  "test_sunlinsol_spfgmr_serial.c\;100 2 100 ${TOL} 0\;passfail"
  )

if (BUILD_FORTRAN_MODULE_INTERFACE)
  set(sunlinsol_spfgmr_fortran_examples
    "test_fsunlinsol_spfgmr_mod_serial\;<none>\;passfail")
endif()

# Dependencies for sunlinsol examples
set(sunlinsol_spfgmr_dependencies
  test_sunlinsol.c
  test_sunlinsol.h
  )
# Adding extra files  
set(sunlinsol_spfgmr_files
  ../../test_sunlinsol.c
  ../../test_sunlinsol.h
  )

# Add source directory to include directories
include_directories(. ../..)

# Specify libraries to link against
set(SUNLINSOL_SPFGMR_LIB sundials_sunlinsolspfgmr)
set(NVECS_LIB sundials_nvecserial)
if (BUILD_FORTRAN_MODULE_INTERFACE)
  list(APPEND NVECS_LIB sundials_fnvecserial_mod)
  list(APPEND SUNLINSOL_SPFGMR_LIB sundials_fsunlinsolspfgmr_mod)
endif()

# Set-up linker flags and link libraries
set(SUNDIALS_LIBS ${SUNLINSOL_SPFGMR_LIB} ${NVECS_LIB} ${EXE_EXTRA_LINK_LIBS})

# Add the build and install targets for each example
foreach(example_tuple ${sunlinsol_spfgmr_examples})

  # parse the example tuple using sundials_parse_tuple macro
  sundials_parse_tuple( sunlinsol_spfgmr_examples 
    FILENAME_VAR example
    PARAMS_VAR example_args
    TYPE_VAR example_type
  )

  # check if this example has already been added, only need to add
  # example source files once for testing with different inputs
  if(NOT TARGET ${example})
    # example source files
    add_executable(${example} ${example})

    # link vector test utilities 
    target_link_libraries(${example} PRIVATE test_sunlinsol_obj)

    # libraries to link against
    target_link_libraries(${example} PRIVATE ${SUNDIALS_LIBS})

    # folder to organize targets in an IDE
    set_target_properties(${example} PROPERTIES FOLDER "Examples")

  endif()

  # check if example args are provided and set the test name
  sundials_create_test_name(${example} test_name
    TEST_ARGS ${example_args})

  # add example to regression tests
  sundials_add_test(${test_name} ${example}
    TEST_ARGS ${example_args}
    EXAMPLE_TYPE ${example_type}
    NODIFF)

endforeach(example_tuple ${sunlinsol_spfgmr_examples})

# Add the build and install targets for each example
foreach(example_tuple ${sunlinsol_spfgmr_fortran_examples})

  # parse the example tuple using sundials_parse_tuple macro
  sundials_parse_tuple( sunlinsol_spfgmr_fortran_examples 
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
      ${SUNDIALS_SOURCE_DIR}/examples/sunlinsol/test_sunlinsol.f90)

    # folder to organize targets in an IDE
    set_target_properties(${example} PROPERTIES FOLDER "Examples")

    # set fortran module directory to avoid name collisions
    set_target_properties(${example} PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

    # libraries to link against
    target_link_libraries(${example} ${SUNDIALS_LIBS})
  endif()

  # check if example args are provided and set the test name
  sundials_create_test_name(${example} test_name
    TEST_ARGS ${example_args})

    # add example to regression tests
    sundials_add_test(${test_name} ${example}
      TEST_ARGS ${example_args}
      EXAMPLE_TYPE ${example_type}
      NODIFF)

endforeach(example_tuple ${sunlinsol_spfgmr_fortran_examples})

# create Makfile and CMakeLists.txt for examples
if(EXAMPLES_INSTALL)

  sundials_install_examples(sunlinsol sunlinsol_spfgmr_examples
    CMAKE_TEMPLATE
      cmakelists_C_ex.in
    MAKE_TEMPLATE
      makefile_serial_C_ex.in
    SUNDIALS_TARGETS
      nvecserial
      sunlinsolspfgmr
    DESTINATION
      sunlinsol/spfgmr_serial
    EXAMPLES_DEPENDENCIES
      ${sunlinsol_spfgmr_dependencies}
    EXTRA_FILES
      #README
      ${sunlinsol_spfgmr_files}
    TEST_INSTALL
      spfgmr_serial
  )

endif()