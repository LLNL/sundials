# ---------------------------------------------------------------
# Programmer(s): Radu Serban @ LLNL
#                David J. Gardner @ LLNL
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
# CMakeLists.txt file for CVODE serial examples
# ---------------------------------------------------------------

# Example lists are tuples "name\;args\;type" where the type is
# 'develop' for examples excluded from 'make test' in releases

# Examples using SUNDIALS linear solvers
set(CVODE_examples
  "cvAdvDiff_bnd.c\;<none>\;develop"
  "cvAnalytic_mels.c\;<none>\;develop"
  "cvDirectDemo_ls.c\;<none>\;develop"
  "cvDisc_dns.c\;<none>\;develop"
  "cvDiurnal_kry_bp.c\;<none>\;develop"
  "cvDiurnal_kry.c\;<none>\;develop"
  "cvKrylovDemo_ls.c\;<none>\;develop"
  "cvKrylovDemo_ls.c\;1\;develop"
  "cvKrylovDemo_ls.c\;2\;develop"
  "cvKrylovDemo_prec.c\;<none>\;develop"
  "cvParticle_dns.c\;<none>\;develop"
  "cvPendulum_dns.c\;<none>\;exclude-single"
  "cvRoberts_dns.c\;<none>\;passfail"
  "cvRoberts_dns_constraints.c\;<none>\;develop"
  "cvRoberts_dns_negsol.c\;<none>\;exclude-single"
  "cvRoberts_dns_uw.c\;<none>\;develop"
  "cvRocket_dns.c\;<none>\;develop"
  )

if(SUNDIALS_BUILD_WITH_MONITORING)
  list(APPEND CVODE_examples "cvKrylovDemo_ls.c\;0 1\;develop")
endif()

# Examples using LAPACK linear solvers
set(CVODE_examples_BL
  "cvAdvDiff_bndL.c\;<none>\;exclude-single"
  "cvRoberts_dnsL.c\;<none>\;exclude-single"
  )

# Examples using KLU linear solver
set(CVODE_examples_KLU
  "cvRoberts_block_klu.c\;<none>\;develop"
  "cvRoberts_klu.c\;<none>\;develop"
  )

# Examples using SuperLU_MT linear solver
set(CVODE_examples_SUPERLUMT
  "cvRoberts_sps.c\;<none>\;develop"
  )

# Auxiliary files to install
set(CVODE_extras
  plot_cvParticle.py
  plot_cvPendulum.py
  cvRoberts_dns_stats.csv
  )

# Specify libraries to link against
set(CVODE_LIB sundials_cvode)
set(NVECS_LIB sundials_nvecserial)

if(SUNDIALS_BUILD_PACKAGE_FUSED_KERNELS)
  list(APPEND CVODE_LIB sundials_cvode_fused_stubs)
endif()

# Set-up linker flags and link libraries
set(SUNDIALS_LIBS ${CVODE_LIB} ${NVECS_LIB} ${EXE_EXTRA_LINK_LIBS})

# Add the build and install targets for each example
foreach(example_tuple ${CVODE_examples})

  # parse the example tuple using sundials_parse_tuple macro
  sundials_parse_tuple( CVODE_examples 
    FILENAME_VAR example
    PARAMS_VAR example_args
    TYPE_VAR example_type
  ) 

  # check if this example has already been added, only need to add
  # example source files once for testing with different inputs
  if(NOT TARGET ${example})
    # example source files
    add_executable(${example} ${example})

    # folder to organize targets in an IDE
    set_target_properties(${example} PROPERTIES FOLDER "Examples")

    # libraries to link against
    target_link_libraries(${example} PRIVATE ${SUNDIALS_LIBS})
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
    ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    ANSWER_FILE ${test_name}.out
    EXAMPLE_TYPE ${example_type})

  # find all .out files for this example
  file(GLOB example_out ${example}*.out)

endforeach(example_tuple ${CVODE_examples})

# Add the build and install targets for each LAPACK example (if needed)
if(BUILD_SUNLINSOL_LAPACKBAND AND BUILD_SUNLINSOL_LAPACKDENSE)

  # Sundials LAPACK linear solver modules
  set(SUNLINSOLLAPACK_LIBS
    sundials_sunlinsollapackband
    sundials_sunlinsollapackdense)

  # LAPACK libraries
  list(APPEND SUNLINSOLLAPACK_LIBS ${LAPACK_LIBRARIES})

  foreach(example_tuple ${CVODE_examples_BL})

    # parse the example tuple using sundials_parse_tuple macro
    sundials_parse_tuple( CVODE_examples 
      FILENAME_VAR example
      PARAMS_VAR example_args
      TYPE_VAR example_type
    ) 

    # check if this example has already been added, only need to add
    # example source files once for testing with different inputs
    if(NOT TARGET ${example})
      # example source files
      add_executable(${example} ${example})

      # folder to organize targets in an IDE
      set_target_properties(${example} PROPERTIES FOLDER "Examples")

      # libraries to link against
      target_link_libraries(${example} ${SUNDIALS_LIBS} ${SUNLINSOLLAPACK_LIBS})
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
      ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
      ANSWER_FILE ${test_name}.out
      EXAMPLE_TYPE ${example_type})

    # find all .out files for this example
    file(GLOB example_out ${example}*.out)

  endforeach(example_tuple ${CVODE_examples_BL})

endif()

# Add the build and install targets for each KLU example (if needed)
if(BUILD_SUNLINSOL_KLU)

  # Sundials KLU linear solver module
  set(SUNLINSOLKLU_LIBS sundials_sunlinsolklu)

  # KLU libraries
  list(APPEND SUNLINSOLKLU_LIBS)

  foreach(example_tuple ${CVODE_examples_KLU})

    # parse the example tuple using sundials_parse_tuple macro
    sundials_parse_tuple( CVODE_examples 
      FILENAME_VAR example
      PARAMS_VAR example_args
      TYPE_VAR example_type
    ) 

    # check if this example has already been added, only need to add
    # example source files once for testing with different inputs
    if(NOT TARGET ${example})
      # add example source files
      add_executable(${example} ${example})

      # folder to organize targets in an IDE
      set_target_properties(${example} PROPERTIES FOLDER "Examples")

      # libraries to link against
      target_link_libraries(${example} ${SUNDIALS_LIBS} ${SUNLINSOLKLU_LIBS})
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
      ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
      ANSWER_FILE ${test_name}.out
      EXAMPLE_TYPE ${example_type})

    # find all .out files for this example
    file(GLOB example_out ${example}*.out)

  endforeach(example_tuple ${CVODE_examples_KLU})

endif()

# Add the build and install targets for each SuperLU_MT example (if needed)
if(BUILD_SUNLINSOL_SUPERLUMT)

  # Sundials SuperLU_MT linear solver module
  set(SUNLINSOLSLUMT_LIBS sundials_sunlinsolsuperlumt)

  # SuperLU_MT libraries
  list(APPEND SUNLINSOLSLUMT_LIBS ${SUPERLUMT_LIBRARIES})

  foreach(example_tuple ${CVODE_examples_SUPERLUMT})

    # parse the example tuple using sundials_parse_tuple macro
    sundials_parse_tuple( CVODE_examples 
      FILENAME_VAR example
      PARAMS_VAR example_args
      TYPE_VAR example_type
    ) 

    # check if this example has already been added, only need to add
    # example source files once for testing with different inputs
    if(NOT TARGET ${example})
      # add example source files
      add_executable(${example} ${example})

      # folder to organize targets in an IDE
      set_target_properties(${example} PROPERTIES FOLDER "Examples")

      # libraries to link against
      target_link_libraries(${example} ${SUNDIALS_LIBS} ${SUNLINSOLSLUMT_LIBS})
    endif()

    # check if example args are provided and set the test name
    if("${example_args}" STREQUAL "")
      set(test_name ${example})
    else()
      string(REGEX REPLACE " " "_" test_name ${example}_${example_args})
    endif()

    # Do not include SuperLUMT examples in testing when the indextype is int64_t.
    # Answer files were generated with int32_t and minor differences in output
    # occur causing a false positive when testing. These tests can be re-enabled
    # when type specific answer files are added.
    if(SUNDIALS_INDEX_SIZE MATCHES "32")
      # add example to regression tests
      sundials_add_test(${test_name} ${example}
        TEST_ARGS ${example_args}
        ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
        ANSWER_FILE ${test_name}.out
        EXAMPLE_TYPE ${example_type})
    endif()

    # find all .out files for this example
    file(GLOB example_out ${example}*.out)

  endforeach(example_tuple ${CVODE_examples_SUPERLUMT})

endif()

# create Makfile and CMakeLists.txt for examples
if(EXAMPLES_INSTALL)

  if(BUILD_SUNLINSOL_SUPERLUMT)
    set(_include_sunlinsolsuperlumt sunlinsolsuperlumt)
    list(APPEND CVODE_examples "cvRoberts_sps.c\;<none>\;develop")
    set(SUPERLUMT_TARGET SUNDIALS::SUPERLUMT)
  endif()

  if(BUILD_SUNLINSOL_KLU)
    set(_include_sunlinsolklu sunlinsolklu)
    list(APPEND CVODE_examples "cvRoberts_block_klu.c\;<none>\;develop" "cvRoberts_klu.c\;<none>\;develop" )
    set(KLU_TARGET SUNDIALS::KLU)
  endif()

  if(BUILD_SUNLINSOL_LAPACKBAND AND BUILD_SUNLINSOL_LAPACKDENSE)
    set(_include_sundials_sunlinsollapackband sunlinsollapackband)
    set(_include_sundials_sunlinsollapackdense sunlinsollapackdense)
    list(APPEND CVODE_examples "cvAdvDiff_bndL.c\;<none>\;exclude-single" "cvRoberts_dnsL.c\;<none>\;exclude-single")
  endif()

  if(SUNDIALS_BUILD_PACKAGE_FUSED_KERNELS)
    set(_include_cvode_fused_stubs cvode_fused_stubs)
  endif()

  sundials_install_examples(cvode CVODE_examples
    CMAKE_TEMPLATE
      cmakelists_C_ex.in
    MAKE_TEMPLATE
      makefile_serial_C_ex.in
    SUNDIALS_TARGETS
      nvecserial
      cvode
      ${_include_sunlinsolsuperlumt}
      ${_include_sunlinsolklu}
      ${_include_sundials_sunlinsollapackband}
      ${_include_sundials_sunlinsollapackdense}
      ${cvode_fused_stubs}
    DESTINATION
      cvode/serial
    OTHER_TARGETS
      ${SUPERLUMT_TARGET}
      ${KLU_TARGET}
    EXTRA_FILES
      README
      ${CVODE_extras}
      ${CVODE_headers}
    TEST_INSTALL
      serial
  )

endif()