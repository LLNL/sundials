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
# CMakeLists.txt file for CVODES serial examples
# ---------------------------------------------------------------

# Example lists are tuples "name\;args\;type" where the type is
# 'develop' for examples excluded from 'make test' in releases

# Examples using SUNDIALS linear solvers
set(CVODES_examples
  "cvsAdvDiff_FSA_non.c\; \"-sensi sim t\" \;exclude-single"
  "cvsAdvDiff_FSA_non.c\; \"-sensi stg t\" \;exclude-single"
  "cvsAdvDiff_bnd.c\; \"\" \;develop"
  "cvsAnalytic_mels.c\; \"\" \;develop"
  "cvsDirectDemo_ls.c\; \"\" \;develop"
  "cvsDiurnal_FSA_kry.c\; \"-sensi sim t\" \;exclude-single"
  "cvsDiurnal_FSA_kry.c\; \"-sensi stg t\" \;exclude-single"
  "cvsDiurnal_kry.c\; \"\" \;develop"
  "cvsDiurnal_kry_bp.c\; \"\" \;develop"
  "cvsFoodWeb_ASAi_kry.c\; \"\" \;develop"
  "cvsFoodWeb_ASAp_kry.c\; \"\" \;develop"
  "cvsHessian_ASA_FSA.c\; \"\" \;exclude-single"
  "cvsKrylovDemo_ls.c\; \"\" \;develop"
  "cvsKrylovDemo_ls.c\; \"1\" \;develop"
  "cvsKrylovDemo_ls.c\; \"2\" \;develop"
  "cvsKrylovDemo_prec.c\; \"\" \;develop"
  "cvsParticle_dns.c\; \"\" \;develop"
  "cvsPendulum_dns.c\; \"\" \;exclude-single"
  "cvsRoberts_ASAi_dns.c\; \"\" \;exclude-single"
  "cvsRoberts_ASAi_dns_constraints.c\; \"\" \;develop"
  "cvsRoberts_FSA_dns.c\; \"-sensi sim t\" \;exclude-single"
  "cvsRoberts_FSA_dns.c\; \"-sensi stg1 t\" \;exclude-single"
  "cvsRoberts_FSA_dns_Switch.c\; \"\" \;exclude-single"
  "cvsRoberts_FSA_dns_constraints.c\; \"-sensi stg1 t\" \;exclude-single"
  "cvsRoberts_dns.c\; \"\" \; \"\" "
  "cvsRoberts_dns_constraints.c\; \"\" \;develop"
  "cvsRoberts_dns_uw.c\; \"\" \;develop"
  )

if(SUNDIALS_BUILD_WITH_MONITORING)
  list(APPEND CVODE_examples "cvsKrylovDemo_ls.c\; \"0 1\" \;develop")
endif()

# Examples using LAPACK linear solvers
set(CVODES_examples_BL
  "cvsAdvDiff_bndL.c\; \"\" \;exclude-single"
  "cvsRoberts_dnsL.c\; \"\" \;exclude-single"
  )

# Examples using KLU linear solver
set(CVODES_examples_KLU
  "cvsRoberts_ASAi_klu.c\; \"\" \;develop"
  "cvsRoberts_FSA_klu.c\; \"-sensi stg1 t\" \;develop"
  "cvsRoberts_klu.c\; \"\" \;develop"
  )

# Examples using SuperLU_MT linear solver
set(CVODES_examples_SUPERLUMT
  "cvsRoberts_ASAi_sps.c\; \"\" \;exclude-single"
  "cvsRoberts_FSA_sps.c\; \"-sensi stg1 t\" \;exclude-single"
  "cvsRoberts_sps.c\; \"\" \;develop"
  )

# Auxiliary files to install
set(CVODES_extras
  plot_cvsParticle.py
  plot_cvsPendulum.py
  cvsRoberts_ASAi_dns_bkw1_stats.csv
  cvsRoberts_ASAi_dns_bkw2_stats.csv
  cvsRoberts_ASAi_dns_fwd_stats.csv
  cvsRoberts_dns_stats.csv
  cvsRoberts_FSA_dns_stats_-sensi_sim_t.csv
  cvsRoberts_FSA_dns_stats_-sensi_stg1_t.csv
  )

# Specify libraries to link against
set(CVODES_LIB sundials_cvodes)
set(NVECS_LIB sundials_nvecserial)

# Set-up linker flags and link libraries
set(SUNDIALS_LIBS ${CVODES_LIB} ${NVECS_LIB} ${EXE_EXTRA_LINK_LIBS})

# Add the build and install targets for each example
foreach(example_tuple ${CVODES_examples})

  # parse the example tuple using sundials_parse_tuple macro
  sundials_parse_tuple( CVODES_examples 
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
    target_link_libraries(${example} ${SUNDIALS_LIBS})
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

endforeach(example_tuple ${CVODES_examples})

# Add the build and install targets for each LAPACK example (if needed)
if(BUILD_SUNLINSOL_LAPACKBAND AND BUILD_SUNLINSOL_LAPACKDENSE)

  # Sundials LAPACK linear solver modules
  set(SUNLINSOLLAPACK_LIBS
    sundials_sunlinsollapackband
    sundials_sunlinsollapackdense)

  # LAPACK libraries
  list(APPEND SUNLINSOLLAPACK_LIBS ${LAPACK_LIBRARIES})

  foreach(example_tuple ${CVODES_examples_BL})

    # parse the example tuple using sundials_parse_tuple macro
    sundials_parse_tuple( CVODES_examples 
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

  endforeach(example_tuple ${CVODES_examples_BL})

endif()

# Add the build and install targets for each KLU example (if needed)
if(BUILD_SUNLINSOL_KLU)

  # Sundials KLU linear solver module
  set(SUNLINSOLKLU_LIBS sundials_sunlinsolklu)

  # KLU libraries
  list(APPEND SUNLINSOLKLU_LIBS)

  foreach(example_tuple ${CVODES_examples_KLU})
    # parse the example tuple using sundials_parse_tuple macro
    sundials_parse_tuple( CVODES_examples 
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

  endforeach(example_tuple ${CVODES_examples_KLU})

endif()

# Add the build and install targets for each SuperLU_MT example (if needed)
if(BUILD_SUNLINSOL_SUPERLUMT)

  # Sundials SuperLU_MT linear solver module
  set(SUNLINSOLSLUMT_LIBS sundials_sunlinsolsuperlumt)

  # SuperLU_MT libraries
  list(APPEND SUNLINSOLSLUMT_LIBS ${SUPERLUMT_LIBRARIES})

  foreach(example_tuple ${CVODES_examples_SUPERLUMT})

    # parse the example tuple using sundials_parse_tuple macro
    sundials_parse_tuple( CVODES_examples 
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

    # Do not include SuperLUMT examples in testing when the indextype is int32_t.
    # Answer files were generated with int64_t and minor differences in output
    # occur causing a false positive when testing. These tests can be re-enabled
    # when type specific answer files are added.
    if(SUNDIALS_INDEX_SIZE MATCHES "64")
      # add example to regression tests
      sundials_add_test(${test_name} ${example}
        TEST_ARGS ${example_args}
        ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
        ANSWER_FILE ${test_name}.out
        EXAMPLE_TYPE ${example_type})
    endif()

    # find all .out files for this example
    file(GLOB example_out ${example}*.out)

  endforeach(example_tuple ${CVODES_examples_SUPERLUMT})

endif()

# create Makfile and CMakeLists.txt for examples
if(EXAMPLES_INSTALL)

  if(BUILD_SUNLINSOL_SUPERLUMT)
    set(_include_sunlinsolsuperlumt sunlinsolsuperlumt)
    list(APPEND CVODES_examples   
      "cvsRoberts_ASAi_sps.c\; \"\" \;exclude-single"
      "cvsRoberts_FSA_sps.c\; \"-sensi stg1 t\" \;exclude-single"
      "cvsRoberts_sps.c\; \"\" \;develop"
    )
    set(SUPERLUMT_TARGET SUNDIALS::SUPERLUMT)
  endif()

  if(BUILD_SUNLINSOL_KLU)
    set(_include_sunlinsolklu sunlinsolklu)
    list(APPEND CVODES_examples
      "cvsRoberts_ASAi_klu.c\; \"\" \;develop"
      "cvsRoberts_FSA_klu.c\; \"-sensi stg1 t\" \;develop"
      "cvsRoberts_klu.c\; \"\" \;develop"
    )
  endif()

  if(BUILD_SUNLINSOL_LAPACKBAND AND BUILD_SUNLINSOL_LAPACKDENSE)
    set(_include_sundials_sunlinsollapackband sunlinsollapackband)
    set(_include_sundials_sunlinsollapackdense sunlinsollapackdense)
    list(APPEND CVODES_examples
      "cvsAdvDiff_bndL.c\; \"\" \;exclude-single"
      "cvsRoberts_dnsL.c\; \"\" \;exclude-single"
    )
  endif()

  sundials_install_examples(cvodes CVODES_examples
    CMAKE_TEMPLATE
      cmakelists_C_ex.in
    MAKE_TEMPLATE
      makefile_serial_C_ex.in
    SUNDIALS_TARGETS
      nvecserial
      cvodes
      ${_include_sunlinsolsuperlumt}
      ${_include_sunlinsolklu}
      ${_include_sundials_sunlinsollapackband}
      ${_include_sundials_sunlinsollapackdense}
    OTHER_TARGETS
      ${SUPERLUMT_TARGET}
    DESTINATION
      cvodes/serial
    EXTRA_FILES
      README
      ${CVODES_extras}
      ${CVODES_headers}
    TEST_INSTALL
      serial
  )

endif()
