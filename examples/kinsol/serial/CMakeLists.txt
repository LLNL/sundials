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
# CMakeLists.txt file for KINSOL serial examples
# ---------------------------------------------------------------

# Example lists are tuples "name\;args\;type" where the type is
# 'develop' for examples excluded from 'make test' in releases

# Examples using SUNDIALS linear solvers
set(KINSOL_examples
  "kinAnalytic_fp.c\;<none>\;passfail"
  "kinAnalytic_fp.c\;--damping_fp 0.5\;develop"
  "kinAnalytic_fp.c\;--m_aa 2\;passfail"
  "kinAnalytic_fp.c\;--m_aa 2 --delay_aa 2\;passfail"
  "kinAnalytic_fp.c\;--m_aa 2 --damping_aa 0.5\;develop"
  "kinAnalytic_fp.c\;--m_aa 2 --orth_aa 1\;passfail"
  "kinAnalytic_fp.c\;--m_aa 2 --orth_aa 2\;passfail"
  "kinAnalytic_fp.c\;--m_aa 2 --orth_aa 3\;passfail"
  "kinFerTron_dns.c\;<none>\;develop"
  "kinFoodWeb_kry.c\;<none>\;exclude-single"
  "kinKrylovDemo_ls.c\;<none>\;exclude-single"
  "kinLaplace_bnd.c\;<none>\;exclude-single"
  "kinLaplace_picard_bnd.c\;<none>\;exclude-single"
  "kinLaplace_picard_kry.c\;<none>\;exclude-single"
  "kinRoberts_fp.c\;<none>\;develop"
  "kinRoboKin_dns.c\;<none>\;exclude-single"
  )

# Examples using LAPACK linear solvers
set(KINSOL_examples_BL
  )

# Examples using KLU linear solver
set(KINSOL_examples_KLU
  "kinFerTron_klu.c\;<none>\;develop"
  )

# Examples using SuperLU_MT linear solver
set(KINSOL_examples_SUPERLUMT
  "kinRoboKin_slu.c\;<none>\;develop"
  )

# Auxiliary files to install
set(KINSOL_extras
  kinRoboKin_dns_stats.csv
  )

# Specify libraries to link against
set(KINSOL_LIB sundials_kinsol)
set(NVECS_LIB sundials_nvecserial)

# Set-up linker flags and link libraries
set(SUNDIALS_LIBS ${KINSOL_LIB} ${NVECS_LIB} ${EXE_EXTRA_LINK_LIBS})

# Add the build and install targets for each example
foreach(example_tuple ${KINSOL_examples})

    # parse the example tuple using sundials_parse_tuple macro
    sundials_parse_tuple( KINSOL_examples 
      FILENAME_VAR example
      PARAMS_VAR example_args
      TYPE_VAR example_type
    )

  # check if this example has already been added, only need to add
  # example source files once for testing with different inputs
  if(NOT TARGET ${example})

    # example source files
    add_executable(${example} ${example})

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

endforeach(example_tuple ${KINSOL_examples})

# Add the build and install targets for each LAPACK example (if needed)
if(BUILD_SUNLINSOL_LAPACKBAND AND BUILD_SUNLINSOL_LAPACKDENSE)

  # Sundials LAPACK linear solver modules
  set(SUNLINSOLLAPACK_LIBS
    sundials_sunlinsollapackband
    sundials_sunlinsollapackdense)

  # LAPACK libraries
  list(APPEND SUNLINSOLLAPACK_LIBS ${LAPACK_LIBRARIES})

  foreach(example_tuple ${KINSOL_examples_BL})
    # parse the example tuple using sundials_parse_tuple macro
    sundials_parse_tuple( KINSOL_examples 
      FILENAME_VAR example
      PARAMS_VAR example_args
      TYPE_VAR example_type
    )
    # example source files
    add_executable(${example} ${example})

    set_target_properties(${example} PROPERTIES FOLDER "Examples")

    # add example to regression tests
    sundials_add_test(${example} ${example}
      ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
      ANSWER_FILE ${example}.out
      EXAMPLE_TYPE ${example_type})

    # libraries to link against
    target_link_libraries(${example} ${SUNDIALS_LIBS} ${SUNLINSOLLAPACK_LIBS})

  endforeach(example_tuple ${KINSOL_examples_BL})

endif()

# Add the build and install targets for each KLU example (if needed)
if(BUILD_SUNLINSOL_KLU)

  # Sundials KLU linear solver module
  set(SUNLINSOLKLU_LIBS sundials_sunlinsolklu)

  # KLU libraries
  list(APPEND SUNLINSOLKLU_LIBS)

  foreach(example_tuple ${KINSOL_examples_KLU})

    # parse the example tuple using sundials_parse_tuple macro
    sundials_parse_tuple( KINSOL_examples 
      FILENAME_VAR example
      PARAMS_VAR example_args
      TYPE_VAR example_type
    )

    # example source files
    add_executable(${example} ${example})

    set_target_properties(${example} PROPERTIES FOLDER "Examples")

    # add example to regression tests
    sundials_add_test(${example} ${example}
      ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
      ANSWER_FILE ${example}.out
      EXAMPLE_TYPE ${example_type})

    # libraries to link against
    target_link_libraries(${example} ${SUNDIALS_LIBS} ${SUNLINSOLKLU_LIBS})

  endforeach(example_tuple ${KINSOL_examples_KLU})

endif()

# Add the build and install targets for each SuperLU_MT example (if needed)
if(BUILD_SUNLINSOL_SUPERLUMT)

  # Sundials SuperLU_MT linear solver module
  set(SUNLINSOLSLUMT_LIBS sundials_sunlinsolsuperlumt)

  # SuperLU_MT libraries
  list(APPEND SUNLINSOLSLUMT_LIBS ${SUPERLUMT_LIBRARIES})

  foreach(example_tuple ${KINSOL_examples_SUPERLUMT})

    # parse the example tuple using sundials_parse_tuple macro
    sundials_parse_tuple( KINSOL_examples 
      FILENAME_VAR example
      PARAMS_VAR example_args
      TYPE_VAR example_type
    )

    # example source files
    add_executable(${example} ${example})

    set_target_properties(${example} PROPERTIES FOLDER "Examples")

    # add example to regression tests
    sundials_add_test(${example} ${example}
      ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
      ANSWER_FILE ${example}.out
      EXAMPLE_TYPE ${example_type})

    # libraries to link against
    target_link_libraries(${example} ${SUNDIALS_LIBS} ${SUNLINSOLSLUMT_LIBS})

  endforeach(example_tuple ${KINSOL_examples_SUPERLUMT})

endif()

# create Makfile and CMakeLists.txt for examples
if(EXAMPLES_INSTALL)

  if(BUILD_SUNLINSOL_SUPERLUMT)
    set(_include_sunlinsolsuperlumt sunlinsolsuperlumt)
    list(APPEND KINSOL_examples 
      "kinRoboKin_slu.c\;<none>\;develop"
    )
    set(SUPERLUMT_TARGET SUNDIALS::SUPERLUMT)
  endif()

  if(BUILD_SUNLINSOL_KLU)
    set(_include_sunlinsolklu sunlinsolklu)
    list(APPEND KINSOL_examples
      "kinFerTron_klu.c\;<none>\;develop"
    )
    set(KLU_TARGET SUNDIALS::KLU)
  endif()

  if(BUILD_SUNLINSOL_LAPACKBAND AND BUILD_SUNLINSOL_LAPACKDENSE)
    set(_include_sundials_sunlinsollapackband sunlinsollapackband)
    set(_include_sundials_sunlinsollapackdense sunlinsollapackdense)
  endif()

  sundials_install_examples(kinsol KINSOL_examples
    CMAKE_TEMPLATE
      cmakelists_C_ex.in
    MAKE_TEMPLATE
      makefile_serial_C_ex.in
    SUNDIALS_TARGETS
      nvecserial
      kinsol
      ${_include_sunlinsolsuperlumt}
      ${_include_sunlinsolklu}
      ${_include_sundials_sunlinsollapackband}
      ${_include_sundials_sunlinsollapackdense}
    OTHER_TARGETS
      ${SUPERLUMT_TARGET}
      ${KLU_TARGET}
    DESTINATION
      kinsol/serial
    EXTRA_FILES
      README
      ${KINSOL_extras}
      ${KINSOL_headers}
    TEST_INSTALL
      serial
  )

endif()