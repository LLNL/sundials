# ---------------------------------------------------------------
# Programmer(s): Daniel R. Reynolds @ SMU
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
# CMakeLists.txt file for ARKODE serial examples
# ---------------------------------------------------------------

# Example lists are tuples "name\;args\;type" where the type is
# 'develop' for examples excluded from 'make test' in releases

# Examples using SUNDIALS linear solvers
set(ARKODE_examples
  # tests that are always run
  "ark_analytic\;\;"
  # develop tests
  "ark_analytic_mels\;\;develop"
  "ark_analytic_nonlin\;\;develop"
  "ark_brusselator_1D_mri\;\;develop"
  "ark_brusselator_fp\;\;exclude-single"
  "ark_brusselator_mri\;\;develop"
  "ark_brusselator\;\;develop"
  "ark_brusselator1D_imexmri\;0 0.001\;exclude-single"
  "ark_brusselator1D_imexmri\;2 0.001\;exclude-single"
  "ark_brusselator1D_imexmri\;3 0.001\;exclude-single"
  "ark_brusselator1D_imexmri\;4 0.001\;exclude-single"
  "ark_brusselator1D_imexmri\;5 0.001\;exclude-single"
  "ark_brusselator1D_imexmri\;6 0.001\;exclude-single"
  "ark_brusselator1D_imexmri\;7 0.001\;exclude-single"
  "ark_brusselator1D\;\;exclude-single"
  "ark_conserved_exp_entropy_ark\;1 0\;develop"
  "ark_conserved_exp_entropy_ark\;1 1\;develop"
  "ark_conserved_exp_entropy_erk\;1\;develop"
  "ark_damped_harmonic_symplectic\;\;exclude-single"
  "ark_dissipated_exp_entropy\;1 0\;develop"
  "ark_dissipated_exp_entropy\;1 1\;develop"
  "ark_harmonic_symplectic\;\;exclude-single"
  "ark_heat1D_adapt\;\;develop"
  "ark_heat1D\;\;develop"
  "ark_kepler\;--stepper ERK --step-mode adapt\;develop"
  "ark_kepler\;--stepper ERK --step-mode fixed --count-orbits\;develop"
  "ark_kepler\;--stepper SPRK --step-mode fixed --count-orbits --use-compensated-sums\;develop"
  "ark_kepler\;--stepper SPRK --step-mode fixed --method ARKODE_SPRK_EULER_1_1 --tf 50 --check-order --nout 1\;exclude-single"
  "ark_kepler\;--stepper SPRK --step-mode fixed --method ARKODE_SPRK_LEAPFROG_2_2 --tf 50 --check-order --nout 1\;exclude-single"
  "ark_kepler\;--stepper SPRK --step-mode fixed --method ARKODE_SPRK_MCLACHLAN_2_2 --tf 50 --check-order --nout 1\;exclude-single"
  "ark_kepler\;--stepper SPRK --step-mode fixed --method ARKODE_SPRK_MCLACHLAN_3_3 --tf 50 --check-order --nout 1\;exclude-single"
  "ark_kepler\;--stepper SPRK --step-mode fixed --method ARKODE_SPRK_MCLACHLAN_4_4 --tf 50 --check-order --nout 1\;exclude-single"
  "ark_kepler\;--stepper SPRK --step-mode fixed --method ARKODE_SPRK_MCLACHLAN_5_6 --tf 50 --check-order --nout 1\;exclude-single"
  "ark_kepler\;--stepper SPRK --step-mode fixed --method ARKODE_SPRK_PSEUDO_LEAPFROG_2_2 --tf 50 --check-order --nout 1\;exclude-single"
  "ark_kepler\;--stepper SPRK --step-mode fixed --method ARKODE_SPRK_RUTH_3_3 --tf 50 --check-order --nout 1\;exclude-single"
  "ark_kepler\;--stepper SPRK --step-mode fixed --method ARKODE_SPRK_YOSHIDA_6_8 --tf 50 --check-order --nout 1\;exclude-single"
  "ark_kepler\;\;develop"
  "ark_kpr_mri\;\;develop"
  "ark_kpr_mri\;0 0.002\;develop"
  "ark_kpr_mri\;1 0.002\;develop"
  "ark_kpr_mri\;2 0.005\;develop"
  "ark_kpr_mri\;3 0.01\;develop"
  "ark_kpr_mri\;4 0.002\;develop"
  "ark_kpr_mri\;5 0.002\;develop"
  "ark_kpr_mri\;6 0.005\;develop"
  "ark_kpr_mri\;7 0.001 -100 100 0.5 1\;exclude-single"
  "ark_kpr_mri\;7 0.001\;exclude-single"
  "ark_kpr_mri\;8 0.001 -100 100 0.5 1\;exclude-single"
  "ark_kpr_mri\;8 0.001\;exclude-single"
  "ark_kpr_mri\;9 0.001 -100 100 0.5 1\;exclude-single"
  "ark_kpr_mri\;9 0.001\;exclude-single"
  "ark_KrylovDemo_prec\;\;exclude-single"
  "ark_KrylovDemo_prec\;1\;exclude-single"
  "ark_KrylovDemo_prec\;2\;exclude-single"
  "ark_onewaycouple_mri\;\;develop"
  "ark_reaction_diffusion_mri\;\;develop"
  "ark_robertson_constraints\;\;exclude-single"
  "ark_robertson_root\;\;exclude-single"
  "ark_robertson\;\;exclude-single"
  "ark_twowaycouple_mri\;\;develop"
  )

if(SUNDIALS_BUILD_WITH_MONITORING)
  list(APPEND ARKODE_examples "ark_brusselator_fp\;1\;exclude-single")
endif()

# Examples using KLU linear solver
set(ARKODE_examples_KLU
  "ark_brusselator1D_klu\;develop"
  )

# Examples using SuperLU_MT linear solver
set(ARKODE_examples_SUPERLUMT
  "ark_brusselator1D_FEM_slu\;exclude-single"
  )

# Auxiliary files to install
set(ARKODE_extras
  ark_analytic_nonlin_stats.csv
  ark_damped_harmonic_symplectic.h
  ark_harmonic_symplectic.h
  ark_kepler_plot.py
  ark_kepler.h
  ark_reaction_diffusion_mri_fast_stats.csv
  ark_reaction_diffusion_mri_slow_stats.csv
  ark_robertson_stats.csv
  plot_brusselator1D_FEM.py
  plot_brusselator1D.py
  plot_heat1D_adapt.py
  plot_heat1D.py
  plot_sol_log.py
  plot_sol.py
  )

# Add the build and install targets for each ARKODE example
foreach(example_tuple ${ARKODE_examples})

  # parse the example tuple
  list(GET example_tuple 0 example)
  list(GET example_tuple 1 example_args)
  list(GET example_tuple 2 example_type)

  if (NOT TARGET ${example})
    # example source files
    add_executable(${example} ${example}.c)

    # folder for IDEs
    set_target_properties(${example} PROPERTIES FOLDER "Examples")

    # libraries to link against
    target_link_libraries(${example}
      sundials_arkode
      sundials_nvecserial
      ${EXE_EXTRA_LINK_LIBS})
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

  # install example source and out files
  if(EXAMPLES_INSTALL)
    install(FILES ${example}.c ${test_name}.out
      DESTINATION ${EXAMPLES_INSTALL_PATH}/arkode/C_serial)
  endif()

endforeach(example_tuple ${ARKODE_examples})


# Add the build and install targets for each KLU example (if needed)
if(BUILD_SUNLINSOL_KLU)

  foreach(example_tuple ${ARKODE_examples_KLU})

    # parse the example tuple
    list(GET example_tuple 0 example)
    list(GET example_tuple 1 example_type)

    # example source files
    add_executable(${example} ${example}.c)

    set_target_properties(${example} PROPERTIES FOLDER "Examples")

    # add example to regression tests
    sundials_add_test(${example} ${example}
      ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
      ANSWER_FILE ${example}.out
      EXAMPLE_TYPE ${example_type})

    # libraries to link against
    target_link_libraries(${example}
      sundials_arkode
      sundials_nvecserial
      sundials_sunlinsolklu
      ${EXE_EXTRA_LINK_LIBS})

    # install example source and out files
    if(EXAMPLES_INSTALL)
      install(FILES ${example}.c ${example}.out
        DESTINATION ${EXAMPLES_INSTALL_PATH}/arkode/C_serial)
    endif()

  endforeach(example_tuple ${ARKODE_examples_KLU})

endif()


# Add the build and install targets for each SuperLU_MT example (if needed)
if(BUILD_SUNLINSOL_SUPERLUMT)

  foreach(example_tuple ${ARKODE_examples_SUPERLUMT})

    # parse the example tuple
    list(GET example_tuple 0 example)
    list(GET example_tuple 1 example_type)

    # example source files
    add_executable(${example} ${example}.c)

    set_target_properties(${example} PROPERTIES FOLDER "Examples")

    # add example to regression tests
    sundials_add_test(${example} ${example}
      ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
      ANSWER_FILE ${example}.out
      EXAMPLE_TYPE ${example_type})

    # libraries to link against
    target_link_libraries(${example}
      sundials_arkode
      sundials_nvecserial
      sundials_sunlinsolsuperlumt
      ${EXE_EXTRA_LINK_LIBS})

    # install example source and out files
    if(EXAMPLES_INSTALL)
      install(FILES ${example}.c ${example}.out
        DESTINATION ${EXAMPLES_INSTALL_PATH}/arkode/C_serial)
    endif()

  endforeach(example_tuple ${ARKODE_examples_SUPERLUMT})

endif()


# create Makfile and CMakeLists.txt for examples
if(EXAMPLES_INSTALL)

  # Install the README file
  install(FILES README DESTINATION ${EXAMPLES_INSTALL_PATH}/arkode/C_serial)

  # Install the extra files
  foreach(extrafile ${ARKODE_extras})
    install(FILES ${extrafile} DESTINATION ${EXAMPLES_INSTALL_PATH}/arkode/C_serial)
  endforeach(extrafile ${ARKODE_extras})

  # Prepare substitution variables for Makefile and/or CMakeLists templates
  set(SOLVER "ARKODE")
  set(SOLVER_LIB "sundials_arkode")

  examples2string(ARKODE_examples EXAMPLES)

  if(BUILD_SUNLINSOL_KLU)
    examples2string(ARKODE_examples_KLU EXAMPLES_KLU)
  else()
    set(EXAMPLES_KLU "")
  endif()

  if(BUILD_SUNLINSOL_SUPERLUMT)
    examples2string(ARKODE_examples_SUPERLUMT EXAMPLES_SLUMT)
    if (SUNDIALS_SUPERLUMT_THREAD_TYPE STREQUAL "PTHREAD")
      set(THREAD_LIBRARY_SLUMT ${CMAKE_THREAD_LIBS_INIT})
    else()
      set(THREAD_LIBRARY_SLUMT "")
    endif()
  else()
    set(EXAMPLES_SLUMT "")
    set(THREAD_LIBRARY_SLUMT "")
  endif()

  # Regardless of the platform we're on, we will generate and install
  # CMakeLists.txt file for building the examples. This file  can then
  # be used as a template for the user's own programs.

  # generate CMakelists.txt in the binary directory
  configure_file(
    ${PROJECT_SOURCE_DIR}/examples/templates/cmakelists_serial_C_ex.in
    ${PROJECT_BINARY_DIR}/examples/arkode/C_serial/CMakeLists.txt
    @ONLY
    )

  # install CMakelists.txt
  install(
    FILES ${PROJECT_BINARY_DIR}/examples/arkode/C_serial/CMakeLists.txt
    DESTINATION ${EXAMPLES_INSTALL_PATH}/arkode/C_serial
    )

  # On UNIX-type platforms, we also  generate and install a makefile for
  # building the examples. This makefile can then be used as a template
  # for the user's own programs.

  if(UNIX)
    # generate Makefile and place it in the binary dir
    configure_file(
      ${PROJECT_SOURCE_DIR}/examples/templates/makefile_serial_C_ex.in
      ${PROJECT_BINARY_DIR}/examples/arkode/C_serial/Makefile_ex
      @ONLY
      )
    # install the configured Makefile_ex as Makefile
    install(
      FILES ${PROJECT_BINARY_DIR}/examples/arkode/C_serial/Makefile_ex
      DESTINATION ${EXAMPLES_INSTALL_PATH}/arkode/C_serial
      RENAME Makefile
      )
  endif()

  # add test_install target
  sundials_add_test_install(arkode C_serial EXECUTABLE ark_analytic)

endif()
