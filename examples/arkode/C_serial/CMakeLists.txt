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
  "ark_analytic.c\;<none>\;passfail"
  # develop tests
  "ark_analytic_mels.c\;<none>\;develop"
  "ark_analytic_nonlin.c\;<none>\;develop"
  "ark_brusselator_1D_mri.c\;<none>\;develop"
  "ark_brusselator_fp.c\;<none>\;exclude-single"
  "ark_brusselator_mri.c\;<none>\;develop"
  "ark_brusselator.c\;<none>\;develop"
  "ark_brusselator1D_imexmri.c\;0 0.001\;exclude-single"
  "ark_brusselator1D_imexmri.c\;2 0.001\;exclude-single"
  "ark_brusselator1D_imexmri.c\;3 0.001\;exclude-single"
  "ark_brusselator1D_imexmri.c\;4 0.001\;exclude-single"
  "ark_brusselator1D_imexmri.c\;5 0.001\;exclude-single"
  "ark_brusselator1D_imexmri.c\;6 0.001\;exclude-single"
  "ark_brusselator1D_imexmri.c\;7 0.001\;exclude-single"
  "ark_brusselator1D.c\;<none>\;exclude-single"
  "ark_conserved_exp_entropy_ark.c\;1 0\;develop"
  "ark_conserved_exp_entropy_ark.c\;1 1\;develop"
  "ark_conserved_exp_entropy_erk.c\;1\;develop"
  "ark_damped_harmonic_symplectic.c\;<none>\;exclude-single"
  "ark_dissipated_exp_entropy.c\;1 0\;develop"
  "ark_dissipated_exp_entropy.c\;1 1\;develop"
  "ark_harmonic_symplectic.c\;<none>\;exclude-single"
  "ark_heat1D_adapt.c\;<none>\;develop"
  "ark_heat1D.c\;<none>\;develop"
  "ark_kepler.c\;--stepper ERK --step-mode adapt\;develop"
  "ark_kepler.c\;--stepper ERK --step-mode fixed --count-orbits\;develop"
  "ark_kepler.c\;--stepper SPRK --step-mode fixed --count-orbits --use-compensated-sums\;develop"
  "ark_kepler.c\;--stepper SPRK --step-mode fixed --method ARKODE_SPRK_EULER_1_1 --tf 50 --check-order --nout 1\;exclude-single"
  "ark_kepler.c\;--stepper SPRK --step-mode fixed --method ARKODE_SPRK_LEAPFROG_2_2 --tf 50 --check-order --nout 1\;exclude-single"
  "ark_kepler.c\;--stepper SPRK --step-mode fixed --method ARKODE_SPRK_MCLACHLAN_2_2 --tf 50 --check-order --nout 1\;exclude-single"
  "ark_kepler.c\;--stepper SPRK --step-mode fixed --method ARKODE_SPRK_MCLACHLAN_3_3 --tf 50 --check-order --nout 1\;exclude-single"
  "ark_kepler.c\;--stepper SPRK --step-mode fixed --method ARKODE_SPRK_MCLACHLAN_4_4 --tf 50 --check-order --nout 1\;exclude-single"
  "ark_kepler.c\;--stepper SPRK --step-mode fixed --method ARKODE_SPRK_MCLACHLAN_5_6 --tf 50 --check-order --nout 1\;exclude-single"
  "ark_kepler.c\;--stepper SPRK --step-mode fixed --method ARKODE_SPRK_PSEUDO_LEAPFROG_2_2 --tf 50 --check-order --nout 1\;exclude-single"
  "ark_kepler.c\;--stepper SPRK --step-mode fixed --method ARKODE_SPRK_RUTH_3_3 --tf 50 --check-order --nout 1\;exclude-single"
  "ark_kepler.c\;--stepper SPRK --step-mode fixed --method ARKODE_SPRK_YOSHIDA_6_8 --tf 50 --check-order --nout 1\;exclude-single"
  "ark_kepler.c\;<none>\;develop"
  "ark_kpr_mri.c\;<none>\;develop"
  "ark_kpr_mri.c\;0 0.002\;develop"
  "ark_kpr_mri.c\;1 0.002\;develop"
  "ark_kpr_mri.c\;2 0.005\;develop"
  "ark_kpr_mri.c\;3 0.01\;develop"
  "ark_kpr_mri.c\;4 0.002\;develop"
  "ark_kpr_mri.c\;5 0.002\;develop"
  "ark_kpr_mri.c\;6 0.005\;develop"
  "ark_kpr_mri.c\;7 0.001 -100 100 0.5 1\;exclude-single"
  "ark_kpr_mri.c\;7 0.001\;exclude-single"
  "ark_kpr_mri.c\;8 0.001 -100 100 0.5 1\;exclude-single"
  "ark_kpr_mri.c\;8 0.001\;exclude-single"
  "ark_kpr_mri.c\;9 0.001 -100 100 0.5 1\;exclude-single"
  "ark_kpr_mri.c\;9 0.001\;exclude-single"
  "ark_KrylovDemo_prec.c\;<none>\;exclude-single"
  "ark_KrylovDemo_prec.c\;1\;exclude-single"
  "ark_KrylovDemo_prec.c\;2\;exclude-single"
  "ark_onewaycouple_mri.c\;<none>\;develop"
  "ark_reaction_diffusion_mri.c\;<none>\;develop"
  "ark_robertson_constraints.c\;<none>\;exclude-single"
  "ark_robertson_root.c\;<none>\;exclude-single"
  "ark_robertson.c\;<none>\;exclude-single"
  "ark_twowaycouple_mri.c\;<none>\;develop"
  )

if(SUNDIALS_BUILD_WITH_MONITORING)
  list(APPEND ARKODE_examples "ark_brusselator_fp.c\;1\;exclude-single")
endif()

# Examples using KLU linear solver
set(ARKODE_examples_KLU
  "ark_brusselator1D_klu.c\;<none>\;develop"
  )

# Examples using SuperLU_MT linear solver
set(ARKODE_examples_SUPERLUMT
  "ark_brusselator1D_FEM_slu.c\;<none>\;exclude-single"
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

  # parse the example tuple using sundials_parse_tuple macro
  sundials_parse_tuple( ARKODE_examples 
    FILENAME_VAR example
    PARAMS_VAR example_args
    TYPE_VAR example_type
  ) 
  if (NOT TARGET ${example})
    # example source files
    add_executable(${example} ${example})

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

endforeach(example_tuple ${ARKODE_examples})

# Add the build and install targets for each KLU example (if needed)
if(BUILD_SUNLINSOL_KLU)

  foreach(example_tuple ${ARKODE_examples_KLU})

    # parse the example tuple using sundials_parse_tuple macro
    sundials_parse_tuple( ARKODE_examples_KLU 
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
    target_link_libraries(${example}
      sundials_arkode
      sundials_nvecserial
      sundials_sunlinsolklu
      ${EXE_EXTRA_LINK_LIBS})

  endforeach(example_tuple ${ARKODE_examples_KLU})

endif()


# Add the build and install targets for each SuperLU_MT example (if needed)
if(BUILD_SUNLINSOL_SUPERLUMT)

  foreach(example_tuple ${ARKODE_examples_SUPERLUMT})

    # parse the example tuple using sundials_parse_tuple macro
    sundials_parse_tuple( ARKODE_examples_SUPERLUMT 
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
    target_link_libraries(${example}
      sundials_arkode
      sundials_nvecserial
      sundials_sunlinsolsuperlumt
      ${EXE_EXTRA_LINK_LIBS})

  endforeach(example_tuple ${ARKODE_examples_SUPERLUMT})

endif()

# create Makfile and CMakeLists.txt for examples
if(EXAMPLES_INSTALL)

  if(BUILD_SUNLINSOL_SUPERLUMT)
    set(_include_sunlinsolsuperlumt sunlinsolsuperlumt)
    list(APPEND ARKODE_examples "ark_brusselator1D_FEM_slu.c\; \"\" \;exclude-single" )
    set(SUPERLUMT_TARGET SUNDIALS::SUPERLUMT)
  endif()

  if(BUILD_SUNLINSOL_KLU)
    set(_include_sunlinsolklu sunlinsolklu)
    list(APPEND ARKODE_examples "ark_brusselator1D_klu.c\; \"\" \;develop" )
    set(KLU_TARGET SUNDIALS::KLU)
  endif()

  sundials_install_examples(arkode ARKODE_examples
    CMAKE_TEMPLATE
      cmakelists_C_ex.in
    MAKE_TEMPLATE
      makefile_serial_C_ex.in
    SUNDIALS_TARGETS
      nvecserial
      arkode
      ${_include_sunlinsolsuperlumt}
      ${_include_sunlinsolklu}
    OTHER_TARGETS
      ${SUPERLUMT_TARGET}
      ${KLU_TARGET}
    DESTINATION
      arkode/C_serial
    EXTRA_FILES
      README
      ${ARKODE_extras}
      ${ARKODE_headers}
    TEST_INSTALL
      C_serial
  )

endif()
