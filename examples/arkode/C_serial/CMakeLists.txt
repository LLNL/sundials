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
  "ark_analytic.c\; \"\" \; \"\" "
  "ark_analytic_mels.c\; \"\" \;develop"
  "ark_analytic_nonlin.c\; \"\" \;develop"
  "ark_brusselator.c\; \"\" \;develop"
  "ark_brusselator_fp.c\; \"\" \;exclude-single"
  "ark_brusselator1D.c\; \"\" \;exclude-single"
  "ark_heat1D.c\; \"\" \;develop"
  "ark_heat1D_adapt.c\; \"\" \;develop"
  "ark_KrylovDemo_prec.c\; \"\" \;exclude-single"
  "ark_KrylovDemo_prec.c\; \"1\" \;exclude-single"
  "ark_KrylovDemo_prec.c\; \"2\" \;exclude-single"
  "ark_robertson.c\; \"\" \;exclude-single"
  "ark_robertson_constraints.c\; \"\" \;exclude-single"
  "ark_robertson_root.c\; \"\" \;exclude-single"
  "ark_brusselator_mri.c\; \"\" \;develop"
  "ark_onewaycouple_mri.c\; \"\" \;develop"
  "ark_twowaycouple_mri.c\; \"\" \;develop"
  "ark_reaction_diffusion_mri.c\; \"\" \;develop"
  "ark_brusselator_1D_mri.c\; \"\" \;develop"
  "ark_kpr_mri.c\; \"\" \;develop"
  "ark_kpr_mri.c\; \"0 0.002\" \;develop"
  "ark_kpr_mri.c\; \"1 0.002\" \;develop"
  "ark_kpr_mri.c\; \"2 0.005\" \;develop"
  "ark_kpr_mri.c\; \"3 0.01\" \;develop"
  "ark_kpr_mri.c\; \"4 0.002\" \;develop"
  "ark_kpr_mri.c\; \"5 0.002\" \;develop"
  "ark_kpr_mri.c\; \"6 0.005\" \;develop"
  "ark_kpr_mri.c\; \"7 0.001\" \;exclude-single"
  "ark_kpr_mri.c\; \"7 0.001 -100 100 0.5 1\" \;exclude-single"
  "ark_kpr_mri.c\; \"8 0.001\" \;exclude-single"
  "ark_kpr_mri.c\; \"8 0.001 -100 100 0.5 1\" \;exclude-single"
  "ark_kpr_mri.c\; \"9 0.001\" \;exclude-single"
  "ark_kpr_mri.c\; \"9 0.001 -100 100 0.5 1\" \;exclude-single"
  "ark_brusselator1D_imexmri.c\; \"0 0.001\" \;exclude-single"
  "ark_brusselator1D_imexmri.c\; \"2 0.001\" \;exclude-single"
  "ark_brusselator1D_imexmri.c\; \"3 0.001\" \;exclude-single"
  "ark_brusselator1D_imexmri.c\; \"4 0.001\" \;exclude-single"
  "ark_brusselator1D_imexmri.c\; \"5 0.001\" \;exclude-single"
  "ark_brusselator1D_imexmri.c\; \"6 0.001\" \;exclude-single"
  "ark_brusselator1D_imexmri.c\; \"7 0.001\" \;exclude-single"
  )

if(SUNDIALS_BUILD_WITH_MONITORING)
  list(APPEND ARKODE_examples "ark_brusselator_fp\; \"1\" \;exclude-single")
endif()

# Examples using KLU linear solver
set(ARKODE_examples_KLU
  "ark_brusselator1D_klu\; \"\" \;develop"
  )

# Examples using SuperLU_MT linear solver
set(ARKODE_examples_SUPERLUMT
  "ark_brusselator1D_FEM_slu\; \"\" \;exclude-single"
  )

# Auxiliary files to install
set(ARKODE_extras
  plot_brusselator1D.py
  plot_brusselator1D_FEM.py
  plot_heat1D.py
  plot_heat1D_adapt.py
  plot_sol.py
  plot_sol_log.py
  ark_analytic_nonlin_stats.csv
  ark_reaction_diffusion_mri_fast_stats.csv
  ark_reaction_diffusion_mri_slow_stats.csv
  ark_robertson_stats.csv
  )

# Specify libraries to link against
set(ARKODE_LIB sundials_arkode)
set(NVECS_LIB sundials_nvecserial)

# Set-up linker flags and link libraries
set(SUNDIALS_LIBS ${ARKODE_LIB} ${NVECS_LIB} ${EXE_EXTRA_LINK_LIBS})

# Add the build and install targets for each ARKODE example
foreach(example_tuple ${ARKODE_examples})

  # parse the example tuple using sundials_parse_tuple macro
  sundials_parse_tuple( CVODES_examples 
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

endforeach(example_tuple ${ARKODE_examples})

# Add the build and install targets for each KLU example (if needed)
if(BUILD_SUNLINSOL_KLU)

  # Sundials KLU linear solver module
  set(SUNLINSOLKLU_LIBS sundials_sunlinsolklu)

  # KLU libraries
  list(APPEND SUNLINSOLKLU_LIBS)

  foreach(example_tuple ${ARKODE_examples_KLU})

    # parse the example tuple using sundials_parse_tuple macro
    sundials_parse_tuple( CVODES_examples 
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

  endforeach(example_tuple ${ARKODE_examples_KLU})

endif()


# Add the build and install targets for each SuperLU_MT example (if needed)
if(BUILD_SUNLINSOL_SUPERLUMT)

  # Sundials SuperLU_MT linear solver module
  set(SUNLINSOLSLUMT_LIBS sundials_sunlinsolsuperlumt)

  # SuperLU_MT libraries
  list(APPEND SUNLINSOLSLUMT_LIBS ${SUPERLUMT_LIBRARIES})

  foreach(example_tuple ${ARKODE_examples_SUPERLUMT})

    # parse the example tuple using sundials_parse_tuple macro
    sundials_parse_tuple( CVODES_examples 
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

  endforeach(example_tuple ${ARKODE_examples_SUPERLUMT})

endif()

# create Makfile and CMakeLists.txt for examples
if(EXAMPLES_INSTALL)

  if(BUILD_SUNLINSOL_SUPERLUMT)
    set(_include_sunlinsolsuperlumt sunlinsolsuperlumt)
  endif()

  if(BUILD_SUNLINSOL_KLU)
    set(_include_sunlinsolklu sunlinsolklu)
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
