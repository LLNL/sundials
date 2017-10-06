# ---------------------------------------------------------------
# Programmer:  Daniel R. Reynolds @ SMU
# ---------------------------------------------------------------
# LLNS/SMU Copyright Start
# Copyright (c) 2015, Southern Methodist University and 
# Lawrence Livermore National Security
#
# This work was performed under the auspices of the U.S. Department 
# of Energy by Southern Methodist University and Lawrence Livermore 
# National Laboratory under Contract DE-AC52-07NA27344.
# Produced at Southern Methodist University and the Lawrence 
# Livermore National Laboratory.
#
# All rights reserved.
# For details, see the LICENSE file.
# LLNS/SMU Copyright End
# Copyright (c) 2013, Southern Methodist University.
# All rights reserved.
# For details, see the LICENSE file.
# -----------------------------------------------------------------
# CMakeLists.txt file for ARKODE serial examples

# Add variable ARKODE_examples with the names of the serial ARKODE examples

SET(ARKODE_examples
  ark_analytic
  ark_analytic_nonlin
  ark_brusselator
  ark_brusselator_fp
  ark_brusselator1D
  ark_heat1D
  ark_heat1D_adapt
  ark_KrylovDemo_prec
  ark_robertson
  ark_robertson_root
  )

# Add variable ARKODE_examples_BL with the names of the serial ARKODE examples
# that use Lapack

SET(ARKODE_examples_BL
  )

# Add variable ARKODE_examples_KLU with the names of the serial ARKODE examples
# that use KLU

SET(ARKODE_examples_KLU
  ark_brusselator1D_klu
  )

# Add variable ARKODE_examples_SUPERLUMT with the names of the serial ARKODE examples
# that use SUPERLUMT

SET(ARKODE_examples_SUPERLUMT
  ark_brusselator1D_FEM_slu
  )

# Add variable ARKODE_extras with the names of auxiliary files to install

SET(ARKODE_extras
  plot_brusselator1D.py
  plot_brusselator1D_FEM.py
  plot_heat1D.py
  plot_heat1D_adapt.py
  plot_sol.py
  plot_sol_log.py
  )

# Specify libraries to link against (through the target that was used to 
# generate them) based on the value of the variable LINK_LIBRARY_TYPE

IF(LINK_LIBRARY_TYPE MATCHES "static")
  SET(ARKODE_LIB sundials_arkode_static)
  SET(NVECS_LIB sundials_nvecserial_static)
ELSE(LINK_LIBRARY_TYPE MATCHES "static")
  SET(ARKODE_LIB sundials_arkode_shared)
  SET(NVECS_LIB sundials_nvecserial_shared)
ENDIF(LINK_LIBRARY_TYPE MATCHES "static")

# Set-up linker flags and link libraries

SET(SUNDIALS_LIBS ${ARKODE_LIB} ${NVECS_LIB} ${EXTRA_LINK_LIBS})
IF(LAPACK_FOUND)
  LIST(APPEND SUNDIALS_LIBS ${LAPACK_LIBRARIES})
ENDIF(LAPACK_FOUND)

IF(KLU_FOUND)
  LIST(APPEND SUNDIALS_LIBS ${KLU_LIBRARIES})
ENDIF(KLU_FOUND)

IF(SUPERLUMT_FOUND)
  LIST(APPEND SUNDIALS_LIBS ${SUPERLUMT_LIBRARIES})
ENDIF(SUPERLUMT_FOUND)

# Add the build and install targets for each ARKODE example

FOREACH(example ${ARKODE_examples})
  ADD_EXECUTABLE(${example} ${example}.c)
  SET_TARGET_PROPERTIES(${example} PROPERTIES FOLDER "Examples")
  SUNDIALS_ADD_TEST(${example} ${example})
  TARGET_LINK_LIBRARIES(${example} ${SUNDIALS_LIBS})
  IF(EXAMPLES_INSTALL)
    INSTALL(FILES ${example}.c ${example}.out DESTINATION ${EXAMPLES_INSTALL_PATH}/arkode/C_serial)
  ENDIF(EXAMPLES_INSTALL)
ENDFOREACH(example ${ARKODE_examples})

# Add the build and install targets for each Lapack ARKODE example (if needed)

IF(LAPACK_FOUND)
  FOREACH(example ${ARKODE_examples_BL})
    ADD_EXECUTABLE(${example} ${example}.c)
    SET_TARGET_PROPERTIES(${example} PROPERTIES FOLDER "Examples")
    SUNDIALS_ADD_TEST(${example} ${example})
    TARGET_LINK_LIBRARIES(${example} ${SUNDIALS_LIBS})
    IF(EXAMPLES_INSTALL)
      INSTALL(FILES ${example}.c ${example}.out DESTINATION ${EXAMPLES_INSTALL_PATH}/arkode/C_serial)
    ENDIF(EXAMPLES_INSTALL)
  ENDFOREACH(example ${ARKODE_examples_BL})
ENDIF(LAPACK_FOUND)

# If KLU support is enabled, add the build and install targets for 
# the examples using KLU
IF(KLU_FOUND)
  FOREACH(example ${ARKODE_examples_KLU})
    ADD_EXECUTABLE(${example} ${example}.c)
    SET_TARGET_PROPERTIES(${example} PROPERTIES FOLDER "Examples")
    SUNDIALS_ADD_TEST(${example} ${example})
    TARGET_LINK_LIBRARIES(${example} ${SUNDIALS_LIBS})
    IF(EXAMPLES_INSTALL)
      INSTALL(FILES ${example}.c ${example}.out DESTINATION ${EXAMPLES_INSTALL_PATH}/arkode/C_serial)
    ENDIF(EXAMPLES_INSTALL)
  ENDFOREACH(example ${ARKODE_examples_KLU})
ENDIF(KLU_FOUND)

# If SuperLUMT support is enabled, add the build and install targets for 
# the examples using SuperLUMT
IF(SUPERLUMT_FOUND)
  FOREACH(example ${ARKODE_examples_SUPERLUMT})
    ADD_EXECUTABLE(${example} ${example}.c)
    SET_TARGET_PROPERTIES(${example} PROPERTIES FOLDER "Examples")
    SUNDIALS_ADD_TEST(${example} ${example})
    TARGET_LINK_LIBRARIES(${example} ${SUNDIALS_LIBS})
    IF(EXAMPLES_INSTALL)
      INSTALL(FILES ${example}.c ${example}.out DESTINATION ${EXAMPLES_INSTALL_PATH}/arkode/C_serial)
    ENDIF(EXAMPLES_INSTALL)
  ENDFOREACH(example ${ARKODE_examples_SUPERLUMT})
ENDIF(SUPERLUMT_FOUND)

IF(EXAMPLES_INSTALL)

  # Install the README file
  INSTALL(FILES README DESTINATION ${EXAMPLES_INSTALL_PATH}/arkode/C_serial)

  # Install the extra files
  FOREACH(extrafile ${ARKODE_extras})
    INSTALL(FILES ${extrafile} DESTINATION ${EXAMPLES_INSTALL_PATH}/arkode/C_serial)
  ENDFOREACH(extrafile ${ARKODE_extras})

  # Prepare substitution variables for Makefile and/or CMakeLists templates
  SET(SOLVER "ARKODE")
  SET(SOLVER_LIB "sundials_arkode")
  LIST2STRING(ARKODE_examples EXAMPLES)

  IF(LAPACK_FOUND)
    LIST2STRING(ARKODE_examples_BL EXAMPLES_BL)
  ELSE(LAPACK_FOUND)
    SET(EXAMPLES_BL "")
  ENDIF(LAPACK_FOUND)

  IF(KLU_FOUND)
    LIST2STRING(ARKODE_examples_KLU EXAMPLES)
  ENDIF(KLU_FOUND)

  IF(SUPERLUMT_FOUND)
    LIST2STRING(ARKODE_examples_SUPERLUMT EXAMPLES)
  ENDIF(SUPERLUMT_FOUND)

  STRING (REPLACE ";" " " TMP_STR ${EXAMPLES})
  SET(EXAMPLES ${TMP_STR})


  # Regardless of the platform we're on, we will generate and install 
  # CMakeLists.txt file for building the examples. This file  can then 
  # be used as a template for the user's own programs.

  # generate CMakelists.txt in the binary directory
  CONFIGURE_FILE(
      ${PROJECT_SOURCE_DIR}/examples/templates/cmakelists_serial_C_ex.in
      ${PROJECT_BINARY_DIR}/examples/arkode/C_serial/CMakeLists.txt
      @ONLY
      )

  # install CMakelists.txt
  INSTALL(
    FILES ${PROJECT_BINARY_DIR}/examples/arkode/C_serial/CMakeLists.txt
    DESTINATION ${EXAMPLES_INSTALL_PATH}/arkode/C_serial 
    )

  # On UNIX-type platforms, we also  generate and install a makefile for 
  # building the examples. This makefile can then be used as a template 
  # for the user's own programs.

  IF(UNIX)
    # generate Makefile and place it in the binary dir
    CONFIGURE_FILE(
      ${PROJECT_SOURCE_DIR}/examples/templates/makefile_serial_C_ex.in
      ${PROJECT_BINARY_DIR}/examples/arkode/C_serial/Makefile_ex
      @ONLY
      )
    # install the configured Makefile_ex as Makefile
    INSTALL(
      FILES ${PROJECT_BINARY_DIR}/examples/arkode/C_serial/Makefile_ex 
      DESTINATION ${EXAMPLES_INSTALL_PATH}/arkode/C_serial 
      RENAME Makefile
      )
  ENDIF(UNIX)

ENDIF(EXAMPLES_INSTALL)
