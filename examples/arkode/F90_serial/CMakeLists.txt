# ---------------------------------------------------------------
# Programmer:  Daniel R. Reynolds
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
# CMakeLists.txt file for the FARKODE serial examples


# Add variable FARKODE_examples with the names of the serial FARKODE examples

SET(FARKODE_examples
  ark_bruss
  )

# Add variable FARKODE_examples_BL with the names of the serial FARKODE examples
# that use Lapack

SET(FARKODE_examples_BL
  )

# Add variable FARKODE_examples_KLU with the names of the serial FARKODE examples
# that use KLU

SET(FARKODE_examples_KLU
  ark_bruss1D_FEM_klu
  )

# Add variable ARKODE_extras with the names of auxiliary files to install

SET(ARKODE_extras
  plot_sol.py
  plot_brusselator1D_FEM.py
  )

# Specify libraries to link against (through the target that was used to 
# generate them) based on the value of the variable LINK_LIBRARY_TYPE

IF(LINK_LIBRARY_TYPE MATCHES "static")
  SET(ARKODE_LIB sundials_arkode_static)
  SET(NVECS_LIB sundials_nvecserial_static)
  SET(FNVECS_LIB sundials_fnvecserial_static)
ELSE(LINK_LIBRARY_TYPE MATCHES "static")
  SET(ARKODE_LIB sundials_arkode_shared)
  SET(NVECS_LIB sundials_nvecserial_shared)
  SET(FNVECS_LIB sundials_fnvecserial_shared)
ENDIF(LINK_LIBRARY_TYPE MATCHES "static")

# Only static FCMIX libraries are available

SET(FARKODE_LIB sundials_farkode_static)

# Set-up linker flags and link libraries

SET(SUNDIALS_LIBS ${FARKODE_LIB} ${ARKODE_LIB} ${FNVECS_LIB} ${NVECS_LIB} ${EXTRA_LINK_LIBS})
IF(LAPACK_FOUND)
  LIST(APPEND SUNDIALS_LIBS ${LAPACK_LIBRARIES})
ENDIF(LAPACK_FOUND)

IF(KLU_FOUND)
  LIST(APPEND SUNDIALS_LIBS ${KLU_LIBRARIES})
ENDIF(KLU_FOUND)

IF(SUPERLUMT_FOUND)
  LIST(APPEND SUNDIALS_LIBS ${SUPERLUMT_LIBRARIES})
ENDIF(SUPERLUMT_FOUND)

# Add the build and install targets for each FARKODE example

FOREACH(example ${FARKODE_examples})
  ADD_EXECUTABLE(${example} ${example}.f90)
  SET_TARGET_PROPERTIES(${example} PROPERTIES FOLDER "Examples")
  SUNDIALS_ADD_TEST(${example} ${example})
  TARGET_LINK_LIBRARIES(${example} ${SUNDIALS_LIBS})
  IF(EXAMPLES_INSTALL)
    INSTALL(FILES ${example}.f90 ${example}.out DESTINATION ${EXAMPLES_INSTALL_PATH}/arkode/F90_serial)
  ENDIF(EXAMPLES_INSTALL)
ENDFOREACH(example ${FARKODE_examples})

# Add the build and install targets for each Lapack FARKODE example (if needed)

IF(LAPACK_FOUND)
  FOREACH(example ${FARKODE_examples_BL})
    ADD_EXECUTABLE(${example} ${example}.f90)
    SET_TARGET_PROPERTIES(${example} PROPERTIES FOLDER "Examples")
    SUNDIALS_ADD_TEST(${example} ${example})
    TARGET_LINK_LIBRARIES(${example} ${SUNDIALS_LIBS})
    IF(EXAMPLES_INSTALL)
      INSTALL(FILES ${example}.f90 ${example}.out DESTINATION ${EXAMPLES_INSTALL_PATH}/arkode/F90_serial)
    ENDIF(EXAMPLES_INSTALL)
  ENDFOREACH(example ${FARKODE_examples_BL})
ENDIF(LAPACK_FOUND)

# If KLU support is enabled, add the build and install targets for 
# the examples using KLU
IF(KLU_FOUND)
  FOREACH(example ${FARKODE_examples_KLU})
    ADD_EXECUTABLE(${example} ${example}.f90)
    SET_TARGET_PROPERTIES(${example} PROPERTIES FOLDER "Examples")
    SUNDIALS_ADD_TEST(${example} ${example})
    TARGET_LINK_LIBRARIES(${example} ${SUNDIALS_LIBS})
    IF(EXAMPLES_INSTALL)
      INSTALL(FILES ${example}.f90 ${example}.out DESTINATION ${EXAMPLES_INSTALL_PATH}/arkode/F90_serial)
    ENDIF(EXAMPLES_INSTALL)
  ENDFOREACH(example ${FARKODE_examples_KLU})
ENDIF(KLU_FOUND)

IF(EXAMPLES_INSTALL)

  # Install the README file
  INSTALL(FILES README DESTINATION ${EXAMPLES_INSTALL_PATH}/arkode/F90_serial)

  # Install the extra files
  FOREACH(extrafile ${ARKODE_extras})
    INSTALL(FILES ${extrafile} DESTINATION ${EXAMPLES_INSTALL_PATH}/arkode/F90_serial)
  ENDFOREACH(extrafile ${ARKODE_extras})

  # Prepare substitution variables for Makefile and/or CMakeLists templates
  SET(SOLVER "ARKODE")
  SET(SOLVER_LIB "sundials_arkode")
  SET(SOLVER_FLIB "sundials_farkode")
  LIST2STRING(FARKODE_examples EXAMPLES)

  IF(LAPACK_FOUND)
    LIST2STRING(FARKODE_examples_BL EXAMPLES_BL)
  ELSE(LAPACK_FOUND)
    SET(EXAMPLES_BL "")
  ENDIF(LAPACK_FOUND)

  IF(KLU_FOUND)
    LIST2STRING(FARKODE_examples_KLU EXAMPLES)
  ENDIF(KLU_FOUND)


  # Regardless of the platform we're on, we will generate and install 
  # CMakeLists.txt file for building the examples. This file  can then 
  # be used as a template for the user's own programs.

  # generate CMakelists.txt in the binary directory
  CONFIGURE_FILE(
      ${PROJECT_SOURCE_DIR}/examples/templates/cmakelists_serial_F90_ex.in
      ${PROJECT_BINARY_DIR}/examples/arkode/F90_serial/CMakeLists.txt
      @ONLY
      )

  # install CMakelists.txt
  INSTALL(
    FILES ${PROJECT_BINARY_DIR}/examples/arkode/F90_serial/CMakeLists.txt
    DESTINATION ${EXAMPLES_INSTALL_PATH}/arkode/F90_serial 
    )

  # On UNIX-type platforms, we also  generate and install a makefile for 
  # building the examples. This makefile can then be used as a template 
  # for the user's own programs.

  IF(UNIX)
    # generate Makefile and place it in the binary dir
    CONFIGURE_FILE(
      ${PROJECT_SOURCE_DIR}/examples/templates/makefile_serial_F90_ex.in
      ${PROJECT_BINARY_DIR}/examples/arkode/F90_serial/Makefile_ex
      @ONLY
      )
    # install the configured Makefile_ex as Makefile
    INSTALL(
      FILES ${PROJECT_BINARY_DIR}/examples/arkode/F90_serial/Makefile_ex 
      DESTINATION ${EXAMPLES_INSTALL_PATH}/arkode/F90_serial 
      RENAME Makefile
      )
  ENDIF(UNIX)

ENDIF(EXAMPLES_INSTALL)
