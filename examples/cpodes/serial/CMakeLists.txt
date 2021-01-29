# ---------------------------------------------------------------
# Programmer(s): Radu Serban @ LLNL
# ---------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2021, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# ---------------------------------------------------------------
# CMakeLists.txt file for CPODES serial examples
# ---------------------------------------------------------------


# Add variable CPODES_examples with the names of the serial CPODES examples

set(CPODES_examples
  cpsAdvDiff_bnd
  cpsAdvDiff_non
  cpsNewtCrd_dns
  cpsPend_dns
  cpsRoberts_dns
  cpsSlCrank_dns
  cpsVanDPol_non
  )

# Add variable CPODES_examples_BL with the names of the serial CPODES examples
# that use Lapack

set(CPODES_examples_BL
  cpsAdvDiff_bndL
  cpsAkzoNob_dnsL
  cpsPend_dnsL
  cpsRoberts_dnsL
  )

# Specify libraries to link against (through the target that was used to
# generate them) based on the value of the variable LINK_LIBRARY_TYPE

if(LINK_LIBRARY_TYPE MATCHES "static")
  set(CPODES_LIB sundials_cpodes${_STATIC_LIB_SUFFIX})
  set(NVECS_LIB sundials_nvecserial${_STATIC_LIB_SUFFIX})
else(LINK_LIBRARY_TYPE MATCHES "static")
  set(CPODES_LIB sundials_cpodes)
  set(NVECS_LIB sundials_nvecserial)
endif(LINK_LIBRARY_TYPE MATCHES "static")

# Set-up linker flags and link libraries

set(SUNDIALS_LIBS ${CPODES_LIB} ${NVECS_LIB} ${EXTRA_LINK_LIBS})
if(LAPACK_FOUND)
  list(APPEND SUNDIALS_LIBS ${LAPACK_LIBRARIES})
endif()

if(BLAS_FOUND)
  list(APPEND SUNDIALS_LIBS ${BLAS_LIBRARIES})
endif()

# Add the build and install targets for each CPODES example

foreach(example ${CPODES_examples})
  add_executable(${example} ${example}.c)
  set_target_properties(${example} PROPERTIES FOLDER "Examples")
  #ADD_TEST(${example} ${TESTRUNNER} ${example})
  target_link_libraries(${example} ${SUNDIALS_LIBS})
  file(GLOB example_out ${example}.out*)
  if(EXAMPLES_INSTALL)
    install(FILES ${example}.c ${example_out} DESTINATION ${EXAMPLES_INSTALL_PATH}/cpodes/serial)
  endif()
endforeach(example ${CPODES_examples})

# If Lapack support is enabled, add the build and install targets for
# the examples using Lapack

if(LAPACK_FOUND)
  foreach(example ${CPODES_examples_BL})
    add_executable(${example} ${example}.c)
    set_target_properties(${example} PROPERTIES FOLDER "Examples")
    sundials_add_test(${example} ${example})
    target_link_libraries(${example} ${SUNDIALS_LIBS})
    file(GLOB example_out ${example}.out*)
    if(EXAMPLES_INSTALL)
      install(FILES ${example}.c ${example_out} DESTINATION ${EXAMPLES_INSTALL_PATH}/cpodes/serial)
    endif()
  endforeach(example ${CPODES_examples_BL})
endif()

if(EXAMPLES_INSTALL)

  # Install the README file
  install(FILES README DESTINATION ${EXAMPLES_INSTALL_PATH}/cpodes/serial)

  # Prepare substitution variables for Makefile and/or CMakeLists templates
  set(SOLVER "CPODES")
  set(SOLVER_LIB "sundials_cpodes")
  list2string(CPODES_examples EXAMPLES)
  if(LAPACK_FOUND)
    list2string(CPODES_examples_BL EXAMPLES_BL)
  else(LAPACK_FOUND)
    set(EXAMPLES_BL "")
  endif()

  # Regardless of the platform we're on, we will generate and install
  # CMakeLists.txt file for building the examples. This file  can then
  # be used as a template for the user's own programs.

  # generate CMakelists.txt in the binary directory
  configure_file(
      ${PROJECT_SOURCE_DIR}/examples/templates/cmakelists_serial_C_ex.in
      ${PROJECT_BINARY_DIR}/examples/cpodes/serial/CMakeLists.txt
      @ONLY
      )

  # install CMakelists.txt
  install(
    FILES ${PROJECT_BINARY_DIR}/examples/cpodes/serial/CMakeLists.txt
    DESTINATION ${EXAMPLES_INSTALL_PATH}/cpodes/serial
    )

  # On UNIX-type platforms, we also  generate and install a makefile for
  # building the examples. This makefile can then be used as a template
  # for the user's own programs.

  if(UNIX)
    # generate Makefile and place it in the binary dir
    configure_file(
      ${PROJECT_SOURCE_DIR}/examples/templates/makefile_serial_C_ex.in
      ${PROJECT_BINARY_DIR}/examples/cpodes/serial/Makefile_ex
      @ONLY
      )
    # install the configured Makefile_ex as Makefile
    install(
      FILES ${PROJECT_BINARY_DIR}/examples/cpodes/serial/Makefile_ex
      DESTINATION ${EXAMPLES_INSTALL_PATH}/cpodes/serial
      RENAME Makefile
      )
  endif()

endif()
