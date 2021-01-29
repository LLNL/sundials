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
# CMakeLists.txt file for CPODES parallel examples
# ---------------------------------------------------------------


# Add variable CPODES_examples with the names of the parallel CPODES examples

set(CPODES_examples
  cpsHeat2D_kry_bbd_p
  )

if(MPI_C_COMPILER)
  # use MPI wrapper as the compiler
  set(CMAKE_C_COMPILER ${MPI_C_COMPILER})
else()
  # add MPI_INCLUDE_PATH to include directories
  include_directories(${MPI_INCLUDE_PATH})
endif()

# Specify libraries to link against (through the target that was used to
# generate them) based on the value of the variable LINK_LIBRARY_TYPE

if(LINK_LIBRARY_TYPE MATCHES "static")
  set(CPODES_LIB sundials_cpodes${_STATIC_LIB_SUFFIX})
  set(NVECP_LIB sundials_nvecparallel${_STATIC_LIB_SUFFIX})
else(LINK_LIBRARY_TYPE MATCHES "static")
  set(CPODES_LIB sundials_cpodes)
  set(NVECP_LIB sundials_nvecparallel)
endif(LINK_LIBRARY_TYPE MATCHES "static")

# Set-up linker flags and link libraries

set(SUNDIALS_LIBS ${CPODES_LIB} ${NVECP_LIB} ${EXTRA_LINK_LIBS})
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
  target_link_libraries(${example} ${SUNDIALS_LIBS})
  if(NOT MPI_C_COMPILER)
    target_link_libraries(${example} ${MPI_LIBRARY} ${MPI_EXTRA_LIBRARIES})
  endif()
  if(EXAMPLES_INSTALL)
    install(FILES ${example}.c ${example}.out DESTINATION ${EXAMPLES_INSTALL_PATH}/cpodes/parallel)
  endif()
endforeach(example ${CPODES_examples})

if(EXAMPLES_INSTALL)

  # Install the README file
  install(FILES README DESTINATION ${EXAMPLES_INSTALL_PATH}/cpodes/parallel)

  # Prepare substitution variables for Makefile and/or CMakeLists templates
  set(SOLVER "CPODES")
  set(SOLVER_LIB "sundials_cpodes")
  list2string(CPODES_examples EXAMPLES)

  # Regardless of the platform we're on, we will generate and install
  # CMakeLists.txt file for building the examples. This file  can then
  # be used as a template for the user's own programs.

  # generate CMakelists.txt in the binary directory
  configure_file(
      ${PROJECT_SOURCE_DIR}/examples/templates/cmakelists_parallel_C_ex.in
      ${PROJECT_BINARY_DIR}/examples/cpodes/parallel/CMakeLists.txt
      @ONLY
      )

  # install CMakelists.txt
  install(
    FILES ${PROJECT_BINARY_DIR}/examples/cpodes/parallel/CMakeLists.txt
    DESTINATION ${EXAMPLES_INSTALL_PATH}/cpodes/parallel
    )

  # On UNIX-type platforms, we also  generate and install a makefile for
  # building the examples. This makefile can then be used as a template
  # for the user's own programs.

  if(UNIX)
    # generate Makefile and place it in the binary dir
    configure_file(
      ${PROJECT_SOURCE_DIR}/examples/templates/makefile_parallel_C_ex.in
      ${PROJECT_BINARY_DIR}/examples/cpodes/parallel/Makefile_ex
      @ONLY
      )
    # install the configured Makefile_ex as Makefile
    install(
      FILES ${PROJECT_BINARY_DIR}/examples/cpodes/parallel/Makefile_ex
      DESTINATION ${EXAMPLES_INSTALL_PATH}/cpodes/parallel
      RENAME Makefile
      )
  endif()

endif()
