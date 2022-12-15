# ---------------------------------------------------------------
# Programmer(s):  Daniel R. Reynolds @ SMU
# ---------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2022, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# ---------------------------------------------------------------

if(BUILD_ARKODE AND BUILD_CVODE AND BUILD_IDA)

  # ----------------------------------------------------------------------------
  # MPI only
  # ----------------------------------------------------------------------------

  add_executable(advection_reaction_3D_kokkos
    advection_reaction_3D.cpp
    arkode_driver.cpp
    cvode_driver.cpp
    ida_driver.cpp
    rhs3D.hpp
    ParallelGrid.hpp
    check_retval.h)

  # ensure the linker language is reset to CXX
  set_target_properties(advection_reaction_3D_kokkos PROPERTIES LINKER_LANGUAGE CXX)

  target_include_directories(advection_reaction_3D_kokkos
    PRIVATE
    ${PROJECT_SOURCE_DIR}/utilities
    ${MPI_CXX_INCLUDE_DIRS})

  target_link_libraries(advection_reaction_3D_kokkos
    PRIVATE
    sundials_arkode
    sundials_cvode
    sundials_ida
    sundials_nvecmpiplusx
    sundials_nveckokkos
    RAJA
    ${MPI_CXX_LIBRARIES}
    ${OTHER_LIBS})

  install(TARGETS advection_reaction_3D_kokkos
    DESTINATION "${CMAKE_INSTALL_BINDIR}/benchmarks/advection_reaction_3D/kokkos")

  install(FILES README.md ../scripts/compare_error.py ../scripts/compute_error.py ../scripts/pickle_solution_output.py
    DESTINATION "${CMAKE_INSTALL_BINDIR}/benchmarks/advection_reaction_3D/kokkos")


endif()
