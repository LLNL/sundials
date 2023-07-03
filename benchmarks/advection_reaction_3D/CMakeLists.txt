# ---------------------------------------------------------------
# Programmer(s): Daniel R. Reynolds @ SMU
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

if(ENABLE_RAJA)
  add_subdirectory(raja)
endif()

  # ----------------------------------------------------------------------------
  # MPI only
  # ----------------------------------------------------------------------------

  add_executable(advection_reaction_3D
    advection_reaction_3D.cpp
    arkode_driver.cpp
    cvode_driver.cpp
    ida_driver.cpp
    rhs3D.hpp
    ParallelGrid.hpp
    backends.hpp)

  # ensure the linker language is reset to CXX
  set_target_properties(advection_reaction_3D PROPERTIES LINKER_LANGUAGE CXX)

  target_include_directories(advection_reaction_3D
    PRIVATE
    ${PROJECT_SOURCE_DIR}/utilities
    ${MPI_CXX_INCLUDE_DIRS})

  target_link_libraries(advection_reaction_3D
    PRIVATE
    sundials_arkode
    sundials_cvode
    sundials_ida
    sundials_nvecmpiplusx
    sundials_nvecserial
    RAJA
    ${MPI_CXX_LIBRARIES}
    ${OTHER_LIBS})

  install(TARGETS advection_reaction_3D
    DESTINATION "${BENCHMARKS_INSTALL_PATH}/advection_reaction_3D")

  install(FILES README.md
    DESTINATION "${BENCHMARKS_INSTALL_PATH}/advection_reaction_3D")

  set(BENCHMARK_VAR
    "--method ARK-IMEX --nls tl-newton --tf 0.01 --dont-save\;arkimex_tlnewton"
    "--method ARK-DIRK --nls newton --tf 0.01 --dont-save\;arkdirk_newton"
    "--method CV-BDF --nls newton --tf 0.01 --dont-save\;cvbdf_newton"
    "--method IDA --nls newton --tf 0.01 --dont-save\;ida_newton"
    )

  # NUM_CORES based on Lassen system for 1 node
  foreach(benchmark_tuple ${BENCHMARK_VAR})
    list(GET benchmark_tuple 0 benchmark_args)
    list(GET benchmark_tuple 1 identifier)

    sundials_add_benchmark(advection_reaction_3D advection_reaction_3D advection_reaction_3D
      NUM_CORES 40
      BENCHMARK_ARGS ${benchmark_args}
      IDENTIFIER ${identifier}
      )
  endforeach()

  # ----------------------------------------------------------------------------
  # MPI + CUDA
  # ----------------------------------------------------------------------------

  if(BUILD_NVECTOR_CUDA)

    set_source_files_properties(advection_reaction_3D.cpp
      PROPERTIES LANGUAGE CUDA)
    set_source_files_properties(arkode_driver.cpp PROPERTIES LANGUAGE CUDA)
    set_source_files_properties(cvode_driver.cpp PROPERTIES LANGUAGE CUDA)
    set_source_files_properties(ida_driver.cpp PROPERTIES LANGUAGE CUDA)

    add_executable(advection_reaction_3D_mpicuda
      advection_reaction_3D.cpp
      arkode_driver.cpp
      cvode_driver.cpp
      ida_driver.cpp
      rhs3D.hpp
      ParallelGrid.hpp
      backends.hpp)

    # ensure the linker language is reset to CXX
    set_target_properties(advection_reaction_3D_mpicuda
      PROPERTIES LINKER_LANGUAGE CXX)

    target_include_directories(advection_reaction_3D_mpicuda
      PRIVATE
      ${PROJECT_SOURCE_DIR}/utilities
      ${MPI_CXX_INCLUDE_DIRS})

    target_link_libraries(advection_reaction_3D_mpicuda
      PRIVATE
      sundials_arkode
      sundials_cvode
      sundials_ida
      sundials_nvecmpiplusx
      sundials_nveccuda
      RAJA
      ${MPI_CXX_LIBRARIES}
      ${OTHER_LIBS})

    target_compile_definitions(advection_reaction_3D_mpicuda PRIVATE USE_CUDA_NVEC)

    install(TARGETS advection_reaction_3D_mpicuda
      DESTINATION "${BENCHMARKS_INSTALL_PATH}/advection_reaction_3D")

    set(BENCHMARK_VAR
        "--method ARK-IMEX --nls tl-newton --tf 0.01 --dont-save\;arkimex_tlnewton_cuda"
        "--method ARK-DIRK --nls newton --tf 0.01 --dont-save\;arkdirk_newton_cuda"
        "--method CV-BDF --nls newton --tf 0.01 --dont-save\;cvbdf_newton_cuda"
        "--method IDA --nls newton --tf 0.01 --dont-save\;ida_newton_cuda"
      )
    # NUM_CORES based on Lassen system for 1 node
    foreach(benchmark_tuple ${BENCHMARK_VAR})
      list(GET benchmark_tuple 0 benchmark_args)
      list(GET benchmark_tuple 1 identifier)

      sundials_add_benchmark(advection_reaction_3D_mpicuda advection_reaction_3D_mpicuda advection_reaction_3D
        ENABLE_GPU
        NUM_CORES 4
        BENCHMARK_ARGS ${benchmark_args}
        IDENTIFIER ${identifier}
        )
    endforeach()

  endif()

  # ----------------------------------------------------------------------------
  # MPI + HIP
  # ----------------------------------------------------------------------------

  if(BUILD_NVECTOR_HIP)
    add_dependencies(benchmark advection_reaction_3D_mpihip)

    add_executable(advection_reaction_3D_mpihip
      advection_reaction_3D.cpp
      arkode_driver.cpp
      cvode_driver.cpp
      ida_driver.cpp
      rhs3D.hpp
      ParallelGrid.hpp
      backends.hpp)

    target_include_directories(advection_reaction_3D_mpihip
      PRIVATE
      ${PROJECT_SOURCE_DIR}/utilities
      ${MPI_CXX_INCLUDE_DIRS})

    target_link_libraries(advection_reaction_3D_mpihip
      PRIVATE
      sundials_arkode
      sundials_cvode
      sundials_ida
      sundials_nvecmpiplusx
      sundials_nvechip
      RAJA
      hip::device
      ${MPI_CXX_LIBRARIES}
      ${OTHER_LIBS})

    target_compile_definitions(advection_reaction_3D_mpihip PRIVATE USE_HIP_NVEC)

    install(TARGETS advection_reaction_3D_mpihip
      DESTINATION "${BENCHMARKS_INSTALL_PATH}/advection_reaction_3D")

    set(BENCHMARK_VAR 
        "--method ARK-IMEX --nls tl-newton --tf 1.0 --dont-save\;arkimex_tlnewton_hip"
        "--method ARK-DIRK --nls newton --tf 1.0 --dont-save\;arkdirk_newton_hip"
        "--method CV-BDF --nls newton --tf 1.0 --dont-save\;cvbdf_newton_hip"
        "--method IDA --nls newton --tf 1.0 --dont-save\;ida_newton_hip"
      )


    foreach(benchmark_tuple ${BENCHMARK_VAR})
      list(GET benchmark_tuple 0 benchmark_args)
      list(GET benchmark_tuple 1 identifier)

      # NUM_CORES based on Corona system for 1 node
      sundials_add_benchmark(advection_reaction_3D_mpihip advection_reaction_3D_mpihip advection_reaction_3D
        ENABLE_GPU
        NUM_CORES 8
        BENCHMARK_ARGS ${benchmark_args}
        IDENTIFIER ${identifier}
        )
    endforeach()
  endif()

endif()
