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

message(STATUS "Added Kokkos NVECTOR benchmark")


set(BENCHMARKS_DIR ${PROJECT_SOURCE_DIR}/benchmarks)

foreach(backend ${KOKKOS_EXAMPLES_BACKENDS})

  sundials_add_nvector_benchmark(test_nvector_performance_kokkos.${backend}
    SOURCES test_nvector_performance_kokkos.cpp
    SUNDIALS_TARGETS sundials_generic sundials_nveckokkos
    )

  target_compile_definitions(test_nvector_performance_kokkos.${backend} PRIVATE USE_${backend})

  install(TARGETS test_nvector_performance_kokkos.${backend}
    DESTINATION "${CMAKE_INSTALL_BINDIR}/benchmarks/nvector/kokkos")

endforeach()
