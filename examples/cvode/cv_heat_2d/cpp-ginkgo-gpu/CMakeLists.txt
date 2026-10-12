# -----------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2025-2026, Lawrence Livermore National Security,
# University of Maryland Baltimore County, and the SUNDIALS contributors.
# Copyright (c) 2013-2025, Lawrence Livermore National Security
# and Southern Methodist University.
# Copyright (c) 2002-2013, Lawrence Livermore National Security.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# -----------------------------------------------------------------
# CVODE Ginkgo heat equation example with GPU backends
set(gpu_examples "../cpp-ginkgo/cv_heat_2d_ginkgo.cpp\;\;develop")
sundials_add_examples_ginkgo(
  gpu_examples
  TARGETS sundials_cvode
  BACKENDS CUDA HIP
  INSTALL_FILES
    ../cpp-ginkgo/cv_heat_2d_ginkgo.hpp cv_heat_2d_ginkgo_gpu.hpp
    "${PROJECT_SOURCE_DIR}/examples/utilities/example_utilities.hpp"
  INCLUDE_DIRECTORIES "${CMAKE_CURRENT_SOURCE_DIR}/../cpp-ginkgo")
