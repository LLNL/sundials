#!/bin/bash
# -------------------------------------------------------------------------------
# Programmer(s): David J. Gardner @ LLNL
# -------------------------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2020, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# -------------------------------------------------------------------------------
# Script that sets up the environment for the suntest.sh when running on Ray.
#
# To use this environment script with suntest.sh, set the environment variable
# SUNDIALS_ENV to this file (e.g., export SUNDIALS_ENV=./env.ray.sh).
# -------------------------------------------------------------------------------

module load xl
module load cmake
module load cuda

# compiler settings
export CC="xlc-gpu"
export CXX="xlC-gpu"
export FC="xlf-gpu"

export CFLAGS="-g -Wall -Wpedantic -Werror"
export CXXFLAGS="-g -Wall -Wpedantic -Werror"
export FFLAGS="-g"

# enable MPI and set compiler wrappers
export MPI_STATUS=ON
export MPICC=$(which mpicc)
export MPICXX=$(which mpicxx)
export MPIFC=$(which mpifort)
export MPIEXEC=$(which jsrun)

# enable CUDA and OpenMP device offloading
export CUDA_STATUS=ON
export OPENMP_STATUS=ON
export OPENMPDEV_STATUS=ON
