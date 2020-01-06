#!/bin/bash
# -------------------------------------------------------------------------------
# Programmer(s): Cody J. Balos @ LLNL
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
# Script that sets up the environment for the suntest_driver.sh when running on
# Quartz.
# -------------------------------------------------------------------------------

# compiler settings
export CC=/usr/tce/packages/gcc/gcc-7.3.0/bin/gcc
export FC=/usr/tce/packages/gcc/gcc-7.3.0/bin/gfortran
export CXX=/usr/tce/packages/gcc/gcc-7.3.0/bin/g++
export COMPILER_SPEC="gcc@7.3.0"
export MPIDIR=/usr/tce/packages/mvapich2/mvapich2-2.3-gcc-7.3.0
export MPIEXEC=$(which srun)

# number of threads in OpenMP examples
export OMP_NUM_THREADS=16

# module loads
module load gcc/7.3.0 mvapich2/2.3
