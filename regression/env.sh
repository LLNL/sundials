#!/bin/bash
# -------------------------------------------------------------------------------
# Programmer(s): Cody J. Balos @ LLNL
# -------------------------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2019, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# -------------------------------------------------------------------------------
# Script that sets up the environment for the suntest_driver.sh when running on
# the regression server.
# -------------------------------------------------------------------------------

# set path shared spack installation
export SPACK_ROOT=/usr/casc/sundials/apps/spack

# setup the spack environment
source ${SPACK_ROOT}/share/spack/setup-env.sh

# compiler settings
export CC=${SPACK_ROOT}/opt/spack/linux-rhel6-x86_64/gcc-4.4.7/gcc-4.9.4-lpphh5vokyvegcvsr5on7unx6pzesbcb/bin/gcc
export FC=${SPACK_ROOT}/opt/spack/linux-rhel6-x86_64/gcc-4.4.7/gcc-4.9.4-lpphh5vokyvegcvsr5on7unx6pzesbcb/bin/gfortran
export CXX=${SPACK_ROOT}/opt/spack/linux-rhel6-x86_64/gcc-4.4.7/gcc-4.9.4-lpphh5vokyvegcvsr5on7unx6pzesbcb/bin/g++
export COMPILER_SPEC="gcc@4.9.4"
export MPIDIR=""
export MPIEXEC=""

# number of threads in OpenMP examples
export OMP_NUM_THREADS=4
