#!/bin/bash
# ------------------------------------------------------------------------------
# Programmer(s): David J. Gardner @ LLNL
# ------------------------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2022, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# ------------------------------------------------------------------------------
# Script that sets up the environment for SUNDIALS testing on Summit.
#
# Usage: source env.summit.sh <compiler spec> <build type>
#
# Optional Inputs:
#   <compiler spec> = Compiler to build sundials with:
#                       e.g., xl@16.1.1-5, gcc@8.1.1, etc.
#   <build type>    = SUNDIALS build type:
#                       dbg : debug build
#                       opt : optimized build
# ------------------------------------------------------------------------------

echo "./summit.sh $*" | tee -a setup_env.log

# set defaults for optional inputs
compiler="xl@16.1.1-5" # compiler spec
bldtype="opt"          # build type dbg = debug or opt = optimized

# set optional inputs if provided
if [ "$#" -gt 0 ]; then
    compiler=$1
fi

if [ "$#" -gt 1 ]; then
    bldtype=$2
fi

# ------------------------------------------------------------------------------
# Setup environment
# ------------------------------------------------------------------------------

# get compiler name and version from spec
compilername="${compiler%%@*}"
compilerversion="${compiler##*@}"

# load default cmake module
module load cmake

# load compiler module
case "$compilername" in
    clang)
        module load clang/${compilerversion}
        if [ $? -ne 0 ]; then return 1; fi
        export CC=$(which clang)
        export CXX=$(which clang++)
        export FC=""
        FORTRAN_STATUS=OFF
        ;;
    gcc)
        module load gcc/${compilerversion}
        if [ $? -ne 0 ]; then return 1; fi
        export CC=$(which gcc)
        export CXX=$(which g++)
        export FC=$(which gfortran)
        FORTRAN_STATUS=ON
        ;;
    xl)
        module load xl/${compilerversion}
        if [ $? -ne 0 ]; then return 1; fi
        export CC=$(which xlc_r)
        export CXX=$(which xlc++_r)
        export FC=$(which xlf2003_r)
        FORTRAN_STATUS=OFF # Build issues with F2003 interface
        ;;
    pgi)
        module load pgi/${compilerversion}
        if [ $? -ne 0 ]; then return 1; fi
        export CC=$(which pgcc)
        export CXX=$(which pgc++)
        export FC=$(which pgfortran)
        FORTRAN_STATUS=ON
        ;;
    *)
        echo "ERROR: Invalid compiler option: $compiler"
        return 1
        ;;
esac

# set compiler flags
if [ "$bldtype" == "dbg" ]; then
    export CFLAGS="-g -O0"
    export CXXFLAGS="-g -O0"
    export FFLAGS="-g -O0"
elif [ "$bldtype" == "opt" ]; then
    export CFLAGS="-g -O3"
    export CXXFLAGS="-g -O3"
    export FFLAGS="-g -O3"
else
    echo "ERROR: Invalid build type: $bldtype"
    return 1
fi

# Fortran settings
if [[ ("$SUNDIALS_PRECISION" == "double") && ("$indexsize" == "64") ]]; then
    export SUNDIALS_FMOD_INTERFACE=${FORTRAN_STATUS}
else
    export SUNDIALS_FMOD_INTERFACE=OFF
fi

# Sundials monitoring
export SUNDIALS_MONITORING=ON

# set MPI compiler wrapper
export SUNDIALS_MPI=ON
export MPICC=$(which mpicc)
export MPICXX=$(which mpicxx)
export MPIFC=$(which mpifort)
export MPIEXEC=$(which srun)

# PThread settings
export SUNDIALS_PTHREAD=ON

# OpenMP settings
export SUNDIALS_OPENMP=ON
export OMP_NUM_THREADS=20

# LAPACK settings
module load essl
export SUNDIALS_LAPACK=ON
export LAPACK_LIBRARIES="${OLCF_ESSL_ROOT}/lib64/libessl.so"

# CUDA settings
module load cuda
export SUNDIALS_CUDA=ON
export CUDAARCHS=70
