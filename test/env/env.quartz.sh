#!/bin/bash
# -------------------------------------------------------------------------------
# Programmer(s): Cody J. Balos @ LLNL
# -------------------------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2021, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# -------------------------------------------------------------------------------
# Script that sets up the environment for SUNDIALS testing on Quartz.
#
# Usage: source env.quartz.sh <real type> <index size> <compiler spec> \
#                             <build type>
#
# Required Inputs:
#   <real type>  = SUNDIALS real type to build/test with:
#                    single   : single (32-bit) precision
#                    double   : double (64-bit) precision
#                    extended : extended (128-bit) precision
#   <index size> = SUNDIALS index size to build/test with:
#                    32       : 32-bit indices
#                    64       : 64-bit indices
#
# Optional Inputs:
#   <compiler spec> = Compiler to build sundials with:
#                       e.g., gcc@8.1.0, intel@19.0.4, clang@9.0.0, etc.
#   <build type>    = SUNDIALS build type:
#                       dbg : debug build
#                       opt : optimized build
# -------------------------------------------------------------------------------

# check number of inputs
if [ "$#" -lt 2 ]; then
    echo "ERROR: TWO (2) inputs required"
    echo "real type  : [single|double|extended]"
    echo "index size : [32|64]"
    return 1
fi

# set required inputs
realtype=$1   # precision for realtypes
indexsize=$2  # integer size for indices

# set defaults for optional inputs
compiler="gcc@8.1.0" # compiler spec
bldtype="dbg"        # build type dbg = debug or opt = optimized

# set optional inputs if provided
if [ "$#" -gt 2 ]; then
    compiler=$3
fi

if [ "$#" -gt 3 ]; then
    bldtype=$4
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
    intel)
        module load intel/${compilerversion}
        if [ $? -ne 0 ]; then return 1; fi
        export CC=$(which icc)
        export CXX=$(which icpc)
        export FC=$(which ifort)
        FORTRAN_STATUS=ON
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
export F77_STATUS=${FORTRAN_STATUS}
if [[ ("$realtype" == "double") && ("$indexsize" == "64") ]]; then
    export F03_STATUS=${FORTRAN_STATUS}
else
    export F03_STATUS=OFF
fi

# set MPI compiler wrapper
export MPI_STATUS=ON
export MPICC=$(which mpicc)
export MPICXX=$(which mpicxx)
export MPIFC=$(which mpifort)
export MPIEXEC=$(which srun)

# PThread settings
export PTHREAD_STATUS=ON

# OpenMP settings
export OPENMP_STATUS=ON
export OMP_NUM_THREADS=16
