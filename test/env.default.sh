#!/bin/bash
# -------------------------------------------------------------------------------
# Programmer(s): Cody J. Balos and David J. Gardner @ LLNL
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
# Script that sets up the default SUNDIALS testing environment.
#
# Usage: source env.default.sh <real type> <index size>
#
# Required Inputs:
#   <real type>  = SUNDIALS real type to build/test with:
#                    single   : single (32-bit) precision
#                    double   : double (64-bit) precision
#                    extended : extended (128-bit) precision
#   <index size> = SUNDIALS index size to build/test with:
#                    32       : 32-bit indices
#                    64       : 64-bit indices
# -------------------------------------------------------------------------------

# check number of inputs
if [ "$#" -lt 2 ]; then
    echo "ERROR: TWO (2) inputs required"
    echo "real type  : [single|double|extended]"
    echo "index size : [32|64]"
    return 1
fi

realtype=$1          # precision for realtypes
indexsize=$2         # integer size for indices
compiler="gcc@4.9.4" # compiler spec for spack

# ------------------------------------------------------------------------------
# Check input values
# ------------------------------------------------------------------------------

# set real types to test
case "$realtype" in
    single|double|extended) ;;
    *)
        echo "ERROR: Unknown real type option: $realtype"
        return 1
        ;;
esac

# set index sizes to test
case "$indexsize" in
    32|64) ;;
    *)
        echo "ERROR: Unknown index size option: $indexsize"
        return 1
        ;;
esac

# ------------------------------------------------------------------------------
# Setup environment
# ------------------------------------------------------------------------------

# set file permissions (rwxrwxr-x)
umask 002

# path to libraries not installed through spack
APPDIR=/usr/casc/sundials/share/rh7/apps/gcc-4.9.4

# compilers
COMPILER_DIR="$(spack location -i "$compiler")"
export CC="${COMPILER_DIR}/bin/gcc"
export CXX="${COMPILER_DIR}/bin/g++"
export FC="${COMPILER_DIR}/bin/gfortran"

# compiler flags (test scripts will append C/C++ standard flags)
export BASE_CFLAGS="-g -Wall -Wpedantic -Werror"
export BASE_CXXFLAGS="-g -Wall -Wpedantic -Werror"
export BASE_FFLAGS="-g -Wall -Wpedantic -ffpe-summary=none"

if [[ "$realtype" == "double" && "$indexsize" == "32" ]]; then
    export BASE_CFLAGS="${BASE_CFLAGS} -Wconversion -Wno-sign-conversion"
    export BASE_CXXFLAGS="${BASE_CXXFLAGS} -Wconversion -Wno-sign-conversion"
fi

# OpenMP settings
export OMP_NUM_THREADS=4

# CUDA settings
if [ "$realtype" == "extended" ]; then
    export CUDASTATUS=OFF
else
    export CUDASTATUS=ON
fi

# MPI
export MPISTATUS=ON
MPIDIR="$(spack location -i openmpi@3.1.4 % "$compiler")"
export MPICC="${MPIDIR}/bin/mpicc"
export MPICXX="${MPIDIR}/bin/mpicxx"
export MPIFC="${MPIDIR}/bin/mpifort"
export MPIEXEC="${MPIDIR}/bin/mpirun"

# LAPACK / BLAS
if [ "$realtype" == "extended" ]; then
    export LAPACKSTATUS=OFF
else
    if [ "$indexsize" == "32" ]; then
        export LAPACKSTATUS=ON
    else
        export LAPACKSTATUS=OFF
    fi
    BLASDIR="$(spack location -i openblas@0.3.7~ilp64 % "$compiler")"
    export BLASLIBS=${BLASDIR}/lib/libopenblas.so
    export LAPACKLIBS=${BLASLIBS}
fi

# PARMETIS
if [ "$indexsize" == "32" ]; then
    PARMETISDIR="$(spack location -i parmetis ^metis~int64~real64 % "$compiler")"
else
    PARMETISDIR="$(spack location -i parmetis ^metis+int64~real64 % "$compiler")"
fi
PARMETISLIB="${PARMETISDIR}/lib/libparmetis.so"

# METIS
if [ "$indexsize" == "32" ]; then
    METISDIR="$(spack location -i metis~int64~real64 % "$compiler")"
else
    METISDIR="$(spack location -i metis+int64~real64 % "$compiler")"
fi
METISLIB="${METISDIR}/lib/libmetis.so"

# KLU
if [ "$realtype" != "double" ]; then
    export KLUSTATUS=OFF
    unset KLUDIR
else
    export KLUSTATUS=ON
    export KLUDIR="$(spack location -i suite-sparse@5.3.0 % "$compiler")"
fi

# SuperLU_MT
if [ "$realtype" == "extended" ]; then
    export SLUMTSTATUS=OFF
    unset SLUMTDIR
else
    export SLUMTSTATUS=ON
    if [ "$indexsize" == "32" ]; then
        export SLUMTDIR="$(spack location -i superlu-mt@3.1~int64~blas % "$compiler")"
    else
        export SLUMTDIR="$(spack location -i superlu-mt@3.1+int64~blas % "$compiler")"
    fi
fi

# SuperLU_DIST
if [ "$realtype" != "double" ]; then
    export SLUDISTSTATUS=OFF
    unset SLUDISTDIR
    unset SLUDISTLIBS
else
    export SLUDISTSTATUS=ON
    if [ "$indexsize" == "32" ]; then
        export SLUDISTDIR="$(spack location -i superlu-dist@6.1.1~int64+openmp % "$compiler")"
    else
        export SLUDISTDIR="$(spack location -i superlu-dist@6.1.1+int64+openmp % "$compiler")"
    fi
    export SLUDISTLIBS="${BLASLIBS};${PARMETISLIB};${METISLIB};${SLUDISTDIR}/lib/libsuperlu_dist.a"
fi

# hypre
if [ "$realtype" != "double" ]; then
    export HYPRESTATUS=OFF
    unset HYPREDIR
else
    export HYPRESTATUS=ON
    if [ "$indexsize" == "32" ]; then
        export HYPREDIR="$(spack location -i hypre@2.16.0~int64 % "$compiler")"
    else
        export HYPREDIR="$(spack location -i hypre@2.16.0+int64 % "$compiler")"
    fi
fi

# petsc
if [ "$realtype" != "double" ]; then
    export PETSCSTATUS=OFF
    unset PETSCDIR
else
    export PETSCSTATUS=ON
    if [ "$indexsize" == "32" ]; then
        export PETSCDIR="$(spack location -i petsc@3.11.3~int64 % $compiler)"
    else
        export PETSCDIR="$(spack location -i petsc@3.11.3+int64 % $compiler)"
    fi
fi

# trilinos
if [ "$realtype" == "double" ] && [ "$indexsize" == "32" ]; then
    TRILINOSSTATUS=ON
    export TRILINOSDIR="$(spack location -i trilinos@12.14.1 % $compiler)"
else
    TRILINOSSTATUS=OFF
fi

# raja
if [ "$realtype" == "double" ]; then
    RAJASTATUS=ON
    export RAJADIR=${APPDIR}/raja-0.9.0
else
    RAJASTATUS=OFF
fi
