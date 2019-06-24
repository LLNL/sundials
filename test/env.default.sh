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

# if SPACK_ROOT is not set, check for the shared spack installation
if [ -z "$SPACK_ROOT" ]; then
    if [ -d "/usr/casc/sundials/apps/spack" ]; then
        export SPACK_ROOT=/usr/casc/sundials/apps/spack
    else
        echo "ERROR: Could not locate spack installation"
        return 1
    fi
fi

# make sure the spack environment is setup
source ${SPACK_ROOT}/share/spack/setup-env.sh

# path to libraries not installed through spack
APPDIR=/usr/casc/sundials/apps/rh6

# compilers
export CC="$(spack location -i "$compiler")/bin/gcc"
export CXX="$(spack location -i "$compiler")/bin/g++"
export FC="$(spack location -i "$compiler")/bin/gfortran"

# compiler flags (test scripts will append C/C++ standard flags)
export CFLAGS="-g -Wall -Wpedantic -Werror"
export CXXFLAGS="-g -Wall -Wpedantic -Werror"
export FFLAGS="-g -Wall -Wpedantic -ffpe-summary=none"

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
MPIDIR="$(spack location -i openmpi@3.1.2 % "$compiler")"
export MPICC="${MPIDIR}/bin/mpicc"
export MPICXX="${MPIDIR}/bin/mpicxx"
export MPIFC="${MPIDIR}/bin/mpifort"
export MPIEXEC="${MPIDIR}/bin/mpirun"

# LAPACK / BLAS
if [ "$realtype" == "extended" ] || [ "$indexsize" == "64" ]; then
    export BLASSTATUS=OFF
    export LAPACKSTATUS=OFF
else
    export BLASSTATUS=ON
    export LAPACKSTATUS=ON
fi
BLASDIR="$(spack location -i openblas@0.3.5~ilp64 % "$compiler")"
export BLASLIBS=${BLASDIR}/lib/libopenblas.so
export LAPACKLIBS=${BLASLIBS}

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
        export SLUMTDIR=${APPDIR}/superlu_mt/SuperLU_MT_3.1_fpic
    else
        export SLUMTDIR=${APPDIR}/superlu_mt/SuperLU_MT_3.1_long_int_fpic
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
        export SLUDISTDIR="$(spack location -i superlu-dist@develop~int64+openmp % "$compiler")"
    else
        export SLUDISTDIR="$(spack location -i superlu-dist@develop+int64+openmp % "$compiler")"
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
        export HYPREDIR="$(spack location -i hypre@2.14.0~int64 % "$compiler")"
    else
        export HYPREDIR="$(spack location -i hypre@2.14.0+int64 % "$compiler")"
    fi
fi

# petsc
if [ "$realtype" != "double" ]; then
    export PETSCSTATUS=OFF
    unset PETSCDIR
else
    export PETSCSTATUS=ON
    if [ "$indexsize" == "32" ]; then
        export PETSCDIR="$(spack location -i petsc@3.10.3~int64 % $compiler)"
    else
        export PETSCDIR="$(spack location -i petsc@3.10.3+int64 % $compiler)"
    fi
fi
