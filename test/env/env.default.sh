#!/bin/bash
# -------------------------------------------------------------------------------
# Programmer(s): Cody J. Balos and David J. Gardner @ LLNL
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
# Script that sets up the default SUNDIALS testing environment.
#
# Usage: source env.default.sh <real type> <index size> <compiler spec> \
#                              <build type>
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
#                       gcc@4.9.4
#   <build type>    = SUNDIALS build type:
#                       dbg : debug build
#                       opt : optimized build
#
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
compiler="gcc@4.9.4" # compiler spec
bldtype="dbg"        # build type dbg = debug or opt = optimized

# set optional inputs if provided
if [ "$#" -gt 2 ]; then
    compiler=$3
fi

if [ "$#" -gt 3 ]; then
    bldtype=$4
fi

# ------------------------------------------------------------------------------
# Check input values
# ------------------------------------------------------------------------------

case "$realtype" in
    single|double|extended) ;;
    *)
        echo "ERROR: Unknown real type option: $realtype"
        return 1
        ;;
esac

case "$indexsize" in
    32|64) ;;
    *)
        echo "ERROR: Unknown index size option: $indexsize"
        return 1
        ;;
esac

case "$compiler" in
    gcc@4.9.4|gcc@8.3.0|gcc@9.2.0) ;;
    clang@9.0.0);;
    *)
        echo "ERROR: Unknown compiler spec: $compiler"
        return 1
        ;;
esac

case "$bldtype" in
    dbg|opt) ;;
    *)
        echo "ERROR: Unknown build type: $bldtype"
        return 1
        ;;
esac

# ------------------------------------------------------------------------------
# Setup environment
# ------------------------------------------------------------------------------

# get compiler name and version from spec
compilername="${compiler%%@*}"
compilerversion="${compiler##*@}"

# set file permissions (rwxrwxr-x)
umask 002

# path to libraries not installed through spack
APPDIR=/usr/casc/sundials/share/sunenv/apps/${compilername}-${compilerversion}

# ------------------------------------------------------------------------------
# Compilers and flags
# ------------------------------------------------------------------------------

if [ "$compilername" == "gcc" ]; then
    COMPILER_DIR="$(spack location -i "$compiler")"
    export CC="${COMPILER_DIR}/bin/gcc"
    export CXX="${COMPILER_DIR}/bin/g++"
    export FC="${COMPILER_DIR}/bin/gfortran"
else
    COMPILER_DIR="$(spack location -i "llvm@$compilerversion")"
    export CC="${COMPILER_DIR}/bin/clang"
    export CXX="${COMPILER_DIR}/bin/clang++"
fi

# compiler flags (test scripts will append C/C++ standard flags)
if [ "$bldtype" == "dbg" ]; then
    export CFLAGS="-g -O0"
    export CXXFLAGS="-g -O0"
    export FFLAGS="-g -O0"
else
    export CFLAGS="-g -O3"
    export CXXFLAGS="-g -O3"
    export FFLAGS="-g -O3"
fi

# append additional compiler flags
export CFLAGS="${CFLAGS} -Wall -Wpedantic -Werror"
export CXXFLAGS="${CXXFLAGS} -Wall -Wpedantic -Werror"
export FFLAGS="${FFLAGS} -Wall -Wpedantic -ffpe-summary=none"

if [[ "$realtype" == "double" && "$indexsize" == "32" ]]; then
    export CFLAGS="${CFLAGS} -Wconversion -Wno-sign-conversion"
    export CXXFLAGS="${CXXFLAGS} -Wconversion -Wno-sign-conversion"
fi

# ------------------------------------------------------------------------------
# SUNDIALS Options
# ------------------------------------------------------------------------------

# Sundials packages
export ARKODE_STATUS=ON
export CVODE_STATUS=ON
export CVODES_STATUS=ON
export IDA_STATUS=ON
export IDAS_STATUS=ON
export KINSOL_STATUS=ON

# Fortran interface status
if [ "$compilername" == "gcc" ]; then
    export F77_STATUS=ON
    export F03_STATUS=ON
else
    export F77_STATUS=OFF
    export F03_STATUS=OFF
fi

# ------------------------------------------------------------------------------
# Third party libraries
# ------------------------------------------------------------------------------

# PThread settings
export PTHREAD_STATUS=ON

# OpenMP settings
export OPENMP_STATUS=ON
export OMP_NUM_THREADS=4

# OpenMP DEV settings
export OPENMPDEV_STATUS=OFF

# CUDA settings
if [[ ("$compiler" != "gcc@9.2.0") && ("$compiler" != "clang@9.0.0") ]]; then
    if [ "$realtype" == "extended" ]; then
        export CUDA_STATUS=OFF
    else
        export CUDA_STATUS=ON
    fi
fi

# MPI
export MPI_STATUS=ON
MPIDIR="$(spack location -i openmpi@3.1.4 % "$compiler")"
export MPICC="${MPIDIR}/bin/mpicc"
export MPICXX="${MPIDIR}/bin/mpicxx"
export MPIFC="${MPIDIR}/bin/mpifort"
export MPIEXEC="${MPIDIR}/bin/mpirun"

# LAPACK / BLAS
if [ "$realtype" == "extended" ]; then
    export LAPACK_STATUS=OFF
    unset LAPACKLIBS
else
    export LAPACK_STATUS=ON
    if [ "$indexsize" == "32" ]; then
        LAPACKDIR="$(spack location -i openblas@0.3.7~ilp64 % "$compiler")"
    else
        LAPACKDIR="$(spack location -i openblas@0.3.7+ilp64 % "$compiler")"
    fi
    export LAPACKLIBS="${LAPACKDIR}/lib/libopenblas.so"
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
    export KLU_STATUS=OFF
    unset KLUDIR
else
    export KLU_STATUS=ON
    export KLUDIR="$(spack location -i suite-sparse@5.3.0 % "$compiler")"
fi

# SuperLU_MT
if [ "$realtype" == "extended" ]; then
    export SLUMT_STATUS=OFF
    unset SLUMTDIR
else
    export SLUMT_STATUS=ON
    if [ "$indexsize" == "32" ]; then
        export SLUMTDIR="$(spack location -i superlu-mt@3.1~int64~blas % "$compiler")"
    else
        export SLUMTDIR="$(spack location -i superlu-mt@3.1+int64~blas % "$compiler")"
    fi
    export SLUMTLIBS="${SLUMTDIR}/lib/libblas_PTHREAD.a"
    export SLUMTTYPE="PTHREAD"
fi

# SuperLU_DIST
if [ "$realtype" != "double" ]; then
    export SLUDIST_STATUS=OFF
    unset SLUDISTDIR
    unset SLUDISTLIBS
else
    export SLUDIST_STATUS=ON
    if [ "$indexsize" == "32" ]; then
        export SLUDISTDIR="$(spack location -i superlu-dist@6.1.1~int64+openmp % "$compiler")"
    else
        export SLUDISTDIR="$(spack location -i superlu-dist@6.1.1+int64+openmp % "$compiler")"
    fi
    # built with 32-bit blas
    BLASDIR="$(spack location -i openblas@0.3.7~ilp64 % "$compiler")"
    BLASLIBS=${BLASDIR}/lib/libopenblas.so
    export SLUDISTLIBS="${BLASLIBS};${PARMETISLIB};${METISLIB};${SLUDISTDIR}/lib/libsuperlu_dist.a"
fi

# hypre
if [ "$realtype" != "double" ]; then
    export HYPRE_STATUS=OFF
    unset HYPREDIR
else
    export HYPRE_STATUS=ON
    if [ "$indexsize" == "32" ]; then
        export HYPREDIR="$(spack location -i hypre@2.18.2~int64 % "$compiler")"
    else
        export HYPREDIR="$(spack location -i hypre@2.18.2+int64 % "$compiler")"
    fi
fi

# petsc
if [ "$realtype" != "double" ]; then
    export PETSC_STATUS=OFF
    unset PETSCDIR
else
    export PETSC_STATUS=ON
    if [ "$indexsize" == "32" ]; then
        export PETSCDIR="$(spack location -i petsc@3.12.1~int64 % $compiler)"
    else
        export PETSCDIR="$(spack location -i petsc@3.12.1+int64 % $compiler)"
    fi
fi

# trilinos
if [ "$realtype" == "double" ] && [ "$indexsize" == "32" ]; then
    export TRILINOS_STATUS=ON
    export TRILINOSDIR="$(spack location -i trilinos@12.14.1 % $compiler)"
else
    export TRILINOS_STATUS=OFF
    unset TRILINOSDIR
fi

# raja
RAJA_STATUS=OFF
if [[ ("$compiler" != "gcc@9.2.0") && ("$compiler" != "clang@9.0.0") ]]; then
    if [ "$realtype" == "double" ]; then
        RAJA_STATUS=ON
        export RAJADIR=${APPDIR}/raja-0.10.0
    fi
fi
