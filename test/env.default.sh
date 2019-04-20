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
# -------------------------------------------------------------------------------

# set the compiler spec for spack
compiler="gcc@4.9.4"

# set file permissions (rwxrwxr-x)
umask 002

# set path shared spack installation
export SPACK_ROOT=/usr/casc/sundials/apps/spack

# setup the spack environment
source ${SPACK_ROOT}/share/spack/setup-env.sh

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

# path to libraries installed without spack
APPDIR=/usr/casc/sundials/apps/rh6

# MPI
export MPIDIR="$(spack location -i openmpi@3.1.2 % "$compiler")"
export MPIEXEC="${MPIDIR}/bin/mpirun"

# LAPACK / BLAS
export BLASDIR="$(spack location -i openblas@0.3.5~ilp64 % "$compiler")"
export BLAS_LIBRARIES=${BLASDIR}/lib/libopenblas.so
export LAPACK_LIBRARIES=${BLAS_LIBRARIES}

# PARMETIS
PARMETISDIR_32="$(spack location -i parmetis ^metis~int64~real64 % "$compiler")"
PARMETISDIR_64="$(spack location -i parmetis ^metis+int64~real64 % "$compiler")"

PARMETISLIB_32="${PARMETISDIR_32}/lib/libparmetis.so"
PARMETISLIB_64="${PARMETISDIR_32}/lib/libparmetis.so"

# METIS
METISDIR_32="$(spack location -i metis~int64~real64 % "$compiler")"
METISDIR_64="$(spack location -i metis+int64~real64 % "$compiler")"

METISLIB_32="${METISDIR_32}/lib/libmetis.so"
METISLIB_64="${METISDIR_32}/lib/libmetis.so"

# KLU
export KLUDIR="$(spack location -i suite-sparse@5.3.0 % "$compiler")"

# SuperLU_MT
export SLUMTDIR_32=${APPDIR}/superlu_mt/SuperLU_MT_3.1_fpic
export SLUMTDIR_64=${APPDIR}/superlu_mt/SuperLU_MT_3.1_long_int_fpic

# SuperLU_DIST
export SLUDISTDIR_32="$(spack location -i superlu-dist@develop~int64+openmp % "$compiler")"
export SLUDISTDIR_64="$(spack location -i superlu-dist@develop+int64+openmp % "$compiler")"

export SLUDISTLIBS_32="${BLAS_LIBRARIES};${PARMETISLIB_32};${METISLIB_32};${SLUDISTDIR_32}/lib/libsuperlu_dist.a"
export SLUDISTLIBS_64="${BLAS_LIBRARIES};${PARMETISLIB_64};${METISLIB_64};${SLUDISTDIR_64}/lib/libsuperlu_dist.a"

# hypre
export HYPREDIR_32="$(spack location -i hypre@2.14.0~int64 % "$compiler")"
export HYPREDIR_64="$(spack location -i hypre@2.14.0+int64 % "$compiler")"

# petsc
export PETSCDIR_32="$(spack location -i petsc@3.10.3~int64 % $compiler)"
export PETSCDIR_64="$(spack location -i petsc@3.10.3+int64 % $compiler)"
