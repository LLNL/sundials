#!/bin/bash
# ------------------------------------------------------------------------------
# Programmer(s): Cody J. Balos and David J. Gardner @ LLNL
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
# Script that sets up the default SUNDIALS testing environment.
#
# Usage: source docker.sh <compiler spec> <build type>
#
# Optional Inputs:
#   <compiler spec> = compiler to build SUNDIALS with
#   <build type>    = SUNDIALS build type:
#                       dbg : debug build
#                       opt : optimized build
# ------------------------------------------------------------------------------

echo "./docker.sh $*" | tee -a configure.log

# set defaults for optional inputs
compiler="gcc" # compiler spec
bldtype="dbg"  # build type dbg = debug or opt = optimized

if [ "$#" -gt 1 ]; then
    bldtype=$2
fi

# get compiler name and version from spec
compilername="${compiler%%@*}"
compilerversion="${compiler##*@}"

# ------------------------------------------------------------------------------
# Check input values
# ------------------------------------------------------------------------------

case "$SUNDIALS_PRECISION" in
    single|double|extended) ;;
    *)
        echo "ERROR: Unknown real type option: $SUNDIALS_PRECISION"
        return 1
        ;;
esac

case "$SUNDIALS_INDEX_SIZE" in
    32|64) ;;
    *)
        echo "ERROR: Unknown index size option: $SUNDIALS_INDEX_SIZE"
        return 1
        ;;
esac

case "$compiler" in
    gcc);;
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

# set file permissions (rwxrwxr-x)
umask 002

# path to shared installs
APPROOT=/opt

# setup the python environment
source ${APPROOT}/python-venv/sundocs/bin/activate

# add CUDA
if [[ ":${PATH}:" != *":/usr/local/cuda-11.5/bin:"* ]]; then
    export PATH="/usr/local/cuda-11.5/bin${PATH:+:${PATH}}"
fi

if [[ ":${LD_LIBRARY_PATH}:" != *":/usr/local/cuda-11.5/lib64:"* ]]; then
    export LD_LIBRARY_PATH="/usr/local/cuda-11.5/lib64${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"
fi

# ------------------------------------------------------------------------------
# Compilers and flags
# ------------------------------------------------------------------------------

if [ "$compilername" == "gcc" ]; then

    export CC=$(which gcc)
    export CXX=$(which g++)
    export FC=$(which gfortran)

    # optimization flags
    if [ "$bldtype" == "dbg" ]; then
        export CFLAGS="-g -O0"
        export CXXFLAGS="-g -O0"
        export FFLAGS="-g -O0"
        export CUDAFLAGS="-g -O0"
    else
        export CFLAGS="-g -O3"
        export CXXFLAGS="-g -O3"
        export FFLAGS="-g -O3"
        export CUDAFLAGS="-g -O3"
    fi

    # append additional warning flags
    if [[ "$SUNDIALS_PRECISION" == "double" && "$SUNDIALS_INDEX_SIZE" == "32" ]]; then
        export CFLAGS="${CFLAGS} -Wconversion -Wno-sign-conversion"
        export CXXFLAGS="${CXXFLAGS} -Wconversion -Wno-sign-conversion"
    fi

fi

# ------------------------------------------------------------------------------
# SUNDIALS Options
# ------------------------------------------------------------------------------

# Verbose build
export CMAKE_VERBOSE_MAKEFILE=ON

# Number of build and test jobs
export SUNDIALS_BUILD_JOBS=4
export SUNDIALS_TEST_JOBS=1

# Sundials packages
export SUNDIALS_ARKODE=ON
export SUNDIALS_CVODE=ON
export SUNDIALS_CVODES=ON
export SUNDIALS_IDA=ON
export SUNDIALS_IDAS=ON
export SUNDIALS_KINSOL=ON

# Fortran interface status
if [ "$compilername" == "gcc" ]; then
    if [[ ("$SUNDIALS_PRECISION" == "double") && ("$SUNDIALS_INDEX_SIZE" == "64") ]]; then
        export SUNDIALS_FMOD_INTERFACE=ON
    else
        export SUNDIALS_FMOD_INTERFACE=OFF
    fi
else
    export SUNDIALS_FMOD_INTERFACE=OFF
fi

# Sundials benchmarks
export SUNDIALS_BENCHMARKS=ON

# Sundials monitoring
export SUNDIALS_MONITORING=ON

# Sundials profiling
export SUNDIALS_PROFILING=ON

# Sundials logging
export SUNDIALS_LOGGING_LEVEL=4
export SUNDIALS_LOGGING_ENABLE_MPI=ON

# Uncomment to override the default output file comparison precisions. The float
# precision is number of digits to compare (0 = all digits) and the integer
# precision is allowed percentage difference (0 = no difference).
export SUNDIALS_TEST_FLOAT_PRECISION=0
export SUNDIALS_TEST_INTEGER_PRECISION=0

# FindMPI fails with this ON
export SUNDIALS_ENABLE_WARNINGS_AS_ERRORS=OFF

# ------------------------------------------------------------------------------
# Third party libraries
# ------------------------------------------------------------------------------

# -------
# PThread
# -------

export SUNDIALS_PTHREAD=ON

# ------
# OpenMP
# ------

export SUNDIALS_OPENMP=ON
export OMP_NUM_THREADS=4

# ---------------------
# OpenMP Device Offload
# ---------------------

export SUNDIALS_OPENMP_OFFLOAD=OFF

# ----
# CUDA
# ----

if [ "$SUNDIALS_PRECISION" != "extended" ]; then
    export SUNDIALS_CUDA=ON
    export CUDAARCHS=60
    export SUNDIALS_FUSED_KERNELS=ON
else
    export SUNDIALS_CUDA=OFF
    export SUNDIALS_FUSED_KERNELS=OFF
fi

# ---
# MPI
# ---

MPI_ROOT=/opt/view

export SUNDIALS_MPI=ON
export MPICC="${MPI_ROOT}/bin/mpicc"
export MPICXX="${MPI_ROOT}/bin/mpicxx"
export MPIFC="${MPI_ROOT}/bin/mpifort"
export MPIEXEC="${MPI_ROOT}/bin/mpirun"

# -------------
# LAPACK / BLAS
# -------------

if [ "$SUNDIALS_PRECISION" != "extended" ]; then
    export SUNDIALS_LAPACK=ON
    export LAPACK_ROOT=/opt/view
    export LAPACK_LIBRARIES="${LAPACK_ROOT}/lib/libopenblas.so"
else
    export SUNDIALS_LAPACK=OFF
    unset LAPACK_LIBRARIES
fi

# -----
# MAGMA
# -----

if [ "$SUNDIALS_PRECISION" != "extended" ] && \
    [ "$SUNDIALS_INDEX_SIZE" == "32" ] && \
    [ "$SUNDIALS_CUDA" == "ON" ]; then
    export SUNDIALS_MAGMA=ON
    export MAGMA_ROOT=/opt/view
    export MAGMA_BACKENDS="CUDA"
else
    export SUNDIALS_MAGMA=OFF
    unset MAGMA_ROOT
    unset MAGMA_BACKENDS
fi

# ------------
# SuperLU_DIST
# ------------

if [ "$SUNDIALS_PRECISION" == "double" ]; then
    export SUNDIALS_SUPERLU_DIST=ON
    export SUPERLU_DIST_ROOT=/opt/view
    export SUPERLU_DIST_INCLUDE_DIR="${SUPERLU_DIST_ROOT}/include"
    export SUPERLU_DIST_LIBRARY_DIR="${SUPERLU_DIST_ROOT}/lib"

    # built with 32-bit blas
    export BLAS_ROOT=/opt/view
    export BLAS_LIBRARIES=${BLAS_ROOT}/lib/libopenblas.so

    # PARMETIS
    export PARMETIS_ROOT=/opt/view
    export PARMETIS_LIBRARIES="${PARMETIS_ROOT}/lib/libparmetis.so"

    # METIS
    export METIS_ROOT=/opt/view
    export METIS_LIBRARIES="${METIS_ROOT}/lib/libmetis.so"

    export SUPERLU_DIST_LIBRARIES="${BLAS_LIBRARIES};${PARMETIS_LIBRARIES};${METIS_LIBRARIES};${SUPERLU_DIST_ROOT}/lib/libsuperlu_dist.a"
    export SUPERLU_DIST_OPENMP=OFF
else
    export SUNDIALS_SUPERLU_DIST=OFF
    unset SUPERLU_DIST_INCLUDE_DIR
    unset SUPERLU_DIST_LIBRARY_DIR
    unset SUPERLU_DIST_LIBRARIES
    unset SUPERLU_DIST_OPENMP
fi

# -----
# hypre
# -----

if [ "$SUNDIALS_PRECISION" == "double" ]; then
    export SUNDIALS_HYPRE=ON
    export HYPRE_ROOT=/opt/view
    export HYPRE_INCLUDE_DIR="${HYPRE_ROOT}/include"
    export HYPRE_LIBRARY_DIR="${HYPRE_ROOT}/lib"
else
    export SUNDIALS_HYPRE=OFF
    unset HYPRE_INCLUDE_DIR
    unset HYPRE_LIBRARY_DIR
fi

# -----
# petsc
# -----

if [ "$SUNDIALS_PRECISION" == "double" ]; then
    export PETSC_ROOT=/opt/view
else
    export SUNDIALS_PETSC=OFF
    unset PETSC_ROOT
fi

# --------
# trilinos
# --------

# if [ "$SUNDIALS_PRECISION" == "double" ] && [ "$SUNDIALS_INDEX_SIZE" == "32" ]; then
#     export SUNDIALS_TRILINOS=ON
#     export TRILINOS_ROOT=/opt/view
# else
#     export SUNDIALS_TRILINOS=OFF
#     unset TRILINOS_ROOT
# fi

# ----
# raja
# ----

if [ "$SUNDIALS_PRECISION" == "double" ]; then
    export SUNDIALS_RAJA=ON
    export RAJA_ROOT="$(spack location -i raja % "$compiler")"
    export RAJA_BACKENDS="CUDA"
else
    export SUNDIALS_RAJA=OFF
    unset RAJA_ROOT
    unset RAJA_BACKENDS
fi

# ------
# xbraid
# ------

# if [ "$SUNDIALS_PRECISION" == "double" ] && [ "$SUNDIALS_INDEX_SIZE" == "32" ]; then
#     export SUNDIALS_XBRAID=ON
#     export XBRAID_ROOT=/opt/view
# else
#     export SUNDIALS_XBRAID=OFF
#     unset XBRAID_ROOT
# fi
