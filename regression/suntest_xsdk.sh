#!/bin/bash
# -------------------------------------------------------------------------------
# Programmer(s): David J. Gardner @ LLNL 
# -------------------------------------------------------------------------------
# LLNS Copyright Start
# Copyright (c) 2014, Lawrence Livermore National Security
# This work was performed under the auspices of the U.S. Department 
# of Energy by Lawrence Livermore National Laboratory in part under 
# Contract W-7405-Eng-48 and in part under Contract DE-AC52-07NA27344.
# Produced at the Lawrence Livermore National Laboratory.
# All rights reserved.
# For details, see the LICENSE file.
# LLNS Copyright End
# -------------------------------------------------------------------------------
# xsdk variant of the SUNDIALS regression testing script with all external 
# libraries enabled
# -------------------------------------------------------------------------------

# check number of inputs
if [ "$#" -lt 2 ]; then
    echo "ERROR: Illegal number of parameters, real and index type required"
    exit 1
fi
realtype=$1     # required, precision for realtypes
indextype=$2    # required, integer type for indices
buildthreads=$3 # optional, number of build threads (if empty will use all threads)

# remove old build and install directories
\rm -rf build_xSDK_${realtype}_${indextype}
\rm -rf install_xSDK_${realtype}_${indextype}

# create new build and install directories
mkdir build_xSDK_${realtype}_${indextype}
mkdir install_xSDK_${realtype}_${indextype}

# move to build directory
cd build_xSDK_${realtype}_${indextype}

# number of threads in OpenMP examples
export OMP_NUM_THREADS=4

# set file permissions (rwxrwxr-x)
umask 002

# -------------------------------------------------------------------------------
# Installed Third Party Libraries
# -------------------------------------------------------------------------------

# path to installed libraries
APPDIR=/usr/casc/sundials/apps/rh6

# MPI
MPIDIR=${APPDIR}/openmpi/1.8.8/bin

# LAPACK / BLAS
BLASSTATUS=ON
BLASDIR=${APPDIR}/lapack/3.6.0/lib64

LAPACKSTATUS=ON
LAPACKDIR=${APPDIR}/lapack/3.6.0/lib64

# LAPACK/BLAS does not support extended precision or 64-bit indices
if [ "$realtype" == "extended" ] || [ "$indextype" == "int64_t" ]; then
    LAPACKSTATUS=OFF
    BLASSTATUS=OFF
fi

# KLU
KLUSTATUS=ON
KLUDIR=${APPDIR}/suitesparse/4.5.3

# KLU does not support single or extended precision
if [ "$realtype" == "single" ] || [ "$realtype" == "extended" ]; then
    KLUSTATUS=OFF
fi

# SuperLU_MT
SUPERLUMTSTATUS=ON

# SuperLU MT index type must be set a build time
if [ "$indextype" == "int32_t" ]; then
    SUPERLUMTDIR=${APPDIR}/superlu_mt/SuperLU_MT_3.1_fpic
else
    SUPERLUMTDIR=${APPDIR}/superlu_mt/SuperLU_MT_3.1_long_int_fpic
fi

# SuperLU MT does not support extended precision
if [ "$realtype" == "extended" ]; then
    SUPERLUMTSTATUS=OFF
fi

# hypre
HYPRESTATUS=ON

# hypre index type must be set a build time
if [ "$indextype" == "int32_t" ]; then
    HYPREDIR=${APPDIR}/hypre/2.11.1_fpic
else
    HYPREDIR=${APPDIR}/hypre/2.11.1_long_int_fpic
fi

# only testing hypre with double precision at this time
if [ "$realtype" != "double" ]; then
    HYPRESTATUS=OFF
fi

# PETSc
PETSCSTATUS=ON

# PETSc index type must be set a build time
if [ "$indextype" == "int32_t" ]; then
    PETSCDIR=${APPDIR}/petsc/3.7.2
else
    PETSCDIR=${APPDIR}/petsc/3.7.2_long_int
fi

# only testing PETSc with double precision at this time
if [ "$realtype" != "double" ]; then
    PETSCSTATUS=OFF
fi

# -------------------------------------------------------------------------------
# Configure SUNDIALS with CMake
#
# NOTE: Helpful options for debugging CMake
#
# The '-LAH' flag lists the non-advanced cached variables (L), the advanced
# variables (A), and help for each variable (H). This will not print any system
# variables.
#
# The CMake option '-D CMAKE_VERBOSE_MAKEFILE=ON' enables additional output during
# compile time which is useful for debugging build issues.
# -------------------------------------------------------------------------------

# set realtype for xSDK build
if [ $realtype == "extended" ]; then
    xsdk_realtype="quad"
else
    xsdk_realtype=$realtype
fi

# set indextype for xSDK build
if [ $indextype == "int32_t" ]; then
    xsdk_indextype="32"
else
    xsdk_indextype="64"
fi

echo "START CMAKE"
cmake \
    -D CMAKE_INSTALL_PREFIX="../install_xSDK_${realtype}_${indextype}" \
    \
    -D XSDK_PRECISION=${xsdk_realtype} \
    -D XSDK_INDEX_SIZE=${xsdk_indextype} \
    \
    -D XSDK_ENABLE_FORTRAN=ON \
    \
    -D EXAMPLES_ENABLE_C=ON \
    -D EXAMPLES_ENABLE_CXX=ON \
    -D EXAMPLES_ENABLE_F77=ON \
    -D EXAMPLES_ENABLE_F90=ON \
    \
    -D OPENMP_ENABLE=ON \
    -D PTHREAD_ENABLE=ON \
    -D CUDA_ENABLE=OFF \
    -D RAJA_ENABLE=OFF \
    \
    -D CMAKE_C_COMPILER="/usr/bin/cc" \
    -D CMAKE_CXX_COMPILER="/usr/bin/c++" \
    -D CMAKE_Fortran_COMPILER="/usr/bin/gfortran" \
    \
    -D CMAKE_C_FLAGS='-g -Wall -std=c99 -pedantic' \
    -D CMAKE_CXX_FLAGS='-g' \
    -D CMAKE_Fortran_FLAGS='-g' \
    \
    -D MPI_ENABLE=ON \
    -D MPI_MPICC="${MPIDIR}/mpicc" \
    -D MPI_MPICXX="${MPIDIR}/mpicxx" \
    -D MPI_MPIF77="${MPIDIR}/mpif77" \
    -D MPI_MPIF90="${MPIDIR}/mpif90" \
    -D MPI_RUN_COMMAND="${MPIDIR}/mpirun" \
    \
    -D TPL_ENABLE_BLAS="${BLASSTATUS}" \
    -D TPL_BLAS_LIBRARIES="${BLASDIR}/libblas.so" \
    \
    -D TPL_ENABLE_LAPACK="${LAPACKSTATUS}" \
    -D TPL_LAPACK_LIBRARIES="${LAPACKDIR}/liblapack.so" \
    \
    -D TPL_ENABLE_KLU="${KLUSTATUS}" \
    -D TPL_KLU_INCLUDE_DIRS="${KLUDIR}/include" \
    -D TPL_KLU_LIBRARIES="${KLUDIR}/lib/libklu.a" \
    \
    -D TPL_ENABLE_HYPRE="${HYPRESTATUS}" \
    -D TPL_HYPRE_INCLUDE_DIRS="${HYPREDIR}/include" \
    -D TPL_HYPRE_LIBRARIES="${HYPREDIR}/lib/libHYPRE.a" \
    \
    -D TPL_ENABLE_PETSC="${PETSCSTATUS}" \
    -D TPL_PETSC_INCLUDE_DIRS="${PETSCDIR}/include" \
    -D TPL_PETSC_LIBRARIES="${PETSCDIR}/lib/libpetsc.so" \
    \
    -D TPL_ENABLE_SUPERLUMT="${SUPERLUMTSTATUS}" \
    -D TPL_SUPERLUMT_INCLUDE_DIRS="${SUPERLUMTDIR}/SRC" \
    -D TPL_SUPERLUMT_LIBRARIES="${SUPERLUMTDIR}/lib/libsuperlu_mt_PTHREAD.a" \
    -D TPL_SUPERLUMT_THREAD_TYPE=Pthread \
    \
    -D SUNDIALS_DEVTESTS=ON \
    ../../. 2>&1 | tee configure.log

# check cmake return code
rc=${PIPESTATUS[0]}
echo -e "\ncmake returned $rc\n" | tee -a configure.log
if [ $rc -ne 0 ]; then exit 1; fi

# -------------------------------------------------------------------------------
# Make SUNDIALS
# -------------------------------------------------------------------------------

echo "START MAKE"
make -j $buildthreads 2>&1 | tee make.log

# check make return code
rc=${PIPESTATUS[0]}
echo -e "\nmake returned $rc\n" | tee -a make.log
if [ $rc -ne 0 ]; then exit 1; fi

# -------------------------------------------------------------------------------
# Test SUNDIALS
# -------------------------------------------------------------------------------

# test sundials
echo "START TEST"
make test 2>&1 | tee test.log

# check make test return code
rc=${PIPESTATUS[0]}
echo -e "\nmake test returned $rc\n" | tee -a test.log
if [ $rc -ne 0 ]; then exit 1; fi

# -------------------------------------------------------------------------------
# Install SUNDIALS
# -------------------------------------------------------------------------------

# install sundials
echo "START INSTALL"
make install 2>&1 | tee install.log

# check make install return code
rc=${PIPESTATUS[0]}
echo -e "\nmake install returned $rc\n" | tee -a install.log
if [ $rc -ne 0 ]; then exit 1; fi

# -------------------------------------------------------------------------------
# Test SUNDIALS Install
# -------------------------------------------------------------------------------

# add make test_install here

# -------------------------------------------------------------------------------
# Return
# -------------------------------------------------------------------------------

# if we make it here all tests have passed
exit 0