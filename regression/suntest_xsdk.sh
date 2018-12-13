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
    echo "ERROR: Illegal number of parameters, real type and index size required"
    exit 1
fi
realtype=$1     # required, precision for realtypes
indexsize=$2    # required, integer size for indices
buildthreads=$3 # optional, number of build threads (if empty will use all threads)

# remove old build and install directories
\rm -rf build_xSDK_${realtype}_${indexsize}
\rm -rf install_xSDK_${realtype}_${indexsize}

# create new build and install directories
mkdir build_xSDK_${realtype}_${indexsize}
mkdir install_xSDK_${realtype}_${indexsize}

# move to build directory
cd build_xSDK_${realtype}_${indexsize}

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
if [ "$realtype" == "extended" ] || [ "$indexsize" == "64" ]; then
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

# SuperLU MT index size must be set a build time
if [ "$indexsize" == "32" ]; then
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

# hypre index size must be set a build time
if [ "$indexsize" == "32" ]; then
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

# PETSc index size must be set a build time
if [ "$indexsize" == "32" ]; then
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
#
# Setting the shared linker flags to
# '-D CMAKE_SHARED_LINKER_FLAGS="-Wl,--no-undefined"'
# is useful for finding undefined references when building shared libraries
# -------------------------------------------------------------------------------

# set realtype for xSDK build
if [ $realtype == "extended" ]; then
    xsdk_realtype="quad"
else
    xsdk_realtype=$realtype
fi

# only run development tests with double precision
if [ "$realtype" != "double" ]; then
    DEVTESTS=OFF
else
    DEVTESTS=ON
fi

echo "START CMAKE"
cmake \
    -D CMAKE_INSTALL_PREFIX="../install_xSDK_${realtype}_${indexsize}" \
    \
    -D BUILD_ARKODE=ON \
    -D BUILD_CVODE=ON \
    -D BUILD_CVODES=ON \
    -D BUILD_IDA=ON \
    -D BUILD_IDAS=ON \
    -D BUILD_KINSOL=ON \
    \
    -D XSDK_PRECISION=${xsdk_realtype} \
    -D XSDK_INDEX_SIZE=${indexsize} \
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
    -D MPI_C_COMPILER="${MPIDIR}/mpicc" \
    -D MPI_CXX_COMPILER="${MPIDIR}/mpicxx" \
    -D MPI_Fortran_COMPILER="${MPIDIR}/mpif90" \
    -D MPIEXEC_EXECUTABLE="${MPIDIR}/mpirun" \
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
    -D SUNDIALS_DEVTESTS="${DEVTESTS}" \
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

# smoke test for installation
echo "START TEST_INSTALL"
make test_install 2>&1 | tee test_install.log

# check make install return code
rc=${PIPESTATUS[0]}
echo -e "\nmake test_install returned $rc\n" | tee -a test_install.log
if [ $rc -ne 0 ]; then exit 1; fi

# -------------------------------------------------------------------------------
# Return
# -------------------------------------------------------------------------------

# if we make it here all tests have passed
exit 0
