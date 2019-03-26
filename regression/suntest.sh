#!/bin/bash
# -------------------------------------------------------------------------------
# Programmer(s): David J. Gardner @ LLNL
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
# SUNDIALS regression testing script with all external libraries enabled
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
\rm -rf build_${realtype}_${indexsize}
\rm -rf install_${realtype}_${indexsize}

# create new build and install directories
mkdir build_${realtype}_${indexsize}
mkdir install_${realtype}_${indexsize}

# move to build directory
cd build_${realtype}_${indexsize}

# set file permissions (rwxrwxr-x)
umask 002

# set compiler spec for spack
# COMPILER_SPEC is an environment variable - we use it if it is not empty
if [[ ! -z $COMPILER_SPEC ]]; then
  compiler="${COMPILER_SPEC}"
else
  compiler="gcc@4.9.4"
fi

# -------------------------------------------------------------------------------
# Installed Third Party Libraries
# -------------------------------------------------------------------------------

# Directory where TPLs not provided by Spack reside
APPDIR=/usr/casc/sundials/apps/rh6

# MPI
# MPIDIR is an environment variable - we use it if it is not empty
if [[ ! -z $MPIEXEC ]]; then
  MPIDIR="${MPIDIR}"
else
  MPIDIR="$(spack location -i openmpi@3.1.2 % "$compiler")"
fi
# MPIEXEC is an environment variable - we use it if it is not empty
if [[ ! -z $MPIEXEC ]]; then
  MPIEXEC="${MPIEXEC}"
else
  MPIEXEC="${MPIDIR}/bin/mpirun"
fi

# LAPACK / BLAS
BLASSTATUS=ON
LAPACKSTATUS=ON
BLASDIR="$(spack location -i openblas@0.3.5~ilp64 % "$compiler")"
BLAS_LIBRARIES=${BLASDIR}/lib/libopenblas.so
LAPACK_LIBRARIES=${BLAS_LIBRARIES}

# LAPACK/BLAS does not support extended precision or 64-bit indices
if [ "$realtype" == "extended" ] || [ "$indexsize" == "64" ]; then
    LAPACKSTATUS=OFF
    BLASSTATUS=OFF
fi

# KLU
KLUSTATUS=ON
KLUDIR="$(spack location -i suite-sparse@5.3.0 % "$compiler")"

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

# SuperLU_DIST
SUPERLUDISTSTATUS=ON

# SuperLU DIST index size must be set a build time
if [ "$indexsize" == "32" ]; then
  SUPERLUDISTDIR="$(spack location -i superlu-dist@develop~int64+openmp % "$compiler")"
else
  SUPERLUDISTDIR="$(spack location -i superlu-dist@develop+int64+openmp % "$compiler")"
fi

# SuperLU DIST only supports double precision
if [ "$realtype" != "double" ]; then
    SUPERLUDISTSTATUS=OFF
fi

# hypre
HYPRESTATUS=ON

# hypre index size must be set a build time
if [ "$indexsize" == "32" ]; then
  HYPREDIR="$(spack location -i hypre@2.14.0~int64 % "$compiler")"
else
  HYPREDIR="$(spack location -i hypre@2.14.0+int64 % "$compiler")"
fi

# only testing hypre with double precision at this time
if [ "$realtype" != "double" ]; then
    HYPRESTATUS=OFF
fi

# PETSc
PETSCSTATUS=ON

# PETSc index size must be set a build time
if [ "$indexsize" == "32" ]; then
  PETSCDIR="$(spack location -i petsc@3.10.3~int64 % $compiler)"
else
  PETSCDIR="$(spack location -i petsc@3.10.3+int64 % $compiler)"
fi

# only testing PETSc with double precision at this time
if [ "$realtype" != "double" ]; then
    PETSCSTATUS=OFF
fi

# ------------------------------------------------------------------------------
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
# ------------------------------------------------------------------------------

# only run development tests with double precision
if [ "$realtype" != "double" ]; then
    DEVTESTS=OFF
else
    DEVTESTS=ON
fi

echo "START CMAKE"
cmake \
    -D CMAKE_INSTALL_PREFIX="../install_${realtype}_${indexsize}" \
    \
    -D BUILD_ARKODE=ON \
    -D BUILD_CVODE=ON \
    -D BUILD_CVODES=ON \
    -D BUILD_IDA=ON \
    -D BUILD_IDAS=ON \
    -D BUILD_KINSOL=ON \
    \
    -D SUNDIALS_PRECISION=$realtype \
    -D SUNDIALS_INDEX_SIZE=$indexsize \
    \
    -D F77_INTERFACE_ENABLE=ON \
    -D F2003_INTERFACE_ENABLE=ON \
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
    -D CMAKE_C_COMPILER=$CC \
    -D CMAKE_CXX_COMPILER=$CXX \
    -D CMAKE_Fortran_COMPILER=$FC \
    \
    -D CMAKE_C_FLAGS="-g -Wall -std=c99 -pedantic" \
    -D CMAKE_CXX_FLAGS="-g -Wall" \
    -D CMAKE_Fortran_FLAGS="-g -ffpe-summary=none" \
    \
    -D MPI_ENABLE=ON \
    -D MPI_C_COMPILER="${MPIDIR}/bin/mpicc" \
    -D MPI_CXX_COMPILER="${MPIDIR}/bin/mpicxx" \
    -D MPI_Fortran_COMPILER="${MPIDIR}/bin/mpif90" \
    -D MPIEXEC_EXECUTABLE="${MPIEXEC}" \
    \
    -D BLAS_ENABLE="${BLASSTATUS}" \
    -D BLAS_LIBRARIES="${BLAS_LIBRARIES}" \
    \
    -D LAPACK_ENABLE="${LAPACKSTATUS}" \
    -D LAPACK_LIBRARIES="${LAPACK_LIBRARIES}" \
    \
    -D KLU_ENABLE="${KLUSTATUS}" \
    -D KLU_INCLUDE_DIR="${KLUDIR}/include" \
    -D KLU_LIBRARY_DIR="${KLUDIR}/lib" \
    \
    -D HYPRE_ENABLE="${HYPRESTATUS}" \
    -D HYPRE_INCLUDE_DIR="${HYPREDIR}/include" \
    -D HYPRE_LIBRARY_DIR="${HYPREDIR}/lib" \
    \
    -D PETSC_ENABLE="${PETSCSTATUS}" \
    -D PETSC_INCLUDE_DIR="${PETSCDIR}/include" \
    -D PETSC_LIBRARY_DIR="${PETSCDIR}/lib" \
    \
    -D SUPERLUMT_ENABLE="${SUPERLUMTSTATUS}" \
    -D SUPERLUMT_INCLUDE_DIR="${SUPERLUMTDIR}/SRC" \
    -D SUPERLUMT_LIBRARY_DIR="${SUPERLUMTDIR}/lib" \
    -D SUPERLUMT_THREAD_TYPE=Pthread \
    \
    -D SUPERLUDIST_ENABLE="${SUPERLUDISTSTATUS}" \
    -D SUPERLUDIST_INCLUDE_DIR="${SUPERLUDISTDIR}/include" \
    -D SUPERLUDIST_LIBRARY_DIR="${SUPERLUDISTDIR}/lib" \
    -D SUPERLUDIST_LIBRARIES="${BLAS_LIBRARIES}" \
    -D SUPERLUDIST_OpenMP=ON \
    -D SKIP_OPENMP_DEVICE_CHECK=ON \
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
