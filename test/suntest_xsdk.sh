#!/bin/bash
# ------------------------------------------------------------------------------
# Programmer(s): David J. Gardner @ LLNL
# ------------------------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2019, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# ------------------------------------------------------------------------------
# SUNDIALS regression testing using xSDK options
#
# Usage: ./suntest_xsdk.sh <real type> <index size> <TPL status> <test type> \
#                          <build threads>
#
# Required Inputs:
#   <real type>  = SUNDIALS real type: single, double or extended
#   <index size> = SUNDIALS index size: 32 or 64
#   <TPL status> = Enable/disable third party libraries: ON or OFF
#   <test type>  = Test type: STD (standard), DEV (development), or
#                  NONE (compile only)
#
# Optional Inputs:
#   <build threads> = number of threads to use in parallel build (default 1)
# ------------------------------------------------------------------------------

# check number of inputs
if [ "$#" -lt 4 ]; then
    echo "ERROR: FOUR (4) inputs required"
    echo "real type  : [single|double|extended]"
    echo "index size : [32|64]"
    echo "TPLs       : [ON|OFF]"
    echo "test type  : [STD|DEV|NONE]"
    exit 1
fi

realtype=$1     # precision for realtypes
indexsize=$2    # integer size for indices
tplstatus=$3    # enable/disable third party libraries
testtype=$4     # run standard tests, dev tests, or no tests (compile only)
buildthreads=1  # default number threads for parallel builds

# check if the number of build threads was set
if [ "$#" -gt 4 ]; then
    buildthreads=$5
fi

# ------------------------------------------------------------------------------
# Check inputs
# ------------------------------------------------------------------------------

# set real types to test
case "$realtype" in
    single|double|extended) ;;
    *)
        echo "ERROR: Unknown real type option: $realtype"
        exit 1
        ;;
esac

# set index sizes to test
case "$indexsize" in
    32|64) ;;
    *)
        echo "ERROR: Unknown index size option: $indexsize"
        exit 1
        ;;
esac

# set TPL status
case "$tplstatus" in
    ON|On|on)
        TPLs=ON
        ;;
    OFF|Off|off)
        TPLs=OFF
        ;;
    *)
        echo "ERROR: Unknown third party library status: $tplstatus"
        exit 1
        ;;
esac

# which tests to run (if any)
case "$testtype" in
    STD|std|Std)
        # only run standard tests
        devtests=OFF
        skiptests=OFF
        ;;
    DEV|dev|Dev)
        # only run development tests (only double precision supported)
        if [ "$realtype" != "double" ]; then
            echo -e "\nWARNING: Development tests only support realtype = double\n"
            devtests=OFF
            skiptests=OFF
        else
            devtests=ON
            skiptests=OFF
        fi
        ;;
    NONE|none|None)
        # only compile sundials, do not test or install
        devtests=OFF
        skiptests=ON
        ;;
    *)
        echo "ERROR: Unknown test option: $testtype"
        exit 1
        ;;
esac

# ------------------------------------------------------------------------------
# Setup test directories
# ------------------------------------------------------------------------------

# build and install directories
if [ "$TPLs" == "ON" ]; then
    builddir=build_xsdk_${realtype}_${indexsize}_tpls
    installdir=install_xsdk_${realtype}_${indexsize}_tpls
else
    builddir=build_xsdk_${realtype}_${indexsize}
    installdir=install_xsdk_${realtype}_${indexsize}
fi

# remove old build and install directories
\rm -rf $builddir
\rm -rf $installdir

# create and move to new build directory
mkdir $builddir
cd $builddir

# ------------------------------------------------------------------------------
# Installed Third Party Libraries
# ------------------------------------------------------------------------------

if [ "$TPLs" == "ON" ]; then

    # C and C++ standard flags to append
    CSTD="-std=c99"
    CXXSTD="-std=c++11"

    # Enable MPI
    MPISTATUS=ON

    # LAPACK/BLAS: Do not currently support extended precision or 64-bit indices
    if [ "$realtype" == "extended" ] || [ "$indexsize" == "64" ]; then
        LAPACKSTATUS=OFF
        BLASSTATUS=OFF
    else
        BLASSTATUS=ON
        LAPACKSTATUS=ON
    fi

    # KLU: Does not support single or extended precision
    if [ "$realtype" == "single" ] || [ "$realtype" == "extended" ]; then
        KLUSTATUS=OFF
    else
        KLUSTATUS=ON
    fi

    # SuperLU_MT: Does not support extended precision
    if [ "$realtype" == "extended" ]; then
        SLUMTSTATUS=OFF
    else
        SLUMTSTATUS=ON
        # SuperLU_MT index size must be set at build time
        if [ "$indexsize" == "32" ]; then
            SLUMTDIR=$SLUMTDIR_32
        else
            SLUMTDIR=$SLUMTDIR_64
        fi
    fi

    # SuperLU_DIST: Only supports double precision
    if [ "$realtype" != "double" ]; then
        SLUDISTSTATUS=OFF
    else
        SLUDISTSTATUS=ON
        # SuperLU DIST index size must be set at build time
        if [ "$indexsize" == "32" ]; then
            SLUDISTDIR=$SLUDISTDIR_32
            SLUDISTLIBS=$SLUDISTLIBS_32
        else
            SLUDISTDIR=$SLUDISTDIR_64
            SLUDISTLIBS=$SLUDISTLIBS_64
        fi
    fi

    # hypre: Only testing hypre with double precision at this time
    if [ "$realtype" != "double" ]; then
        HYPRESTATUS=OFF
    else
        HYPRESTATUS=ON
        # hypre index size must be set at build time
        if [ "$indexsize" == "32" ]; then
            HYPREDIR=$HYPREDIR_32
        else
            HYPREDIR=$HYPREDIR_64
        fi
    fi

    # PETSc: Only testing PETSc with double precision at this time
    if [ "$realtype" != "double" ]; then
        PETSCSTATUS=OFF
    else
        PETSCSTATUS=ON
        # PETSc index size must be set at build time
        if [ "$indexsize" == "32" ]; then
            PETSCDIR=$PETSCDIR_32
        else
            PETSCDIR=$PETSCDIR_64
        fi
    fi

    # CUDA does not support extended precision
    if [ "$realtype" == "extended" ]; then
        CUDASTATUS=OFF
    else
        CUDASTATUS=ON
    fi

else

    # C and C++ standard flags to append
    CSTD="-std=c90"
    CXXSTD="-std=c++11"

    # disable all TPLs
    MPISTATUS=OFF
    LAPACKSTATUS=OFF
    BLASSTATUS=OFF
    KLUSTATUS=OFF
    SLUMTSTATUS=OFF
    SLUDISTSTATUS=OFF
    HYPRESTATUS=OFF
    PETSCSTATUS=OFF
    CUDASTATUS=OFF

fi

# ------------------------------------------------------------------------------
# Configure SUNDIALS with CMake
# ------------------------------------------------------------------------------

# set realtype for xSDK build
if [ $realtype == "extended" ]; then
    xsdk_realtype="quad"
else
    xsdk_realtype=$realtype
fi

echo "START CMAKE"
cmake \
    -D USE_XSDK_DEFAULTS=ON \
    \
    -D CMAKE_INSTALL_PREFIX="../$installdir" \
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
    -D EXAMPLES_ENABLE_CUDA=${CUDASTATUS} \
    \
    -D CMAKE_C_COMPILER=$CC \
    -D CMAKE_CXX_COMPILER=$CXX \
    -D CMAKE_Fortran_COMPILER=$FC \
    \
    -D CMAKE_C_FLAGS="${CFLAGS} ${CSTD}" \
    -D CMAKE_CXX_FLAGS="${CXXFLAGS} ${CXXSTD}" \
    -D CMAKE_Fortran_FLAGS="${FFLAGS}" \
    -D CUDA_NVCC_FLAGS="--compiler-options;-Wall;--compiler-options;-Werror" \
    -D CUDA_PROPAGATE_HOST_FLAGS=OFF \
    \
    -D OPENMP_ENABLE=ON \
    -D PTHREAD_ENABLE=ON \
    -D XSDK_ENABLE_CUDA=${CUDASTATUS} \
    -D RAJA_ENABLE=OFF \
    \
    -D MPI_ENABLE="${MPISTATUS}" \
    -D MPI_C_COMPILER="${MPIDIR}/bin/mpicc" \
    -D MPI_CXX_COMPILER="${MPIDIR}/bin/mpicxx" \
    -D MPI_Fortran_COMPILER="${MPIDIR}/bin/mpif90" \
    -D MPIEXEC_EXECUTABLE="${MPIEXEC}" \
    \
    -D TPL_ENABLE_BLAS="${BLASSTATUS}" \
    -D TPL_BLAS_LIBRARIES="${BLAS_LIBRARIES}" \
    \
    -D TPL_ENABLE_LAPACK="${LAPACKSTATUS}" \
    -D TPL_LAPACK_LIBRARIES="${LAPACK_LIBRARIES}" \
    \
    -D TPL_ENABLE_KLU="${KLUSTATUS}" \
    -D TPL_KLU_INCLUDE_DIRS="${KLUDIR}/include" \
    -D TPL_KLU_LIBRARIES="${KLUDIR}/lib/libklu.so" \
    \
    -D TPL_ENABLE_HYPRE="${HYPRESTATUS}" \
    -D TPL_HYPRE_INCLUDE_DIRS="${HYPREDIR}/include" \
    -D TPL_HYPRE_LIBRARIES="${HYPREDIR}/lib/libHYPRE.so" \
    \
    -D TPL_ENABLE_PETSC="${PETSCSTATUS}" \
    -D TPL_PETSC_INCLUDE_DIRS="${PETSCDIR}/include" \
    -D TPL_PETSC_LIBRARIES="${PETSCDIR}/lib/libpetsc.so" \
    \
    -D TPL_ENABLE_SUPERLUMT="${SLUMTSTATUS}" \
    -D TPL_SUPERLUMT_INCLUDE_DIRS="${SLUMTDIR}/SRC" \
    -D TPL_SUPERLUMT_LIBRARIES="${SLUMTDIR}/lib/libsuperlu_mt_PTHREAD.a" \
    -D TPL_SUPERLUMT_THREAD_TYPE=Pthread \
    \
    -D TPL_ENABLE_SUPERLUDIST="${SLUDISTSTATUS}" \
    -D TPL_SUPERLUDIST_INCLUDE_DIRS="${SLUDISTDIR}/include" \
    -D TPL_SUPERLUDIST_LIBRARIES="${SLUDISTLIBS}" \
    -D TPL_SUPERLUDIST_OpenMP=ON \
    -D SKIP_OPENMP_DEVICE_CHECK=ON \
    \
    -D SUNDIALS_DEVTESTS="${devtests}" \
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

# check if tests should be skipped (compile check only)
if [ "$skiptests" = "ON" ]; then exit 0; fi

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
