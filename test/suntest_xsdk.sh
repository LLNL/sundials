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
# Usage: ./suntest_xsdk.sh <real type> <index size> <library type> <TPL status>
#                          <test type> <build threads>
#
# Required Inputs:
#   <real type>  = SUNDIALS real type to build/test with:
#                    single   : single (32-bit) precision
#                    double   : double (64-bit) precision
#                    extended : extended (128-bit) precision
#   <index size> = SUNDIALS index size to build/test with:
#                    32       : 32-bit indices
#                    64       : 64-bit indices
#   <lib type>   = Which library type to test:
#                    static   : only build static libraries
#                    shared   : only build shared libraries
#                    both     : build static and shared simultaneously
#   <TPL status> = Enable/disable third party libraries:
#                    ON       : All possible TPLs enabled
#                    OFF      : No TPLs enabled
#   <test type>  = Test type to run:
#                    STD      : standard tests
#                    DEV      : development tests
#                    NONE     : no test, configure and compile only
#
# Optional Inputs:
#   <build threads> = number of threads to use in parallel build (default 1)
# ------------------------------------------------------------------------------

# check number of inputs
if [ "$#" -lt 5 ]; then
    echo "ERROR: FIVE (5) inputs required"
    echo "real type    : [single|double|extended]"
    echo "index size   : [32|64]"
    echo "library type : [static|shared|both]"
    echo "TPLs         : [ON|OFF]"
    echo "test type    : [STD|DEV|NONE]"
    exit 1
fi

realtype=$1     # precision for realtypes
indexsize=$2    # integer size for indices
libtype=$3      # library type to build
tplstatus=$4    # enable/disable third party libraries
testtype=$5     # run standard tests, dev tests, or no tests (compile only)
buildthreads=1  # default number threads for parallel builds

# check if the number of build threads was set
if [ "$#" -gt 5 ]; then
    buildthreads=$6
fi

# ------------------------------------------------------------------------------
# Check inputs
# ------------------------------------------------------------------------------

# build and install directory names
builddir=build_xsdk
installdir=install_xsdk

# set real types to test
case "$realtype" in
    single|double|extended) ;;
    *)
        echo "ERROR: Unknown real type option: $realtype"
        exit 1
        ;;
esac
builddir=${builddir}_${realtype}
installdir=${installdir}_${realtype}

# set index sizes to test
case "$indexsize" in
    32|64) ;;
    *)
        echo "ERROR: Unknown index size option: $indexsize"
        exit 1
        ;;
esac
builddir=${builddir}_${indexsize}
installdir=${installdir}_${indexsize}

# set library types
case "$libtype" in
    STATIC|Static|static)
        STATIC=ON
        SHARED=OFF
        ;;
    SHARED|Shared|shared)
        STATIC=OFF
        SHARED=ON
        ;;
    BOTH|Both|both)
        STATIC=ON
        SHARED=ON
        ;;
    *)
        echo "ERROR: Unknown library type: $libtype"
        exit 1
        ;;
esac
builddir=${builddir}_${libtype}
installdir=${installdir}_${libtype}


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
builddir=${builddir}_${tplstatus}
installdir=${installdir}_${tplstatus}

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
builddir=${builddir}_${testtype}
installdir=${installdir}_${testtype}

# ------------------------------------------------------------------------------
# Setup the test environment
# 1. User defined environment script
# 2. User's local environment script
# 3. User's global environment script
# 4. Sundials default environment script
# ------------------------------------------------------------------------------

if [ ! -z "$SUNDIALS_ENV" ]; then
    echo "Setting up environment with $SUNDIALS_ENV"
    source $SUNDIALS_ENV $realtype $indexsize
elif [ -f env.sh ]; then
    echo "Setting up environment with ./env.sh"
    source env.sh $realtype $indexsize
elif [ -f ~/.sundials_config/env.sh ]; then
    echo "Setting up environment with ~/.sundials_config/env.sh"
    source ~/.sundials_config/env.sh $realtype $indexsize
else
    echo "Setting up environment with ./env.default.sh"
    source env.default.sh $realtype $indexsize
fi

# check return value
if [ $? -ne 0 ]; then
    echo "environment setup failed"
    exit 1;
fi

# ------------------------------------------------------------------------------
# Check third party library settings
# ------------------------------------------------------------------------------

if [ "$TPLs" == "ON" ]; then

    # C and C++ standard flags to append
    CSTD="-std=c99"
    CXXSTD="-std=c++11"

    # CUDA
    CUDASTATUS=${CUDASTATUS:-"OFF"}

    # MPI
    MPISTATUS=${MPISTATUS:-"OFF"}
    if [ "$MPISTATUS" == "ON" ] && [ -z "$MPICC" ]; then
        echo "ERROR: MPISTATUS = ON but MPICC is not set"
        exit 1
    fi

    # BLAS
    BLASSTATUS=${BLASSTATUS:-"OFF"}
    if [ "$BLASSTATUS" == "ON" ] && [ -z "$BLASLIBS" ]; then
        echo "ERROR: BLASSTATUS = ON but BLASLIBS is not set"
        exit 1
    fi

    # LAPACK
    LAPACKSTATUS=${LAPACKSTATUS:-"OFF"}
    if [ "$LAPACKSTATUS" == "ON" ] && [ -z "$LAPACKLIBS" ]; then
        echo "ERROR: LAPACKSTATUS = ON but LAPACKLIBS is not set"
        exit 1
    fi

    # KLU
    KLUSTATUS=${KLUSTATUS:-"OFF"}
    if [ "$KLUSTATUS" == "ON" ] && [ -z "$KLUDIR" ]; then
        echo "ERROR: KLUSTATUS = ON but KLUDIR is not set"
        exit 1
    fi

    # SuperLU_MT
    SLUMTSTATUS=${SLUMTSTATUS:-"OFF"}
    if [ "$SLUMTSTATUS" == "ON" ] && [ -z "$SLUMTDIR" ]; then
        echo "ERROR: SLUMTSTATUS = ON but SLUMTDIR is not set"
        exit 1
    fi

    # SuperLU_DIST
    SLUDISTSTATUS=${SLUDISTSTATUS:-"OFF"}
    if [ "$SLUDISTSTATUS" == "ON" ] && [ -z "$SLUDISTDIR" ]; then
        echo "ERROR: SLUDISTSTATUS = ON but SLUDISTDIR is not set"
        exit 1
    fi

    # hypre
    HYPRESTATUS=${HYPRESTATUS:-"OFF"}
    if [ "$HYPRESTATUS" == "ON" ] && [ -z "$HYPREDIR" ]; then
        echo "ERROR: HYPRESTATUS = ON but HYPREDIR is not set"
        exit 1
    fi

    # PETSc
    PETSCSTATUS=${PETSCSTATUS:-"OFF"}
    if [ "$PETSCSTATUS" == "ON" ] && [ -z "$PETSCDIR" ]; then
        echo "ERROR: PETSCSTATUS = ON but PETSCDIR is not set"
        exit 1
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
# Setup test directories
# ------------------------------------------------------------------------------

# remove old build and install directories
\rm -rf $builddir
\rm -rf $installdir

# create and move to new build directory
mkdir $builddir
cd $builddir

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
    -D BUILD_STATIC_LIBS="${STATIC}" \
    -D BUILD_SHARED_LIBS="${SHARED}" \
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
    -D MPI_C_COMPILER="${MPICC}" \
    -D MPI_CXX_COMPILER="${MPICXX}" \
    -D MPI_Fortran_COMPILER="${MPIFC}" \
    -D MPIEXEC_EXECUTABLE="${MPIEXEC}" \
    \
    -D TPL_ENABLE_BLAS="${BLASSTATUS}" \
    -D TPL_BLAS_LIBRARIES="${BLASLIBS}" \
    \
    -D TPL_ENABLE_LAPACK="${LAPACKSTATUS}" \
    -D TPL_LAPACK_LIBRARIES="${LAPACKLIBS}" \
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
if [ $rc -ne 0 ]; then cd ..; exit 1; fi

# -------------------------------------------------------------------------------
# Make SUNDIALS
# -------------------------------------------------------------------------------

echo "START MAKE"
make -j $buildthreads 2>&1 | tee make.log

# check make return code
rc=${PIPESTATUS[0]}
echo -e "\nmake returned $rc\n" | tee -a make.log
if [ $rc -ne 0 ]; then cd ..; exit 1; fi

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
if [ $rc -ne 0 ]; then cd ..; exit 1; fi

# -------------------------------------------------------------------------------
# Install SUNDIALS
# -------------------------------------------------------------------------------

# install sundials
echo "START INSTALL"
make -j $buildthread install 2>&1 | tee install.log

# check make install return code
rc=${PIPESTATUS[0]}
echo -e "\nmake install returned $rc\n" | tee -a install.log
if [ $rc -ne 0 ]; then cd ..; exit 1; fi

# -------------------------------------------------------------------------------
# Test SUNDIALS Install
# -------------------------------------------------------------------------------

# smoke test for installation
echo "START TEST_INSTALL"
make test_install 2>&1 | tee test_install.log

# check make install return code
rc=${PIPESTATUS[0]}
echo -e "\nmake test_install returned $rc\n" | tee -a test_install.log
if [ $rc -ne 0 ]; then cd ..; exit 1; fi

# -------------------------------------------------------------------------------
# Return
# -------------------------------------------------------------------------------

# if we make it here all tests have passed
cd ..
exit 0
