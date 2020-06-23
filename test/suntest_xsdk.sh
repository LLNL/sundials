#!/bin/bash
# ------------------------------------------------------------------------------
# Programmer(s): David J. Gardner @ LLNL
# ------------------------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2020, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# ------------------------------------------------------------------------------
# SUNDIALS regression testing script using xSDK options
#
# Usage: ./suntest_xsdk.sh <real type> <index size> <library type> <TPL status>
#                          <test type> <build threads> <compiler spec>
#                          <build type>
#
# Required Inputs:
#   <real type>  = SUNDIALS real type to build/test with:
#                    single   : single (32-bit) precision
#                    double   : double (64-bit) precision
#                    extended : extended (80-bit) precision
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
#   <compiler spec> = compiler spec (compiler-name@compiler-version)
#   <build type>    = debug (dbg) or optimized (opt)
# ------------------------------------------------------------------------------

# check number of inputs
if [ "$#" -lt 5 ]; then
    echo "ERROR: FIVE (5) inputs required"
    echo "  real type     : [single|double|extended]"
    echo "  index size    : [32|64]"
    echo "  library type  : [static|shared|both]"
    echo "  TPLs          : [ON|OFF]"
    echo "  test type     : [STD|DEV|NONE]"
    echo "Optional inputs"
    echo "  build threads : [number of build threads e.g., 8]"
    echo "  compiler spec : [compiler spec e.g., gcc@4.9.4]"
    echo "  build type    : [dbg|opt]"
    exit 1
fi

realtype=$1     # precision for realtypes
indexsize=$2    # integer size for indices
libtype=$3      # library type to build
tplstatus=$4    # enable/disable third party libraries
testtype=$5     # run standard tests, dev tests, or no tests (compile only)

# set defaults for optional inputs
buildthreads=1  # number threads for parallel builds
compiler=""     # compiler spec
bldtype=""      # build type

# set optional inputs if provided
if [ "$#" -gt 5 ]; then
    buildthreads=$6
fi

if [ "$#" -gt 6 ]; then
    compiler=$7
fi

if [ "$#" -gt 7 ]; then
    bldtype=$8
fi

# ------------------------------------------------------------------------------
# Check inputs
# ------------------------------------------------------------------------------

# build and install directory names
builddir=build_xsdk
installdir=install_xsdk

# add compiler spec to directory names
if [ "$compiler" != "" ]; then
    # replace @ with -
    compilername=${compiler/@/-}
    builddir=${builddir}_${compilername}
    installdir=${installdir}_${compilername}
fi

# add build type to directory names
if [ "$bldtype" != "" ]; then
    builddir=${builddir}_${bldtype}
    installdir=${installdir}_${bldtype}
fi

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
# ------------------------------------------------------------------------------

if [ -n "$SUNDIALS_ENV" ]; then
    echo "Setting up environment with $SUNDIALS_ENV"
    time source $SUNDIALS_ENV $realtype $indexsize $compiler $bldtype
elif [ -f env.sh ]; then
    echo "Setting up environment with ./env.sh"
    time source env.sh $realtype $indexsize $compiler $bldtype
elif [ -f ~/.sundials_config/env.sh ]; then
    echo "Setting up environment with ~/.sundials_config/env.sh"
    time source ~/.sundials_config/env.sh $realtype $indexsize $compiler $bldtype
elif [ -f env/env.${HOSTNAME}.sh ]; then
    echo "Setting up environment with env/env.${HOSTNAME}.sh"
    time source env/env.${HOSTNAME}.sh $realtype $indexsize $compiler $bldtype
elif [ -f env/env.${HOST}.sh ]; then
    echo "Setting up environment with env/env.${HOST}.sh"
    time source env/env.${HOST}.sh $realtype $indexsize $compiler $bldtype
elif [ -f env/env.default.sh ]; then
    echo "Setting up environment with ./env.default.sh"
    time source env/env.default.sh $realtype $indexsize $compiler $bldtype
else
    echo "WARNING: No environment setup script found"
fi

# check return value
if [ $? -ne 0 ]; then
    echo "environment setup failed"
    exit 1;
fi

# ------------------------------------------------------------------------------
# SUNDIALS test settings
# ------------------------------------------------------------------------------

# Check that only one of SUNDIALS_TEST_OUTPUT_DIR and SUNDIALS_TEST_ANSWER_DIR
# are set to ensure tests do not pass erronously

if [ -n "${SUNDIALS_TEST_OUTPUT_DIR}" ] && [ -n "${SUNDIALS_TEST_ANSWER_DIR}" ]
then
    echo "ERROR: Both SUNDIALS_TEST_OUTPUT_DIR and SUNDIALS_TEST_ANSWER_DIR are set"
    echo "SUNDIALS_TEST_OUTPUT_DIR = ${SUNDIALS_TEST_OUTPUT_DIR}"
    echo "SUNDIALS_TEST_ANSWER_DIR = ${SUNDIALS_TEST_ANSWER_DIR}"
    exit 1
fi

# ------------------------------------------------------------------------------
# Check third party library settings
# ------------------------------------------------------------------------------

if [ "$TPLs" == "ON" ]; then

    # C and C++ standard flags to append
    CSTD="-std=c99"
    CXXSTD="-std=c++11"
    C90MATH=OFF

    # CUDA
    CUDA_STATUS=${CUDA_STATUS:-"OFF"}
    if [ "$CUDA_STATUS" == "ON" ] && [ -z "$CUDA_ARCH" ]; then
        export CUDA_ARCH=sm_30
    fi

    # MPI
    MPI_STATUS=${MPI_STATUS:-"OFF"}
    if [ "$MPI_STATUS" == "ON" ] && [ -z "$MPICC" ]; then
        echo "ERROR: MPI_STATUS = ON but MPICC is not set"
        exit 1
    fi

    # LAPACK
    LAPACK_STATUS=${LAPACK_STATUS:-"OFF"}
    if [ "$LAPACK_STATUS" == "ON" ] && [ -z "$LAPACKLIBS" ]; then
        echo "ERROR: LAPACK_STATUS = ON but LAPACKLIBS is not set"
        exit 1
    fi

    # KLU
    KLU_STATUS=${KLU_STATUS:-"OFF"}
    if [ "$KLU_STATUS" == "ON" ] && [ -z "$KLUDIR" ]; then
        echo "ERROR: KLU_STATUS = ON but KLUDIR is not set"
        exit 1
    fi

    # SuperLU_MT
    SLUMT_STATUS=${SLUMT_STATUS:-"OFF"}
    if [ "$SLUMT_STATUS" == "ON" ] && [ -z "$SLUMTDIR" ]; then
        echo "ERROR: SLUMT_STATUS = ON but SLUMTDIR is not set"
        exit 1
    fi

    # SuperLU_DIST
    SLUDIST_STATUS=${SLUDIST_STATUS:-"OFF"}
    if [ "$SLUDIST_STATUS" == "ON" ] && [ -z "$SLUDISTDIR" ]; then
        echo "ERROR: SLUDIST_STATUS = ON but SLUDISTDIR is not set"
        exit 1
    fi

    # hypre
    HYPRE_STATUS=${HYPRE_STATUS:-"OFF"}
    if [ "$HYPRE_STATUS" == "ON" ] && [ -z "$HYPREDIR" ]; then
        echo "ERROR: HYPRE_STATUS = ON but HYPREDIR is not set"
        exit 1
    fi

    # PETSc
    PETSC_STATUS=${PETSC_STATUS:-"OFF"}
    if [ "$PETSC_STATUS" == "ON" ] && [ -z "$PETSCDIR" ]; then
        echo "ERROR: PETSC_STATUS = ON but PETSCDIR is not set"
        exit 1
    fi

    # Trilinos
    TRILINOS_STATUS=${TRILINOS_STATUS:-"OFF"}
    if [ "$TRILINOS_STATUS" == "ON" ] && [ -z "$TRILINOSDIR" ]; then
        echo "ERROR: TRILINOS_STATUS = ON but TRILINOSDIR is not set"
        exit 1
    fi

    # RAJA
    RAJA_STATUS=${RAJA_STATUS:-"OFF"}
    if [ "$RAJA_STATUS" == "ON" ] && [ -z "$RAJADIR" ]; then
        echo "ERROR: RAJA_STATUS = ON but RAJADIR is not set"
        exit 1
    fi

else

    # C and C++ standard flags to append
    if [ "$realtype" != "double" ]; then
        CSTD="-std=c99"
        C90MATH=OFF
    else
        CSTD="-std=c90"
        C90MATH=ON
    fi
    CXXSTD="-std=c++11"

    # disable all TPLs
    MPI_STATUS=OFF
    LAPACK_STATUS=OFF
    KLU_STATUS=OFF
    SLUMT_STATUS=OFF
    SLUDIST_STATUS=OFF
    HYPRE_STATUS=OFF
    PETSC_STATUS=OFF
    CUDA_STATUS=OFF
    TRILINOS_STATUS=OFF
    RAJA_STATUS=OFF

fi

# Ensure OpenMP and PThread options are set (default to OFF)
OPENMP_STATUS=${OPENMP_STATUS:-"OFF"}
OPENMPDEV_STATUS=${OPENMPDEV_STATUS:-"OFF"}
PTHREAD_STATUS=${PTHREAD_STATUS:-"OFF"}

# Ensure Fortran interface options are set (default to OFF)
F77_STATUS=${F77_STATUS:-"OFF"}
F03_STATUS=${F03_STATUS:-"OFF"}

# Ensure SUNDIALS package options are set (default to ON)
ARKODE_STATUS=${ARKODE_STATUS:-"ON"}
CVODE_STATUS=${CVODE_STATUS:-"ON"}
CVODES_STATUS=${CVODES_STATUS:-"ON"}
IDA_STATUS=${IDA_STATUS:-"ON"}
IDAS_STATUS=${IDAS_STATUS:-"ON"}
KINSOL_STATUS=${KINSOL_STATUS:-"ON"}

# Ensure monitoring is set (default to OFF)
MONITOR_STATUS=${MONITOR_STATUS:-"OFF"}

# Ensure fused kernel status is set (default is OFF)
FUSED_STATUS=${FUSED_STATUS:-"OFF"}

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

# set Fortran status for xSDK build
if [[ "${F77_STATUS}" == "ON" || "${F03_STATUS}" == "ON" ]]; then
    FORTRAN_STATUS=ON
else
    FORTRAN_STATUS=OFF
fi

echo "START CMAKE"
time cmake \
    -D USE_XSDK_DEFAULTS=ON \
    \
    -D CMAKE_INSTALL_PREFIX="../$installdir" \
    \
    -D BUILD_STATIC_LIBS="${STATIC}" \
    -D BUILD_SHARED_LIBS="${SHARED}" \
    \
    -D BUILD_ARKODE="${ARKODE_STATUS}" \
    -D BUILD_CVODE="${CVODE_STATUS}" \
    -D BUILD_CVODES="${CVODES_STATUS}" \
    -D BUILD_IDA="${IDA_STATUS}" \
    -D BUILD_IDAS="${IDAS_STATUS}" \
    -D BUILD_KINSOL="${KINSOL_STATUS}" \
    \
    -D SUNDIALS_BUILD_WITH_MONITORING=${MONITOR_STATUS} \
    \
    -D XSDK_PRECISION=${xsdk_realtype} \
    -D XSDK_INDEX_SIZE=${indexsize} \
    \
    -D XSDK_ENABLE_FORTRAN="${FORTRAN_STATUS}" \
    \
    -D EXAMPLES_ENABLE_C=ON \
    -D EXAMPLES_ENABLE_CXX=ON \
    -D EXAMPLES_ENABLE_F77="${FORTRAN_STATUS}" \
    -D EXAMPLES_ENABLE_F90="${FORTRAN_STATUS}" \
    -D EXAMPLES_ENABLE_CUDA="${CUDA_STATUS}" \
    \
    -D CMAKE_C_COMPILER=$CC \
    -D CMAKE_CXX_COMPILER=$CXX \
    -D CMAKE_Fortran_COMPILER=$FC \
    \
    -D CMAKE_C_FLAGS="${CFLAGS} ${CSTD}" \
    -D CMAKE_CXX_FLAGS="${CXXFLAGS} ${CXXSTD}" \
    -D CMAKE_CUDA_FLAGS="${CUDAFLAGS}" \
    -D CMAKE_Fortran_FLAGS="${FFLAGS}" \
    \
    -D OPENMP_ENABLE="${OPENMP_STATUS}" \
    -D PTHREAD_ENABLE="${PTHREAD_STATUS}" \
    \
    -D XSDK_ENABLE_CUDA="${CUDA_STATUS}" \
    -D CUDA_ARCH="${CUDA_ARCH}" \
    \
    -D OPENMP_DEVICE_ENABLE="${OPENMPDEV_STATUS}" \
    -D SKIP_OPENMP_DEVICE_CHECK=TRUE \
    \
    -D MPI_ENABLE="${MPI_STATUS}" \
    -D MPI_C_COMPILER="${MPICC}" \
    -D MPI_CXX_COMPILER="${MPICXX}" \
    -D MPI_Fortran_COMPILER="${MPIFC}" \
    -D MPIEXEC_EXECUTABLE="${MPIEXEC}" \
    \
    -D TPL_ENABLE_LAPACK="${LAPACK_STATUS}" \
    -D TPL_LAPACK_LIBRARIES="${LAPACKLIBS}" \
    \
    -D TPL_ENABLE_KLU="${KLU_STATUS}" \
    -D TPL_KLU_INCLUDE_DIRS="${KLUDIR}/include" \
    -D TPL_KLU_LIBRARIES="${KLUDIR}/lib/libklu.so" \
    \
    -D TPL_ENABLE_HYPRE="${HYPRE_STATUS}" \
    -D TPL_HYPRE_INCLUDE_DIRS="${HYPREDIR}/include" \
    -D TPL_HYPRE_LIBRARIES="${HYPREDIR}/lib/libHYPRE.so" \
    \
    -D TPL_ENABLE_PETSC="${PETSC_STATUS}" \
    -D TPL_PETSC_DIR="${PETSCDIR}" \
    \
    -D TPL_ENABLE_SUPERLUMT="${SLUMT_STATUS}" \
    -D TPL_SUPERLUMT_INCLUDE_DIRS="${SLUMTDIR}/include" \
    -D TPL_SUPERLUMT_LIBRARIES="${SLUMTLIBS};${SLUMTDIR}/lib/libsuperlu_mt_PTHREAD.a" \
    -D TPL_SUPERLUMT_THREAD_TYPE="${SLUMTTYPE}" \
    \
    -D TPL_ENABLE_SUPERLUDIST="${SLUDIST_STATUS}" \
    -D TPL_SUPERLUDIST_INCLUDE_DIRS="${SLUDISTDIR}/include" \
    -D TPL_SUPERLUDIST_LIBRARIES="${SLUDISTLIBS}" \
    -D TPL_SUPERLUDIST_OpenMP=ON \
    -D SKIP_OPENMP_DEVICE_CHECK=ON \
    \
    -D TPL_ENABLE_TRILINOS="${TRILINOS_STATUS}" \
    -D Trilinos_DIR="${TRILINOSDIR}" \
    \
    -D TPL_ENABLE_RAJA="${RAJA_STATUS}" \
    -D RAJA_DIR="${RAJADIR}" \
    \
    -D USE_GENERIC_MATH="${C90MATH}" \
    \
    -D SUNDIALS_BUILD_PACKAGE_FUSED_KERNELS="${FUSED_STATUS}" \
    \
    -D SUNDIALS_TEST_DEVTESTS="${devtests}" \
    -D SUNDIALS_TEST_UNITTESTS=ON \
    -D SUNDIALS_TEST_OUTPUT_DIR="${SUNDIALS_TEST_OUTPUT_DIR}" \
    -D SUNDIALS_TEST_ANSWER_DIR="${SUNDIALS_TEST_ANSWER_DIR}" \
    -D SUNDIALS_TEST_FLOAT_PRECISION="${SUNDIALS_TEST_FLOAT_PRECISION}" \
    -D SUNDIALS_TEST_INTEGER_PRECISION="${SUNDIALS_TEST_INTEGER_PRECISION}" \
    \
    -D CMAKE_VERBOSE_MAKEFILE=OFF \
    \
    ../../. 2>&1 | tee configure.log

# check cmake return code
rc=${PIPESTATUS[0]}
echo -e "\ncmake returned $rc\n" | tee -a configure.log
if [ $rc -ne 0 ]; then cd ..; exit 1; fi

# -------------------------------------------------------------------------------
# Make SUNDIALS
# -------------------------------------------------------------------------------

echo "START MAKE"
time make -j $buildthreads 2>&1 | tee make.log

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
time ctest -j $buildthreads test 2>&1 | tee test.log

# check make test return code
rc=${PIPESTATUS[0]}
echo -e "\nmake test returned $rc\n" | tee -a test.log
if [ $rc -ne 0 ]; then cd ..; exit 1; fi

# -------------------------------------------------------------------------------
# Install SUNDIALS
# -------------------------------------------------------------------------------

# install sundials
echo "START INSTALL"
time make -j $buildthreads install 2>&1 | tee install.log

# check make install return code
rc=${PIPESTATUS[0]}
echo -e "\nmake install returned $rc\n" | tee -a install.log
if [ $rc -ne 0 ]; then cd ..; exit 1; fi

# -------------------------------------------------------------------------------
# Test SUNDIALS Install
# -------------------------------------------------------------------------------

# smoke test for installation
echo "START TEST_INSTALL"
time make test_install 2>&1 | tee test_install.log

# check make install return code
rc=${PIPESTATUS[0]}
echo -e "\nmake test_install returned $rc\n" | tee -a test_install.log
if [ $rc -ne 0 ]; then cd ..; exit 1; fi

# -------------------------------------------------------------------------------
# Test SUNDIALS Install All
# -------------------------------------------------------------------------------

# smoke test for installation
echo "START TEST_INSTALL_ALL"
time make test_install_all 2>&1 | tee test_install_all.log

# check make install all return code
rc=${PIPESTATUS[0]}
echo -e "\nmake test_install_all returned $rc\n" | tee -a test_install_all.log
if [ $rc -ne 0 ]; then cd ..; exit 1; fi

# -------------------------------------------------------------------------------
# Return
# -------------------------------------------------------------------------------

# if we make it here all tests have passed
cd ..
exit 0
