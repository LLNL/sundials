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
# SUNDIALS tarball regression testing script
#
# Usage: ./suntest_tarscript.sh <package> <lib type> <real type> <index size> \
#                               <TPL status> <test type> <build threads>
#
# Required Inputs:
#   <package>    = Which tarball to make and test:
#                    arkode   : create ARKode tarball only
#                    cvode    : create CVODE tarball only
#                    cvodes   : create CVODES tarball only
#                    ida      : create IDA tarball only
#                    idas     : create IDAS tarball only
#                    kinsol   : create KINSOL tarball only
#                    sundials : create sundials tarball containing all packages
#                    all      : all of the above options
#   <lib type>   = Which library type to test:
#                    static : only build static libraries
#                    shared : only build shared libraries
#                    each   : build static and shared separately
#                    both   : build static and shared simultaneously
#   <real type>  = SUNDIALS real type to build/test with:
#                    single   : single precision only
#                    double   : double precision only
#                    extended : extended (quad) precision only
#                    all      : all of the above options
#   <index size> = SUNDIALS index size to build/test with:
#                    32   : 32-bit indices only
#                    64   : 64-bit indices only
#                    both : both of the above options
#   <TPL status> = Enable/disable third party libraries:
#                    ON  : All possible TPLs enabled
#                    OFF : No TPLs enabled
#   <test type>  = Test type to run:
#                    STD  : standard tests
#                    DEV  : development tests
#                    NONE : no test, configure and compile only
#
# Optional Inputs:
#   <build threads> = number of threads to use in parallel build (default 1)
# ------------------------------------------------------------------------------

# check number of inputs
if [ "$#" -lt 6 ]; then
    echo "ERROR: SIX (6) inputs required"
    echo "package      : [arkode|cvode|cvodes|ida|idas|kinsol|sundials|all]"
    echo "library type : [static|shared|each|both]"
    echo "real type    : [single|double|extended|all]"
    echo "index size   : [32|64|both]"
    echo "TPLs         : [ON|OFF]"
    echo "test type    : [STD|DEV|NONE]"
    exit 1
fi

package=$1        # sundials package to test
tmplibtype=$2     # library type to build
tmprealtype=$3    # precision for realtypes
tmpindexsize=$4   # integer size for indices
tplstatus=$5      # enable/disable third party libraries
testtype=$6       # run standard tests, dev tests, or no tests (compile only)
buildthreads=1    # default number threads for parallel builds

# check if the number of build threads was set
if [ "$#" -gt 6 ]; then
    buildthreads=$7
fi

# ------------------------------------------------------------------------------
# Check inputs
# ------------------------------------------------------------------------------

# check package option
case $package in
    arkode|cvode|cvodes|ida|idas|kinsol|sundials|all) ;;
    *)
        echo "ERROR: Unknown package option: $package"
        exit 1
        ;;
esac

# set library types to test
case "$tmplibtype" in
    static) libtype=( "static" );;
    shared) libtype=( "shared" );;
    each)   libtype=( "static" "shared") ;;
    both)   libtype=( "both" ) ;;
    *)
        echo "ERROR: Unknown library type option: $tmplibtype"
        exit 1
        ;;
esac

# set real types to test
case "$tmprealtype" in
    single)   realtype=( "single" );;
    double)   realtype=( "double" );;
    extended) realtype=( "extended" );;
    all)      realtype=( "single" "double" "extended" );;
    *)
        echo "ERROR: Unknown real type option: $tmprealtype"
        exit 1
        ;;
esac

# set index sizes to test
case "$tmpindexsize" in
    32)   indexsize=( "32" );;
    64)   indexsize=( "64" );;
    both) indexsize=( "32" "64" );;
    *)
        echo "ERROR: Unknown index size option: $tmpindexsize"
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
        # only run development tests
        devtests=ON
        skiptests=OFF
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
# Create tarballs
# ------------------------------------------------------------------------------

# location of testing directory
testdir=`pwd`

# remove old tarball directory and create new directory
\rm -rf tarballs || exit 1
mkdir tarballs   || exit 1

# run tarscript to create tarballs
cd ../scripts || exit 1

echo "START TARSCRIPT"
./tarscript -s $package | tee -a tar.log

# check tarscript return code
rc=${PIPESTATUS[0]}
echo -e "\ntarscript returned $rc\n" | tee -a tar.log
if [ $rc -ne 0 ]; then
    # remove temporary file created by tarscript and exit with error
    \rm -rf ../../tmp_dir.*
    exit 1;
fi

# relocate tarballs
mv tar.log $testdir/tarballs/.

# move tarballs to tarball directory
case $package in
    arkode)
        mv ../../arkode-*.tar.gz $testdir/tarballs/.   || exit 1
        ;;
    cvode)
        mv ../../cvode-*.tar.gz $testdir/tarballs/.    || exit 1
        ;;
    cvodes)
        mv ../../cvodes-*.tar.gz $testdir/tarballs/.   || exit 1
        ;;
    ida)
        mv ../../ida-*.tar.gz $testdir/tarballs/.      || exit 1
        ;;
    idas)
        mv ../../idas-*.tar.gz $testdir/tarballs/.     || exit 1
        ;;
    kinsol)
        mv ../../kinsol-*.tar.gz $testdir/tarballs/.   || exit 1
        ;;
    sundials)
        mv ../../sundials-*.tar.gz $testdir/tarballs/. || exit 1
        ;;
    all)
        mv ../../sundials-*.tar.gz $testdir/tarballs/. || exit 1
        mv ../../arkode-*.tar.gz $testdir/tarballs/.   || exit 1
        mv ../../cvode-*.tar.gz $testdir/tarballs/.    || exit 1
        mv ../../cvodes-*.tar.gz $testdir/tarballs/.   || exit 1
        mv ../../ida-*.tar.gz $testdir/tarballs/.      || exit 1
        mv ../../idas-*.tar.gz $testdir/tarballs/.     || exit 1
        mv ../../kinsol-*.tar.gz $testdir/tarballs/.   || exit 1
        ;;
esac

# ------------------------------------------------------------------------------
# Test tarballs
# ------------------------------------------------------------------------------

# move to tarball directory
cd $testdir/tarballs/. || exit 1

# loop over tarballs and test each one
for tarball in *.tar.gz; do

    # get package name
    package=${tarball%.tar.gz}

    # --------------------------------------------------------------------------
    # Uncompress tarball and setup build
    # --------------------------------------------------------------------------

    echo "START UNTAR"
    tar -xvzf $tarball 2>&1 | tee -a tar.log

    # check tar return code
    rc=${PIPESTATUS[0]}
    echo -e "\ntar -xzvf returned $rc\n" | tee -a tar.log
    if [ $rc -ne 0 ]; then exit 1; fi

    # move log to package directory
    mv tar.log $package/. || exit 1

    # move to package directory
    cd $package || exit 1

    # loop over build options
    for lt in "${libtype[@]}"; do
        for rt in "${realtype[@]}"; do
            for is in "${indexsize[@]}"; do

                # print test label for Jenkins section collapsing
                echo "TEST: $lt $rt $is"

                # build and install directories
                if [ "$TPLs" == "ON" ]; then
                    builddir=build_${lt}_${rt}_${is}_tpls
                    installdir=install_${lt}_${rt}_${is}_tpls
                else
                    builddir=build_${lt}_${rt}_${is}
                    installdir=install_${lt}_${rt}_${is}
                fi

                # remove old build and install directories
                \rm -rf $builddir   || exit 1
                \rm -rf $installdir || exit 1

                # create and move to new build directory
                mkdir $builddir || exit 1
                cd $builddir    || exit 1

                # set library type to build
                if [ "${lt}" == "static" ]; then
                    STATIC=ON
                    SHARED=OFF
                elif [ "${lt}" == "shared" ]; then
                    STATIC=OFF
                    SHARED=ON
                else
                    STATIC=ON
                    SHARED=ON
                fi

                # --------------------------------------------------------------
                # Installed Third Party Libraries
                # --------------------------------------------------------------

                if [ "$TPLs" == "ON" ]; then

                    # C and C++ standard flags to append
                    CSTD="-std=c99"
                    CXXSTD="-std=c++11"

                    # Enable MPI
                    MPISTATUS=ON

                    # LAPACK/BLAS: Do currently support extended precision or 64-bit indices
                    if [ "$rt" == "extended" ] || [ "$is" == "64" ]; then
                        LAPACKSTATUS=OFF
                        BLASSTATUS=OFF
                    else
                        BLASSTATUS=ON
                        LAPACKSTATUS=ON
                    fi

                    # KLU: Does not support single or extended precision
                    if [ "$rt" == "single" ] || [ "$rt" == "extended" ]; then
                        KLUSTATUS=OFF
                    else
                        KLUSTATUS=ON
                    fi

                    # SuperLU_MT: Does not support extended precision
                    if [ "$rt" == "extended" ]; then
                        SLUMTSTATUS=OFF
                    else
                        SLUMTSTATUS=ON
                        # SuperLU_MT index size must be set at build time
                        if [ "$is" == "32" ]; then
                            SLUMTDIR=$SLUMTDIR_32
                        else
                            SLUMTDIR=$SLUMTDIR_64
                        fi
                    fi

                    # SuperLU_DIST: Only supports double precision
                    if [ "$rt" != "double" ]; then
                        SLUDISTSTATUS=OFF
                    else
                        SLUDISTSTATUS=ON
                        # SuperLU DIST index size must be set at build time
                        if [ "$is" == "32" ]; then
                            SLUDISTDIR=$SLUDISTDIR_32
                        else
                            SLUDISTDIR=$SLUDISTDIR_64
                        fi
                    fi

                    # hypre: Only testing hypre with double precision at this time
                    if [ "$rt" != "double" ]; then
                        HYPRESTATUS=OFF
                    else
                        HYPRESTATUS=ON
                        # hypre index size must be set at build time
                        if [ "$is" == "32" ]; then
                            HYPREDIR=$HYPREDIR_32
                        else
                            HYPREDIR=$HYPREDIR_64
                        fi
                    fi

                    # PETSc: Only testing PETSc with double precision at this time
                    if [ "$rt" != "double" ]; then
                        PETSCSTATUS=OFF
                    else
                        PETSCSTATUS=ON
                        # PETSc index size must be set at build time
                        if [ "$is" == "32" ]; then
                            PETSCDIR=$PETSCDIR_32
                        else
                            PETSCDIR=$PETSCDIR_64
                        fi
                    fi

                    # CUDA does not support extended precision
                    if [ "$rt" == "extended" ]; then
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

                # -------------------------------------------------------------------
                # Configure SUNDIALS with CMake
                # -------------------------------------------------------------------

                # only run development tests with double precision
                if [ "$rt" != "double" ] && [ "$devtests" == "ON" ]; then
                    echo -e "\nWARNING: Development tests only support realtype = double\n"
                    dt=OFF
                else
                    dt=$devtests
                fi

                echo "START CMAKE"
                cmake \
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
                    -D SUNDIALS_PRECISION=$rt \
                    -D SUNDIALS_INDEX_SIZE=$is \
                    \
                    -D F77_INTERFACE_ENABLE=ON \
                    -D F2003_INTERFACE_ENABLE=ON \
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
                    -D CUDA_ENABLE=${CUDASTATUS} \
                    -D RAJA_ENABLE=OFF \
                    \
                    -D MPI_ENABLE="${MPISTATUS}" \
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
                    -D SUPERLUMT_ENABLE="${SLUMTSTATUS}" \
                    -D SUPERLUMT_INCLUDE_DIR="${SLUMTDIR}/SRC" \
                    -D SUPERLUMT_LIBRARY_DIR="${SLUMTDIR}/lib" \
                    -D SUPERLUMT_THREAD_TYPE=Pthread \
                    \
                    -D SUPERLUDIST_ENABLE="${SLUDISTSTATUS}" \
                    -D SUPERLUDIST_INCLUDE_DIR="${SLUDISTDIR}/include" \
                    -D SUPERLUDIST_LIBRARY_DIR="${SLUDISTDIR}/lib" \
                    -D SUPERLUDIST_LIBRARIES="${BLAS_LIBRARIES}" \
                    -D SUPERLUDIST_OpenMP=ON \
                    -D SKIP_OPENMP_DEVICE_CHECK=ON \
                    \
                    -D SUNDIALS_DEVTESTS=$dt \
                    ../. 2>&1 | tee configure.log

                # check cmake return code
                rc=${PIPESTATUS[0]}
                echo -e "\ncmake returned $rc\n" | tee -a configure.log
                if [ $rc -ne 0 ]; then exit 1; fi

                # -------------------------------------------------------------------
                # Make SUNDIALS
                # -------------------------------------------------------------------

                echo "START MAKE"
                make -j $buildthreads 2>&1 | tee make.log

                # check make return code
                rc=${PIPESTATUS[0]}
                echo -e "\nmake returned $rc\n" | tee -a make.log
                if [ $rc -ne 0 ]; then exit 1; fi

                # check if tests should be skipped (compile check only)
                if [ "$skiptests" = "ON" ]; then cd ..; continue; fi

                # -------------------------------------------------------------------
                # Test SUNDIALS
                # -------------------------------------------------------------------

                # test sundials
                echo "START TEST"
                make test 2>&1 | tee test.log

                # check make test return code
                rc=${PIPESTATUS[0]}
                echo -e "\nmake test returned $rc\n" | tee -a test.log
                if [ $rc -ne 0 ]; then exit 1; fi

                # -------------------------------------------------------------------
                # Install SUNDIALS
                # -------------------------------------------------------------------

                # install sundials
                echo "START INSTALL"
                make install 2>&1 | tee install.log

                # check make install return code
                rc=${PIPESTATUS[0]}
                echo -e "\nmake install returned $rc\n" | tee -a install.log
                if [ $rc -ne 0 ]; then exit 1; fi

                # -------------------------------------------------------------------
                # Test SUNDIALS Install
                # -------------------------------------------------------------------

                # smoke test for installation
                echo "START TEST_INSTALL"
                make test_install 2>&1 | tee test_install.log

                # check make install return code
                rc=${PIPESTATUS[0]}
                echo -e "\nmake test_install returned $rc\n" | tee -a test_install.log
                if [ $rc -ne 0 ]; then exit 1; fi

                # return to package directory
                echo "PASSED"
                cd .. || exit 1

            done
        done
    done

    # return to tarball directory
    cd .. || exit 1

done

# -------------------------------------------------------------------------------
# Return
# -------------------------------------------------------------------------------

# if we make it here all tests have passed
exit 0
