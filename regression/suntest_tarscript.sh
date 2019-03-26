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
# Testing script that creates a SUNDIALS release tarballs, builds and installs
# SUNDIALS packages from the tarballs, and tests the installed packages with all
# external libraries enabled.
# ------------------------------------------------------------------------------

# check if a module was specified, otherwise default to all modules
if [ "$#" -gt 0 ]; then
    module=$1
else
    module="all"
fi

# common environment setup
source env.sh

# set file permissions (rwxrwxr-x)
umask 002

# set compiler spec for spack
# COMPILER_SPEC is an environment variable - we use it if it is not empty
if [[ ! -z $COMPILER_SPEC ]]; then
  compiler="${COMPILER_SPEC}"
else
  compiler="gcc@4.9.4"
fi

# library types, real types, and index sizes to test
# NOTE: may need to create answer files for different realtypes
libtype=( "static" "shared" )
realtype=( "single" "double" "extended" )
indexsize=( "32" "64" )

# ------------------------------------------------------------------------------
# Installed third party libraries
# ------------------------------------------------------------------------------

# path to installed libraries
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
BLASDIR="$(spack location -i openblas@0.3.5~ilp64 % "$compiler")"
BLAS_LIBRARIES=${BLASDIR}/lib/libopenblas.so
LAPACK_LIBRARIES=${BLAS_LIBRARIES}

# KLU
KLUDIR="$(spack location -i suite-sparse@5.3.0 % "$compiler")"

# SuperLU_MT
SLUMTDIR_32=${APPDIR}/superlu_mt/SuperLU_MT_3.1_fpic
SLUMTDIR_64=${APPDIR}/superlu_mt/SuperLU_MT_3.1_long_int_fpic

# SuperLU_DIST
SLUDISTDIR_32="$(spack location -i superlu-dist@develop~int64+openmp % "$compiler")"
SLUDISTDIR_64="$(spack location -i superlu-dist@develop+int64+openmp % "$compiler")"

# hypre
HYPREDIR_32="$(spack location -i hypre@2.14.0~int64 % "$compiler")"
HYPREDIR_64="$(spack location -i hypre@2.14.0+int64 % "$compiler")"

# petsc
PETSCDIR_32="$(spack location -i petsc@3.10.3~int64 % $compiler)"
PETSCDIR_64="$(spack location -i petsc@3.10.3+int64 % $compiler)"

# number of threads in OpenMP examples
export OMP_NUM_THREADS=4

# ------------------------------------------------------------------------------
# Create tarballs
# ------------------------------------------------------------------------------

# location of testing directory
testdir=`pwd`

# remove old tarball directory and create new directory
\rm -rf tarballs
mkdir tarballs

# run tarscript to create tarballs
cd ../scripts
./tarscript $module

# move tarballs to tarball directory
case $module in
    all)
        mv ../../sundials-*.tar.gz $testdir/tarballs/.
        mv ../../arkode-*.tar.gz $testdir/tarballs/.
        mv ../../cvode-*.tar.gz $testdir/tarballs/.
        mv ../../cvodes-*.tar.gz $testdir/tarballs/.
        mv ../../ida-*.tar.gz $testdir/tarballs/.
        mv ../../idas-*.tar.gz $testdir/tarballs/.
        mv ../../kinsol-*.tar.gz $testdir/tarballs/.
        ;;
    sundials)
        mv ../../sundials-*.tar.gz $testdir/tarballs/.
        ;;
    arkode)
        mv ../../arkode-*.tar.gz $testdir/tarballs/.
        ;;
    cvode)
        mv ../../cvode-*.tar.gz $testdir/tarballs/.
        ;;
    cvodes)
        mv ../../cvodes-*.tar.gz $testdir/tarballs/.
        ;;
    ida)
        mv ../../ida-*.tar.gz $testdir/tarballs/.
        ;;
    idas)
        mv ../../idas-*.tar.gz $testdir/tarballs/.
        ;;
    kinsol)
        mv ../../kinsol-*.tar.gz $testdir/tarballs/.
        ;;
    *)
        echo "Invalid module $module"
        exit 1
        ;;
esac

cd $testdir/tarballs/.

for tarball in *.tar.gz; do

    # get package name
    package=${tarball%.tar.gz}

    # --------------------------------------------------------------------------
    # Uncompress package and setup build
    # --------------------------------------------------------------------------

    tar -xvzf $tarball 2>&1 | tee -a tar.log

    rc=${PIPESTATUS[0]}
    echo -e "\ntar -xzvf returned $rc\n" | tee -a tar.log
    if [ $rc -ne 0 ]; then exit 1; fi

    # move log to package directory
    mv tar.log $package/.

    # move to package directory
    cd $package

    # loop over build options
    for lt in "${libtype[@]}"; do
        for rt in "${realtype[@]}"; do
            for it in "${indexsize[@]}"; do

                # remove old build and install directories
                \rm -rf build_${rt}_${it}_${lt}
                \rm -rf install_${rt}_${it}_${lt}

                # create new build and install directories
                mkdir build_${rt}_${it}_${lt}
                mkdir install_${rt}_${it}_${lt}

                # move to build directory
                cd build_${rt}_${it}_${lt}

                # set library type to build
                if [ "${lt}" == "static" ]; then
                    STATIC=ON
                    SHARED=OFF
                else
                    STATIC=OFF
                    SHARED=ON
                fi

                # -------------------------------------------------------------------
                # Enable third party library settings based on precision and index size
                # -------------------------------------------------------------------

                # LAPACK/BLAS does not support extended precision or 64-bit indices
                if [ "$rt" == "extended" ] || [ "$it" == "64" ]; then
                    LAPACKSTATUS=OFF
                    BLASSTATUS=OFF
                else
                    LAPACKSTATUS=ON
                    BLASSTATUS=ON
                fi

                # KLU does not support single or extended precision
                if [ "$rt" == "single" ] || [ "$rt" == "extended" ]; then
                    KLUSTATUS=OFF
                else
                    KLUSTATUS=ON
                fi

                # SuperLU MT index size must be set a build time
                if [ "$it" == "32" ]; then
                    SLUMTDIR=$SLUMTDIR_32
                else
                    SLUMTDIR=$SLUMTDIR_64
                fi

                # SuperLU MT does not support extended precision
                if [ "$rt" == "extended" ]; then
                    SLUMTSTATUS=OFF
                else
                    SLUMTSTATUS=ON
                fi

                # SuperLU MT index size must be set a build time
                if [ "$it" == "32" ]; then
                    SLUDISTDIR=$SLUDISTDIR_32
                else
                    SLUDISTDIR=$SLUDISTDIR_64
                fi

                # SuperLU DIST does not support single or extended precision
                if [ "$rt" == "double" ]; then
                    SLUDISTSTATUS=ON
                else
                    SLUDISTSTATUS=OFF
                fi

                # hypre index size must be set a build time
                if [ "$it" == "32" ]; then
                    HYPREDIR=$HYPREDIR_32
                else
                    HYPREDIR=$HYPREDIR_64
                fi

                # only testing hypre with double precision at this time
                if [ "$rt" != "double" ]; then
                    HYPRESTATUS=OFF
                else
                    HYPRESTATUS=ON
                fi

                # PETSc index size must be set a build time
                if [ "$it" == "32" ]; then
                    PETSCDIR=$PETSCDIR_32
                else
                    PETSCDIR=$PETSCDIR_64
                fi

                # only testing PETSc with double precision at this time
                if [ "$rt" != "double" ]; then
                    PETSCSTATUS=OFF
                else
                    PETSCSTATUS=ON
                fi

                # -------------------------------------------------------------------
                # Configure SUNDIALS with CMake
                #
                # NOTE: Helpful options for debugging CMake
                #
                # The '-LAH' flag lists the non-advanced cached variables (L), the
                # advanced variables (A), and help for each variable (H). This will
                # not print any system variables.
                #
                # The CMake option '-D CMAKE_VERBOSE_MAKEFILE=ON' enables additional
                # output during compile time which is useful for debugging build
                # issues.
                #
                # Setting the shared linker flags to
                # '-D CMAKE_SHARED_LINKER_FLAGS="-Wl,--no-undefined"'
                # is useful for finding undefined references when building shared
                # libraries
                # -------------------------------------------------------------------

                # only run development tests with double precision
                if [ "$rt" != "double" ]; then
                    DEVTESTS=OFF
                else
                    DEVTESTS=ON
                fi

                echo "START CMAKE"
                cmake \
                    -D CMAKE_INSTALL_PREFIX="../install_${rt}_${it}_${lt}" \
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
                    -D SUNDIALS_INDEX_SIZE=$it \
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
                    -D CMAKE_C_FLAGS='-g -Wall -std=c99 -pedantic' \
                    -D CMAKE_CXX_FLAGS='-g' \
                    -D CMAKE_Fortran_FLAGS='-g -ffpe-summary=none' \
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
                    -D SUPERLUMT_ENABLE="${SLUMTSTATUS}" \
                    -D SUPERLUMT_INCLUDE_DIR="${SLUMTDIR}/SRC" \
                    -D SUPERLUMT_LIBRARY_DIR="${SLUMTDIR}/lib" \
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
                    ../. 2>&1 | tee configure.log

                # check cmake return code
                rc=${PIPESTATUS[0]}
                echo -e "\ncmake returned $rc\n" | tee -a configure.log
                if [ $rc -ne 0 ]; then exit 1; fi

                # -------------------------------------------------------------------
                # Make SUNDIALS
                # -------------------------------------------------------------------

                echo "START MAKE"
                make -j 4 2>&1 | tee make.log

                # check make return code
                rc=${PIPESTATUS[0]}
                echo -e "\nmake returned $rc\n" | tee -a make.log
                if [ $rc -ne 0 ]; then exit 1; fi

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
                cd ..

            done
        done
    done

    # return to tarball directory
    cd ..

done

# if we make it here all packages and tests have passed
exit 0
