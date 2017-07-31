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
# SUNDIALS regression testing script with no external libraries enabled
# ------------------------------------------------------------------------------

# check number of inputs
if [ "$#" -lt 2 ]; then
    echo "ERROR: Illegal number of parameters, real and index type required"
    exit 1
fi
realtype=$1   # required, precision for realtypes
indextype=$2  # required, integer type for indices
nbt=$3        # optional, number of build threads

# add newer python install to path
export PATH=/usr/apps/python/latest/bin:$PATH

# number of threads in OpenMP examples
export OMP_NUM_THREADS=4

# set file permissions (rwxrwxr-x)
umask 002

# path to cmake
MYCMAKE=/usr/casc/sundials/apps/rh6/cmake/cmake-2.8.10.2/bin/cmake

# path to installed libraries
APPDIR=/usr/casc/sundials/apps/rh6
MPIDIR=${APPDIR}/openmpi/1.8.8/bin

# create build directory
\rm -rf suntest_minimal_${realtype}_${indextype}
mkdir suntest_minimal_${realtype}_${indextype}
cd suntest_minimal_${realtype}_${indextype}

# configure sundials with CMake
# Note the -LAH flag lists the non-advanced cached variables (L), 
# the dvanced variables (A), and help for each variable (H). This
# will not print any system variables.
$MYCMAKE \
    -D SUNDIALS_PRECISION=$realtype \
    -D SUNDIALS_INDEX_TYPE=$indextype \
    \
    -D CXX_ENABLE=ON \
    -D FCMIX_ENABLE=ON \
    -D F90_ENABLE=ON \
    \
    -D OPENMP_ENABLE=ON \
    -D PTHREAD_ENABLE=ON \
    \
    -D CMAKE_C_COMPILER="/usr/bin/cc" \
    -D CMAKE_CXX_COMPILER="/usr/bin/c++" \
    -D CMAKE_Fortran_COMPILER="/usr/bin/gfortran" \
    \
    -D CMAKE_C_FLAGS='-Wall -std=c99 -pedantic' \
    \
    -D MPI_ENABLE=ON \
    -D MPI_MPICC="${MPIDIR}/mpicc" \
    -D MPI_MPICXX="${MPIDIR}/mpicxx" \
    -D MPI_MPIF77="${MPIDIR}/mpif77" \
    -D MPI_MPIF90="${MPIDIR}/mpif90" \
    -D MPI_RUN_COMMAND="${MPIDIR}/mpirun" \
    \
    -LAH \
    ../../. 2>&1 | tee configure.log

# check cmake return code
rc=${PIPESTATUS[0]}
echo -e "\ncmake returned $rc\n" | tee -a configure.log
if [ $rc -ne 0 ]; then
    exit 1
fi

# build sundials
make -j $nbt 2>&1 | tee make.log

# check make return code
rc=${PIPESTATUS[0]}
echo -e "\nmake returned $rc\n" | tee -a make.log
if [ $rc -ne 0 ]; then
    exit 1
fi

# test sundials
make test 2>&1 | tee test.log

# check make test return code
rc=${PIPESTATUS[0]}
echo -e "\nmake test returned $rc\n" | tee -a test.log
if [ $rc -ne 0 ]; then
    exit 1
fi

# if we make it here all tests have passed
exit 0