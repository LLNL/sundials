#!/bin/bash
# ------------------------------------------------------------------------------
# SUNDIALS Minimal Regression Test
# ------------------------------------------------------------------------------
# configure sundials with CMake
#
# set compiler flags to check for non-standard code
# -ansi  OR  -std-c89  OR  -std=c99  OR  -std=c11
# NOTE: PETSC requires -std=c99 or newer
# 
# enable mpi
# enable openmp
# enable pthreads
#
# enable FCMIX
# enable F90
# enable C++
# ------------------------------------------------------------------------------

# check number of inputs
if [ "$#" -lt 2 ]; then
    echo "ERROR: Illegal number of parameters, real and index type required"
    exit 1
fi
realtype=$1
indextype=$2

# optional input, number of build threads
nbt=$3

# Python
export PATH=/usr/apps/python/latest/bin:$PATH

# cmake
source /usr/casc/sundials/apps/rh6/cmake/cmake-2.8.10.2/setup.sh

# openmpi
source /usr/casc/sundials/apps/rh6/openmpi/1.8.8/setup.sh

# number of threads in OpenMP examples
export OMP_NUM_THREADS=4

umask 002

# create build directory
\rm -rf suntest_minimal_${realtype}_${indextype}
mkdir suntest_minimal_${realtype}_${indextype}
cd suntest_minimal_${realtype}_${indextype}

# configure sundials with CMake
# Note the -LAH flag lists the non-advanced cached variables (L), 
# the dvanced variables (A), and help for each variable (H). This
# will not print any system variables.
cmake \
    -D SUNDIALS_PRECISION=$realtype \
    -D SUNDIALS_INDEX_TYPE=$indextype \
    \
    -D CMAKE_C_FLAGS='-Wall -std=c99 -pedantic' \
    \
    -D MPI_ENABLE=ON \
    -D OPENMP_ENABLE=ON \
    -D PTHREAD_ENABLE=ON \
    \
    -D FCMIX_ENABLE=TRUE \
    -D F90_ENABLE=TRUE \
    -D CXX_ENABLE=TRUE \
    \
    -LAH \
    ../../. 2>&1 | tee configure.log

# check return code
rc=${PIPESTATUS[0]}
echo -e "\ncmake returned $rc\n" | tee -a configure.log
if [ $rc -ne 0 ]; then
    exit 1
fi

# build sundials
make -j $nbt 2>&1 | tee make.log

# check return code
rc=${PIPESTATUS[0]}
echo -e "\nmake returned $rc\n" | tee -a make.log
if [ $rc -ne 0 ]; then
    exit 1
fi

# test sundials
make test 2>&1 | tee test.log

# check return code
rc=${PIPESTATUS[0]}
echo -e "\nmake test returned $rc\n" | tee -a test.log
if [ $rc -ne 0 ]; then
    exit 1
fi

# if we make it here all tests have passed
exit 0