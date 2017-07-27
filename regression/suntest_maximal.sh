#!/bin/bash
# ------------------------------------------------------------------------------
# SUNDIALS Maximal Regression Test
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
#
# enable lapack
# enable klu
# enable hypre
# enable PETSc
# enable SUPERLU_MT
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
\rm -rf suntest_maximal_${realtype}_${indextype}
mkdir suntest_maximal_${realtype}_${indextype}
cd suntest_maximal_${realtype}_${indextype}

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
    -D LAPACK_ENABLE=ON \
    -D LAPACK_LIBRARIES="/usr/casc/sundials/apps/rh6/lapack/3.6.0/lib64/liblapack.so;/usr/casc/sundials/apps/rh6/lapack/3.6.0/lib64/libblas.so;" \
    \
    -D KLU_ENABLE=ON \
    -D KLU_INCLUDE_DIR="/usr/casc/sundials/apps/rh6/suitesparse/4.5.3/include" \
    -D KLU_LIBRARY_DIR="/usr/casc/sundials/apps/rh6/suitesparse/4.5.3/lib" \
    \
    -D HYPRE_ENABLE=ON \
    -D HYPRE_INCLUDE_DIR="/usr/casc/sundials/apps/rh6/hypre/2.11.1/include" \
    -D HYPRE_LIBRARY="/usr/casc/sundials/apps/rh6/hypre/2.11.1/lib/libHYPRE.a" \
    \
    -D PETSC_ENABLE=ON \
    -D PETSC_INCLUDE_DIR="/usr/casc/sundials/apps/rh6/petsc/3.7.2/include" \
    -D PETSC_LIBRARY_DIR="/usr/casc/sundials/apps/rh6/petsc/3.7.2/lib" \
    \
    -D SUPERLUMT_ENABLE=ON \
    -D SUPERLUMT_INCLUDE_DIR="/usr/casc/sundials/apps/rh6/superlu_mt/SuperLU_MT_3.1/SRC" \
    -D SUPERLUMT_LIBRARY_DIR="/usr/casc/sundials/apps/rh6/superlu_mt/SuperLU_MT_3.1/lib" \
    -D SUPERLUMT_THREAD_TYPE=Pthread \
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