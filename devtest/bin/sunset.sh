#!/bin/bash

# set path for SUNDIALS build
export PATH=/usr/casc/sundials/devtest/bin:$PATH

# Python
export PATH=/usr/apps/python/latest/bin:$PATH

# git
source /usr/apps/git/2.3.1/setup.sh

# cmake
source /usr/casc/sundials/apps/rh6/cmake/3.0.2/setup.sh

# lapack
source /usr/casc/sundials/apps/rh6/lapack/3.6.0/setup.sh

# suitesparse klu
source /usr/casc/sundials/apps/rh6/suitesparse/4.5.3/setup.sh

# SuperLU_MT
source /usr/casc/sundials/apps/rh6/superlu_mt/SuperLU_MT_3.1_long_int_fpic/setup.sh

# openmpi
source /usr/casc/sundials/apps/rh6/openmpi/1.8.8/setup.sh

# hypre
source /usr/casc/sundials/apps/rh6/hypre/2.11.1_long_int_fpic/setup.sh

# petsc
source /usr/casc/sundials/apps/rh6/petsc/3.7.2_long_int/setup.sh

# variable used by Arkode OpenMP examples
export OMP_NUM_THREADS=4

umask 002
