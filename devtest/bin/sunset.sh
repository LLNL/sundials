#!/bin/bash

# set path for SUNDIALS build
export PATH=/usr/casc/sundials/devtest/bin:$PATH

# Python
export PATH=/usr/apps/python/latest/bin:$PATH

# cmake
source /usr/casc/sundials/apps/rh6/cmake/cmake-2.8.10.2/setup.sh

# lapack
source /usr/casc/sundials/apps/rh6/lapack/3.5.0/setup.sh

# openmpi
source /usr/casc/sundials/apps/rh6/openmpi/1.4.5/setup.sh

# SuperLU_MT
source /usr/casc/sundials/apps/rh6/superlu_mt/SuperLU_MT_2.4/setup.sh

# variable used by Arkode OpenMP examples
export OMP_NUM_THREADS=4

umask 002
