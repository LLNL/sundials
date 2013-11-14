#!/bin/bash

# set path for SUNDIALS build
export PATH=/usr/casc/sundials/nightly/bin:$PATH

# Python
source /usr/casc/apps/python/2.7.3/setup.sh

# cmake
export CMAKE_HOME=/usr/casc/sundials/nightly/bin/cmake-2.8.10.2
export PATH=$CMAKE_HOME/bin:$PATH

# openmpi
source /usr/casc/sundials/apps/rh5/openmpi/1.4.5/setup.sh
#export PATH=/usr/casc/sundials/apps/rh5/openmpi/1.4.5/bin:$PATH
umask 002
