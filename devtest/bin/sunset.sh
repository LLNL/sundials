#!/bin/bash

# set path for SUNDIALS build
export PATH=/usr/casc/sundials/devtest/bin:$PATH

# Python
source /usr/apps/python/2.7.3/setup.sh

# cmake
export CMAKE_HOME=/usr/casc/sundials/apps/rh6/cmake-2.8.10.2
export PATH=$CMAKE_HOME/bin:$PATH

# openmpi
source /usr/casc/sundials/apps/rh6/openmpi/1.4.5/setup.sh
#export PATH=/usr/casc/sundials/apps/rh6/openmpi/1.4.5/bin:$PATH
umask 002
