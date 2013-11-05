#!/bin/bash

# set path for SUNDIALS build
export PATH=/usr/casc/sundials/nightly/bin:$PATH

# cmake
export CMAKE_HOME=/usr/casc/sundials/nightly/bin/cmake-2.8.10.2
export PATH=$CMAKE_HOME/bin:$PATH

umask 002
