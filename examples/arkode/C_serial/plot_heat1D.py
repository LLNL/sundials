#!/usr/bin/env python
# ----------------------------------------------------------------
# Programmer(s): Daniel R. Reynolds @ SMU
# ----------------------------------------------------------------
# LLNS/SMU Copyright Start
# Copyright (c) 2015, Southern Methodist University and 
# Lawrence Livermore National Security
#
# This work was performed under the auspices of the U.S. Department 
# of Energy by Southern Methodist University and Lawrence Livermore 
# National Laboratory under Contract DE-AC52-07NA27344.
# Produced at Southern Methodist University and the Lawrence 
# Livermore National Laboratory.
#
# All rights reserved.
# For details, see the LICENSE file.
# LLNS/SMU Copyright End
# Copyright (c) 2013, Southern Methodist University.
# All rights reserved.
# For details, see the LICENSE file.
# ----------------------------------------------------------------
# matplotlib-based plotting script for heat1D.c example

# imports
import sys
import pylab as plt
import numpy as np

# load mesh data file
mesh = np.loadtxt('heat_mesh.txt', dtype=np.double)

# load solution data file
data = np.loadtxt('heat1D.txt', dtype=np.double)

# determine number of time steps, mesh size
nt,nx = np.shape(data)

# determine maximum temperature
maxtemp = 1.1*data.max()

# generate plots of results
for tstep in range(nt):

    # set string constants for output plots, current time, mesh size
    pname = 'heat1d.' + repr(tstep).zfill(3) + '.png'
    tstr  = repr(tstep)
    nxstr = repr(nx)

    # plot current solution and save to disk
    plt.figure(1)
    plt.plot(mesh,data[tstep,:])
    plt.xlabel('x')
    plt.ylabel('solution')
    plt.title('u(x) at output ' + tstr + ', mesh = ' + nxstr)
    plt.axis((0.0, 1.0, 0.0, maxtemp))
    plt.grid()
    plt.savefig(pname)
    plt.close()


##### end of script #####
