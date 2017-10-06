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
# matplotlib-based plotting script for brusselator1D.c example

# imports
import sys
import pylab as plt
import numpy as np

# load mesh data file
mesh = np.loadtxt('bruss_mesh.txt', dtype=np.double)

# load solution data files
udata = np.loadtxt('bruss_u.txt', dtype=np.double)
vdata = np.loadtxt('bruss_v.txt', dtype=np.double)
wdata = np.loadtxt('bruss_w.txt', dtype=np.double)

# determine number of time steps, mesh size
nt,nx = np.shape(udata)

# determine min/max values
umin = 0.9*udata.min()
umax = 1.1*udata.max()
vmin = 0.9*vdata.min()
vmax = 1.1*vdata.max()
wmin = 0.9*wdata.min()
wmax = 1.1*wdata.max()
minval = np.array([umin, vmin, wmin]).min()
maxval = np.array([umax, vmax, wmax]).max()

# generate plots of results
for tstep in range(nt):

    # set string constants for output plots, current time, mesh size
    pname = 'brusselator1D.' + repr(tstep).zfill(3) + '.png'
    tstr  = repr(tstep)
    nxstr = repr(nx)

    # plot current solution and save to disk
    plt.figure(1)
    plt.plot(mesh,udata[tstep,:],label="u")
    plt.plot(mesh,vdata[tstep,:],label="v")
    plt.plot(mesh,wdata[tstep,:],label="w")
    plt.xlabel('x')
    plt.ylabel('solution')
    plt.title('Solutions at output ' + tstr + ', mesh = ' + nxstr)
    plt.axis((0.0, 1.0, minval, maxval))
    plt.grid()
    plt.legend(loc='upper right', shadow=True)
    plt.savefig(pname)
    plt.close()


##### end of script #####
