#!/usr/bin/env python
# ----------------------------------------------------------------
# Programmer(s):  Daniel R. Reynolds @ SMU
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
# matplotlib-based plotting script for ODE examples

# imports
import sys
import pylab as plt
import numpy as np

# load solution data file
data = np.loadtxt('solution.txt', dtype=np.double)

# determine number of time steps, number of fields
nt,nv = np.shape(data)

# extract time array
times = data[:,0]

# parse comment line to determine solution names
f = open('solution.txt', 'r')
commentline = f.readline()
commentsplit = commentline.split()
names = commentsplit[2:]

# create plot
plt.figure()

# add curves to figure
for i in range(nv-1):
    plt.plot(times,data[:,i+1],label=names[i])
plt.xlabel('t')
if (nv > 2):
    plt.ylabel('solutions')
else:
    plt.ylabel('solution')
plt.legend(loc='upper right', shadow=True)
plt.grid()
plt.savefig('solution.png')




##### end of script #####
