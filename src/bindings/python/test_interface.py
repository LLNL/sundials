#!/bin/python

import numpy as np
from pysundials.core import *
from pysundials.arkode import *

sunctx = SUNContextView()
nv = NVectorView(N_VNew_Serial(1, sunctx.get()))

# Get the array and change a value in it
arr = N_VGetArrayPointer(nv.get()) # Option 1: have to call get when passing the NVectorView
# arr = nv.GetArrayPointer() # Option 2: wrap N_V calls as NVectorView class methods
arr[0] = 0.0

def rhs(t, yvec, ydotvec, user_data):
  y = N_VGetArrayPointer(yvec)
  ydot = N_VGetArrayPointer(ydotvec)
  ydot[0] = 1.0 * y[0] + 1.0 / (1.0 + t * t) - 1.0 * np.atan(t)
  return 0

erk = ARKodeView(ERKStepCreate(rhs, 0, nv.get(), sunctx.get()))
ARKodeSStolerances(erk.get(), 1e-6, 1e-6)

tout, tret = 10.0, 0.0
status = ARKodeEvolve(erk.get(), tout, nv.get(), tret, ARK_NORMAL)
print(arr)

# Try again
erk = ARKodeView(ERKStepCreate(rhs, 0, nv.get(), sunctx.get()))
ARKodeSStolerances(erk.get(), 1e-6, 1e-6)

tret = 0.0
status = ARKodeEvolve(erk.get(), tout, nv.get(), tret, ARK_NORMAL)
