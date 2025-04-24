#!/bin/python

import numpy as np
from pysundials.core import *
from pysundials.arkode import *



ode_problem = MyODEProblem2()

sunctx = SUNContextView()
nv = NVectorView(N_VNew_Serial(1, sunctx.get()))

# Get the array and change a value in it
arr = N_VGetArrayPointer(nv.get())

# set initial condition
arr[0] = 0.0

erk = ARKodeView(
    ARKStepCreate(
        lambda t, y, ydot, _: ode_problem.fe(t, y, ydot),
        lambda t, y, ydot, _: ode_problem.fi(t, y, ydot),
        0,
        nv.get(),
        sunctx.get(),
    )
)
ARKodeSStolerances(erk.get(), 1e-6, 1e-6)
tret = 0.0
status = ARKodeEvolve(erk.get(), tout, nv.get(), tret, ARK_NORMAL)
print(arr)
