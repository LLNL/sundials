#!/bin/python

import numpy as np
from pysundials.core import *
from pysundials.arkode import *
from analytic_ode_problem import *


def test_erkstep():
    sunctx = SUNContextView.Create()
    nv = NVectorView.Create(N_VNew_Serial(1, sunctx.get()))

    # Get the array and change a value in it
    arr = N_VGetArrayPointer(nv.get())  # Option 1: have to call get when passing the NVectorView
    # arr = nv.GetArrayPointer() # Option 2: wrap N_V calls as NVectorView class methods
    arr[0] = 0.0  # set initial condition

    ode_problem = AnalyticODEProblem()
    erk = ARKodeView.Create(
        ERKStepCreate(lambda t, y, ydot, _: ode_problem.rhs(t, y, ydot), 0, nv.get(), sunctx.get())
    )
    status = ARKodeSStolerances(erk.get(), 1e-6, 1e-6)
    tout, tret = 10.0, 0.0
    status = ARKodeEvolve(erk.get(), tout, nv.get(), tret, ARK_NORMAL)
    print(f"status={status}, ans={arr}")

if __name__ == "__main__":
    test_erkstep()
