#!/bin/python

import sys
import numpy as np
from pysundials.core import *
from pysundials.arkode import *
from analytic_ode_problem import *

def test_sprkstep():
    sunctx = SUNContextView.Create()
    nv = NVectorView.Create(N_VNew_Serial(1, sunctx.get()))

    arr = N_VGetArrayPointer(nv.get())
    arr[0] = 0.0

    # TODO(CJB): not the right problem for this
    ode_problem = AnalyticIMEXODEProblem()

    def f1(t, y, ydot, _):
        return ode_problem.fe(t, y, ydot)

    def f2(t, y, ydot, _):
        return ode_problem.fi(t, y, ydot)

    sprk = ARKodeView.Create(SPRKStepCreate(f1, f2, 0, nv.get(), sunctx.get()))
    status = ARKodeSStolerances(sprk.get(), 1e-6, 1e-6)
    tout, tret = 10.0, 0.0
    status = ARKodeEvolve(sprk.get(), tout, nv.get(), tret, ARK_NORMAL)
    print(f"status={status}, ans={arr}")

if __name__ == "__main__":
    test_sprkstep()
