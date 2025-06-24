#!/bin/python

import sys
import numpy as np
from pysundials.core import *
from pysundials.arkode import *
from analytic_ode_problem import *

def test_lsrkstep():
    sunctx = SUNContextView.Create()
    nv = NVectorView.Create(N_VNew_Serial(1, sunctx.get()))

    arr = N_VGetArrayPointer(nv.get())
    arr[0] = 0.0

    ode_problem = AnalyticODEProblem()

    def rhs(t, y, ydot, _):
        return ode_problem.rhs(t, y, ydot)

    lsrk = ARKodeView.Create(LSRKStepCreateSTS(rhs, 0, nv.get(), sunctx.get()))
    status = ARKodeSStolerances(lsrk.get(), 1e-6, 1e-6)
    tout, tret = 10.0, 0.0
    status = ARKodeEvolve(lsrk.get(), tout, nv.get(), tret, ARK_NORMAL)
    print(f"status={status}, ans={arr}")

if __name__ == "__main__":
    test_lsrkstep()
