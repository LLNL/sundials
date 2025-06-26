#!/bin/python

import numpy as np
from pysundials.core import *
from pysundials.cvodes import *
from ode_problems import AnalyticODE

def test_implicit():
    print("  testing implicit")
    sunctx = SUNContextView.Create()
    nv = NVectorView.Create(N_VNew_Serial(1, sunctx.get()))
    ls = SUNLinearSolverView.Create(SUNLinSol_SPGMR(nv.get(), 0, 0, sunctx.get()))

    # Get the array and change a value in it
    arr = N_VGetArrayPointer(nv.get())
    arr[0] = 0.0  # set initial condition

    ode_problem = AnalyticODE()
    ark = CVodeView.Create(
        CVodeCreate(
            lambda t, y, ydot, _: ode_problem.fe(t, y, ydot), 0, nv.get(), sunctx.get()
        )
    )
    status = CVodeSStolerances(ark.get(), 1e-6, 1e-6)
    status = CVodeSetLinearSolver(ark.get(), ls.get(), None)

    tout, tret = 10.0, 0.0
    status = CVodeEvolve(ark.get(), tout, nv.get(), tret, CV_NORMAL)
    print(f"status={status}, ans={arr}")


if __name__ == "__main__":
    test_implicit()
