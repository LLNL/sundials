#!/bin/python

import numpy as np
from pysundials.core import *
from pysundials.cvodes import *
from ode_problems import AnalyticODE

def test_bdf():
    sunctx = SUNContextView.Create()
    nv = NVectorView.Create(N_VNew_Serial(1, sunctx.get()))
    ls = SUNLinearSolverView.Create(SUNLinSol_SPGMR(nv.get(), 0, 0, sunctx.get()))

    # Get the array and change a value in it
    arr = N_VGetArrayPointer(nv.get())
    arr[:] = 0.0  # set initial condition

    ode_problem = AnalyticODE()

    solver = CVodeView.Create(
        CVodeCreate(
            CV_BDF, sunctx.get()
        )
    )

    status = CVodeInit(solver.get(),
        lambda t, y, ydot, _: ode_problem.f(t, y, ydot), 0, nv.get(),
    )
    status = CVodeSStolerances(solver.get(), 1e-6, 1e-6)
    status = CVodeSetLinearSolver(solver.get(), ls.get(), None)

    tout, tret = 10.0, 0.0
    status = CVode(solver.get(), tout, nv.get(), tret, CV_NORMAL)
    print(f"status={status}, ans={arr}")

    nv2 = NVectorView.Create(N_VNew_Serial(1, sunctx.get()))
    arr2 = N_VGetArrayPointer(nv2.get())
    arr2[:] = 1.0 

    status = CVodeGetCurrentState(solver.get(), [nv2.get()])
    print(arr2)

    status, num_steps = CVodeGetNumSteps(solver.get(), 0)
    print(f"num_steps={num_steps}")


if __name__ == "__main__":
    test_bdf()
