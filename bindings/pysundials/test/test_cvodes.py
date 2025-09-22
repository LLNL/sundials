#!/bin/python
# -----------------------------------------------------------------
# Programmer(s): Cody J. Balos @ LLNL
# -----------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2025, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# -----------------------------------------------------------------

import numpy as np
from pysundials.core import *
from pysundials.cvodes import *
from problems import AnalyticODE


def test_bdf():
    sunctx = SUNContextView.Create()
    nv = NVectorView.Create(N_VNew_Serial(1, sunctx.get()))
    ls = SUNLinearSolverView.Create(SUNLinSol_SPGMR(nv.get(), 0, 0, sunctx.get()))

    # Get the array and change a value in it
    arr = N_VGetArrayPointer(nv.get())
    arr[:] = 0.0  # set initial condition

    ode_problem = AnalyticODE()

    solver = CVodeView.Create(CVodeCreate(CV_BDF, sunctx.get()))

    status = CVodeInit(solver.get(), lambda t, y, ydot, _: ode_problem.f(t, y, ydot), 0, nv.get())
    status = CVodeSStolerances(solver.get(), 1e-6, 1e-6)
    status = CVodeSetLinearSolver(solver.get(), ls.get(), None)

    tout, tret = 10.0, 0.0
    status, tret = CVode(solver.get(), tout, nv.get(), tret, CV_NORMAL)
    print(f"status={status}, tret={tret}, ans={arr}")

    status, num_steps = CVodeGetNumSteps(solver.get(), 0)
    print(f"num_steps={num_steps}")


if __name__ == "__main__":
    test_bdf()
