#!/bin/python
# -----------------------------------------------------------------
# Programmer(s): Cody J. Balos @ LLNL
# -----------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2025-2025, Lawrence Livermore National Security,
# University of Maryland Baltimore County, and the SUNDIALS contributors.
# Copyright (c) 2013-2025, Lawrence Livermore National Security
# and Southern Methodist University.
# Copyright (c) 2002-2013, Lawrence Livermore National Security.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# -----------------------------------------------------------------

import numpy as np
from pysundials.core import *
from pysundials.arkode import *
from problems import AnalyticIMEXODE


def test_explicit():
    print("  testing explicit")
    sunctx = SUNContextView.Create()
    nv = NVectorView.Create(N_VNew_Serial(1, sunctx.get()))

    # Get the array and change a value in it
    arr = N_VGetArrayPointer(nv.get())
    arr[0] = 0.0  # set initial condition

    ode_problem = AnalyticIMEXODE()
    ark = ARKodeView.Create(
        ARKStepCreate(
            lambda t, y, ydot, _: ode_problem.fe(t, y, ydot), None, 0, nv.get(), sunctx.get()
        )
    )
    status = ARKodeSStolerances(ark.get(), 1e-6, 1e-6)
    tout, tret = 10.0, 0.0
    status = ARKodeEvolve(ark.get(), tout, nv.get(), tret, ARK_NORMAL)
    print(f"status={status}, ans={arr}")


def test_implicit():
    print("  testing implicit")
    sunctx = SUNContextView.Create()
    nv = NVectorView.Create(N_VNew_Serial(1, sunctx.get()))
    ls = SUNLinearSolverView.Create(SUNLinSol_SPGMR(nv.get(), 0, 0, sunctx.get()))

    # Get the array and change a value in it
    arr = N_VGetArrayPointer(nv.get())
    arr[0] = 0.0  # set initial condition

    ode_problem = AnalyticIMEXODE()
    ark = ARKodeView.Create(
        ARKStepCreate(
            None, lambda t, y, ydot, _: ode_problem.fe(t, y, ydot), 0, nv.get(), sunctx.get()
        )
    )
    status = ARKodeSStolerances(ark.get(), 1e-6, 1e-6)
    status = ARKodeSetLinearSolver(ark.get(), ls.get(), None)

    tout, tret = 10.0, 0.0
    status = ARKodeEvolve(ark.get(), tout, nv.get(), tret, ARK_NORMAL)
    print(f"status={status}, ans={arr}")


def test_imex():
    print("  testing imex")
    sunctx = SUNContextView.Create()
    nv = NVectorView.Create(N_VNew_Serial(1, sunctx.get()))
    ls = SUNLinearSolverView.Create(SUNLinSol_SPGMR(nv.get(), 0, 0, sunctx.get()))

    # Get the array and change a value in it
    arr = N_VGetArrayPointer(nv.get())
    arr[0] = 0.0  # set initial condition

    ode_problem = AnalyticIMEXODE()
    ark = ARKodeView.Create(
        ARKStepCreate(
            lambda t, y, ydot, _: ode_problem.fe(t, y, ydot),
            lambda t, y, ydot, _: ode_problem.fi(t, y, ydot),
            0,
            nv.get(),
            sunctx.get(),
        )
    )
    status = ARKodeSStolerances(ark.get(), 1e-6, 1e-6)
    status = ARKodeSetLinearSolver(ark.get(), ls.get(), None)

    tout, tret = 10.0, 0.0
    status = ARKodeEvolve(ark.get(), tout, nv.get(), tret, ARK_NORMAL)
    print(f"status={status}, ans={arr}")


if __name__ == "__main__":
    test_explicit()
    test_implicit()
    test_imex()
