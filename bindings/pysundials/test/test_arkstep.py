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

import pytest
import numpy as np
from fixtures import *
from pysundials.core import *
from pysundials.arkode import *
from problems import AnalyticODE, AnalyticMultiscaleODE


def test_explicit(sunctx):
    y = NVectorView.Create(N_VNew_Serial(1, sunctx.get()))

    ode_problem = AnalyticODE()

    ode_problem.set_init_cond(y.get())

    ark = ARKodeView.Create(
        ARKStepCreate(
            lambda t, y, ydot, _: ode_problem.f(t, y, ydot), None, 0, y.get(), sunctx.get()
        )
    )

    status = ARKodeSStolerances(ark.get(), 1e-10, 1e-10)
    assert status == ARK_SUCCESS

    tout = 10.0
    status, tret = ARKodeEvolve(ark.get(), tout, y.get(), ARK_NORMAL)
    assert status == ARK_SUCCESS

    sol = NVectorView.Create(N_VClone(y.get()))
    ode_problem.solution(y.get(), sol.get(), tret)

    assert np.allclose(N_VGetArrayPointer(sol.get()), N_VGetArrayPointer(y.get()), atol=1e-2)


def test_implicit(sunctx):
    y = NVectorView.Create(N_VNew_Serial(1, sunctx.get()))
    ls = SUNLinearSolverView.Create(SUNLinSol_SPGMR(y.get(), 0, 0, sunctx.get()))

    ode_problem = AnalyticODE()

    ode_problem.set_init_cond(y.get())

    ark = ARKodeView.Create(
        ARKStepCreate(
            None, lambda t, y, ydot, _: ode_problem.f(t, y, ydot), 0, y.get(), sunctx.get()
        )
    )

    status = ARKodeSStolerances(ark.get(), 1e-10, 1e-10)
    assert status == ARK_SUCCESS

    status = ARKodeSetLinearSolver(ark.get(), ls.get(), None)
    assert status == ARK_SUCCESS

    tout = 10.0
    status, tret = ARKodeEvolve(ark.get(), tout, y.get(), ARK_NORMAL)
    assert status == ARK_SUCCESS

    sol = NVectorView.Create(N_VClone(y.get()))
    ode_problem.solution(y.get(), sol.get(), tret)

    assert np.allclose(N_VGetArrayPointer(sol.get()), N_VGetArrayPointer(y.get()), atol=1e-2)


def test_imex(sunctx):
    sunctx = SUNContextView.Create()
    y = NVectorView.Create(N_VNew_Serial(1, sunctx.get()))
    ls = SUNLinearSolverView.Create(SUNLinSol_SPGMR(y.get(), 0, 0, sunctx.get()))

    ode_problem = AnalyticMultiscaleODE()

    ode_problem.set_init_cond(y.get())

    ark = ARKodeView.Create(
        ARKStepCreate(
            lambda t, y, ydot, _: ode_problem.f_nonlinear(t, y, ydot),
            lambda t, y, ydot, _: ode_problem.f_linear(t, y, ydot),
            0,
            y.get(),
            sunctx.get(),
        )
    )

    status = ARKodeSStolerances(ark.get(), 1e-10, 1e-10)
    assert status == ARK_SUCCESS

    status = ARKodeSetLinearSolver(ark.get(), ls.get(), None)
    assert status == ARK_SUCCESS

    tout = 10.0
    status, tret = ARKodeEvolve(ark.get(), tout, y.get(), ARK_NORMAL)
    assert status == 0

    sol = NVectorView.Create(N_VClone(y.get()))
    ode_problem.solution(y.get(), sol.get(), tret)

    assert np.allclose(N_VGetArrayPointer(sol.get()), N_VGetArrayPointer(y.get()), atol=1e-2)
