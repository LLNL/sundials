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
from pysundials.cvodes import *
from problems import AnalyticODE


def test_bdf(sunctx):
    y = NVectorView.Create(N_VNew_Serial(1, sunctx.get()))
    ls = SUNLinearSolverView.Create(SUNLinSol_SPGMR(y.get(), 0, 0, sunctx.get()))

    ode_problem = AnalyticODE()
    ode_problem.set_init_cond(y.get())

    solver = CVodeView.Create(CVodeCreate(CV_BDF, sunctx.get()))

    status = CVodeInit(solver.get(), lambda t, y, ydot, _: ode_problem.f(t, y, ydot), 0, y.get())
    assert status == CV_SUCCESS

    status = CVodeSStolerances(solver.get(), 1e-10, 1e-10)
    assert status == CV_SUCCESS

    status = CVodeSetLinearSolver(solver.get(), ls.get(), None)
    assert status == CV_SUCCESS

    tout = 10.0
    status, tret = CVode(solver.get(), tout, y.get(), CV_NORMAL)
    assert status == CV_SUCCESS

    sol = NVectorView.Create(N_VClone(y.get()))
    ode_problem.solution(y.get(), sol.get(), tret)
    assert np.allclose(N_VGetArrayPointer(sol.get()), N_VGetArrayPointer(y.get()), atol=1e-2)

    status, num_steps = CVodeGetNumSteps(solver.get())
    assert status == CV_SUCCESS
    assert num_steps > 0


def test_cvodes_fsa(sunctx):
    # TODO(CJB): implement
    pass


def test_cvodes_adjoint(sunctx):
    # TODO(CJB): implement
    pass

