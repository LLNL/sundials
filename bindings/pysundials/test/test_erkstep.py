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
from problems import AnalyticODE


def test_erkstep_with_postprocess(sunctx):
    y = NVectorView.Create(N_VNew_Serial(1, sunctx.get()))
    ode_problem = AnalyticODE()
    ode_problem.set_init_cond(y.get())

    def rhs(t, y, ydot, _):
        return ode_problem.f(t, y, ydot)

    postprocess_called = {"count": 0}

    def postprocess_fn(t, y, _):
        postprocess_called["count"] += 1
        return 0  # success

    erk = ARKodeView.Create(ERKStepCreate(rhs, 0, y.get(), sunctx.get()))
    ARKodeSetPostprocessStepFn(erk.get(), postprocess_fn)
    status = ARKodeSStolerances(erk.get(), 1e-10, 1e-10)
    assert status == ARK_SUCCESS

    tout, tret = 10.0, 0.0
    status, tret = ARKodeEvolve(erk.get(), tout, y.get(), tret, ARK_NORMAL)
    assert status == ARK_SUCCESS
    assert postprocess_called["count"] > 0

    sol = NVectorView.Create(N_VClone(y.get()))
    ode_problem.solution(y.get(), sol.get(), tret)
    assert np.allclose(N_VGetArrayPointer(sol.get()), N_VGetArrayPointer(y.get()), atol=1e-2)


def test_erkstep(sunctx):
    y = NVectorView.Create(N_VNew_Serial(1, sunctx.get()))
    ode_problem = AnalyticODE()
    ode_problem.set_init_cond(y.get())

    def rhs(t, y, ydot, _):
        return ode_problem.f(t, y, ydot)

    erk = ARKodeView.Create(ERKStepCreate(rhs, 0, y.get(), sunctx.get()))
    status = ARKodeSStolerances(erk.get(), 1e-10, 1e-10)
    assert status == ARK_SUCCESS

    tout, tret = 10.0, 0.0
    status, tret = ARKodeEvolve(erk.get(), tout, y.get(), tret, ARK_NORMAL)
    assert status == ARK_SUCCESS

    sol = NVectorView.Create(N_VClone(y.get()))
    ode_problem.solution(y.get(), sol.get(), tret)
    assert np.allclose(N_VGetArrayPointer(sol.get()), N_VGetArrayPointer(y.get()), atol=1e-2)
