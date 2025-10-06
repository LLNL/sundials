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


def test_lsrkstep(sunctx):
    y = NVectorView.Create(N_VNew_Serial(1, sunctx.get()))
    ode_problem = AnalyticODE()
    ode_problem.set_init_cond(y.get())

    def rhs(t, y, ydot, _):
        return ode_problem.f(t, y, ydot)

    def dom_eig(t, yvec, fnvec, lambdaR, lambdaI, _, tempv1, tempv2, tempv3):
        return ode_problem.dom_eig(t, yvec, fnvec, lambdaR, lambdaI, tempv1, tempv2, tempv3)

    lsrk = ARKodeView.Create(LSRKStepCreateSTS(rhs, 0, y.get(), sunctx.get()))
    status = LSRKStepSetDomEigFn(lsrk.get(), dom_eig)
    assert status == 0

    status = ARKodeSStolerances(lsrk.get(), 1e-10, 1e-10)
    assert status == ARK_SUCCESS

    status = ARKodeSetMaxNumSteps(lsrk.get(), 100000)

    tout = 10.0
    status, tret = ARKodeEvolve(lsrk.get(), tout, y.get(), ARK_NORMAL)
    assert status == ARK_SUCCESS

    sol = NVectorView.Create(N_VClone(y.get()))
    ode_problem.solution(y.get(), sol.get(), tret)
    assert np.allclose(N_VGetArrayPointer(sol.get()), N_VGetArrayPointer(y.get()), atol=1e-2)
