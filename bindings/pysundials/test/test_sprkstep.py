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
from problems import HarmonicOscillatorODE


def test_sprkstep(sunctx):
    tout = 2 * np.pi
    dt = 0.01
    y = NVectorView.Create(N_VNew_Serial(2, sunctx.get()))
    ode_problem = HarmonicOscillatorODE()

    def f1(t, y, ydot, _):
        return ode_problem.xdot(t, y, ydot)

    def f2(t, y, ydot, _):
        return ode_problem.vdot(t, y, ydot)

    ode_problem.set_init_cond(y.get())

    sprk = ARKodeView.Create(SPRKStepCreate(f1, f2, 0, y.get(), sunctx.get()))

    status = ARKodeSetFixedStep(sprk.get(), dt)
    assert status == ARK_SUCCESS

    status = ARKodeSetMaxNumSteps(sprk.get(), int(np.ceil(tout / dt)))
    assert status == ARK_SUCCESS

    status, tret = ARKodeEvolve(sprk.get(), tout, y.get(), ARK_NORMAL)
    assert status == ARK_SUCCESS

    sol = NVectorView.Create(N_VClone(y.get()))
    ode_problem.solution(y.get(), sol.get(), tret)
    assert np.allclose(N_VGetArrayPointer(sol.get()), N_VGetArrayPointer(y.get()), atol=1e-2)
