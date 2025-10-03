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
    tout, tret = 2 * np.pi, 0.0
    dt = 0.01
    y = NVectorView.Create(N_VNew_Serial(2, sunctx.get()))
    ode_problem = HarmonicOscillatorODE()

    def f1(t, y, ydot, _):
        return ode_problem.xdot(t, y, ydot)

    def f2(t, y, ydot, _):
        return ode_problem.vdot(t, y, ydot)

    arr = N_VGetArrayPointer(y.get())
    ode_problem.initial_conditions(arr)

    sprk = ARKodeView.Create(SPRKStepCreate(f1, f2, 0, y.get(), sunctx.get()))

    status = ARKodeSetFixedStep(sprk.get(), dt)
    assert status == 0

    status = ARKodeSetMaxNumSteps(sprk.get(), int(np.ceil(tout / dt)))
    assert status == 0

    status, tret = ARKodeEvolve(sprk.get(), tout, y.get(), tret, ARK_NORMAL)
    assert status == 0
