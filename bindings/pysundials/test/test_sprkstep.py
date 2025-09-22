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

import sys
import numpy as np
from pysundials.core import *
from pysundials.arkode import *
from problems import HarmonicOscillatorODE


def test_sprkstep():
    tout, tret = 2 * np.pi, 0.0
    dt = 0.01

    sunctx = SUNContextView.Create()

    nv = NVectorView.Create(N_VNew_Serial(2, sunctx.get()))
    arr = N_VGetArrayPointer(nv.get())

    ode_problem = HarmonicOscillatorODE()

    def f1(t, y, ydot, _):
        return ode_problem.xdot(t, y, ydot)

    def f2(t, y, ydot, _):
        return ode_problem.vdot(t, y, ydot)

    ode_problem.initial_conditions(arr)

    sprk = ARKodeView.Create(SPRKStepCreate(f1, f2, 0, nv.get(), sunctx.get()))

    status = ARKodeSetFixedStep(sprk.get(), dt)
    status = ARKodeSetMaxNumSteps(sprk.get(), int(np.ceil(tout / dt)))

    status = ARKodeEvolve(sprk.get(), tout, nv.get(), tret, ARK_NORMAL)
    print(f"status={status}, ans={arr}")


if __name__ == "__main__":
    test_sprkstep()
