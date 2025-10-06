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
from problems import AnalyticMultiscaleODE


def test_multirate():
    print("  testing multirate")

    sunctx = SUNContextView.Create()

    ode_problem = AnalyticMultiscaleODE()

    t0, tf = AnalyticMultiscaleODE.T0, 0.01

    def fslow(t, y, ydot, _):
        return ode_problem.f_linear(t, y, ydot)

    def ffast(t, y, ydot, _):
        return ode_problem.f_nonlinear(t, y, ydot)

    yview = NVectorView.Create(N_VNew_Serial(1, sunctx.get()))
    y0view = NVectorView.Create(N_VNew_Serial(1, sunctx.get()))

    ode_problem.set_init_cond(yview.get())
    ode_problem.set_init_cond(y0view.get())

    # create fast integrator
    inner_ark = ARKodeView.Create(ERKStepCreate(ffast, t0, yview.get(), sunctx.get()))
    status = ARKodeSetFixedStep(inner_ark.get(), 5e-3)

    status, inner_stepper = ARKodeCreateMRIStepInnerStepper(inner_ark.get())
    inner_stepper = MRIStepInnerStepperView.Create(inner_stepper)

    # create slow integrator
    ark = ARKodeView.Create(
        MRIStepCreate(fslow, None, t0, yview.get(), inner_stepper.get(), sunctx.get())
    )
    ARKodeSetFixedStep(ark.get(), 1e-3)

    tout = tf
    status = ARKodeEvolve(ark.get(), tout, yview.get(), ARK_NORMAL)

    ode_problem.solution(y0view.get(), yview.get(), tf)


if __name__ == "__main__":
    test_multirate()
