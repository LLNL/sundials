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
from problems import AnalyticMultiscaleODE


def test_forcingstep(sunctx):
    ode_problem = AnalyticMultiscaleODE()
    t0, tf = AnalyticMultiscaleODE.T0, 0.01

    def f_linear(t, y, ydot, _):
        return ode_problem.f_linear(t, y, ydot)

    def f_nonlinear(t, y, ydot, _):
        return ode_problem.f_nonlinear(t, y, ydot)

    y = NVectorView.Create(N_VNew_Serial(1, sunctx.get()))
    y0 = NVectorView.Create(N_VNew_Serial(1, sunctx.get()))
    ode_problem.set_init_cond(y.get())
    ode_problem.set_init_cond(y0.get())

    linear_ark = ARKodeView.Create(ERKStepCreate(f_linear, t0, y.get(), sunctx.get()))
    status = ARKodeSetFixedStep(linear_ark.get(), 5e-3)
    assert status == 0

    nonlinear_ark = ARKodeView.Create(
        ARKStepCreate(f_nonlinear, None, t0, y.get(), sunctx.get())
    )
    status = ARKodeSetFixedStep(nonlinear_ark.get(), 1e-3)
    assert status == 0

    status, linear_stepper = ARKodeCreateSUNStepper(linear_ark.get())
    linear_stepper = SUNStepperView.Create(linear_stepper)
    status, nonlinear_stepper = ARKodeCreateSUNStepper(nonlinear_ark.get())
    nonlinear_stepper = SUNStepperView.Create(nonlinear_stepper)

    ark = ARKodeView.Create(
        ForcingStepCreate(
            linear_stepper.get(), nonlinear_stepper.get(), t0, y.get(), sunctx.get()
        )
    )
    status = ARKodeSetFixedStep(ark.get(), 1e-2)
    assert status == 0

    tout = tf
    status, tret = ARKodeEvolve(ark.get(), tout, y.get(), ARK_NORMAL)
    assert status == 0

    sol = NVectorView.Create(N_VClone(y.get()))
    ode_problem.solution(y0.get(), sol.get(), tf)
    assert np.allclose(N_VGetArrayPointer(sol.get()), N_VGetArrayPointer(y.get()), atol=1e-2)
