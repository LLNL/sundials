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

import numpy as np
from pysundials.core import *
from pysundials.arkode import *
from problems import AnalyticMultiscaleODE


def test_splittingstep():
    print("  testing SplittingStep")

    sunctx = SUNContextView.Create()

    ode_problem = AnalyticMultiscaleODE()

    t0, tf = AnalyticMultiscaleODE.T0, 0.01

    def f_linear(t, y, ydot, _):
        return ode_problem.f_linear(t, y, ydot)

    def f_nonlinear(t, y, ydot, _):
        return ode_problem.f_nonlinear(t, y, ydot)

    yview = NVectorView.Create(N_VNew_Serial(1, sunctx.get()))
    y0view = NVectorView.Create(N_VNew_Serial(1, sunctx.get()))

    ode_problem.set_init_cond(yview.get())
    ode_problem.set_init_cond(y0view.get())

    # create linear partition integrator
    linear_ark = ARKodeView.Create(ERKStepCreate(f_linear, t0, yview.get(), sunctx.get()))
    status = ARKodeSetFixedStep(linear_ark.get(), 5e-3)

    # create nonlinear partition integrator
    nonlinear_ark = ARKodeView.Create(
        ARKStepCreate(f_nonlinear, None, t0, yview.get(), sunctx.get())
    )
    status = ARKodeSetFixedStep(nonlinear_ark.get(), 1e-3)

    # create SUNStepper for each partition
    status, linear_stepper = ARKodeCreateSUNStepper(linear_ark.get())
    linear_stepper = SUNStepperView.Create(linear_stepper)
    status, nonlinear_stepper = ARKodeCreateSUNStepper(nonlinear_ark.get())
    nonlinear_stepper = SUNStepperView.Create(nonlinear_stepper)

    # create the outer integrator
    steppers = [linear_stepper.get(), nonlinear_stepper.get()]
    ark = ARKodeView.Create(
        SplittingStepCreate(steppers, len(steppers), t0, yview.get(), sunctx.get())
    )

    ARKodeSetFixedStep(ark.get(), 1e-2)

    tret, tout = t0, tf
    status = ARKodeEvolve(ark.get(), tout, yview.get(), tret, ARK_NORMAL)
    N_VPrint(yview.get())

    ode_problem.solution(y0view.get(), yview.get(), tf)
    N_VPrint(yview.get())


if __name__ == "__main__":
    test_splittingstep()
