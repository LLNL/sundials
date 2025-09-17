#!/bin/python

import numpy as np
from pysundials.core import *
from pysundials.arkode import *
from problems import AnalyticMultiscaleODE


def test_forcingstep():
    print("  testing ForcingStep")

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
    linear_stepper = SUNStepperView.Create(
        ARKodeCreateSUNStepper(linear_ark.get())
    )
    nonlinear_stepper = SUNStepperView.Create(
        ARKodeCreateSUNStepper(nonlinear_ark.get())
    )

    # create the outer integrator
    ark = ARKodeView.Create(
        ForcingStepCreate(
            linear_stepper.get(), nonlinear_stepper.get(), t0, yview.get(), sunctx.get()
        )
    )

    ARKodeSetFixedStep(ark.get(), 1e-2)

    tret, tout = t0, tf
    status = ARKodeEvolve(ark.get(), tout, yview.get(), tret, ARK_NORMAL)
    N_VPrint(yview.get())

    ode_problem.solution(y0view.get(), yview.get(), tf)
    N_VPrint(yview.get())


if __name__ == "__main__":
    test_forcingstep()
