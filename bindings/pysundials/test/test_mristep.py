#!/bin/python

import numpy as np
from pysundials.core import *
from pysundials.arkode import *
from ode_problems import AnalyticMultiscaleODE


def test_multirate():
    print("  testing multirate")

    sunctx = SUNContextView.Create()

    ode_problem = AnalyticMultiscaleODE()

    t0, tf = AnalyticMultiscaleODE.T0, AnalyticMultiscaleODE.TF

    def fslow(t, y, ydot, _):
        return ode_problem.f_linear(t, y, ydot)

    def ffast(t, y, ydot, _):
        return ode_problem.f_nonlinear(t, y, ydot)

    yview = NVectorView.Create(N_VNew_Serial(1, sunctx.get()))

    ode_problem.set_init_cond(yview.get())

    # create fast integrator
    inner_ark = ARKodeView.Create(ERKStepCreate(ffast, t0, yview.get(), sunctx.get()))
    status = ARKodeSetFixedStep(inner_ark.get(), 5e-3)

    inner_stepper = MRIStepInnerStepperView.Create(
        ARKodeCreateMRIStepInnerStepper(inner_ark.get())
    )

    # # create slow integrator
    # ark = ARKodeView.Create(MRIStepCreate(fslow, None, t0, yview.get(), inner_stepper.get(), sunctx.get()))
    # ARKodeSetFixedStep(ark.get(), 1e-3)

    tret, tout = t0, tf
    status = ARKodeEvolve(inner_ark.get(), tout, yview.get(), tret, ARK_NORMAL)
    N_VPrint(yview.get())


if __name__ == "__main__":
    test_multirate()
