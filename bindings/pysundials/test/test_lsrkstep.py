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
from problems import AnalyticODE


def test_lsrkstep():
    sunctx = SUNContextView.Create()
    nv = NVectorView.Create(N_VNew_Serial(1, sunctx.get()))

    arr = N_VGetArrayPointer(nv.get())
    arr[0] = 0.0

    ode_problem = AnalyticODE()

    def rhs(t, y, ydot, _):
        return ode_problem.f(t, y, ydot)

    def dom_eig(t, yvec, fnvec, lambdaR, lambdaI, _, tempv1, tempv2, tempv3):
        return ode_problem.dom_eig(t, yvec, fnvec, lambdaR, lambdaI, tempv1, tempv2, tempv3)

    lsrk = ARKodeView.Create(LSRKStepCreateSTS(rhs, 0, nv.get(), sunctx.get()))

    status = LSRKStepSetDomEigFn(lsrk.get(), dom_eig)
    status = ARKodeSStolerances(lsrk.get(), 1e-6, 1e-6)

    tout, tret = 10.0, 0.0
    status = ARKodeEvolve(lsrk.get(), tout, nv.get(), tret, ARK_NORMAL)
    print(f"status={status}, ans={arr}")


if __name__ == "__main__":
    test_lsrkstep()
