#!/bin/python

import sys
import numpy as np
from pysundials.core import *
from pysundials.arkode import *
from analytic_ode_problem import *


def test_erkstep_with_postprocess():
    sunctx = SUNContextView.Create()
    nv = NVectorView.Create(N_VNew_Serial(1, sunctx.get()))

    arr = N_VGetArrayPointer(nv.get())
    arr[0] = 0.0

    ode_problem = AnalyticODEProblem()

    def rhs(t, y, ydot, _):
        return ode_problem.rhs(t, y, ydot)

    postprocess_called = {"count": 0}

    def postprocess_fn(t, y, _):
        postprocess_called["count"] += 1
        return 0  # success

    erk = ARKodeView.Create(ERKStepCreate(rhs, 0, nv.get(), sunctx.get()))
    ARKodeSetPostprocessStepFn(erk.get(), postprocess_fn)
    status = ARKodeSStolerances(erk.get(), 1e-6, 1e-6)
    tout, tret = 10.0, 0.0
    status = ARKodeEvolve(erk.get(), tout, nv.get(), tret, ARK_NORMAL)
    print(f"status={status}, ans={arr}, postprocess_called={postprocess_called['count']}")

    # TODO(CJB): interface functions which take FILE* args
    # ARKodePrintAllStats(erk.get(), sys.stdout, SUN_OUTPUTFORMAT_TABLE)


def test_erkstep():
    sunctx = SUNContextView.Create()
    nv = NVectorView.Create(N_VNew_Serial(1, sunctx.get()))

    arr = N_VGetArrayPointer(nv.get())
    arr[0] = 0.0

    ode_problem = AnalyticODEProblem()

    def rhs(t, y, ydot, _):
        return ode_problem.rhs(t, y, ydot)

    erk = ARKodeView.Create(ERKStepCreate(rhs, 0, nv.get(), sunctx.get()))
    status = ARKodeSStolerances(erk.get(), 1e-6, 1e-6)
    tout, tret = 10.0, 0.0
    status = ARKodeEvolve(erk.get(), tout, nv.get(), tret, ARK_NORMAL)
    print(f"status={status}, ans={arr}")

    # TODO(CJB): interface functions which take FILE* args
    # ARKodePrintAllStats(erk.get(), sys.stdout, SUN_OUTPUTFORMAT_TABLE)


if __name__ == "__main__":
    test_erkstep()
    test_erkstep_with_postprocess()
