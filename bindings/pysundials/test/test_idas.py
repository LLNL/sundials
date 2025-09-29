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
from pysundials.idas import *
from problems import AnalyticDAE


def test_bdf_idas():
    sunctx = SUNContextView.Create()

    ode_problem = AnalyticDAE()

    solver = IDAView.Create(IDACreate(sunctx.get()))
    yy = NVectorView.Create(N_VNew_Serial(2, sunctx.get()))
    yp = NVectorView.Create(N_VNew_Serial(2, sunctx.get()))

    # y and y' initial conditions
    ode_problem.set_init_cond(yy.get(), yp.get(), ode_problem.T0)

    ls = SUNLinearSolverView.Create(SUNLinSol_SPGMR(yy.get(), SUN_PREC_LEFT, 0, sunctx.get()))

    def resfn(t, y, yp, rr, _):
        return ode_problem.res(t, y, yp, rr)

    def psetup(t, y, yp, rr, cj, _):
        print("setup")
        return 0

    def psolve(t, y, yp, rr, rvec, zvec, cj, delta, _):
        print("solve")
        return ode_problem.psolve(t, y, yp, rr, rvec, zvec, cj, delta)

    status = IDAInit(solver.get(), resfn, 0.0, yy.get(), yp.get())
    status = IDASStolerances(solver.get(), 1e-4, 1e-9)
    status = IDASetLinearSolver(solver.get(), ls.get(), None)
    status = IDASetPreconditioner(solver.get(), None, psolve)

    tout, tret = ode_problem.TF, ode_problem.T0
    status, tret = IDASolve(solver.get(), tout, tret, yy.get(), yp.get(), IDA_NORMAL)
    print(f"status={status}, tret={tret}")
    N_VPrint(yy.get())
    N_VPrint(yp.get())

    status, num_steps = IDAGetNumSteps(solver.get(), 0)
    print(f"num_steps={num_steps}")


if __name__ == "__main__":
    test_bdf_idas()
