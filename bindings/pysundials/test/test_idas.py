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
from pysundials.core import *
from pysundials.idas import *
from problems import AnalyticDAE
from fixtures import *

def test_bdf_idas(sunctx):
    ode_problem = AnalyticDAE()

    solver = IDAView.Create(IDACreate(sunctx.get()))
    yy = NVectorView.Create(N_VNew_Serial(2, sunctx.get()))
    yp = NVectorView.Create(N_VNew_Serial(2, sunctx.get()))

    # y and y' initial conditions
    ode_problem.set_init_cond(yy.get(), yp.get(), ode_problem.T0)

    ls = SUNLinearSolverView.Create(SUNLinSol_SPGMR(yy.get(), SUN_PREC_LEFT, 0, sunctx.get()))

    def resfn(t, yy, yp, rr, _):
        return ode_problem.res(t, yy, yp, rr)

    def psolve(t, yy, yp, rr, r, z, cj, delta, _):
        return ode_problem.psolve(t, yy, yp, rr, r, z, cj, delta)

    status = IDAInit(solver.get(), resfn, 0.0, yy.get(), yp.get())
    assert status == IDA_SUCCESS

    status = IDASStolerances(solver.get(), 1e-4, 1e-9)
    assert status == IDA_SUCCESS

    status = IDASetLinearSolver(solver.get(), ls.get(), None)
    assert status == IDA_SUCCESS

    status = IDASetPreconditioner(solver.get(), None, psolve)
    assert status == IDA_SUCCESS

    tout = ode_problem.TF
    status, tret = IDASolve(solver.get(), tout, yy.get(), yp.get(), IDA_NORMAL)
    assert status == IDA_SUCCESS

    status, num_steps = IDAGetNumSteps(solver.get())
    assert status == IDA_SUCCESS
    print("Number of steps: ", num_steps)

    sol_yy = NVectorView.Create(N_VClone(yy.get()))
    sol_yp = NVectorView.Create(N_VClone(yp.get()))

    ode_problem.solution(sol_yy.get(), sol_yp.get(), tret)
    assert np.allclose(N_VGetArrayPointer(sol_yy.get()), N_VGetArrayPointer(yy.get()), rtol=1e-2)
    assert np.allclose(N_VGetArrayPointer(sol_yp.get()), N_VGetArrayPointer(yp.get()), rtol=1e-2)


def test_idas_fsa(sunctx):
    # TODO(CJB): implement
    pass


def test_idas_adjoint(sunctx):
    # TODO(CJB): implement
    pass
