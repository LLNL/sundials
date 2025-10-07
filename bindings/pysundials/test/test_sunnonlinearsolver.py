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
# Unit/smoke tests for SUNNonlinearSolver module
# -----------------------------------------------------------------

import pytest
import numpy as np
from fixtures import *
from pysundials.core import *


@pytest.mark.parametrize("solver_type", ["fixedpoint", "newton"])
def test_create_solver(solver_type, sunctx, nvec):
    if solver_type == "fixedpoint":
        nls = SUNNonlinearSolverView.Create(SUNNonlinSol_FixedPoint(nvec.get(), 5, sunctx.get()))
    elif solver_type == "newton":
        nls = SUNNonlinearSolverView.Create(SUNNonlinSol_Newton(nvec.get(), sunctx.get()))
    else:
        raise ValueError("Unknown solver type")
    assert nls is not None


@pytest.mark.parametrize("solver_type", ["fixedpoint", "newton"])
def make_solver(solver_type, sunctx, nvec):
    if solver_type == "fixedpoint":
        return SUNNonlinearSolverView.Create(SUNNonlinSol_FixedPoint(nvec.get(), 5, sunctx.get()))
    elif solver_type == "newton":
        return SUNNonlinearSolverView.Create(SUNNonlinSol_Newton(nvec.get(), sunctx.get()))
    else:
        raise ValueError("Unknown solver type")


@pytest.mark.parametrize("solver_type, expected_type", [
    ("newton", SUNNONLINEARSOLVER_ROOTFIND),
    ("fixedpoint", SUNNONLINEARSOLVER_FIXEDPOINT)
])
def test_gettype(solver_type, expected_type, sunctx, nvec):
    nls = make_solver(solver_type, sunctx, nvec)
    typ = SUNNonlinSolGetType(nls.get())
    assert typ is expected_type


@pytest.mark.parametrize("solver_type", ["fixedpoint", "newton"])
def test_initialize(solver_type, sunctx, nvec):
    nls = make_solver(solver_type, sunctx, nvec)
    ret = SUNNonlinSolInitialize(nls.get())
    assert ret == 0


@pytest.mark.parametrize("solver_type,max_iters", [
    ("newton", 5),
    ("fixedpoint", 10)
])
def test_set_max_iters_and_get_num_iters(solver_type, max_iters, sunctx, nvec):
    nls = make_solver(solver_type, sunctx, nvec)
    ret = SUNNonlinSolSetMaxIters(nls.get(), max_iters)
    assert ret == 0
    err, niters = SUNNonlinSolGetNumIters(nls.get())
    assert err == 0
    assert isinstance(niters, int)


@pytest.mark.parametrize("solver_type", ["fixedpoint", "newton"])
def test_get_cur_iter(solver_type, sunctx, nvec):
    nls = make_solver(solver_type, sunctx, nvec)
    err, cur_iter = SUNNonlinSolGetCurIter(nls.get())
    assert err == 0
    assert isinstance(cur_iter, int)


@pytest.mark.parametrize("solver_type", ["fixedpoint", "newton"])
def test_get_num_conv_fails(solver_type, sunctx, nvec):
    nls = make_solver(solver_type, sunctx, nvec)
    err, nconvfails = SUNNonlinSolGetNumConvFails(nls.get())
    assert err == 0
    assert isinstance(nconvfails, int)


def test_fixedpoint_setup_and_solve(sunctx):
    from problems import AnalyticNonlinearSys

    NEQ = AnalyticNonlinearSys.NEQ
    ucor = NVectorView.Create(N_VNew_Serial(NEQ, sunctx.get()))
    u0 = NVectorView.Create(N_VNew_Serial(NEQ, sunctx.get()))
    w = NVectorView.Create(N_VNew_Serial(NEQ, sunctx.get()))
    ucur = NVectorView.Create(N_VNew_Serial(NEQ, sunctx.get()))

    # Initial guess
    udata = N_VGetArrayPointer(u0.get())
    udata[:] = [0.1, 0.1, -0.1]

    # Initial correction
    N_VConst(0.0, ucor.get())

    # Set the weights
    N_VConst(1.0, w.get())

    # Create the problem
    with AnalyticNonlinearSys(u0.get()) as problem:

        # Create the solver
        nls = SUNNonlinearSolverView.Create(SUNNonlinSol_FixedPoint(u0.get(), 2, sunctx.get()))

        # System function
        def g_fn(u, g, _):
            return problem.corrector_fp_fn(u, g)

        # Convergence test
        def conv_test(nls, u, delta, tol, ewt, _):
            return problem.conv_test(nls, u, delta, tol, ewt)

        ret = SUNNonlinSolSetSysFn(nls.get(), g_fn)
        assert ret == 0

        ret = SUNNonlinSolSetConvTestFn(nls.get(), conv_test)
        assert ret == 0

        ret = SUNNonlinSolSetMaxIters(nls.get(), 50)

        ret = SUNNonlinSolSetup(nls.get(), u0.get())
        assert ret == 0

        tol = 1e-10
        ret = SUNNonlinSolSolve(nls.get(), u0.get(), ucor.get(), w.get(), tol, 0)
        assert ret == 0

        # Update the initial guess with the correction
        N_VLinearSum(1.0, u0.get(), 1.0, ucor.get(), ucur.get())

        # Compare to analytic solution
        utrue = NVectorView.Create(N_VNew_Serial(NEQ, sunctx.get()))
        problem.solution(utrue.get())
        assert np.allclose(
            N_VGetArrayPointer(ucur.get()), N_VGetArrayPointer(utrue.get()), atol=1e-2
        )


if __name__ == "__main__":
    test_fixedpoint_setup_and_solve(SUNContextView.Create())
