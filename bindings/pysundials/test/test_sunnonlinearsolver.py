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
# Unit tests for SUNNonlinearSolver module
# -----------------------------------------------------------------

import pytest
import numpy as np
from fixtures import *
from pysundials.core import *


def test_create_fixedpoint(sunctx, nvec):
    nls = SUNNonlinearSolverView.Create(SUNNonlinSol_FixedPoint(nvec.get(), 5, sunctx.get()))
    assert nls is not None


def test_create_newton(sunctx, nvec):
    nls = SUNNonlinearSolverView.Create(SUNNonlinSol_Newton(nvec.get(), sunctx.get()))
    assert nls is not None


@pytest.fixture
def nls_fixedpoint(sunctx, nvec):
    nls = SUNNonlinearSolverView.Create(SUNNonlinSol_FixedPoint(nvec.get(), 5, sunctx.get()))
    yield nls


@pytest.fixture
def nls_newton(sunctx, nvec):
    nls = SUNNonlinearSolverView.Create(SUNNonlinSol_Newton(nvec.get(), sunctx.get()))
    yield nls


def test_gettype(sunctx, nvec, nls_newton, nls_fixedpoint):
    typ = SUNNonlinSolGetType(nls_newton.get())
    assert typ is SUNNONLINEARSOLVER_ROOTFIND

    typ = SUNNonlinSolGetType(nls_fixedpoint.get())
    assert typ is SUNNONLINEARSOLVER_FIXEDPOINT


def test_initialize(sunctx, nvec, nls_newton, nls_fixedpoint):
    # Test initialization of Newton solver
    ret = SUNNonlinSolInitialize(nls_newton.get())
    assert ret == 0  # 0 usually means success

    # Test initialization of FixedPoint solver
    ret = SUNNonlinSolInitialize(nls_fixedpoint.get())
    assert ret == 0


def test_set_max_iters_and_get_num_iters(sunctx, nvec, nls_newton, nls_fixedpoint):
    # Set max iters for Newton solver
    ret = SUNNonlinSolSetMaxIters(nls_newton.get(), 15)
    assert ret == 0

    # Set max iters for FixedPoint solver
    ret = SUNNonlinSolSetMaxIters(nls_fixedpoint.get(), 20)
    assert ret == 0

    # Get num iters for Newton solver
    err, niters = SUNNonlinSolGetNumIters(nls_newton.get(), 0)
    assert err == 0
    assert isinstance(niters, int)

    # Get num iters for FixedPoint solver
    err, niters = SUNNonlinSolGetNumIters(nls_fixedpoint.get(), 0)
    assert err == 0
    assert isinstance(niters, int)


def test_get_cur_iter(sunctx, nvec, nls_newton, nls_fixedpoint):
    # Get current iteration for Newton solver
    err, cur_iter = SUNNonlinSolGetCurIter(nls_newton.get(), 0)
    assert err == 0
    assert isinstance(cur_iter, int)

    # Get current iteration for FixedPoint solver
    err, cur_iter = SUNNonlinSolGetCurIter(nls_fixedpoint.get(), 0)
    assert err == 0
    assert isinstance(cur_iter, int)


def test_get_num_conv_fails(sunctx, nvec, nls_newton, nls_fixedpoint):
    # Get number of convergence failures for Newton solver
    err, nconvfails = SUNNonlinSolGetNumConvFails(nls_newton.get(), 0)
    assert err == 0
    assert isinstance(nconvfails, int)

    # Get number of convergence failures for FixedPoint solver
    err, nconvfails = SUNNonlinSolGetNumConvFails(nls_fixedpoint.get(), 0)
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
            return problem.fixed_point_fn(u, g)

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
