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

def test_gettype(sunctx, nvec, nls_newton):
    typ = SUNNonlinSolGetType(nls_newton.get())
    assert type is SUNNONLINEARSOLVER_ROOTFIND

# def test_sunnonlinearsolver_api(sunctx):
#     nls = SUNNonlinSol_FixedPoint(nvec.get(), 5, sunctx.get())
#     # Type
#     typ = SUNNonlinSolGetType(nls)
#     assert typ is not None
#     # Initialize
#     ret = SUNNonlinSolInitialize(nls)
#     # Setup (dummy args)
#     try:
#         SUNNonlinSolSetup(nls, nvec.get(), None)
#     except Exception:
#         pass  # Accept failure if memory arg is required
#     # Set max iters
#     ret = SUNNonlinSolSetMaxIters(nls, 10)
#     # Get num iters
#     err, niters = SUNNonlinSolGetNumIters(nls, 0)
#     assert isinstance(niters, int)
#     # Get cur iter
#     err, iter = SUNNonlinSolGetCurIter(nls, 0)
#     assert isinstance(iter, int)
#     # Get num conv fails
#     err, nconvfails = SUNNonlinSolGetNumConvFails(nls, 0)
#     assert isinstance(nconvfails, int)

