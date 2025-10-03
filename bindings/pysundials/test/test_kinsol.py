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
from fixtures import *
from pysundials.core import *
from pysundials.kinsol import *
from problems import AnalyticNonlinearSys


def test_kinsol(sunctx):
    NEQ = 3
    m_aa = 2
    tol = 1e-4
    kin_view = KINView.Create(KINCreate(sunctx.get()))
    problem = AnalyticNonlinearSys(None)
    u = NVectorView.Create(N_VNew_Serial(NEQ, sunctx.get()))

    def fp_function(u, g, _):
        return problem.fixed_point_fn(u, g)

    kin_status = KINSetMAA(kin_view.get(), m_aa)
    assert kin_status == 0
    kin_status = KINInit(kin_view.get(), fp_function, u.get())
    assert kin_status == 0
    kin_status = KINSetFuncNormTol(kin_view.get(), tol)
    assert kin_status == 0

    # initial guess
    u_data = N_VGetArrayPointer(u.get())
    u_data[:] = [0.1, 0.1, -0.1]

    # no scaling used
    scale = NVectorView.Create(N_VNew_Serial(NEQ, sunctx.get()))
    N_VConst(1.0, scale.get())

    kin_status = KINSol(kin_view.get(), u.get(), KIN_FP, scale.get(), scale.get())
    assert kin_status == 0

    u_expected = NVectorView.Create(N_VNew_Serial(NEQ, sunctx.get()))
    u_expected_data = N_VGetArrayPointer(u_expected.get())
    problem.solution(u_expected.get())
    assert np.allclose(u_data, u_expected_data, atol=1e-6)

    kin_status, nni = KINGetNumNonlinSolvIters(kin_view.get(), 0)
    kin_status, nfe = KINGetNumFuncEvals(kin_view.get(), 0)
    assert kin_status == 0
