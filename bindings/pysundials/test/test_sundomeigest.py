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
# Unit/smoke tests for SUNDomEigEstimator module
# -----------------------------------------------------------------

import pytest
from fixtures import *
from pysundials.core import *


def make_estimator(estimator_type, sunctx):
	if estimator_type == "power":
		nv_view = NVectorView.Create(N_VNew_Serial(5, sunctx.get()))
		e = SUNDomEigEstimator_Power(nv_view.get(), 10, 1.0, sunctx.get())
		return SUNDomEigEstimatorView.Create(e), nv_view
	else:
		raise ValueError("Unknown estimator type")

@pytest.mark.parametrize("estimator_type", ["power"])
def test_create_estimator(estimator_type, sunctx):
	e_view, nv_view = make_estimator(estimator_type, sunctx)
	assert e_view is not None

@pytest.mark.parametrize("estimator_type", ["power"])
def test_set_max_iters(estimator_type, sunctx):
	e_view, nv_view = make_estimator(estimator_type, sunctx)
	status = SUNDomEigEstimator_SetMaxIters(e_view.get(), 10)
	assert status == 0

@pytest.mark.parametrize("estimator_type", ["power"])
def test_set_num_preprocess_iters(estimator_type, sunctx):
	e_view, nv_view = make_estimator(estimator_type, sunctx)
	status = SUNDomEigEstimator_SetNumPreprocessIters(e_view.get(), 2)
	assert status == 0

@pytest.mark.parametrize("estimator_type", ["power"])
def test_set_reltol(estimator_type, sunctx):
	e_view, nv_view = make_estimator(estimator_type, sunctx)
	status = SUNDomEigEstimator_SetRelTol(e_view.get(), 1e-6)
	assert status == 0

@pytest.mark.parametrize("estimator_type", ["power"])
def test_set_initial_guess(estimator_type, sunctx):
	e_view, nv_view = make_estimator(estimator_type, sunctx)
	status = SUNDomEigEstimator_SetInitialGuess(e_view.get(), nv_view.get())
	assert status == 0

@pytest.mark.parametrize("estimator_type", ["power"])
def test_initialize(estimator_type, sunctx):
	e_view, nv_view = make_estimator(estimator_type, sunctx)
	status = SUNDomEigEstimator_Initialize(e_view.get())
	assert status == 0

# @pytest.mark.parametrize("estimator_type", ["power"])
# def test_estimate(estimator_type, sunctx):
# 	e_view, nv_view = make_estimator(estimator_type, sunctx)
# 	err, lambdaR, lambdaI = SUNDomEigEstimator_Estimate(e_view.get(), 1.0, 0.0)
# 	assert isinstance(err, int)
# 	assert isinstance(lambdaR, float)
# 	assert isinstance(lambdaI, float)

@pytest.mark.parametrize("estimator_type", ["power"])
def test_get_res(estimator_type, sunctx):
	e_view, nv_view = make_estimator(estimator_type, sunctx)
	err, res = SUNDomEigEstimator_GetRes(e_view.get())
	assert isinstance(err, int)
	assert isinstance(res, float)

@pytest.mark.parametrize("estimator_type", ["power"])
def test_get_num_iters(estimator_type, sunctx):
	e_view, nv_view = make_estimator(estimator_type, sunctx)
	err, num_iters = SUNDomEigEstimator_GetNumIters(e_view.get())
	assert isinstance(err, int)
	assert isinstance(num_iters, int)

@pytest.mark.parametrize("estimator_type", ["power"])
def test_get_num_atimes_calls(estimator_type, sunctx):
	e_view, nv_view = make_estimator(estimator_type, sunctx)
	err, num_atimes = SUNDomEigEstimator_GetNumATimesCalls(e_view.get())
	assert isinstance(err, int)
	assert isinstance(num_atimes, int)
