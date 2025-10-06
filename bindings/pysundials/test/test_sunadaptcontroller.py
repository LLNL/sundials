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
# Unit/smoke tests for SUNAdaptController module
# -----------------------------------------------------------------

import pytest
from fixtures import *
from pysundials.core import *

def make_controller(controller_type, sunctx):
	if controller_type == "soderlind":
		c = SUNAdaptController_Soderlind(sunctx.get())
		return SUNAdaptControllerView.Create(c), None, None
	elif controller_type == "imexgus":
		c = SUNAdaptController_ImExGus(sunctx.get())
		return SUNAdaptControllerView.Create(c), None, None
	elif controller_type == "mrihtol":
		c1 = SUNAdaptController_ImExGus(sunctx.get())
		c2 = SUNAdaptController_Soderlind(sunctx.get())
		c1_view = SUNAdaptControllerView.Create(c1)
		c2_view = SUNAdaptControllerView.Create(c2)
		c = SUNAdaptController_MRIHTol(c1_view.get(), c2_view.get(), sunctx.get())
		return SUNAdaptControllerView.Create(c), c1_view, c2_view
	else:
		raise ValueError("Unknown controller type")

@pytest.mark.parametrize("controller_type", ["soderlind", "imexgus", "mrihtol"])
def test_create_controller(controller_type, sunctx):
	c_view, c1_view, c2_view = make_controller(controller_type, sunctx)
	assert c_view is not None

@pytest.mark.parametrize("controller_type", ["soderlind", "imexgus", "mrihtol"])
def test_get_type(controller_type, sunctx):
	c_view, c1_view, c2_view = make_controller(controller_type, sunctx)
	t = SUNAdaptController_GetType(c_view.get())
	assert isinstance(t, int)

@pytest.mark.parametrize("controller_type", ["soderlind", "imexgus", "mrihtol"])
def test_estimate_step(controller_type, sunctx):
	c_view, c1_view, c2_view = make_controller(controller_type, sunctx)
	err, hnew = SUNAdaptController_EstimateStep(c_view.get(), 1.0, 1, 0.1, 1.0)
	assert isinstance(err, int)
	assert isinstance(hnew, float)

@pytest.mark.parametrize("controller_type", ["soderlind", "imexgus", "mrihtol"])
def test_reset(controller_type, sunctx):
    c_view, c1_view, c2_view = make_controller(controller_type, sunctx)
    status = SUNAdaptController_Reset(c_view.get())
    assert status == 0

@pytest.mark.parametrize("controller_type", ["soderlind", "imexgus", "mrihtol"])
def test_set_defaults(controller_type, sunctx):
    c_view, c1_view, c2_view = make_controller(controller_type, sunctx)
    status = SUNAdaptController_SetDefaults(c_view.get())
    assert status == 0

@pytest.mark.parametrize("controller_type", ["soderlind", "imexgus", "mrihtol"])
def test_set_error_bias(controller_type, sunctx):
    c_view, c1_view, c2_view = make_controller(controller_type, sunctx)
    status = SUNAdaptController_SetErrorBias(c_view.get(), 1.0)
    assert status == 0

@pytest.mark.parametrize("controller_type", ["soderlind", "imexgus", "mrihtol"])
def test_update_h(controller_type, sunctx):
    c_view, c1_view, c2_view = make_controller(controller_type, sunctx)
    status = SUNAdaptController_UpdateH(c_view.get(), 1.0, 0.1)
    assert status == 0

@pytest.mark.parametrize("controller_type", ["soderlind", "imexgus", "mrihtol"])
def test_estimate_step_tol(controller_type, sunctx):
    c_view, c1_view, c2_view = make_controller(controller_type, sunctx)
    err, hnew, tolfacnew = SUNAdaptController_EstimateStepTol(c_view.get(), 1.0, 1.0, 1, 0.1, 0.1, 1.0, 1.0)
    assert isinstance(err, int)
    assert isinstance(hnew, float)
    assert isinstance(tolfacnew, float)

@pytest.mark.parametrize("controller_type", ["mrihtol"])
def test_update_mrihtol(controller_type, sunctx):
    c_view, c1_view, c2_view = make_controller(controller_type, sunctx)
    status = SUNAdaptController_UpdateMRIHTol(c_view.get(), 1.0, 1.0, 0.1, 0.1)
    assert status == 0