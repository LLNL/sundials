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
# Unit/smoke tests for SUNStepper module
# -----------------------------------------------------------------

import pytest
from fixtures import *
from pysundials.core import *

def make_stepper(sunctx):
	# Create am empty stepper
	status, s = SUNStepper_Create(sunctx.get())
	return SUNStepperView.Create(s)

def test_create_stepper(sunctx):
	s_view = make_stepper(sunctx)
	assert s_view.get() is not None

def test_stepper_evolve(sunctx, nvec):
	s_view = make_stepper(sunctx)
	vret = nvec.get()
	err, tret = SUNStepper_Evolve(s_view.get(), 1.0, vret)
	assert isinstance(err, int)
	assert isinstance(tret, float)

def test_stepper_one_step(sunctx, nvec):
	s_view = make_stepper(sunctx)
	vret = nvec.get()
	err, tret = SUNStepper_OneStep(s_view.get(), 1.0, vret)
	assert isinstance(err, int)
	assert isinstance(tret, float)

def test_stepper_reset(sunctx, nvec):
	s_view = make_stepper(sunctx)
	err = SUNStepper_Reset(s_view.get(), 0.0, nvec.get())
	assert isinstance(err, int)

def test_stepper_set_evolve_fn(sunctx, nvec):
	s_view = make_stepper(sunctx)
	called = {"flag": False}
	def evolve_fn(stepper, tout, vret, tret):
		called["flag"] = True
		return 0
	err = SUNStepper_SetEvolveFn(s_view.get(), evolve_fn)
	assert err == 0
	# Call evolve to trigger callback
	SUNStepper_Evolve(s_view.get(), 1.0, nvec.get())
	assert called["flag"]

def test_stepper_set_one_step_fn(sunctx, nvec):
	s_view = make_stepper(sunctx)
	called = {"flag": False}
	def one_step_fn(stepper, tout, vret, tret):
		called["flag"] = True
		return 0
	err = SUNStepper_SetOneStepFn(s_view.get(), one_step_fn)
	assert err == 0
	SUNStepper_OneStep(s_view.get(), 1.0, nvec.get())
	assert called["flag"]

def test_stepper_set_full_rhs_fn(sunctx, nvec):
	s_view = make_stepper(sunctx)
	called = {"flag": False}
	def full_rhs_fn(stepper, t, v, f, mode):
		called["flag"] = True
		return 0
	err = SUNStepper_SetFullRhsFn(s_view.get(), full_rhs_fn)
	assert err == 0
	# Call with dummy args
	SUNStepper_FullRhs(s_view.get(), 0.0, nvec.get(), nvec.get(), 0)
	assert called["flag"]

def test_stepper_set_reinit_fn(sunctx, nvec):
	s_view = make_stepper(sunctx)
	called = {"flag": False}
	def reinit_fn(stepper, t, y):
		called["flag"] = True
		return 0
	err = SUNStepper_SetReInitFn(s_view.get(), reinit_fn)
	assert err == 0
	SUNStepper_ReInit(s_view.get(), 0.0, nvec.get())
	assert called["flag"]

def test_stepper_set_reset_fn(sunctx, nvec):
	s_view = make_stepper(sunctx)
	called = {"flag": False}
	def reset_fn(stepper, t, y):
		called["flag"] = True
		return 0
	err = SUNStepper_SetResetFn(s_view.get(), reset_fn)
	assert err == 0
	SUNStepper_Reset(s_view.get(), 0.0, nvec.get())
	assert called["flag"]

def test_stepper_set_reset_ckpt_idx_fn(sunctx):
	s_view = make_stepper(sunctx)
	called = {"flag": False}
	def reset_ckpt_idx_fn(stepper, idx):
		called["flag"] = True
		return 0
	err = SUNStepper_SetResetCheckpointIndexFn(s_view.get(), reset_ckpt_idx_fn)
	assert err == 0
	SUNStepper_ResetCheckpointIndex(s_view.get(), 1)
	assert called["flag"]

def test_stepper_set_stop_time_fn(sunctx):
	s_view = make_stepper(sunctx)
	called = {"flag": False}
	def stop_time_fn(stepper, tstop):
		called["flag"] = True
		return 0
	err = SUNStepper_SetStopTimeFn(s_view.get(), stop_time_fn)
	assert err == 0
	SUNStepper_SetStopTime(s_view.get(), 2.0)
	assert called["flag"]

def test_stepper_set_step_direction_fn(sunctx):
	s_view = make_stepper(sunctx)
	called = {"flag": False}
	def step_direction_fn(stepper, direction):
		called["flag"] = True
		return 0
	err = SUNStepper_SetStepDirectionFn(s_view.get(), step_direction_fn)
	assert err == 0
	SUNStepper_SetStepDirection(s_view.get(), 1)
	assert called["flag"]

def test_stepper_set_forcing_fn(sunctx, nvec):
	s_view = make_stepper(sunctx)
	called = {"flag": False}
	def forcing_fn(stepper, tshift, tscale, forcing, nforcing):
		called["flag"] = True
		return 0
	err = SUNStepper_SetForcingFn(s_view.get(), forcing_fn)
	assert err == 0
	SUNStepper_SetForcing(s_view.get(), 0.0, 1.0, [nvec.get()], 1)
	assert called["flag"]

def test_stepper_set_get_num_steps_fn(sunctx):
	s_view = make_stepper(sunctx)
	called = {"flag": False}
	def get_num_steps_fn(stepper, nst):
		called["flag"] = True
		return 0
	err = SUNStepper_SetGetNumStepsFn(s_view.get(), get_num_steps_fn)
	assert err == 0
	status, nst = SUNStepper_GetNumSteps(s_view.get())
	assert called["flag"]
	assert isinstance(nst, int)
