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
from fixtures import *
from pysundials.core import *

def test_smoke_create_dense(sunctx, nvec):
	A = SUNMatrixView.Create(SUNDenseMatrix(2, 2, sunctx.get()))
	LS = SUNLinearSolverView.Create(SUNLinSol_Dense(nvec.get(), A.get(), sunctx.get()))
	assert LS is not None

def test_smoke_create_band(sunctx, nvec):
	A = SUNMatrixView.Create(SUNBandMatrix(2, 1, 1, sunctx.get()))
	LS = SUNLinearSolverView.Create(SUNLinSol_Band(nvec.get(), A.get(), sunctx.get()))
	assert LS is not None

def test_smoke_create_spgmr(sunctx, nvec):
	LS = SUNLinearSolverView.Create(SUNLinSol_SPGMR(nvec.get(), SUN_PREC_NONE, 0, sunctx.get()))
	assert LS is not None

def test_smoke_create_pcg(sunctx, nvec):
	LS = SUNLinearSolverView.Create(SUNLinSol_PCG(nvec.get(), SUN_PREC_NONE, 0, sunctx.get()))
	assert LS is not None

def test_smoke_create_spbcgs(sunctx, nvec):
	LS = SUNLinearSolverView.Create(SUNLinSol_SPBCGS(nvec.get(), SUN_PREC_NONE, 0, sunctx.get()))
	assert LS is not None

def test_smoke_create_sptfqmr(sunctx, nvec):
	LS = SUNLinearSolverView.Create(SUNLinSol_SPTFQMR(nvec.get(), SUN_PREC_NONE, 0, sunctx.get()))
	assert LS is not None

def test_smoke_get_type_and_id(sunctx, nvec):
	A = SUNMatrixView.Create(SUNDenseMatrix(2, 2, sunctx.get()))
	LS = SUNLinearSolverView.Create(SUNLinSol_Dense(nvec.get(), A.get(), sunctx.get()))
	typ = SUNLinSolGetType(LS.get())
	id_ = SUNLinSolGetID(LS.get())
	assert isinstance(typ, int)
	assert isinstance(id_, int)

def test_smoke_initialize_setup(sunctx, nvec):
	A = SUNMatrixView.Create(SUNDenseMatrix(2, 2, sunctx.get()))
	LS = SUNLinearSolverView.Create(SUNLinSol_Dense(nvec.get(), A.get(), sunctx.get()))
	ret_init = SUNLinSolInitialize(LS.get())
	ret_setup = SUNLinSolSetup(LS.get(), A.get())
	assert isinstance(ret_init, int)
	assert isinstance(ret_setup, int)

def test_smoke_num_iters_resnorm_lastflag(sunctx, nvec):
    LS = SUNLinearSolverView.Create(SUNLinSol_SPGMR(nvec.get(), 0, 0, sunctx.get()))
    niters = SUNLinSolNumIters(LS.get())
    resnorm = SUNLinSolResNorm(LS.get())
    lastflag = SUNLinSolLastFlag(LS.get())
    assert isinstance(niters, int)
    assert isinstance(resnorm, float)
    assert isinstance(lastflag, int)
