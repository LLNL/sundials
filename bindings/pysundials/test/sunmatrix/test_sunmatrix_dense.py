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

def test_create_dense_matrix(sunctx):
	rows, cols = 3, 2
	A = SUNMatrixView.Create(SUNDenseMatrix(rows, cols, sunctx.get()))
	assert A is not None
	assert SUNMatGetID(A.get()) == SUNMATRIX_DENSE
	# Ensure the shape is being translated correctly
	dataA_shape = np.shape(SUNDenseMatrix_Data(A.get()))
	assert dataA_shape[0] == rows and dataA_shape[1] == cols

def test_clone_matrix(sunctx):
	A = SUNMatrixView.Create(SUNDenseMatrix(2, 2, sunctx.get()))
	B = SUNMatClone(A.get())
	assert B is not None

def test_zero_matrix(sunctx):
	A = SUNMatrixView.Create(SUNDenseMatrix(2, 2, sunctx.get()))
	ret = SUNMatZero(A.get())
	assert ret == 0

def test_copy_matrix(sunctx):
	A = SUNMatrixView.Create(SUNDenseMatrix(2, 2, sunctx.get()))
	B = SUNMatrixView.Create(SUNDenseMatrix(2, 2, sunctx.get()))
	# Set some values in A
	dataA = SUNDenseMatrix_Data(A.get())
	dataA[0, 0] = 1.0
	ret = SUNMatCopy(A.get(), B.get())
	assert ret == 0
	dataB = SUNDenseMatrix_Data(B.get())
	assert dataB[0, 0] == 1.0

def test_scale_add_matrix(sunctx):
	A = SUNMatrixView.Create(SUNDenseMatrix(2, 2, sunctx.get()))
	B = SUNMatrixView.Create(SUNDenseMatrix(2, 2, sunctx.get()))
	dataA = SUNDenseMatrix_Data(A.get())
	dataB = SUNDenseMatrix_Data(B.get())
	dataA[:, :] = 1.0
	dataB[:, :] = 2.0
	ret = SUNMatScaleAdd(3.0, A.get(), B.get())
	assert ret == 0
	# A should now be 3*A + B = 3*1 + 2 = 5
	assert np.allclose(dataA, 5.0)

def test_scale_add_identity(sunctx):
	A = SUNMatrixView.Create(SUNDenseMatrix(2, 2, sunctx.get()))
	dataA = SUNDenseMatrix_Data(A.get())
	dataA[:, :] = 2.0
	ret = SUNMatScaleAddI(3.0, A.get())
	assert ret == 0
	# A should now be 3*A + I
	expected = np.eye(2) + 6.0
	assert np.allclose(dataA, expected)

def test_matvec(sunctx, nvec):
	A = SUNMatrixView.Create(SUNDenseMatrix(2, 2, sunctx.get()))
	dataA = SUNDenseMatrix_Data(A.get())
	dataA[:, :] = [[1.0, 2.0], [3.0, 4.0]]
	x = NVectorView.Create(N_VNew_Serial(2, sunctx.get()))
	y = NVectorView.Create(N_VNew_Serial(2, sunctx.get()))
	N_VConst(1.0, x.get())
	ret = SUNMatMatvec(A.get(), x.get(), y.get())
	assert ret == 0
	arr = N_VGetArrayPointer(y.get())
	assert np.allclose(arr, [3.0, 7.0])
