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

def test_create_sparse_matrix(sunctx):
    rows, cols, nnz = 3, 3, 4
    A = SUNMatrixView.Create(SUNSparseMatrix(rows, cols, nnz, SUN_CSR_MAT, sunctx.get()))
    assert A is not None
    assert SUNMatGetID(A.get()) == SUNMATRIX_SPARSE
    # Check shape of data
    dataA = SUNSparseMatrix_Data(A.get())
    assert dataA.shape[0] == nnz

def test_clone_matrix(sunctx):
    rows, cols, nnz = 3, 3, 4
    A = SUNMatrixView.Create(SUNSparseMatrix(rows, cols, nnz, SUN_CSR_MAT, sunctx.get()))
    B = SUNMatClone(A.get())
    assert B is not None

def test_zero_matrix(sunctx):
    rows, cols, nnz = 3, 3, 4
    A = SUNMatrixView.Create(SUNSparseMatrix(rows, cols, nnz, SUN_CSR_MAT, sunctx.get()))
    ret = SUNMatZero(A.get())
    assert ret == 0

def test_copy_matrix(sunctx):
    rows, cols, nnz = 3, 3, 4
    A = SUNMatrixView.Create(SUNSparseMatrix(rows, cols, nnz, SUN_CSR_MAT, sunctx.get()))
    B = SUNMatrixView.Create(SUNSparseMatrix(rows, cols, nnz, SUN_CSR_MAT, sunctx.get()))
    dataA = SUNSparseMatrix_Data(A.get())
    idx_vals = SUNSparseMatrix_IndexValues(A.get())
    idx_ptrs = SUNSparseMatrix_IndexPointers(A.get())
    # CSR: row 0: col 0 (1.0), row 1: col 1 (1.0), row 2: col 2 (1.0), row 2: col 0 (2.0)
    dataA[:] = [1.0, 1.0, 1.0, 2.0]
    idx_vals[:] = [0, 1, 2, 0]
    idx_ptrs[:] = [0, 1, 2, 4]
    ret = SUNMatCopy(A.get(), B.get())
    assert ret == 0
    dataB = SUNSparseMatrix_Data(B.get())
    idx_valsB = SUNSparseMatrix_IndexValues(B.get())
    idx_ptrsB = SUNSparseMatrix_IndexPointers(B.get())
    assert np.allclose(dataB, dataA)
    assert np.allclose(idx_valsB, idx_vals)
    assert np.allclose(idx_ptrsB, idx_ptrs)

def test_scale_add_matrix(sunctx):
    rows, cols, nnz = 3, 3, 4
    A = SUNMatrixView.Create(SUNSparseMatrix(rows, cols, nnz, SUN_CSR_MAT, sunctx.get()))
    B = SUNMatrixView.Create(SUNSparseMatrix(rows, cols, nnz, SUN_CSR_MAT, sunctx.get()))
    dataA = SUNSparseMatrix_Data(A.get())
    dataB = SUNSparseMatrix_Data(B.get())
    idx_vals = SUNSparseMatrix_IndexValues(A.get())
    idx_ptrs = SUNSparseMatrix_IndexPointers(A.get())
    idx_valsB = SUNSparseMatrix_IndexValues(B.get())
    idx_ptrsB = SUNSparseMatrix_IndexPointers(B.get())
    # CSR: row 0: col 0 (1.0), row 1: col 1 (1.0), row 2: col 2 (1.0), row 2: col 0 (2.0)
    dataA[:] = [1.0, 1.0, 1.0, 2.0]
    dataB[:] = [2.0, 2.0, 2.0, 4.0]
    idx_vals[:] = [0, 1, 2, 0]
    idx_ptrs[:] = [0, 1, 2, 4]
    idx_valsB[:] = [0, 1, 2, 0]
    idx_ptrsB[:] = [0, 1, 2, 4]
    ret = SUNMatScaleAdd(3.0, A.get(), B.get())
    assert ret == 0
    # 3*A + B = [3+2, 3+2, 3+2, 6+4] = [5, 5, 5, 10]
    assert np.allclose(dataA, [5.0, 5.0, 5.0, 10.0])

def test_scale_add_identity(sunctx):
    rows, cols, nnz = 3, 3, 4
    A = SUNMatrixView.Create(SUNSparseMatrix(rows, cols, nnz, SUN_CSR_MAT, sunctx.get()))
    dataA = SUNSparseMatrix_Data(A.get())
    idx_vals = SUNSparseMatrix_IndexValues(A.get())
    idx_ptrs = SUNSparseMatrix_IndexPointers(A.get())
    # CSR: row 0: col 0 (2.0), row 1: col 1 (2.0), row 2: col 2 (2.0), row 2: col 0 (2.0)
    dataA[:] = [2.0, 2.0, 2.0, 2.0]
    idx_vals[:] = [0, 1, 2, 0]
    idx_ptrs[:] = [0, 1, 2, 4]
    ret = SUNMatScaleAddI(3.0, A.get())
    assert ret == 0
    # Diagonal elements should be 3*2+1=7, off-diagonal unchanged
    # So dataA = [7.0, 7.0, 7.0, 2.0] (assuming diagonal at idx_vals[0:3])
    # But since idx_vals = [0,1,2,0], diagonal is at positions 0,1,2
    assert np.allclose(dataA[:3], [7.0, 7.0, 7.0])

def test_matvec(sunctx):
    rows, cols, nnz = 3, 3, 4
    A = SUNMatrixView.Create(SUNSparseMatrix(rows, cols, nnz, SUN_CSR_MAT, sunctx.get()))
    x = NVectorView.Create(N_VNew_Serial(cols, sunctx.get()))
    y = NVectorView.Create(N_VNew_Serial(rows, sunctx.get()))
    N_VConst(1.0, x.get())
    # Fill a simple sparse matrix: 3x3 identity with one extra off-diagonal
    dataA = SUNSparseMatrix_Data(A.get())
    idx_vals = SUNSparseMatrix_IndexValues(A.get())
    idx_ptrs = SUNSparseMatrix_IndexPointers(A.get())
    # CSR: row 0: col 0 (1.0), row 1: col 1 (1.0), row 2: col 2 (1.0), row 2: col 0 (2.0)
    dataA[:] = [1.0, 1.0, 1.0, 2.0]
    idx_vals[:] = [0, 1, 2, 0]
    idx_ptrs[:] = [0, 1, 2, 4]
    ret = SUNMatMatvec(A.get(), x.get(), y.get())
    assert ret == 0
    # y[0] = 1.0*x[0] = 1.0
    # y[1] = 1.0*x[1] = 1.0
    # y[2] = 1.0*x[2] + 2.0*x[0] = 1.0 + 2.0 = 3.0
    assert np.allclose(N_VGetArrayPointer(y.get()), [1.0, 1.0, 3.0])
