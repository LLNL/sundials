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

def test_create_band_matrix(sunctx):
    rows, mu, ml = 4, 1, 1
    A = SUNMatrixView.Create(SUNBandMatrix(rows, mu, ml, sunctx.get()))
    assert A is not None
    assert SUNMatGetID(A.get()) == SUNMATRIX_BAND
    # Ensure the shape is being translated correctly
    dataA = SUNBandMatrix_Data(A.get())
    ldata = SUNBandMatrix_LData(A.get())
    assert dataA.shape[0] == ldata

def test_clone_matrix(sunctx):
    rows, mu, ml = 4, 1, 1
    A = SUNMatrixView.Create(SUNBandMatrix(rows, mu, ml, sunctx.get()))
    B = SUNMatClone(A.get())
    assert B is not None

def test_zero_matrix(sunctx):
    rows, mu, ml = 4, 1, 1
    A = SUNMatrixView.Create(SUNBandMatrix(rows, mu, ml, sunctx.get()))
    ret = SUNMatZero(A.get())
    assert ret == 0

def test_copy_matrix(sunctx):
    rows, mu, ml = 4, 1, 1
    A = SUNMatrixView.Create(SUNBandMatrix(rows, mu, ml, sunctx.get()))
    B = SUNMatrixView.Create(SUNBandMatrix(rows, mu, ml, sunctx.get()))
    smu = SUNBandMatrix_StoredUpperBandwidth(A.get())
    dataA = SUNBandMatrix_Data(A.get())
    dataA[smu-mu] = 1.0
    ret = SUNMatCopy(A.get(), B.get())
    assert ret == 0
    dataB = SUNBandMatrix_Data(B.get())
    assert dataB[smu-mu] == 1.0

def test_scale_add_matrix(sunctx):
    rows, mu, ml = 4, 1, 1
    A = SUNMatrixView.Create(SUNBandMatrix(rows, mu, ml, sunctx.get()))
    B = SUNMatrixView.Create(SUNBandMatrix(rows, mu, ml, sunctx.get()))
    smu = SUNBandMatrix_StoredUpperBandwidth(A.get())
    dataA = SUNBandMatrix_Data(A.get())
    dataB = SUNBandMatrix_Data(B.get())
    dataA[smu-mu:smu+ml] = 1.0 # column 0 set to 1.0
    dataB[smu-mu:smu+ml] = 2.0
    ret = SUNMatScaleAdd(3.0, A.get(), B.get())
    assert ret == 0
    # A should now be 3*A + B = 3*1 + 2 = 5
    assert np.allclose(dataA[smu-mu:smu+ml], 5.0)

def test_scale_add_identity(sunctx):
    rows, mu, ml = 4, 1, 1
    A = SUNMatrixView.Create(SUNBandMatrix(rows, mu, ml, sunctx.get()))
    ldim = SUNBandMatrix_LDim(A.get())
    smu = SUNBandMatrix_StoredUpperBandwidth(A.get())
    dataA = SUNBandMatrix_Data(A.get())
    ret = SUNMatScaleAddI(0.0, A.get())
    assert ret == 0
    # A should now be I
    diag = np.array([dataA[smu + i * ldim] for i in range(rows)])
    assert np.allclose(diag, 1.0)

def test_matvec(sunctx):
    rows, mu, ml = 4, 1, 1
    A = SUNMatrixView.Create(SUNBandMatrix(rows, mu, ml, sunctx.get()))
    x = NVectorView.Create(N_VNew_Serial(rows, sunctx.get()))
    y = NVectorView.Create(N_VNew_Serial(rows, sunctx.get()))

    N_VConst(1.0, x.get())

    # Fill band matrix data for a simple 4x4 banded matrix
    # [3 2 0 0]
    # [1 3 2 0]
    # [0 1 3 2]
    # [0 0 1 3]
    dataA = SUNBandMatrix_Data(A.get())
    ldim = SUNBandMatrix_LDim(A.get())
    smu = SUNBandMatrix_StoredUpperBandwidth(A.get())
    for j in range(rows):
        # Diagonal
        dataA[smu + j * ldim] = 3.0
        # Lower diagonal
        if j > 0:
            dataA[smu - 1 + j * ldim] = 2.0
        # Upper diagonal
        if j < rows - 1:
            dataA[smu + 1 + j * ldim] = 1.0
    
    ret = SUNMatMatvec(A.get(), x.get(), y.get())
    assert ret == 0

    assert np.allclose(N_VGetArrayPointer(y.get()), [5.0, 6.0, 6.0, 4.0])  
