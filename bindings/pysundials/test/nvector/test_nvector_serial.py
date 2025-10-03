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

def test_create_serial_nvector(sunctx):
	nvec = NVectorView.Create(N_VNew_Serial(5, sunctx.get()))
	assert nvec is not None
	arr = N_VGetArrayPointer(nvec.get())
	assert arr.shape[0] == 5

def test_nvconst(sunctx):
	nvec = NVectorView.Create(N_VNew_Serial(5, sunctx.get()))
	N_VConst(2.0, nvec.get())
	arr = N_VGetArrayPointer(nvec.get())
	assert np.allclose(arr, 2.0)

def test_setarraypointer(sunctx):
	nvec = NVectorView.Create(N_VNew_Serial(5, sunctx.get()))
	arr = np.array([10, 20, 30, 40, 50], dtype=float)
	N_VSetArrayPointer(arr, nvec.get())
	# Perform a C-side operation
	N_VConst(7.0, nvec.get())
	# Python should see the change
	arr2 = N_VGetArrayPointer(nvec.get())
	assert np.allclose(arr2, 7.0)
	# The original numpy array should also reflect the change
	assert np.allclose(arr, 7.0)

# Test an operation that involves vector arrays
def test_nvlinearcombination(sunctx):

    # Create two serial nvectors
    nvec1 = NVectorView.Create(N_VNew_Serial(5, sunctx.get()))
    nvec2 = NVectorView.Create(N_VNew_Serial(5, sunctx.get()))

    # Set their values
    arr1 = N_VGetArrayPointer(nvec1.get())
    arr2 = N_VGetArrayPointer(nvec2.get())
    arr1[:] = np.array([1, 2, 3, 4, 5], dtype=float)
    arr2[:] = np.array([10, 20, 30, 40, 50], dtype=float)

    # Prepare coefficients and vectors
    c = np.array([0.5, 2.0], dtype=float)
    X = [nvec1.get(), nvec2.get()]
    z = NVectorView.Create(N_VNew_Serial(5, sunctx.get()))

    # Perform linear combination: z = 0.5*x1 + 2.0*x2
    N_VLinearCombination(2, c, X, z.get())

    arr_z = N_VGetArrayPointer(z.get())
    expected = 0.5 * arr1 + 2.0 * arr2
    assert np.allclose(arr_z, expected)

# Test an operation that involves arrays of vector arrays
def test_nvlinearcombinationvectorarray(sunctx):
    # Create arrays of nvectors
    nvecs1 = [NVectorView.Create(N_VNew_Serial(5, sunctx.get())) for _ in range(3)]
    nvecs2 = [NVectorView.Create(N_VNew_Serial(5, sunctx.get())) for _ in range(3)]

    # Set their values
    for i in range(3):
        arr1 = N_VGetArrayPointer(nvecs1[i].get())
        arr2 = N_VGetArrayPointer(nvecs2[i].get())
        arr1[:] = np.arange(1, 6, dtype=float) + i
        arr2[:] = np.arange(10, 15, dtype=float) + i

    # Prepare coefficients: shape (num_vecs, num_terms)
    c = np.array([[1.0, -0.5], [2.0, 0.5], [-1.0, 1.0]], dtype=float)
    # Prepare X: shape (num_terms, num_vecs)
    X = np.array([[nvecs1[i].get() for i in range(3)],
                  [nvecs2[i].get() for i in range(3)]])
    # Prepare Z: shape (num_vecs,)
    Z = [NVectorView.Create(N_VNew_Serial(5, sunctx.get())) for _ in range(3)]

    # Perform linear combination vector array: for each i, Z[i] = c[i,0]*X[0,i] + c[i,1]*X[1,i]
    N_VLinearCombinationVectorArray(2, 3, c, X, [z.get() for z in Z])

    # Check results
    for i in range(3):
        arr1 = N_VGetArrayPointer(nvecs1[i].get())
        arr2 = N_VGetArrayPointer(nvecs2[i].get())
        arr_z = N_VGetArrayPointer(Z[i].get())
        expected = c[i,0] * arr1 + c[i,1] * arr2
        assert np.allclose(arr_z, expected)
        
