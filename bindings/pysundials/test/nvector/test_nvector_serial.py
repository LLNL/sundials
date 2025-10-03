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
	nvec = NVectorView.Create(N_VNew_Serial(2, sunctx.get()))
	arr = np.array([10, 10], dtype=np.float64)
	N_VSetArrayPointer(arr, nvec.get())
	# Perform a C-side operation
	N_VScale(2.0, nvec.get(), nvec.get())
	# Python should see the change
	arr2 = N_VGetArrayPointer(nvec.get())
	assert np.allclose(arr2, 20.0)
	# The original numpy array should also reflect the change
	assert np.allclose(arr, 20.0)

# # Test an operation that involves vector arrays
# def test_nvlinearcombination(sunctx):

#     # Create two serial nvectors
#     nvec1 = NVectorView.Create(N_VNew_Serial(5, sunctx.get()))
#     nvec2 = NVectorView.Create(N_VNew_Serial(5, sunctx.get()))

#     # Set their values
#     arr1 = N_VGetArrayPointer(nvec1.get())
#     arr1[:] = [1.0, 2.0, 3.0, 4.0, 5.0]

#     arr2 = N_VGetArrayPointer(nvec2.get())
#     arr2[:] = [10.0, 20.0, 30.0, 40.0, 50.0]

#     # Prepare coefficients and vectors
#     c = np.array([1.0, 0.1], dtype=np.float64)
#     X = [nvec1.get(), nvec2.get()]

#     z = NVectorView.Create(N_VNew_Serial(5, sunctx.get()))
#     N_VConst(0.0, z.get())

#     # Perform linear combination: z = 1.0*x1 + 0.1*x2
#     N_VLinearCombination(2, c, X, z.get())

#     assert np.allclose(N_VGetArrayPointer(z.get()), [2.0, 4.0, 6.0, 8.0, 10.0])