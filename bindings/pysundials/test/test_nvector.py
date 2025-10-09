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


@pytest.mark.parametrize("vector_type", ["serial"])
def test_create_nvector(vector_type, sunctx):
	if vector_type == "serial":
		nvec = NVectorView.Create(N_VNew_Serial(5, sunctx.get()))
	else:
		raise ValueError("Unknown vector type")
	assert nvec.get() is not None

	arr = N_VGetArrayPointer(nvec.get())
	assert arr.shape[0] == 5

	arr[:] = np.array([5.0, 4.0, 3.0, 2.0, 1.0], dtype=np.float64)
	assert np.allclose(N_VGetArrayPointer(nvec.get()), [5.0, 4.0, 3.0, 2.0, 1.0])

	N_VConst(2.0, nvec.get())
	assert np.allclose(arr, 2.0)


@pytest.mark.parametrize("vector_type", ["serial"])
def test_make_nvector(vector_type, sunctx):
	arr = np.array([1.0, 2.0, 3.0, 4.0, 5.0], dtype=np.float64)
	if vector_type == "serial":
		nvec = NVectorView.Create(N_VMake_Serial(5, arr, sunctx.get()))
	else:
		raise ValueError("Unknown vector type")
	assert nvec.get() is not None

	N_VConst(2.0, nvec.get())
	assert np.allclose(arr, 2.0)

	arr[:] = np.array([5.0, 4.0, 3.0, 2.0, 1.0], dtype=np.float64)
	assert np.allclose(N_VGetArrayPointer(nvec.get()), [5.0, 4.0, 3.0, 2.0, 1.0])



@pytest.mark.parametrize("vector_type", ["serial"])
def test_setarraypointer(vector_type, sunctx):
	if vector_type == "serial":
		nvec = NVectorView.Create(N_VNew_Serial(5, sunctx.get()))
	else:
		raise ValueError("Unknown vector type")
	assert nvec.get() is not None

	arr = np.array([1.0, 2.0, 3.0, 4.0, 5.0], dtype=np.float64)
	N_VSetArrayPointer(arr, nvec.get())

	assert np.allclose(N_VGetArrayPointer(nvec.get()), arr)

	N_VScale(2.0, nvec.get(), nvec.get())

	assert np.allclose(arr, [2.0, 4.0, 6.0, 8.0, 10.0])


# Test an operation that involves vector arrays
@pytest.mark.parametrize("vector_type", ["serial"])
def test_nvlinearcombination(vector_type, sunctx):
	if vector_type == "serial":
		nvec1 = NVectorView.Create(N_VNew_Serial(5, sunctx.get()))
		nvec2 = NVectorView.Create(N_VNew_Serial(5, sunctx.get()))
	else:
		raise ValueError("Unknown vector type")

	arr1 = N_VGetArrayPointer(nvec1.get())
	arr1[:] = [1.0, 2.0, 3.0, 4.0, 5.0]

	arr2 = N_VGetArrayPointer(nvec2.get())
	arr2[:] = [10.0, 20.0, 30.0, 40.0, 50.0]

	c = np.array([1.0, 0.1], dtype=np.float64)
	X = [nvec1.get(), nvec2.get()]

	z = NVectorView.Create(N_VNew_Serial(5, sunctx.get()))
	N_VConst(0.0, z.get())

	N_VLinearCombination(2, c, X, z.get())

	assert np.allclose(N_VGetArrayPointer(z.get()), [2.0, 4.0, 6.0, 8.0, 10.0])

