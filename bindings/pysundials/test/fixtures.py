# -----------------------------------------------------------------
# Programmer(s): Cody J. Balos @ LLNL
# -----------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2025, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# -----------------------------------------------------------------

import pytest

from pysundials.core import *


@pytest.fixture
def sunctx():
    sunctx = SUNContextView.Create()
    yield sunctx


@pytest.fixture
def nvec(sunctx):
    nvec = NVectorView.Create(N_VNew_Serial(10, sunctx.get()))
    yield nvec
