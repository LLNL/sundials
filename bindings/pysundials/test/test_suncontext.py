#!/bin/python

import pytest
import numpy as np
from pysundials.core import *


def test_suncontext_wo_comm():
    # Create without comm
    sunctx = SUNContextView.Create()

    # Try calling a SUNContext_ function
    last_err = SUNContext_GetLastError(sunctx.get())

    assert last_err == SUN_SUCCESS

def test_with_null_comm():
    # Create a new context with a null comm
    sunctx = SUNContextView.Create(SUN_COMM_NULL)
