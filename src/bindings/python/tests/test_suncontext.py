#!/bin/python

import numpy as np
from pysundials.core import *


def test_suncontext():
    print("  testing SUNContextView")

    # Create without comm
    sunctx = SUNContextView.Create()

    # Try calling a SUNContext_ function
    last_err = SUNContext_GetLastError(sunctx.get())

    # Create a new context with a null comm
    sunctx = SUNContextView.Create(SUN_COMM_NULL)

if __name__ == "__main__":
    test_suncontext()
