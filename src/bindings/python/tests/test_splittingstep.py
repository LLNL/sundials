#!/bin/python

import numpy as np
from pysundials.core import *
from pysundials.arkode import *
from analytic_ode_problem import AnalyticIMEXODEProblem

def test_splittingstep():
  print('  testing implicit')

  sunctx = SUNContextView()
  stepper = SUNStepperView.Create(sunctx.get())

if __name__ == "__main__":
  test_splittingstep()
