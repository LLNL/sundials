#!/usr/bin/env python3

"""
 -----------------------------------------------------------------
 Programmer(s): Cody J. Balos @ LLNL
 -----------------------------------------------------------------
 Based on the example kinRoberts_fp.c by Carol Woodward @ LLNL
 -----------------------------------------------------------------
 Example problem:

 The following is a simple example problem, with the coding
 needed for its solution by the accelerated fixed point solver in
 KINSOL.
 The problem is from chemical kinetics, and consists of solving
 the first time step in a Backward Euler solution for the
 following three rate equations:
    dy1/dt = -.04*y1 + 1.e4*y2*y3
    dy2/dt = .04*y1 - 1.e4*y2*y3 - 3.e2*(y2)^2
    dy3/dt = 3.e2*(y2)^2
 on the interval from t = 0.0 to t = 0.1, with initial
 conditions: y1 = 1.0, y2 = y3 = 0. The problem is stiff.
 Run statistics (optional outputs) are printed at the end.
 -----------------------------------------------------------------
"""

import kinsol as kin
import numpy as np


class Problem:

  # constants
  NEQ    = 3     # number of equations
  Y10    = 1.0   # initial y components
  Y20    = 0.0
  Y30    = 0.0
  TOL    = 1e-10 # function tolerance
  DSTEP  = 0.1   # size of the single time step used
  PRIORS = 2

  # system defining function
  def funcRoberts(sunvec_y, sunvec_g, user_data):
    y = kin.N_VGetData(sunvec_y)
    g = kin.N_VGetData(sunvec_g)

    y1 = y[0]
    y2 = y[1]
    y3 = y[2]

    yd1 = Problem.DSTEP * (-0.04*y1 + 1.0e4*y2*y3)
    yd3 = Problem.DSTEP * 3.0e2*y2*y2

    g[0] = yd1 + Problem.Y10
    g[1] = -yd1 - yd3 + Problem.Y20
    g[2] = yd3 + Problem.Y30

    return 0


def solve():
  print("Example problem from chemical kinetics solving")
  print("the first time step in a Backward Euler solution for the")
  print("following three rate equations:")
  print("    dy1/dt = -.04*y1 + 1.e4*y2*y3")
  print("    dy2/dt = .04*y1 - 1.e4*y2*y3 - 3.e2*(y2)^2")
  print("    dy3/dt = 3.e2*(y2)^2")
  print("on the interval from t = 0.0 to t = 0.1, with initial")
  print("conditions: y1 = 1.0, y2 = y3 = 0.")
  print("Solution method: Anderson accelerated fixed point iteration.\n")

  # --------------------------------------
  # Create vectors for solution and scales
  # --------------------------------------

  # data arrays
  y = np.zeros(Problem.NEQ)
  scale = np.zeros(Problem.NEQ)

  # create N_Vector objects
  sunvec_y = kin.N_VMake_Serial(y)
  sunvec_scale = kin.N_VMake_Serial(scale)

  # -----------------------------------------
  # Initialize and allocate memory for KINSOL
  # -----------------------------------------

  kmem = kin.KINCreate()

  # set number of prior residuals used in Anderson acceleration
  flag = kin.KINSetMAA(kmem, Problem.PRIORS)
  if flag < 0: raise RuntimeError(f'KINSetMAA returned {flag}')

  # have to wrap the python system function so that it is callable from C
  sysfn = kin.WrapPythonSysFn(Problem.funcRoberts)

  # must call KINInitPy instead of KINInit
  flag = kin.KINInitPy(kmem, sysfn, sunvec_y)
  if flag < 0: raise RuntimeError(f'KINInitPy returned {flag}')

  # -------------------
  # Set optional inputs
  # -------------------

  # specify stopping tolerance based on residual
  fnormtol = Problem.TOL
  flag = kin.KINSetFuncNormTol(kmem, fnormtol)
  if flag < 0: raise RuntimeError(f'KINSetFuncNormTol returned {flag}')

  flag = kin.KINSetErrFilename(kmem, "error.log")
  if flag < 0: raise RuntimeError(f'KINSetErrFilename returned {flag}')


  # -------------
  # Initial guess
  # -------------

  y[0] = 1.0

  # ----------------------------
  # Call KINSOL to solve problem
  # ----------------------------

  # no scaling used
  kin.N_VConst(1.0, sunvec_scale)

  # call main solver
  flag = kin.KINSol(kmem,          # KINSOL memory block
                    sunvec_y,      # initial guess on input; solution vector
                    kin.KIN_FP,    # global strategy choice
                    sunvec_scale,  # scaling vector for the variable cc
                    sunvec_scale)  # scaling vector for function values fval
  if flag < 0: raise RuntimeError(f'KINSol returned {flag}')

  # ------------------------------------
  # Print solution and solver statistics
  # ------------------------------------
  flag, fnorm = kin.KINGetFuncNorm(kmem)
  if flag < 0: raise RuntimeError(f'KINGetFuncNorm returned {flag}')

  print('Computed solution (||F|| = %Lg):\n' % fnorm)
  print('y = %14.6e  %14.6e  %14.6e' % (y[0], y[1], y[2]))
  PrintFinalStats(kmem)

  # -----------
  # Free memory
  # -----------

  # kin.KINFree(kmem)
  kin.N_VDestroy(sunvec_y)
  kin.N_VDestroy(sunvec_scale)


def PrintFinalStats(kmem):
  flag, nni = kin.KINGetNumNonlinSolvIters(kmem)
  if flag < 0: raise RuntimeError(f'KINGetNumNonlinSolvIters returned {flag}')
  flag, nfe = kin.KINGetNumFuncEvals(kmem)
  if flag < 0: raise RuntimeError(f'KINGetNumFuncEvals returned {flag}')

  print('\nFinal Statistics.. \n')
  print('nni      = %6ld     nfe     = %6ld' % (nni, nfe))


if __name__ == '__main__':
  solve()
