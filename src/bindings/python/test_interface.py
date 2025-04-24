#!/bin/python

import numpy as np
# import jax
# import jax.numpy as jnp
# from numba import jit
from pysundials.core import *
from pysundials.arkode import *

class MyODEProblem:
  lamb = 1.0

  def rhs(t, yvec, ydotvec, user_data):
    # NOTE: user_data should not be messed with since we use it for the callback wrapper
    y = N_VGetArrayPointer(yvec)
    ydot = N_VGetArrayPointer(ydotvec)
    ydot[0] = MyODEProblem.lamb * y[0] + 1.0 / (1.0 + t * t) - 1.0 * np.arctan(t)
    return 0

class MyODEProblem2:
  def __init__(self, lamb = 1.1):
    self.lamb = lamb

  def rhs(self, t, yvec, ydotvec):
    y = N_VGetArrayPointer(yvec)
    ydot = N_VGetArrayPointer(ydotvec)
    ydot[0] = self.lamb * y[0] + 1.0 / (1.0 + t * t) - 1.0 * np.arctan(t)
    return 0

# class MyODEProblem3:
#   lamb = 1.2

#   def rhs(t, yvec, ydotvec, user_data):
#     y = N_VGetArrayPointer(yvec)
#     ydot = N_VGetArrayPointer(ydotvec)
#     ydot[:] = MyODEProblem3.rhs_jit(t, y, ydot)
#     return 0

#   @jax.jit
#   def rhs_jit(t, y, ydot):
#     return ydot.at[0].set(MyODEProblem3.lamb * y[0] + 1.0 / (1.0 + t * t) - 1.0 * jnp.atan(t))

# class MyODEProblem4:
#   lamb = 1.0

#   @jit(nopython=True)
#   def rhs(t, yvec, ydotvec, user_data):
#     # NOTE: user_data should not be messed with since we use it for the callback wrapper
#     y = N_VGetArrayPointer(yvec)
#     ydot = N_VGetArrayPointer(ydotvec)
#     ydot[0] = MyODEProblem.lamb * y[0] + 1.0 / (1.0 + t * t) - 1.0 * np.atan(t)
#     return 0

def main():
  sunctx = SUNContextView()
  nv = NVectorView(N_VNew_Serial(1, sunctx.get()))

  # Get the array and change a value in it
  arr = N_VGetArrayPointer(nv.get()) # Option 1: have to call get when passing the NVectorView
  # # arr = nv.GetArrayPointer() # Option 2: wrap N_V calls as NVectorView class methods
  arr[0] = 0.0 # set initial condition

  print('Try with static class')
  erk = ARKodeView(ERKStepCreate(MyODEProblem.rhs, 0, nv.get(), sunctx.get()))
  status = ARKodeSStolerances(erk.get(), 1e-6, 1e-6)
  tout, tret = 10.0, 0.0
  status = ARKodeEvolve(erk.get(), tout, nv.get(), tret, ARK_NORMAL)
  print(arr)

  print('Try again with lambda rhs')  
  ode_problem = MyODEProblem2()
  arr[0] = 0.0 # reset initial condition
  erk = ARKodeView(ERKStepCreate(lambda t, y, ydot, user_data: ode_problem.rhs(t, y, ydot), 0, nv.get(), sunctx.get()))
  status = ARKodeSStolerances(erk.get(), 1e-6, 1e-6)
  tret = 0.0
  status = ARKodeEvolve(erk.get(), tout, nv.get(), tret, ARK_NORMAL)
  print(arr)

  # # Try again with jax
  # arr[0] = 0.0 # reset initial condition
  # erk = ARKodeView(ERKStepCreate(MyODEProblem3.rhs, 0, nv.get(), sunctx.get()))
  # ARKodeSStolerances(erk.get(), 1e-6, 1e-6)
  # tret = 0.0
  # status = ARKodeEvolve(erk.get(), tout, nv.get(), tret, ARK_NORMAL)
  # print(arr)

if __name__ == "__main__":
  main()
  