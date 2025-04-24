import numpy as np
# import jax
# import jax.numpy as jnp
# from numba import jit
from pysundials.core import *

class AnalyticODEProblem:
  def __init__(self, lamb = 1.0):
    self.lamb = lamb

  def rhs(self, t, yvec, ydotvec):
    y = N_VGetArrayPointer(yvec)
    ydot = N_VGetArrayPointer(ydotvec)
    ydot[0] = self.lamb * y[0] + 1.0 / (1.0 + t * t) - 1.0 * np.arctan(t)
    return 0

# class AnalyticODEProblemJaxJit:
#   lamb = 1.0

#   def rhs(t, yvec, ydotvec, user_data):
#     y = N_VGetArrayPointer(yvec)
#     ydot = N_VGetArrayPointer(ydotvec)
#     ydot[:] = AnalyticODEProblemJaxJit.rhs_jit(t, y, ydot)
#     return 0

#   @jax.jit
#   def rhs_jit(t, y, ydot):
#     return ydot.at[0].set(AnalyticODEProblemJaxJit.lamb * y[0] + 1.0 / (1.0 + t * t) - 1.0 * jnp.arctan(t))

# class AnalyticODEProblemNumbaJit:
#   lamb = 1.0

#   @jit(nopython=True)
#   def rhs(t, yvec, ydotvec, user_data):
#     # NOTE: user_data should not be messed with since we use it for the callback wrapper
#     y = N_VGetArrayPointer(yvec)
#     ydot = N_VGetArrayPointer(ydotvec)
#     ydot[0] = AnalyticODEProblemNumbaJit.lamb * y[0] + 1.0 / (1.0 + t * t) - 1.0 * np.arctan(t)
#     return 0

class MyIMEXODEProblem:
    def __init__(self, lamb=1.0):
        self.lamb = lamb

    def fe(self, t, yvec, ydotvec):
        y = N_VGetArrayPointer(yvec)
        ydot = N_VGetArrayPointer(ydotvec)
        ydot[0] = self.lamb * y[0] + 1.0 / (1.0 + t * t) - 1.0 * np.arctan(t)
        return 0

    def fi(self, t, yvec, ydotvec):
        y = N_VGetArrayPointer(yvec)
        ydot = N_VGetArrayPointer(ydotvec)
        ydot[0] = self.lamb * y[0] + 1.0 / (1.0 + t * t) - 1.0 * np.arctan(t)
        return 0
