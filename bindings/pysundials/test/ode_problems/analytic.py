import numpy as np

# import jax
# import jax.numpy as jnp
# from numba import jit
from pysundials.core import *


class AnalyticODE:
    def __init__(self, lamb=1.0):
        self.lamb = lamb

    def f(self, t, yvec, ydotvec):
        y = N_VGetArrayPointer(yvec)
        ydot = N_VGetArrayPointer(ydotvec)
        ydot[0] = self.lamb * y[0] + 1.0 / (1.0 + t * t) - 1.0 * np.arctan(t)
        return 0

    def dom_eig(self, t, yvec, fnvec, lambdaR, lambdaI, tempv1, tempv2, tempv3):
        lamdbaR = self.lamb
        lamdbaI = 0.0
        return 0


# class AnalyticODEJaxJit:
#   lamb = 1.0

#   def f(t, yvec, ydotvec, user_data):
#     y = N_VGetArrayPointer(yvec)
#     ydot = N_VGetArrayPointer(ydotvec)
#     ydot[:] = AnalyticODEJaxJit.f_jit(t, y, ydot)
#     return 0

#   @jax.jit
#   def f_jit(t, y, ydot):
#     return ydot.at[0].set(AnalyticODEJaxJit.lamb * y[0] + 1.0 / (1.0 + t * t) - 1.0 * jnp.arctan(t))

# class AnalyticODENumbaJit:
#   lamb = 1.0

#   @jit(nopython=True)
#   def f(t, yvec, ydotvec, user_data):
#     # NOTE: user_data should not be messed with since we use it for the callback wrapper
#     y = N_VGetArrayPointer(yvec)
#     ydot = N_VGetArrayPointer(ydotvec)
#     ydot[0] = AnalyticODENumbaJit.lamb * y[0] + 1.0 / (1.0 + t * t) - 1.0 * np.arctan(t)
#     return 0


class AnalyticIMEXODE:
    def __init__(self, lamb=1.0):
        self.lamb = lamb

    def fe(self, t, yvec, ydotvec):
        y = N_VGetArrayPointer(yvec)
        ydot = N_VGetArrayPointer(ydotvec)
        ydot[0] = self.lamb * y[0] + 1.0 / (1.0 + t * t) - 1.0 * np.arctan(t)
        return 0

    def fi(self, t, yvec, ydotvec):
        # y = N_VGetArrayPointer(yvec)
        # ydot = N_VGetArrayPointer(ydotvec)
        # ydot[0] = self.lamb * y[0] + 1.0 / (1.0 + t * t) - 1.0 * np.arctan(t)
        return 0
