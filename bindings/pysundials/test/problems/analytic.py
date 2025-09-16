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


class AnalyticMultiscaleODE:
    """
    * We consider the initial value problem
    *    y' + lambda*y = y^2, y(0) = 1
    * proposed in
    *
    * Estep, D., et al. "An a posteriori–a priori analysis of multiscale operator
    * splitting." SIAM Journal on Numerical Analysis 46.3 (2008): 1116-1146.
    *
    * The parameter lambda is positive, t is in [0, 1], and the exact solution is
    *
    *    y(t) = lambda*y / (y(0) - (y(0) - lambda)*exp(lambda*t))
    *
    """

    T0 = 0.0
    TF = 1.0

    def __init__(self, lamb=2.0):
        self.lamb = lamb

    def f_linear(self, t, yvec, ydotvec):
        y = N_VGetArrayPointer(yvec)
        ydot = N_VGetArrayPointer(ydotvec)
        ydot[:] = -self.lamb * y
        return 0

    def f_nonlinear(self, t, yvec, ydotvec):
        y = N_VGetArrayPointer(yvec)
        ydot = N_VGetArrayPointer(ydotvec)
        ydot[:] = y * y
        return 0

    def solution(self, y0vec, yvec, tf):
        y0 = N_VGetArrayPointer(y0vec)[0]
        y = N_VGetArrayPointer(yvec)
        y[0] = self.lamb * y0 / (y0 - (y0 - self.lamb) * np.exp(self.lamb * tf))
        return 0

    def set_init_cond(self, y0vec):
        y0 = N_VGetArrayPointer(y0vec)
        y0[0] = 1.0
        return 0


class AnalyticDAE:
    """
    * The following is a simple example problem with analytical
    * solution adapted from example 10.2 of Ascher & Petzold, "Computer
    * Methods for Ordinary Differential Equations and
    * Differential-Algebraic Equations," SIAM, 1998, page 267:
    *    x1'(t) = (1-alpha)/(t-2)*x1 - x1 + (alpha-1)*x2 + 2*exp(t)
    *         0 = (t+2)*x1 - (t+2)*exp(t)
    * for t in the interval [0.0, 1.0], with initial condition:
    *    x1(0) = 1   and   x2(0) = -1/2.
    * The problem has true solution
    *    x1(t) = exp(t)  and  x2(t) = exp(t)/(t-2)
    """

    T0 = 0.0
    TF = 1.0

    def __init__(self, alpha=10.0):
        self.alpha = alpha

    def res(self, t, yyvec, ypvec, resvec):
        yy = N_VGetArrayPointer(yyvec)
        yp = N_VGetArrayPointer(ypvec)
        res = N_VGetArrayPointer(resvec)
        alpha = self.alpha

        # x1'(t) = (1-alpha)/(t-2)*x1 - x1 + (alpha-1)*x2 + 2*exp(t)
        res[0] = yp[0] - (
            (1 - alpha) / (t - 2.0) * yy[0] - yy[0] + (alpha - 1) * yy[1] + 2 * np.exp(t)
        )
        # 0 = (t+2)*x1 - (t+2)*exp(t)
        res[1] = (t + 2.0) * yy[0] - (t + 2.0) * np.exp(t)
        return 0

    def solution(self, yyvec, ypvec, t):
        yy = N_VGetArrayPointer(yyvec)
        yp = N_VGetArrayPointer(ypvec)
        yy[0] = np.exp(t)
        yy[1] = np.exp(t) / (t - 2.0)
        yp[0] = np.exp(t)
        yp[1] = np.exp(t) / (t - 2.0) - np.exp(t) / (t - 2.0) / (t - 2.0)
        return 0

    def set_init_cond(self, yyvec, ypvec, t0):
        return self.solution(yyvec, ypvec, t0)

    def psolve(self, t, yyvec, ypvec, rrvec, rvec, zvec, cj, delta):
        """
        Exact solution as preconditioner
        A = df/dy + cj*df/dyp
                =>
                    A = [ - cj - (alpha - 1)/(t - 2) - 1, alpha - 1]
                        [                          t + 2,         0]

        */
        """
        yy = N_VGetArrayPointer(yyvec)
        yp = N_VGetArrayPointer(ypvec)
        r = N_VGetArrayPointer(rvec)
        z = N_VGetArrayPointer(zvec)
        alpha = self.alpha
        A = np.array([[-cj - (alpha - 1.0) / (t - 2.0) - 1.0, alpha - 1.0], [t + 2.0, 0.0]])
        z[0] = r[1] / A[1, 0]
        z[1] = -(A[0, 0] * r[1] - A[1, 0] * r[0]) / (A[0, 1] * A[1, 0])
        print(z)
        return 0
