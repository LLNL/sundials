#
# Python implementation of the harmonic oscillator ODE from ark_harmonic_symplectic.c
#

import numpy as np
from pysundials.core import *


class HarmonicOscillatorODE:
    def __init__(self, A=10.0, phi=0.0, omega=1.0):
        self.A = A
        self.phi = phi
        self.omega = omega

    def xdot(self, t, yvec, ydotvec):
        y = N_VGetArrayPointer(yvec)
        ydot = N_VGetArrayPointer(ydotvec)
        ydot[0] = y[1]
        return 0

    def vdot(self, t, yvec, ydotvec):
        y = N_VGetArrayPointer(yvec)
        ydot = N_VGetArrayPointer(ydotvec)
        ydot[0] = -self.omega * self.omega * y[0]
        return 0

    def initial_conditions(self, y):
        y[0] = self.A * np.cos(self.phi)
        y[1] = -self.A * self.omega * np.sin(self.phi)

    def energy(self, t, y):
        return 0.5 * (y[0] * y[0] + y[1] * y[1]) + 0.5 * self.omega * self.omega * y[0] * y[0]

    def solution(self, t):
        return self.A * np.cos(self.omega * t + self.phi), -self.A * self.omega * np.sin(
            self.omega * t + self.phi
        )
