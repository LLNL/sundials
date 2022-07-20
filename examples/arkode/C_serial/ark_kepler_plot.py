#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

def Hamiltonian(y):
  q1 = y[0]
  q2 = y[1]
  p1 = y[2]
  p2 = y[3]

  h = 0.5 * ( p1**2 + p2**2 ) - 1.0 / np.sqrt ( q1**2 + q2**2 ) \
    - 0.005 / ( np.sqrt ( q1**2 + q2**2 ) )**3 / 2.0

  return h

def AngularMomentum(y):
  q1 = y[0]
  q2 = y[1]
  p1 = y[2]
  p2 = y[3]

  L = q1*p2 - q2*p1

  return L


t = np.loadtxt('ark_kepler_times.txt', dtype=np.float64)
y = np.loadtxt('ark_kepler_solution.txt', dtype=np.float64)
y = np.reshape(y, (y.shape[0]//4, 4))

plt.figure(dpi=200)
plt.plot(y[:,0], y[:,1])
plt.savefig('ark_kepler_phase_plot.png')

conserved = np.loadtxt('ark_kepler_conserved.txt', delimiter=',', dtype=np.float64)
energy = conserved[:,0]
energy_0 = Hamiltonian(y[0])
L = conserved[:,1]
L_0 = AngularMomentum(y[0])

plt.figure(dpi=200)
plt.title('Error in Hamiltonian')
plt.plot(t, np.abs(energy-energy_0))
plt.ylabel('| error |')
plt.xlabel('time')
plt.savefig('ark_kepler_energy.png')

plt.figure(dpi=200)
plt.title('Error in Momentum')
plt.plot(t, np.abs(L-L_0))
plt.ylabel('| error |')
plt.xlabel('time')
plt.savefig('ark_kepler_momentum.png')