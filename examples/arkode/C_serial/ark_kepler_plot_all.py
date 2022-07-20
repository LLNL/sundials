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


def load_results(case):
   t = np.loadtxt('ark_kepler_times' + case + '.txt', dtype=np.float64)
   y = np.loadtxt('ark_kepler_solution' + case + '.txt', dtype=np.float64)
   y = np.reshape(y, (y.shape[0]//4, 4))
   conserved = np.loadtxt('ark_kepler_conserved' + case + '.txt', delimiter=',', dtype=np.float64)
   return t, y, conserved


t, y_erk_2, consv_erk_2 = load_results('_erk-2')
t, y_sprk_2, consv_sprk_2 = load_results('_sprk-2')

plt.figure(dpi=200)
plt.plot(y_erk_2[:,0], y_erk_2[:,1])
# plt.plot(y_sprk_2[:,0], y_sprk_2[:,1])
plt.savefig('ark_kepler_phase_plot_compare.png')

energy = np.array([consv_erk_2[:,0], consv_sprk_2[:,0]])
energy_0 = np.array([Hamiltonian(y_erk_2[0]), Hamiltonian(y_sprk_2[0])])
energy_diff = np.abs(energy.T - energy_0)
L = np.array([consv_erk_2[:,1], consv_sprk_2[:,1]])
L_0 = np.array([AngularMomentum(y_erk_2[0]), AngularMomentum(y_sprk_2[0])])
L_diff = np.abs(L.T - L_0)

plt.figure(dpi=200)
plt.title('Error in Hamiltonian, h = 0.01')
plt.plot(t, energy_diff)
plt.ylabel('| error |')
plt.xlabel('t')
plt.legend({r'$O(h^2)$ ERK', r'$O(h^2)$ SPRK'})
plt.savefig('ark_kepler_energy_compare.png')

plt.figure(dpi=200)
plt.title('Error in Momentum, h = 0.01')
plt.plot(t, L_diff)
plt.ylabel('| error |')
plt.xlabel('t')
plt.legend({r'$O(h^2)$ ERK', r'$O(h^2)$ SPRK'})
plt.savefig('ark_kepler_momentum_compare.png')