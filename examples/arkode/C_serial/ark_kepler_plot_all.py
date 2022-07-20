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

def Momentum(y):
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


t, y_sprk_1, consv_sprk_1 = load_results('_sprk-1')
t, y_erk_2, consv_erk_2 = load_results('_erk-2')
t, y_sprk_2, consv_sprk_2 = load_results('_sprk-2')
t, y_erk_3, consv_erk_3 = load_results('_erk-3')
t, y_sprk_3, consv_sprk_3 = load_results('_sprk-3')
t, y_erk_4, consv_erk_4 = load_results('_erk-4')

h = 0.01
trend_line_erk_2 = t*h*h

energy = []
energy.append(consv_sprk_1[:,0])
energy.append(consv_erk_2[:,0])
# energy.append(trend_line_erk_2)
energy.append(consv_sprk_2[:,0])
energy.append(consv_erk_3[:,0])
energy.append(consv_sprk_3[:,0])
energy.append(consv_erk_4[:,0])
energy = np.array(energy)

energy_0 = []
energy_0.append(Hamiltonian(y_sprk_1[0]))
energy_0.append(Hamiltonian(y_erk_2[0]))
# energy_0.append(0.0)
energy_0.append(Hamiltonian(y_sprk_2[0]))
energy_0.append(Hamiltonian(y_erk_3[0]))
energy_0.append(Hamiltonian(y_sprk_3[0]))
energy_0.append(Hamiltonian(y_erk_4[0]))
energy_0 = np.array(energy_0)

energy_diff = np.abs(energy.T - energy_0)

# L = np.array([consv_sprk_1[:,1],
#               consv_erk_2[:,1], consv_sprk_2[:,1],
#               consv_erk_3[:,1], consv_sprk_3[:,1],
#               consv_erk_4[:,1]])
# L_0 = np.array([
#   Momentum(y_sprk_1[0]),
#   Momentum(y_erk_2[0]), Momentum(y_sprk_2[0]),
#   Momentum(y_erk_3[0]), Momentum(y_sprk_3[0]),
#   Momentum(y_erk_4[0])
# ])
# L_diff = np.abs(L.T - L_0)

plt.figure(dpi=200)
plt.title('Error in Hamiltonian, h = 0.01')
plt.plot(t, energy_diff)
plt.ylabel('| error |')
plt.xlabel('t')
plt.legend([r'$O(h^1)$ SPRK',
            r'$O(h^2)$ ERK', r'$O(h^2)$ SPRK',
            r'$O(h^3)$ ERK', r'$O(h^3)$ SPRK',
            r'$O(h^4)$ ERK'])
plt.savefig('ark_kepler_energy_compare.png')

# plt.figure(dpi=200)
# plt.title('Error in Momentum, h = 0.01')
# plt.plot(t, L_diff)
# plt.ylabel('| error |')
# plt.xlabel('t')
# plt.legend({r'$O(h^1)$ SPRK',
#             r'$O(h^2)$ ERK', r'$O(h^2)$ SPRK',
#             r'$O(h^3)$ ERK', r'$O(h^3)$ SPRK',
#             r'$O(h^4)$ ERK'})
# plt.savefig('ark_kepler_momentum_compare.png')