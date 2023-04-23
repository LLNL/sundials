#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as lines
import seaborn as sns

def load_results(case):
  t = np.loadtxt('ark_kepler_times' + case + '.txt', dtype=np.float64)
  y = np.loadtxt('ark_kepler_solution' + case + '.txt', dtype=np.float64)
  y = np.reshape(y, (y.shape[0]//4, 4))
  conserved = np.loadtxt('ark_kepler_conserved' + case + '.txt', delimiter=',', dtype=np.float64)
  return t, y, conserved

def compute_energy_err(case):
  _, sol, consv = case
  energy = consv[:,0]
  energy_0 = consv[0,0]
  return np.max(np.abs(energy - energy_0))

def compute_momentum_err(case):
  _, sol, consv = case
  momentum = consv[:,1]
  momentum_0 = consv[0,1]
  return np.max(np.abs(momentum - momentum_0))

sprk = []
erk = []
# sprk.append((1000010909,load_results('_sprk-1-dt-1.00e-05')))
# sprk.append((100006814, load_results('_sprk-1-dt-1.00e-04')))
# sprk.append((10004035, load_results('_sprk-1-dt-1.00e-03')))
# sprk.append((1003175, load_results('_sprk-1-dt-1.00e-02')))
# sprk.append((200010444,load_results('_sprk-2-dt-1.00e-04')))
# sprk.append((20004886, load_results('_sprk-2-dt-1.00e-03')))
# sprk.append((2003183, load_results('_sprk-2-dt-1.00e-02')))
# sprk.append((203098, load_results('_sprk-2-dt-1.00e-01')))
sprk.append((400017704,load_results('_sprk-4-dt-1.00e-04')))
sprk.append((40006588, load_results('_sprk-4-dt-1.00e-03')))
sprk.append((4003184, load_results('_sprk-4-dt-1.00e-02')))
sprk.append((403163, load_results('_sprk-4-dt-1.00e-01')))
erk.append((581420890, load_results('_erk-4-dt-1.00e-16')))
erk.append((58012802, load_results('_erk-4-dt-1.00e-12')))
erk.append((5823147, load_results('_erk-4-dt-1.00e-08')))
erk.append((633933, load_results('_erk-4-dt-1.00e-04')))

legend = []
legend.append(r'$O(h^4) ERK$')
legend.append(r'$O(h^4) SPRK$')

blue_star = lines.Line2D([], [], color='blue', marker='*', linestyle='None',
                          markersize=10, label='Blue stars')
red_square = lines.Line2D([], [], color='red', marker='s', linestyle='None',
                          markersize=10, label='Red squares')

plt.figure(dpi=200)
plt.title('Work vs Precision')
for work, method in erk:
  plt.scatter(compute_energy_err(method), work, color='b', marker='*', label=r'$O(h^4) ERK$')
  # plt.scatter(compute_momentum_err(method), work, color='b', marker='o', label=r'$O(h^4) ERK$')
for work, method in sprk:
  plt.scatter(compute_energy_err(method), work, color='r', marker='s', label=r'$O(h^4) SPRK$')
  # plt.scatter(compute_momentum_err(method), work, color='r', marker='^', label=r'$O(h^4) SPRK$')
plt.xscale('log')
plt.yscale('log')
plt.ylabel('Number of function evals.')
plt.xlabel('Error in Energy')
plt.legend(labels=legend, handles=[blue_star, red_square])
plt.savefig('ark_kepler_work_precision.png')
