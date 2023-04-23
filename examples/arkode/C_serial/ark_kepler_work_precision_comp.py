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
sprkinc = []
sprk.append((105.695, load_results('_sprk-6-dt-1.00e-03')))
sprk.append((12.773, load_results('_sprk-6-dt-1.00e-02')))
sprk.append((3.706, load_results('_sprk-6-dt-1.00e-01')))
sprk.append((66.015, load_results('_sprk-4-dt-1.00e-03')))
sprk.append((6.772, load_results('_sprk-4-dt-1.00e-02')))
sprk.append((3.646, load_results('_sprk-4-dt-1.00e-01')))
sprk.append((45.116, load_results('_sprk-2-dt-1.00e-03')))
sprk.append((6.707, load_results('_sprk-2-dt-1.00e-02')))
sprk.append((3.439, load_results('_sprk-2-dt-1.00e-01')))
sprkinc.append((137.000, load_results('_sprkinc-6-dt-1.00e-03')))
sprkinc.append((15.720, load_results('_sprkinc-6-dt-1.00e-02')))
sprkinc.append((4.165, load_results('_sprkinc-6-dt-1.00e-01')))
sprkinc.append((84.957, load_results('_sprkinc-4-dt-1.00e-03')))
sprkinc.append((10.859, load_results('_sprkinc-4-dt-1.00e-02')))
sprkinc.append((3.597, load_results('_sprkinc-4-dt-1.00e-01')))
sprkinc.append((59.029, load_results('_sprkinc-2-dt-1.00e-03')))
sprkinc.append((8.100, load_results('_sprkinc-2-dt-1.00e-02')))
sprkinc.append((3.514, load_results('_sprkinc-2-dt-1.00e-01')))

legend = []
legend.append(r'$O(h^4)$ SPRK Compensated')
legend.append(r'$O(h^4) SPRK$')

blue_star = lines.Line2D([], [], color='blue', marker='*', linestyle='None',
                          markersize=10, label='Blue stars')
red_square = lines.Line2D([], [], color='red', marker='s', linestyle='None',
                          markersize=10, label='Red squares')

plt.figure(dpi=200)
plt.title('Work vs Precision')
for work, method in sprkinc:
  plt.scatter(compute_energy_err(method), work, color='b', marker='*', label=r'$O(h^4) sprkinc$')
  # plt.scatter(compute_momentum_err(method), work, color='b', marker='o', label=r'$O(h^4) sprkinc$')
for work, method in sprk:
  plt.scatter(compute_energy_err(method), work, color='r', marker='s', label=r'$O(h^4) SPRK$')
  # plt.scatter(compute_momentum_err(method), work, color='r', marker='^', label=r'$O(h^4) SPRK$')
plt.xscale('log')
plt.yscale('log')
plt.ylabel('Number of function evals.')
plt.xlabel('Error in Energy')
plt.legend(labels=legend, handles=[blue_star, red_square])
plt.savefig('ark_kepler_work_precision.png')
