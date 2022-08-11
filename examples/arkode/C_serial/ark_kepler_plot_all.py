#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

def load_results(case):
   t = np.loadtxt('ark_kepler_times' + case + '.txt', dtype=np.float64)
   y = np.loadtxt('ark_kepler_solution' + case + '.txt', dtype=np.float64)
   y = np.reshape(y, (y.shape[0]//4, 4))
   conserved = np.loadtxt('ark_kepler_conserved' + case + '.txt', delimiter=',', dtype=np.float64)
   return t, y, conserved


sprk1 = load_results('_sprk-1')
sprk2 = load_results('_sprk-2')
sprk3 = load_results('_sprk-3')
sprk4 = load_results('_sprk-4')
erk2 = load_results('_erk-2')
erk3 = load_results('_erk-3')
# sprkinc1 = load_results('_sprkinc-1')
# sprkinc2 = load_results('_sprkinc-2')
# sprkinc3 = load_results('_sprkinc-3')

all_to_compare = []
all_to_compare.append(sprk1)
all_to_compare.append(sprk2)
all_to_compare.append(sprk3)
all_to_compare.append(sprk4)
all_to_compare.append(erk2)
all_to_compare.append(erk3)
# all_to_compare.append(sprkinc1)
# all_to_compare.append(sprkinc2)
# all_to_compare.append(sprkinc3)

h = 0.01
t = sprk1[0]
trend_line_erk_2 = t*h*h

energy = []
energy_0 = []
for _, sol, consv in all_to_compare:
  energy.append(consv[:,0])
  energy_0.append(consv[0,0])
energy = np.array(energy)
energy_0 = np.array(energy_0)
energy_diff = np.abs(energy.T - energy_0)

plt.figure(dpi=200)
plt.title('Error in Hamiltonian, h = %g' % h)
plt.plot(t, energy_diff)
plt.ylabel('| error |')
plt.xlabel('<---  t  --->')
plt.xscale('log')
plt.yscale('log')
legend = []
legend.append(r'$O(h^1)$ SPRK')
legend.append(r'$O(h^2)$ SPRK')
legend.append(r'$O(h^3)$ SPRK')
legend.append(r'$O(h^4)$ SPRK')
legend.append(r'$O(h^2)$ ERK')
legend.append(r'$O(h^3)$ ERK')
# legend.append(r'$O(h^1)$ SPRK (inc)')
# legend.append(r'$O(h^2)$ SPRK (inc)')
# legend.append(r'$O(h^3)$ SPRK (inc)')
plt.legend(legend)
plt.savefig('ark_kepler_energy_compare.png')

momentum = []
momentum_0 = []
for _, sol, consv in all_to_compare:
  momentum.append(consv[:,1])
  momentum_0.append(consv[0,1])
momentum = np.array(momentum)
momentum_0 = np.array(momentum_0)
momentum_diff = np.abs(momentum.T - momentum_0)

plt.figure(dpi=200)
plt.title('Error in Momentum, h = %g' % h)
plt.plot(t, momentum_diff)
plt.ylabel('| error |')
plt.xlabel('<---  t  --->')
plt.xscale('log')
plt.yscale('log')
legend = []
legend.append(r'$O(h^1)$ SPRK')
legend.append(r'$O(h^2)$ SPRK')
legend.append(r'$O(h^3)$ SPRK')
legend.append(r'$O(h^4)$ SPRK')
legend.append(r'$O(h^2)$ ERK')
legend.append(r'$O(h^3)$ ERK')
# legend.append(r'$O(h^1)$ SPRK (inc)')
# legend.append(r'$O(h^2)$ SPRK (inc)')
# legend.append(r'$O(h^3)$ SPRK (inc)')
plt.legend(legend)
plt.savefig('ark_kepler_momentum_compare.png')

# _, sol1, _ = sprk1
# _, sol2, _ = sprkinc1
# phase_diff = np.abs(sol1 - sol2)
# print(phase_diff)

# plt.figure(dpi=200)
# plt.title('Error in Phase, h = %g' % h)
# plt.plot(t, momentum_diff)
# plt.ylabel('| error |')
# plt.xlabel('t')
# plt.yscale('log')
# legend = []
# legend.append(r'$O(h^1)$ SPRK')
# # legend.append(r'$O(h^2)$ SPRK')
# # legend.append(r'$O(h^3)$ SPRK')
# legend.append(r'$O(h^1)$ SPRK (inc)')
# # legend.append(r'$O(h^2)$ SPRK (inc)')
# # legend.append(r'$O(h^3)$ SPRK (inc)')
# plt.legend(legend)
# plt.savefig('ark_kepler_phase_compare.png')