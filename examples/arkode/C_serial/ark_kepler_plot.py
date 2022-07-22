#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

t = np.loadtxt('ark_kepler_times_sprk-1.txt', dtype=np.float64)
y = np.loadtxt('ark_kepler_solution_sprk-1.txt', dtype=np.float64)
y = np.reshape(y, (y.shape[0]//4, 4))

plt.figure(dpi=200)
plt.plot(y[:,0], y[:,1])
plt.savefig('ark_kepler_phase_plot.png')

conserved = np.loadtxt('ark_kepler_conserved_sprk-1.txt', delimiter=',', dtype=np.float64)
energy = conserved[:,0]
energy_0 = conserved[0,0]
L = conserved[:,1]
L_0 = conserved[0,1]

plt.figure(dpi=200)
plt.title('Error in Hamiltonian')
plt.plot(t, np.abs(energy-energy_0))
plt.ylabel('| error |')
plt.xlabel('time')
plt.savefig('ark_kepler_energy.png')
plt.close()

plt.figure(dpi=200)
plt.title('Error in Momentum')
plt.plot(t, np.abs(L-L_0))
plt.ylabel('| error |')
plt.xlabel('time')
plt.savefig('ark_kepler_momentum.png')
plt.close()

#
#  Time plot.
#
plt.figure(dpi=200)
plt.plot(t, y[:,0], linewidth = 2)
plt.plot(t, y[:,1], linewidth = 2)
plt.plot(t, y[:,2], linewidth = 2)
plt.plot(t, y[:,3], linewidth = 2)
plt.grid(True)
plt.xlabel('<---  t  --->')
plt.ylabel('<---  y(1:4)  --->')
plt.title('ark_kepler: Time Plot')
plt.savefig('ark_kepler_plot.png')
plt.close()

#
#  Phase plot.
#
plt.figure(dpi=200)
plt.plot(y[:,0], y[:,1], linewidth = 2)
plt.grid(True)
plt.xlabel('<---  y1  --->')
plt.ylabel('<---  y2  --->')
plt.title('ark_kepler: Phase Plot')
plt.savefig('ark_kepler_phase.png')
plt.show(block = False)
plt.close()