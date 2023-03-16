#!/usr/bin/env python
# ----------------------------------------------------------------
# Programmer(s): Cody J. Balos @ LLNL
# ----------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2022, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# ----------------------------------------------------------------
# matplotlib-based plotting script for ark_kepler.c
# ----------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt

t = np.loadtxt('ark_kepler_times_sprk-2.txt', dtype=np.float64)
y = np.loadtxt('ark_kepler_solution_sprk-2.txt', dtype=np.float64)
y = np.reshape(y, (y.shape[0]//4, 4))

plt.figure(dpi=200)
plt.plot(y[:,0], y[:,1])
plt.savefig('ark_kepler_phase.png')
plt.close()

conserved = np.loadtxt('ark_kepler_conserved_sprk-2.txt', delimiter=',', dtype=np.float64)
energy = conserved[:,0]
energy_0 = conserved[0,0]
L = conserved[:,1]
L_0 = conserved[0,1]

plt.figure(dpi=200)
plt.title('Error in Hamiltonian')
plt.plot(t, np.abs(energy-energy_0))
plt.ylabel('| error |')
plt.xlabel('<---  t  --->')
plt.xscale('log')
plt.yscale('log')
plt.savefig('ark_kepler_energy.png')
plt.close()

plt.figure(dpi=200)
plt.title('Error in Momentum')
plt.plot(t, np.abs(L-L_0))
plt.ylabel('| error |')
plt.xlabel('<---  t  --->')
plt.xscale('log')
plt.yscale('log')
plt.savefig('ark_kepler_momentum.png')
plt.close()

# # Step history
# hhist = np.loadtxt('ark_kepler_hhist_sprk-2.txt', dtype=np.float64)
# plt.figure(dpi=200)
# plt.title('Step Size History')
# plt.plot(t[:-1], hhist)
# plt.xlabel('<---  t  --->')
# plt.yscale('log')
# plt.savefig('ark_kepler_hhist.png')
# plt.close()

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
plt.close()
