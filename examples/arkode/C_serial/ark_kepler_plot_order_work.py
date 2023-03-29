#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

def load_results(case):
   t = np.loadtxt('ark_kepler_times' + case + '.txt', dtype=np.float64)
   y = np.loadtxt('ark_kepler_solution' + case + '.txt', dtype=np.float64)
   y = np.reshape(y, (y.shape[0]//4, 4))
   conserved = np.loadtxt('ark_kepler_conserved' + case + '.txt', delimiter=',', dtype=np.float64)
   return t, y, conserved

sprk1 = [] # Euler
sprk2 = [] # Leapforg
sprk22 = [] # Pseudo Leapfrog
sprk222 = [] # McLachlan 2nd order
sprk3 = [] # Ruth 3rd order
sprk33 = [] # McLachlan 3rd order
sprk4 = [] # Candy Rozmus 4th order
sprk44 = [] # McLachlan 4th order
sprk5 = [] # McLachlan 5th order
sprk6 = [] # Yoshida 6th order
sprk8 = [] # McLachlan 8th order
sprk10 = [] # Sofroniou 10th order

# step_sizes = [0.000010, 0.000100, 0.001000, 0.010000, 0.100000]
step_sizes = [0.000100, 0.001000, 0.010000, 0.100000]
for dt in step_sizes:
  sprk1.append({
     'method': 'Symplectic Euler',
     'order': 1,
     'data': load_results('_sprk-1-dt-%.6f' % dt),
     'dt': dt
  })
  sprk2.append({
     'method': 'Leapfrog',
     'order': 2,
     'data': load_results('_sprk-2-dt-%.6f' % dt),
     'dt': dt
  })
  sprk22.append({
     'method': 'Pseudo Leapfrog',
     'order': 2,
     'data': load_results('_sprk-22-dt-%.6f' % dt),
     'dt': dt
  })
  sprk222.append({
     'method': 'McLachlan2',
     'order': 2,
     'data': load_results('_sprk-222-dt-%.6f' % dt),
     'dt': dt
  })
  sprk3.append({
     'method': 'Ruth3',
     'order': 3,
     'data': load_results('_sprk-3-dt-%.6f' % dt),
     'dt': dt
  })
  sprk33.append({
     'method': 'McLachlan3',
     'order': 3,
     'data': load_results('_sprk-33-dt-%.6f' % dt),
     'dt': dt
  })
  sprk4.append({
     'method': 'CandyRozmus4',
     'order': 4,
     'data': load_results('_sprk-4-dt-%.6f' % dt),
     'dt': dt
  })
  sprk44.append({
     'method': 'McLachlan4',
     'order': 4,
     'data': load_results('_sprk-44-dt-%.6f' % dt),
     'dt': dt
  })
  sprk5.append({
     'method': 'McLachlan5',
     'order': 5,
     'data': load_results('_sprk-5-dt-%.6f' % dt),
     'dt': dt
  })
  sprk6.append({
     'method': 'Yoshida6',
     'order': 6,
     'data': load_results('_sprk-6-dt-%.6f' % dt),
     'dt': dt
  })
  sprk8.append({
     'method': 'McLachlan8',
     'order': 8,
     'data': load_results('_sprk-8-dt-%.6f' % dt),
     'dt': dt
  })
  sprk10.append({
     'method': 'Sofroniou10',
     'order': 10,
     'data': load_results('_sprk-10-dt-%.6f' % dt),
     'dt': dt
  })

orders = [1, 2, 3, 4, 5, 6, 8, 10]
all_methods = {}
if 1 in orders:
   all_methods[1] = []
   all_methods[1].append(sprk1)
if 2 in orders:
   all_methods[2] = []
   all_methods[2].append(sprk2)
   all_methods[2].append(sprk22)
   all_methods[2].append(sprk222)
if 3 in orders:
   all_methods[3] = []
   all_methods[3].append(sprk3)
   all_methods[3].append(sprk33)
if 4 in orders:
   all_methods[4] = []
   all_methods[4].append(sprk4)
   all_methods[4].append(sprk44)
if 5 in orders:
   all_methods[5] = []
   all_methods[5].append(sprk5)
if 6 in orders:
   all_methods[6] = []
   all_methods[6].append(sprk6)
if 8 in orders:
   all_methods[8] = []
   all_methods[8].append(sprk8)
if 10 in orders:
   all_methods[10] = []
   all_methods[10].append(sprk10)

# Reference solution
_, y_ref, _ = load_results('_erk-8')

#
# Solution error plot
#
for order in orders:
  plt.figure(dpi=200)
  legend = []
  legend.append('$O(h^{%d})$' % order)
  order_line = np.power(step_sizes, order)
  plt.plot(np.array(step_sizes), np.array(order_line).T, 'k--')
  for method in all_methods[order]:
    errs = []
    dts = []
    for d in method:
      _, y, _ = d['data']
      err = np.linalg.norm(y - y_ref, np.inf) / np.linalg.norm(y_ref, np.inf)
      print('method: %15s,  dt = %.6f, err = %g' % (d['method'], d['dt'], err))
      if err >= 10.:
        continue
      else:
         dts.append(d['dt'])
         errs.append(err)
    legend.append(method[0]['method'])
    plt.plot(dts, errs)
  plt.xscale('log')
  plt.yscale('log')
  plt.xlabel('h')
  plt.ylabel('error')
  plt.title('Order plot for $O(h^{%d})$ methods' % order)
  plt.legend(legend)
  plt.savefig('ark_kepler_sol_order%d.png' % order)
  plt.close()

#
# Energy error plot
#
for order in orders:
  plt.figure(dpi=200)
  legend = []
  legend.append('$O(h^{%d})$' % order)
  order_line = np.power(step_sizes, order)
  plt.plot(np.array(step_sizes), np.array(order_line).T, 'k--')
  for method in all_methods[order]:
    errs = []
    dts = []
    for d in method:
      _, y, conserved = d['data']
      energy_0 = conserved[0,0]
      energy = conserved[:,0]
      err = np.linalg.norm(energy-energy_0, np.inf)
      print('method: %15s,  dt = %.6f, energy err = %g' % (d['method'], d['dt'], err))
      if err >= 10.:
        continue
      else:
         dts.append(d['dt'])
         errs.append(err)
    legend.append(method[0]['method'])
    plt.plot(dts, errs)
  plt.xscale('log')
  plt.yscale('log')
  plt.xlabel('h')
  plt.ylabel('error in energy')
  plt.title('Order plot for $O(h^{%d})$ methods' % order)
  plt.legend(legend)
  plt.savefig('ark_kepler_energy_order%d.png' % order)
  plt.close()
