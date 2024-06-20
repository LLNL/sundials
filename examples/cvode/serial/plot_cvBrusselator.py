#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2024, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# -----------------------------------------------------------------------------

# imports
import argparse
import numpy as np
import matplotlib.pyplot as plt
import sys
import pprint

sys.path.append("/Volumes/Workspace/SUNDIALS/repos/feature/nls-switching-v7/scripts")
from sundialsdev.logs import *

# command line options
parser = argparse.ArgumentParser(description='Plots cvBrusselator output')
parser.add_argument('sfile', type=str,
                    help='solution output file to read')
parser.add_argument('dfile', type=str,
                    help='debug output file to read')
parser.add_argument('elefile', type=str,
                    help='est local error output file to read')

# parse inputs
args = parser.parse_args()

# read solution output file
data = np.loadtxt(args.sfile, dtype=np.double)

# extract times and individual components
t = data[:, 0]
uvw = data[:, 1:]

# extract the times when the nls switch was made
steps = cvode_debug_file_to_list(args.dfile)
switched_to_fp = []
switched_to_newton = []
for step in steps:
  t_n = np.double(step[0]['payload']['t_n'])
  for entry in step[1:]:
    if 'switch to fixed point' in entry['payload']:
      switched_to_fp.append(t_n)
    elif 'switch to Newton' in entry['payload']:
      switched_to_newton.append(t_n)

# for plotting purposes, find nearest output time to switch
def find_nearest(arr, value):
    arr = np.asarray(arr)
    idx = (np.abs(arr - value)).argmin()
    return idx

nearest_fp_idxs = [find_nearest(t, x) for x in switched_to_fp]
nearest_newt_idxs = [find_nearest(t, x) for x in switched_to_newton]

fpt = t[nearest_fp_idxs]
newtt = t[nearest_newt_idxs]

# plot solution against time
plt.plot(t, uvw)
plt.legend(['u', 'v', 'w'])

# plot when switches were made
for t in fpt:
  plt.axvline(x=t, color='r', linewidth=0.2)
for t in newtt:
  plt.axvline(x=t, color='b', linewidth=0.2)

plt.title('Plot of Brusselator Solution\n blue => switched to newton\n red => switched to fixedpoint')
plt.savefig('cvBrusselator.png', dpi=300)
# plt.show()
plt.close()


# read est local error output file
data = np.loadtxt(args.elefile, dtype=np.double)

t = data[:, 0]
ele = data[:, 1:]

# extract the times when halpha was used
used_halpha = []
for step in steps:
  t_n = np.double(step[0]['payload']['t_n'])
  for entry in step[1:]:
    if 'consider-halpha' in entry['label']:
      halpha = np.double(entry['payload']['halpha'])
      hprime = np.double(entry['payload']['hprime'])
      if halpha < hprime:
        used_halpha.append(t_n)

plt.plot(t, ele)
plt.legend(['u', 'v', 'w'])
# plot when halpha was used
for t in used_halpha:
  plt.axvline(x=t, color='r', linewidth=0.2)
plt.title("Plot of Est. Local Errors\n red => halpha used for next step")
plt.savefig('cvBrusselator_ele.png', dpi=300)
# plt.show()
plt.close()

##### end of script #####
