#!/usr/bin/env python3
# ------------------------------------------------------------------------------
# Programmer(s):  Daniel R. Reynolds @ SMU
# ------------------------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2024, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# ------------------------------------------------------------------------------
# matplotlib-based plotting script for the serial ark_sod_lsrk example
# ------------------------------------------------------------------------------

# imports
import sys, os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

# data file name
datafile = "sod.out"

# return with an error if the file does not exist
if not os.path.isfile(datafile):
    msg = "Error: file " + datafile + " does not exist"
    sys.exit(msg)

# read solution file, storing each line as a string in a list
with open(datafile, "r") as file:
    lines = file.readlines()

    # extract header information
    title = lines.pop(0)
    nvar = int((lines.pop(0).split())[2])
    varnames = lines.pop(0)
    nt = int((lines.pop(0).split())[2])
    nx = int((lines.pop(0).split())[2])
    xl = float((lines.pop(0).split())[2])
    xr = float((lines.pop(0).split())[2])

    # allocate solution data as 2D Python arrays
    t = np.zeros((nt), dtype=float)
    rho = np.zeros((nt, nx), dtype=float)
    mx = np.zeros((nt, nx), dtype=float)
    my = np.zeros((nt, nx), dtype=float)
    mz = np.zeros((nt, nx), dtype=float)
    et = np.zeros((nt, nx), dtype=float)

    # store remaining data into numpy arrays
    for it in range(nt):
        line = (lines.pop(0)).split()
        t[it] = line.pop(0)
        for ix in range(nx):
            rho[it, ix] = line.pop(0)
            mx[it, ix] = line.pop(0)
            my[it, ix] = line.pop(0)
            mz[it, ix] = line.pop(0)
            et[it, ix] = line.pop(0)

    gamma = 1.4
    u = mx / rho
    p = (gamma - 1.0) * (et - (mx * mx + my * my + mz * mz) / (2.0 * rho))

# generate plots
x = np.linspace(xl, xr, nx)

#   plot defaults: increase default font size, increase plot width, enable LaTeX rendering
plt.rc("font", size=15)
plt.rcParams["figure.figsize"] = [7.2, 4.8]
plt.rcParams["text.usetex"] = True
plt.rcParams["figure.constrained_layout.use"] = True

#   subplots with time snapshots of the density, x-velocity, and pressure
fig = plt.figure(figsize=(10, 5))
gs = GridSpec(3, 3, figure=fig)
ax00 = fig.add_subplot(gs[0, 0])  # left column
ax10 = fig.add_subplot(gs[1, 0])
ax20 = fig.add_subplot(gs[2, 0])
ax01 = fig.add_subplot(gs[0, 1])  # middle column
ax11 = fig.add_subplot(gs[1, 1])
ax21 = fig.add_subplot(gs[2, 1])
ax02 = fig.add_subplot(gs[0, 2])  # right column
ax12 = fig.add_subplot(gs[1, 2])
ax22 = fig.add_subplot(gs[2, 2])
it = 0
ax00.plot(x, rho[it, :])
ax10.plot(x, u[it, :])
ax20.plot(x, p[it, :])
ax00.set_title(r"$t =$ " + repr(t[it]).zfill(3))
ax00.set_ylabel(r"$\rho$")
ax10.set_ylabel(r"$v_x$")
ax20.set_ylabel(r"$p$")
ax20.set_xlabel(r"$x$")
it = nt // 2
ax01.plot(x, rho[it, :])
ax11.plot(x, u[it, :])
ax21.plot(x, p[it, :])
ax01.set_title(r"$t =$ " + repr(t[it]).zfill(3))
ax21.set_xlabel(r"$x$")
it = nt - 1
ax02.plot(x, rho[it, :])
ax12.plot(x, u[it, :])
ax22.plot(x, p[it, :])
ax02.set_title(r"$t =$ " + repr(t[it]).zfill(3))
ax22.set_xlabel(r"$x$")
plt.savefig("sod_frames.png")

plt.show()

##### end of script #####
