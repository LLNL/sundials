#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# Programmer(s): David J. Gardner @ LLNL
# -----------------------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2025, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# -----------------------------------------------------------------------------

import argparse
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.animation as animation

from suntools import logs as sunlog

parser = argparse.ArgumentParser(description="Plots")

parser.add_argument("outfile", type=str, help="Output file to plot")
parser.add_argument("logfile", type=str, help="Log file to plot")

parser.add_argument(
    "--save",
    type=str,
    nargs="?",
    const="ani.gif",
    default=None,
    metavar="FILE_NAME",
    help="Save figure to file",
)

# parse command line args
args = parser.parse_args()

# ----------------------
# Solution data and plot
# ----------------------

fig, axes = plt.subplots(2, sharex=True)

sol_data = np.loadtxt(args.outfile)
sol_x = sol_data[:,0]
sol_y = sol_data[:,1:]

sol_plts = []
for i in range(sol_y.shape[1]):
    ln, = axes[0].plot([], [], label=f"Species {i+1}")
    sol_plts.append(ln)

# Pre-set x/y limits
axes[0].set_xlim(sol_x.min(), sol_x.max())
min_y, max_y = sol_y.min(), sol_y.max()
margin = (max_y - min_y) * 0.1
axes[0].set_ylim(min_y - margin, max_y + margin)

# -----------------
# Log data and plot
# -----------------

# parse log file
log = sunlog.log_file_to_list(args.logfile)

# get successful step data
_, pass_x, pass_y = sunlog.get_history(log, "h", "success")
pass_x = np.array(pass_x)
pass_y = np.array(pass_y)

pass_plt = axes[1].scatter([], [], marker='.', color="green", label="Pass")

# get data for error test failures
_, fail_x, fail_y = sunlog.get_history(log, "h", "failed")
fail_x = np.array(fail_x)
fail_y = np.array(fail_y)

fail_plt = axes[1].scatter([], [], marker="X", color="red", label="Fail")

# Pre-set x/y limits
axes[1].set_xlim(sol_x.min(), sol_x.max())
min_y, max_y = np.min(pass_y), np.max(fail_y)
margin = (max_y - min_y) * 0.1
axes[1].set_ylim(min_y - margin, max_y + margin)

# ------------
# Plot options
# ------------

axes[0].grid(alpha=0.3, linestyle="--")
axes[1].grid(alpha=0.3, linestyle="--")

axes[1].set_xlabel("time")
axes[0].set_ylabel("concentration")
axes[1].set_ylabel("step size")
axes[0].set_title("Brusselator")

axes[0].legend(loc="upper right")
axes[1].legend(loc="upper right")

# ------------------
# Animation function
# ------------------

# print(f"sol pts  = {len(sol_x)}")
# print(f"pass pts = {len(pass_x)}")
# print(f"fail pts = {len(fail_x)}")

def update(frame):
    # Current time for this animation frame
    current_time = sol_x[frame]

    # Clear prior plots when repeating animation
    if frame == 0:
        pass_plt.set_offsets(np.empty((0, 2)))
        fail_plt.set_offsets(np.empty((0, 2)))

    # Plot solution
    for idx, ln in enumerate(sol_plts):
        ln.set_data(sol_x[:frame], sol_y[:frame, idx])

    # Plot solution up to current time
    mask_p = pass_x <= current_time
    if np.any(mask_p):
        sdata = np.stack([pass_x[mask_p], pass_y[mask_p]]).T
        pass_plt.set_offsets(sdata)

    # Failed: plot all up to current_time
    mask_f = fail_x <= current_time
    if np.any(mask_f):
        fdata = np.stack([fail_x[mask_f], fail_y[mask_f]]).T
        fail_plt.set_offsets(fdata)

    # print(f"Frame {frame}: current_time={current_time:.3f}, "
    #       f"solution_pts={frame}, success_pts={mask_p.sum()}, failed_pts={mask_f.sum()}")

    return sol_plts, pass_plt, fail_plt

# ----------------
# Create animation
# ----------------

ani = animation.FuncAnimation(fig=fig, func=update, frames=len(sol_x),
                              interval=60, repeat=True)

if args.save:
    print("saving animation...", end="")
    ani.save("ani.gif")
    print("done.")
else:
    plt.show()


# Make the full static plot
fig, axes = plt.subplots(2, sharex=True)
sol_x = sol_data[:,0]
sol_y = sol_data[:,1:]

for i in range(sol_y.shape[1]):
    axes[0].plot(sol_x, sol_y[:,i], label=f"Species {i+1}")

axes[1].scatter(pass_x, pass_y, marker='.', color="green", label="Pass")
axes[1].scatter(fail_x, fail_y, marker="X", color="red", label="Fail")

# Pre-set x/y limits
axes[0].set_xlim(sol_x.min(), sol_x.max())
min_y, max_y = sol_y.min(), sol_y.max()
margin = (max_y - min_y) * 0.1
axes[0].set_ylim(min_y - margin, max_y + margin)

# Pre-set x/y limits
axes[1].set_xlim(sol_x.min(), sol_x.max())
min_y, max_y = np.min(pass_y), np.max(fail_y)
margin = (max_y - min_y) * 0.1
axes[1].set_ylim(min_y - margin, max_y + margin)

axes[0].grid(alpha=0.3, linestyle="--")
axes[1].grid(alpha=0.3, linestyle="--")

axes[1].set_xlabel("time")
axes[0].set_ylabel("concentration")
axes[1].set_ylabel("step size")
axes[0].set_title("Brusselator")

axes[0].legend(loc="upper right")
axes[1].legend(loc="upper right")

if args.save:
    print("saving plot...", end="")
    plt.savefig("fig.pdf", bbox_inches="tight")
    print("done.")
else:
    plt.show()
