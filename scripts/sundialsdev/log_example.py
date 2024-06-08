#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# Programmer(s): David J. Gardner @ LLNL
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
# Example script demonstrating how to use Python functions to extract and plot
# logs produced by SUNLogger.
# -----------------------------------------------------------------------------

def main():

    import argparse
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib.ticker as tik

    import logs as sunlog

    parser = argparse.ArgumentParser(description='Plots')

    parser.add_argument('logfile', type=str,
                        help='Log file to plot')

    parser.add_argument('--val', type=str, default='h',
                        choices=['h', 'q', 'dsm'],
                        help='Value to plot')

    parser.add_argument('--step-number', action='store_true',
                        help='Plot value vs step number')

    parser.add_argument('--step-range', type=int, nargs=2,
                        default=None, help='Step range to plot')

    parser.add_argument('--time-range', type=float, nargs=2,
                        default=None, help='Time range to plot')

    parser.add_argument('--save', type=str, nargs='?',
                        const='fig.pdf', default=None,
                        help='''Save figure to file''')

    # parse command line args
    args = parser.parse_args()

    # parse log file
    log = sunlog.cvode_debug_file_to_list(args.logfile)

    # plot log data
    steps_s, times_s, vals_s = sunlog.get_history(log, args.val,
                                                  step_range=args.step_range,
                                                  time_range=args.time_range)

    steps_f, times_f, vals_f = sunlog.get_history(log, args.val, 'failed',
                                                  step_range=args.step_range,
                                                  time_range=args.time_range)

    if args.step_number:
        x_s = steps_s
        x_f = steps_f
    else:
        x_s = times_s
        x_f = times_f

    fig, ax = plt.subplots()

    ax.scatter(x_s, vals_s, color='green', marker='o', label='successful',
               zorder=0.1)

    ax.scatter(x_f, vals_f, color='red', marker='x', label='failed',
               zorder=0.2)

    if args.step_number:
        ax.set_xlabel("step")
    else:
        ax.set_xlabel("time")

    if args.val == 'h':
        ax.set_ylabel("step size")
    elif args.val == 'q':
        ax.set_ylabel("order")
        ax.yaxis.set_major_locator(tik.MaxNLocator(integer=True))
    elif args.val == 'dsm':
        ax.set_ylabel("LTE estimate")

    ax.legend(loc='best')

    ax.grid(alpha=0.3, linestyle='--')

    if args.save:
        plt.savefig(args.save, bbox_inches='tight')
    else:
        plt.show()

# -----------------------------------------------------------------------------
# run the main routine
# -----------------------------------------------------------------------------

if __name__ == '__main__':
    import sys
    sys.exit(main())
