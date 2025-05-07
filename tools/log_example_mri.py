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
# Example script demonstrating how to use Python functions to extract and plot
# logs produced by the SUNLogger with an MRI method.
# -----------------------------------------------------------------------------


def main():

    import argparse
    import matplotlib.pyplot as plt

    from suntools import logs as sunlog

    parser = argparse.ArgumentParser(description="Plots")

    parser.add_argument("logfile", type=str, help="Log file to plot")

    parser.add_argument("--step-number", action="store_true", help="Plot value vs step number")

    parser.add_argument(
        "--step-range",
        type=int,
        nargs=2,
        default=None,
        metavar=("LOWER_BOUND", "UPPER_BOUND"),
        help="Step range to plot",
    )

    parser.add_argument(
        "--time-range",
        type=float,
        nargs=2,
        default=None,
        metavar=("LOWER_BOUND", "UPPER_BOUND"),
        help="Time range to plot",
    )

    parser.add_argument(
        "--save",
        type=str,
        nargs="?",
        const="fig.pdf",
        default=None,
        metavar="FILE_NAME",
        help="""Save figure to file""",
    )

    # parse command line args
    args = parser.parse_args()

    # parse log file
    log = sunlog.log_file_to_list(args.logfile)

    # plot log data
    steps_s, times_s, vals_s = sunlog.get_history(
        log, "h", step_range=args.step_range, time_range=args.time_range, group_by_level=True
    )

    # get data for error test failures
    steps_etf, times_etf, vals_etf = sunlog.get_history(
        log, "h", "failed error test", step_range=args.step_range, time_range=args.time_range
    )

    # get data for solver failures
    steps_sf, times_sf, vals_sf = sunlog.get_history(
        log, "h", "failed solve", step_range=args.step_range, time_range=args.time_range
    )

    if args.step_number:
        x_s = steps_s
        x_etf = steps_etf
        x_sf = steps_sf
    else:
        x_s = times_s
        x_etf = times_etf
        x_sf = times_sf

    fig, ax = plt.subplots()
    colors = plt.get_cmap("tab10")

    for idx, level in enumerate(x_s):
        ax.plot(
            x_s[level],
            vals_s[level],
            color=colors(idx),
            marker=".",
            zorder=0.1,
            label=f"level {level} successful",
        )

    ax.scatter(x_etf, vals_etf, color="red", marker="x", zorder=0.2, label="error test failed")
    ax.scatter(x_sf, vals_sf, color="darkorange", marker="d", zorder=0.2, label="solver failed")

    if args.step_number:
        ax.set_xlabel("step")
    else:
        ax.set_xlabel("time")

    ax.set_ylabel("step size")
    ax.legend(loc="best")
    ax.grid(alpha=0.3, linestyle="--")

    if args.save:
        plt.savefig(args.save, bbox_inches="tight")
    else:
        plt.show()


# run the main routine
if __name__ == "__main__":
    import sys

    sys.exit(main())
