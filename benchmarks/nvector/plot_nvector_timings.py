#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# Programmer(s): David J. Gardner @ LLNL
# -----------------------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2022, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# -----------------------------------------------------------------------------
# This script plot data from the vector performance tests
# -----------------------------------------------------------------------------

def main():

    import argparse
    import os
    import sys
    import shlex
    import matplotlib.pyplot as plt

    parser = argparse.ArgumentParser(
        description='Plot data from NVector performance tests')

    parser.add_argument('operation', type=str,
                        choices=['N_VLinearSum-1a',
                                 'N_VLinearSum-1b',
                                 'N_VLinearSum-1c',
                                 'N_VLinearSum-2a',
                                 'N_VLinearSum-2b',
                                 'N_VLinearSum-2c',
                                 'N_VLinearSum-3',
                                 'N_VLinearSum-4a',
                                 'N_VLinearSum-4b',
                                 'N_VLinearSum-5a',
                                 'N_VLinearSum-5b',
                                 'N_VLinearSum-6a',
                                 'N_VLinearSum-6b',
                                 'N_VLinearSum-7',
                                 'N_VLinearSum-8',
                                 'N_VLinearSum-9',
                                 'N_VConst',
                                 'N_VProd',
                                 'N_VDiv',
                                 'N_VScale-1',
                                 'N_VScale-2',
                                 'N_VScale-3',
                                 'N_VScale-4',
                                 'N_VAbs',
                                 'N_VInv',
                                 'N_VAddConst',
                                 'N_VDotProd',
                                 'N_VMaxNorm',
                                 'N_VWrmsNorm',
                                 'N_VWrmsNormMask',
                                 'N_VMin',
                                 'N_VWL2Norm',
                                 'N_VL1Norm',
                                 'N_VCompare',
                                 'N_VInvTest',
                                 'N_VConstrMask',
                                 'N_VMinQuotient',
                                 'N_VLinearCombination-1',
                                 'N_VLinearCombination-2',
                                 'N_VLinearCombination-3',
                                 'N_VScaleAddMulti-1',
                                 'N_VScaleAddMulti-2',
                                 'N_VDotProdMulti',
                                 'N_VLinearSumVectorArray',
                                 'N_VScaleVectorArray-1',
                                 'N_VScaleVectorArray-2',
                                 'N_VConstVectorArray',
                                 'N_VWrmsNormVectorArray',
                                 'N_VWrmsNormMaskVectorArray',
                                 'N_VScaleAddMultiVectorArray-1',
                                 'N_VScaleAddMultiVectorArray-2',
                                 'N_VLinearCombinationVectorArray-1',
                                 'N_VLinearCombinationVectorArray-2',
                                 'N_VLinearCombinationVectorArray-3'],
                        help='Which NVector operation to plot')

    parser.add_argument('datafiles', type=str, nargs='+',
                        help='Data files to plot')

    parser.add_argument('-d', '--debug', action='count', default=0,
                        help='Enable debugging output')

    # parse command line args
    args = parser.parse_args()

    if args.debug > 0:
        print(args)

    # check that data files exist
    for f in args.datafiles:
        if not os.path.isfile(f):
            print(f"ERROR: {f} does not exist")
            sys.exit()

    # list of dictionaries with test parameters and data
    tests = [{} for x in range(len(args.datafiles))]

    # dictionary keys
    keys = ['vector', 'length', 'vectors', 'sums', 'tests', 'avg', 'sdev',
            'min', 'max']

    for d in tests:
        for k in keys:
            d[k] = None

    # parse file names to get input parameters
    for i, f in enumerate(args.datafiles):

        if args.debug > 1:
            print(f"Reading file {i}: {f}")

        with open(f) as fout:
            for line in fout:

                # skip empty lines
                if not line.strip():
                    continue

                # split line into list
                split_line = shlex.split(line)

                if args.debug > 1:
                    print(line)

                # get test parameters
                if "Vector Name" in line:
                    tests[i]["vector"] = split_line[-1]

                if "vector length" in line:
                    tests[i]["length"] = int(split_line[-1])

                if "max number of vectors" in line:
                    tests[i]["vectors"] = int(split_line[-1])

                if "max number of sums" in line:
                    tests[i]["sums"] = int(split_line[-1])

                if "number of tests" in line:
                    tests[i]["tests"] = int(split_line[-1])

                # extract timing data for the desired operation
                if args.operation == split_line[0]:
                    tests[i]["avg"] = float(split_line[1])
                    tests[i]["sdev"] = float(split_line[2])
                    tests[i]["min"] = float(split_line[3])
                    tests[i]["max"] = float(split_line[4])

                    if args.debug > 1:
                        print(tests[i])

                    break

    if args.debug > 0:
        for d in tests:
            print(d)

    # get vector types to plot
    vectors = []
    for d in tests:
        if d["vector"] not in vectors:
            vectors.append(d["vector"])

    if args.debug > 0:
        print(vectors)

    # plot data for each vector type
    fig, ax = plt.subplots()

    for v in vectors:

        # extract vector data
        x = []
        y = []
        for d in tests:
            if v == d["vector"]:
                x.append(d["length"])
                y.append(d["avg"])

        # sort data by x values
        x, y = (list(i) for i in zip(*sorted(zip(x, y))))

        if args.debug > 0:
            print(x)
            print(y)

        ax.plot(x, y, label=v)

    plt.show()


# =============================================================================


if __name__ == "__main__":
    main()
