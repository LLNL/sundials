#!/usr/bin/env python
# ------------------------------------------------------------------------------
# Programmer(s): David J. Gardner @ LLNL
# ------------------------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2021, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# ------------------------------------------------------------------------------
# Script to update example output files for failed tests
#
# Example usage:
#   $ ./updateOutFiles.py sundials/build sundials/examples
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# main routine
# ------------------------------------------------------------------------------
def main():

    import argparse
    import sys, os, shutil

    parser = argparse.ArgumentParser(description='Update output files')

    parser.add_argument('build', type=str,
                        help='Full path to build directory to read from')
    parser.add_argument('examples', type=str,
                        help='Full path to examples directory to write to')
    parser.add_argument('--all','-a', action='store_true',
                        help='Update all output files')
    parser.add_argument('--verbose','-v', action='store_true',
                        help='Enable verbose output')


    # parse command line args
    args = parser.parse_args()

    # check inputs
    if (not os.path.isdir(args.examples)):
        print("Error: could not find {}".format(args.examples))
        return -1

    if (not os.path.isdir(args.build)):
        print("Error: could not find {}".format(args.build))
        return -1

    # check that the output directory exists
    output = os.path.join(args.build, "Testing", "output")
    if (not os.path.isdir(output)):
        print("Error: could not find {}".format(output))
        return -1

    # create a list of all test run or just the failed tests
    tests = []

    if (args.all):
        # get names of all .out files
        for f in os.listdir(output):
            if f.endswith(".out"):
                tests.append(f)
    else:
        failed = os.path.join(args.build, "Testing", "Temporary", "LastTestsFailed.log")

        # extract test names from list and append .out
        with open(failed, 'r') as f:
            for line in f:
                tests.append(line.split(':')[1].rstrip() + ".out")

    # if failed tests were found update the output files
    if tests:
        for t in tests:
            found = False
            for root, dirs, files in os.walk(args.examples):
                if t in files:
                    if (args.verbose):
                        print("Updating: {}".format(t))
                    shutil.copy(os.path.join(output, t), os.path.join(root, t))
                    found = True
                    break
            if not found:
                print("Warning: did not find {}".format(t))

# ------------------------------------------------------------------------------
# run the main routine
# ------------------------------------------------------------------------------
if __name__ == '__main__':
    import sys
    sys.exit(main())
