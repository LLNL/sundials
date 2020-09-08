#!/usr/bin/env python
# ------------------------------------------------------------------------------
# Programmer(s): David J. Gardner @ LLNL
# ------------------------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2020, Lawrence Livermore National Security
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
#   $ ./updateOutFiles.py sundials/examples sundials/build/Testing
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# main routine
# ------------------------------------------------------------------------------
def main():

    import argparse
    import sys, os, shutil

    parser = argparse.ArgumentParser(
        description='Update output files',
        formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('examples',type=str,
                        help='Full path to SUNDIALS examples directory')
    parser.add_argument('testing',type=str,
                        help='Full path to CTest Testing directory')

    # parse command line args
    args = parser.parse_args()

    # check inputs
    if (not os.path.isdir(args.examples)):
        print("Error: could not find {}".format(args.examples))
        return -1

    if (not os.path.isdir(args.testing)):
        print("Error: could not find {}".format(args.testing))
        return -1

    # check that the output directory exists
    output = os.path.join(args.testing, "output")
    if (not os.path.isdir(output)):
        print("Error: could not find {}".format(output))
        return -1

    # get list of failed tests
    failed = os.path.join(args.testing, "Temporary", "LastTestsFailed.log")

    # extract test names from list and append .out
    tests = []
    with open(failed, 'r') as f:
        for line in f:
            tests.append(line.split(':')[1].rstrip() + ".out")

    # if failed tests were found update the output files
    if tests:
        for t in tests:
            found = False
            for root, dirs, files in os.walk(args.examples):
                if t in files:
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
