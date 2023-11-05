#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# Programmer(s): David J. Gardner @ LLNL
# -----------------------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2023, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# -----------------------------------------------------------------------------
# Script to update example output files for failed tests
#
# Example usage:
#   $ ./updateOutFiles.py ../build ../.
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# main routine
# -----------------------------------------------------------------------------
def main():

    import argparse
    import sys, os, shutil

    parser = argparse.ArgumentParser(description='Update output files')

    parser.add_argument('source', type=str,
                        help='Full path of build directory to read files from')
    parser.add_argument('destination', type=str,
                        help='Full path of sundials location to write files to')
    parser.add_argument('--all','-a', action='store_true',
                        help='Update all output files')
    parser.add_argument('--copy','-c', action='store_true',
                        help='Copy file to destination if not found')
    parser.add_argument('--verbose','-v', action='count', default=0,
                        help='Enable verbose output')


    # parse command line args
    args = parser.parse_args()

    # check inputs
    if (not os.path.isdir(args.source)):
        print(f"Error: could not find {args.source}")
        return -1

    if (not os.path.isdir(args.destination)):
        print(f"Error: could not find {args.destination}")
        return -1

    # check that the output directory exists
    output = os.path.join(args.source, "Testing", "output")
    if (not os.path.isdir(output)):
        print(f"Error: could not find {output}")
        return -1

    # create a list of all test run or just the failed tests
    tests = []

    if (args.all):
        # get names of all .out files
        for f in os.listdir(output):
            if f.endswith(".out"):
                tests.append(f)
    else:
        failed = os.path.join(args.source, "Testing", "Temporary", "LastTestsFailed.log")

        # extract test names from list and append .out
        with open(failed, 'r') as f:
            for line in f:
                tests.append(line.split(':')[1].rstrip() + ".out")

    # if failed tests were found update the output files
    if tests:
        paths = [os.path.join(args.destination, "examples"),
                 os.path.join(args.destination, "test", "unit_tests"),
                 args.destination]
        for t in tests:
            if (args.verbose > 0):
                print(f"Searching for {t}")
            found = False
            for p in paths:
                if (args.verbose == 2):
                    print(f"  Looking in {p}")
                for root, dirs, files in os.walk(p):
                    if (args.verbose == 3):
                        print(f"  Looking in {root}")
                    if t in files:
                        if (args.verbose == 1):
                            print(f"  Found file in {root}")
                        if (args.verbose > 1):
                            print("  Found file")
                        shutil.copy(os.path.join(output, t), os.path.join(root, t))
                        found = True
                        break
                if found:
                    break
            if not found:
                if args.copy:
                    print(f"Warning: did not find {t}, copying to {args.destination}")
                    shutil.copy(os.path.join(output, t), os.path.join(args.destination, t))
                else:
                    print(f"Warning: did not find {t}")

# -----------------------------------------------------------------------------
# run the main routine
# -----------------------------------------------------------------------------
if __name__ == '__main__':
    import sys
    sys.exit(main())
