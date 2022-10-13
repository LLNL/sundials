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
# Test all installed Makefiles or CMakeLists.txt files under a given directory
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# main routine
# -----------------------------------------------------------------------------

def main():

    import argparse
    import fnmatch
    import os
    import re
    import subprocess

    parser = argparse.ArgumentParser(
        description='Find and test installed build files in a given directory',
        formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('directory', type=str,
                        help='Directory to search for build files')

    parser.add_argument('--cmake', action='store_true',
                        help='Test CMake build')

    parser.add_argument('--regex', type=str,
                        help='Regular expression for filtering example directories')

    parser.add_argument('-v', '--verbose', action='count', default=0,
                        help='Verbose output')

    parser.add_argument('--failfast', action='store_true',
                        help='Stop on first failure')

    # parse command line args
    args = parser.parse_args()

    # get current directory
    cwd = os.getcwd()

    # which build file to search for
    if args.cmake:
        filename = "CMakeLists.txt"
    else:
        filename = "Makefile"

    # get list of build files
    buildfiles = []
    for root, dirs, files in os.walk(args.directory):
        for fn in fnmatch.filter(files, filename):
            buildfiles.append(os.path.join(root, fn))

    # output list of build files
    if args.verbose > 0:
        print(f"Total files: {len(buildfiles)}")
    if args.verbose > 1:
        for bf in buildfiles:
            print(bf)

    # filter files
    if args.regex:
        regex = re.compile(args.regex)
        buildfiles = [ bf for bf in buildfiles if re.search(regex, bf) ]
        if args.verbose > 0:
            print(f"Total files (filtered): {len(buildfiles)}")
        if args.verbose > 1:
            for bf in buildfiles:
                print(bf)

    # total files to test
    total = len(buildfiles)

    # build examples
    buildfail = False
    errors = []
    for i, bf in enumerate(buildfiles):

        print(f"[{i+1} of {total}] Building: {bf}")

        # move to example directory
        os.chdir(os.path.dirname(bf))

        # confgure cmake if necessary
        configfail=False
        if args.cmake:
            ret = subprocess.call('cmake -DCMAKE_VERBOSE_MAKEFILE=ON .',
                                  shell=True, stdout=subprocess.DEVNULL,
                                  stderr=subprocess.DEVNULL)
            if ret != 0:
                errors.append(os.path.dirname(bf))
                configfail = True
                buildfail = True

        # make examples
        if not configfail:
            ret = subprocess.call('make', shell=True, stdout=subprocess.DEVNULL,
                                  stderr=subprocess.DEVNULL)
            if ret != 0:
                errors.append(os.path.dirname(bf))
                buildfail = True

        # return to original directory
        os.chdir(cwd)

        # stop on first error
        if args.failfast and buildfail:
            break

    # print list of failed builds
    if buildfail:
        print("The following builds failed:")
        for err in errors:
            print(err)
        sys.exit(1)
    else:
        print("All builds successful.")

# -----------------------------------------------------------------------------
# run the main routine
# -----------------------------------------------------------------------------

if __name__ == '__main__':
    import sys
    sys.exit(main())
