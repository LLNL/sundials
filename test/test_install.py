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
# Test all installed Makefiles or CMakeLists.txt files in a given directory
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# main routine
# ------------------------------------------------------------------------------
def main():

    import argparse
    import fnmatch
    import os
    from subprocess import call

    parser = argparse.ArgumentParser(
        description='Find and test installed build files in a given directory',
        formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('directory', type=str,
                        help='Directory to search for makefiles')

    parser.add_argument('--cmake', action='store_true',
                        help='Test CMake build')

    parser.add_argument('--debug', action='store_true',
                        help='Enable debugging output')

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

    # output list of makefiles
    if args.debug:
        for bf in buildfiles:
            print bf

    # run make for each makefile
    buildfail = False
    errors = []
    for bf in buildfiles:

        # move to example directory
        os.chdir(os.path.dirname(bf))

        # confgure cmake if necessary
        if args.cmake:
            ret = call('cmake -DCMAKE_VERBOSE_MAKEFILE=ON .', shell=True)
            if ret != 0:
                errors.append(os.path.dirname(bf))
                buildfail = True
                os.chdir(cwd)
                continue

        # make examples
        ret = call('make', shell=True)
        if ret != 0:
            errors.append(os.path.dirname(bf))
            buildfail = True

        # return to original directory
        os.chdir(cwd)

        # stop on first error when debugging
        if args.debug:
            if buildfail:
                break

    # print list of failed builds
    if buildfail:
        print "The following builds failed:"
        for err in errors:
            print err
        sys.exit(1)
    else:
        print "All builds successful."

# ------------------------------------------------------------------------------
# run the main routine
# ------------------------------------------------------------------------------
if __name__ == '__main__':
    import sys
    sys.exit(main())
