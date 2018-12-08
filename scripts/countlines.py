#!/usr/bin/env python
# ------------------------------------------------------------------------------
# Programmer(s): David J. Gardner @ LLNL
# ------------------------------------------------------------------------------
# LLNS Copyright Start
# Copyright (c) 2014, Lawrence Livermore National Security
# This work was performed under the auspices of the U.S. Department
# of Energy by Lawrence Livermore National Laboratory in part under
# Contract W-7405-Eng-48 and in part under Contract DE-AC52-07NA27344.
# Produced at the Lawrence Livermore National Laboratory.
# All rights reserved.
# For details, see the LICENSE file.
# LLNS Copyright End
# ------------------------------------------------------------------------------
# Find unique lines from a file containing a specific string and write them to
# a separate file.
#
# git diff --numstat | awk '{ added += $1; removed += $2 } END { print "+" added " -" removed }'
#
# Example usage:
#   Find unique compiler warnings in make.log and write them to warnings.txt:
#   $ ./findlines.py make.log warning warnings.txt
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# main routine
# ------------------------------------------------------------------------------
def main():

    import argparse
    import sys, os
    import shlex
    from subprocess import call

    parser = argparse.ArgumentParser(
        description='Count the number lines changed in a certain path',
        formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('key1',type=str,
                        help='First git identifier')

    parser.add_argument('key2',type=str,
                        help='Second git identifier')

    parser.add_argument('patterns',type=str,nargs='+',
                        help='File path patterns to look for in diff')


    # parse command line args
    args = parser.parse_args()

    # git stats from git
    cmd = 'git diff --numstat '+args.key1+' '+args.key2+' > tmp.txt'
    print cmd
    call(cmd, shell=True)

    totaladded   = 0
    totalremoved = 0
    netchange    = 0

    with open('tmp.txt',"r") as fr:
        for line in fr:
            split_line = shlex.split(line)
            added    = split_line[0]
            removed  = split_line[1]
            filename = split_line[2]
            if (any(x in filename for x in args.patterns)):
                if ((added == '-') or (removed == '-')):
                    continue
                else:
                    print added, removed, filename
                    totaladded   += int(added)
                    totalremoved += int(removed)

    netchange = totaladded - totalremoved
    print totaladded, totalremoved, netchange

# ------------------------------------------------------------------------------
# run the main routine
# ------------------------------------------------------------------------------
if __name__ == '__main__':
    import sys
    sys.exit(main())
