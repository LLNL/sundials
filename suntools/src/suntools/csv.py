#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# Programmer(s): David J. Gardner @ LLNL
# -----------------------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2025-2026, Lawrence Livermore National Security,
# University of Maryland Baltimore County, and the SUNDIALS contributors.
# Copyright (c) 2013-2025, Lawrence Livermore National Security
# and Southern Methodist University.
# Copyright (c) 2002-2013, Lawrence Livermore National Security.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# -----------------------------------------------------------------------------
# Functions to parse SUNDIALS CSV output files
# -----------------------------------------------------------------------------

"""Read SUNDIALS statistics written in CSV format.

SUNDIALS CSV output stores a key and value in alternating columns.  The
functions in this module return a dictionary whose values are lists, which is
convenient for comparing several output rows.
"""

from .utils import str2num


def keys(filename):
    """Extract the keys from a SUNDIALS CSV file.

    :param str filename: Path to the SUNDIALS CSV file.
    :returns: The dictionary keys in the first row of the file.
    :rtype: list[str]

    The file is expected to contain a header row in which the key columns are
    interleaved with value columns.  Only the key columns are returned.
    """

    # Get keys from the first row
    with open(filename, mode="r") as csvfile:
        keys = csvfile.readline().split(",")[::2]

    return keys


def read(filename):
    """Read a SUNDIALS CSV file into a dictionary.

    :param str filename: Path to the SUNDIALS CSV file.
    :returns: A dictionary mapping each key to a list of values. Numeric
              values are converted to :class:`int` or :class:`float` when
              possible.
    :rtype: dict[str, list]

    The output has one list per key.  Values are read from alternating columns
    and converted with :func:`suntools.utils.str2num`.
    """

    import csv

    # Get dictionary keys
    fields = keys(filename)

    # Initialize dictionary
    csv_dict = {}
    for k in fields:
        csv_dict[k] = []

    # Get values from each row
    with open(filename, mode="r") as csvfile:
        reader = csv.reader(csvfile)
        for row in reader:
            values = row[1::2]
            for k, v in zip(fields, values):
                csv_dict[k].append(str2num(v))

    return csv_dict


def write(filename):
    """Print a SUNDIALS CSV file as key-value lists.

    :param str filename: Path to the SUNDIALS CSV file.

    This convenience function reads the file with :func:`read` and writes a
    human-readable representation to standard output.  It does not write a
    new CSV file and returns ``None``.
    """

    csv_dict = read(filename)

    for key in csv_dict.keys():
        print(f"{key:29} = {csv_dict[key]}")
    print()
