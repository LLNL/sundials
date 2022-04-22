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
# Function to parse SUNDIALS CSV output files
# -----------------------------------------------------------------------------

def num(s):
    """Try to convert a string to an int or float"""

    try:
        return int(s)
    except ValueError:
        try:
            return float(s)
        except ValueError:
            return s


def keys(filename):
    """Extracts keys from a SUNDIALS CSV file

    Parameters
    ----------
    filename : str
        The file location of the SUNDIALS CSV file

    Returns
    -------
    list
        A list of dictionary keys
    """

    # Get keys from the first row
    with open(filename, mode='r') as csvfile:
        keys = csvfile.readline().split(",")[::2]

    return keys


def read(filename):
    """Reads a SUNDIALS CSV file

    Parameters
    ----------
    filename : str
        The file location of the SUNDIALS CSV file

    Returns
    -------
    dict
        A dictionary containing the CSV keys and values
    """

    import csv

    # Get dictionary keys
    fields = keys(filename)

    # Initialize dictionary
    csv_dict = {}
    for k in fields:
        csv_dict[k] = []

    # Get values from each row
    with open(filename, mode='r') as csvfile:
        reader = csv.reader(csvfile)
        for row in reader:
            values = row[1::2]
            for k, v in zip(fields, values):
                csv_dict[k].append(num(v))

    return csv_dict


def write(filename):
    """Prints a SUNDIALS CSV file

    Parameters
    ----------
    filename : str
        The file location of the SUNDIALS CSV file
    """

    csv_dict = read(filename)

    for key in csv_dict.keys():
        print(f"{key:29} = {csv_dict[key]}")
    print()
