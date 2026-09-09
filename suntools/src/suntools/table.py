#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# Programmer(s): Cody J. Balos @ LLNL
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
# Function to parse SUNDIALS statistics in table form
# -----------------------------------------------------------------------------

"""Read SUNDIALS statistics written in table format."""

import re
from .utils import str2num


def parse_table(data):
    """Parse table-formatted SUNDIALS statistics.

    :param str data: Multiline text containing one or more key-value tables.
    :returns: A dictionary mapping each key to a list of values. Numeric values
              are converted to :class:`int` or :class:`float` when possible.
    :rtype: dict[str, list]

    Each matching line must have the form ``key = value``. Repeated keys are
    accumulated in insertion order, giving the same dictionary shape as
    :func:`suntools.csv.read`.
    """

    stats = {}
    lines = data.splitlines()
    regexp = re.compile(r"^(.*?)\s+=\s+(.*)$")
    for line in lines:
        # Extract key-value pairs
        match = regexp.match(line)
        if match:
            key = match.group(1).strip()
            value = match.group(2).strip()

            # Attempt to convert value to a float or int if numeric
            value = str2num(value)

            if key in stats:
                stats[key].append(value)
            else:
                stats[key] = [value]

    return stats
