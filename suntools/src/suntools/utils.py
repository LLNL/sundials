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
# Shared suntools utilities.
# -----------------------------------------------------------------------------

"""Small helpers shared by the :mod:`suntools` parsers."""


def str2num(s):
    """Convert a numeric string to an integer or floating-point value.

    :param str s: String to convert.
    :returns: An :class:`int` if ``s`` is an integer, a :class:`float` if it is
              a floating-point value, or the original string otherwise.
    :rtype: int, float, or str
    """

    try:
        return int(s)
    except ValueError:
        try:
            return float(s)
        except ValueError:
            return s
