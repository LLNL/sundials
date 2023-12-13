#!/bin/bash
# ---------------------------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2023, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# ---------------------------------------------------------------------------------
# This script will use clang-tidy and clang-format to format code.
#
# Usage:
#    ./format.sh <path to directory to format>
#
# We require clang-format 17.0.4. Other versions may produce different styles!
# ---------------------------------------------------------------------------------

if [ "$#" -ne 1 ]; then
    echo "ERROR: One input required"
    echo "1) Path to directory to format"
    exit 1
fi

find $1 -iname '*.h' -o -iname '*.hpp' -o \
  -iname '*.c' -o -iname '*.cpp' -o \
  -iname '*.cuh' -o -iname '*.cu' | grep -v fmod | xargs clang-format -i
