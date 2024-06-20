#!/bin/bash
# ---------------------------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2024, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# ---------------------------------------------------------------------------------
# This script will use clang-tidy and clang-format to format C/C++ code and
# fprettify for Fortran code.
#
# Usage:
#    ./format.sh <paths to directories or files to format>
#
# We require clang-format 17.0.4. Other versions may produce different styles!
# ---------------------------------------------------------------------------------

if [ $# -lt 1 ]; then
    echo "ERROR: At least one path to format required"
    exit 1
fi

paths=( "$@" )

find "${paths[@]}" -iname '*.h' -o -iname '*.hpp' -o \
  -iname '*.c' -o -iname '*.cpp' -o \
  -iname '*.cuh' -o -iname '*.cu' | grep -v fmod | xargs clang-format -i

find "${paths[@]}" -iname '*.f90' | grep -v fmod | xargs fprettify --indent 2 --enable-replacements --c-relations
