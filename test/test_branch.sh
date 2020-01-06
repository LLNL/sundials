#!/bin/bash
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
# SUNDIALS branch testing driver script
#
# Usage: ./test_branch.sh <num build threads>
#
# Inputs (all optional):
#   <num build threads> = number of threads in parallel build
# ------------------------------------------------------------------------------

# Number of threads for parallel build (default 1)
bt=1
if [ "$#" -ge 1 ]; then
    bt=$1
fi

# ======================================================================
# Compile tests (check for C90 compliance and compiler warnings)
# ======================================================================

realtype=( "single" "double" "extended" )
indexsize=( "32" "64" )

for rt in "${realtype[@]}"; do
    for is in "${indexsize[@]}"; do

        # print test header for Jenkins section collapsing
        echo "TEST: ./suntest.sh $rt $is both OFF NONE $bt"

        # run tests
        ./suntest.sh $rt $is both OFF NONE $bt

        # check return flag
        if [ $? -ne 0 ]; then
            echo "FAILED: ./suntest.sh $rt $is both OFF NONE $bt" | tee -a suntest.log
            exit 1
        else
            echo "PASSED"
        fi

    done
done

# ======================================================================
# Test with TPLs and DEVTESTS ON
# ======================================================================

realtype=( "double" )
indexsize=( "32" "64" )

for rt in "${realtype[@]}"; do
    for is in "${indexsize[@]}"; do

        # print test header for Jenkins section collapsing
        echo "TEST: ./suntest.sh $rt $is both ON DEV $bt"

        # run tests
        ./suntest.sh $rt $is both ON DEV $bt

        # check return flag
        if [ $? -ne 0 ]; then
            echo "FAILED: ./suntest.sh $rt $is both ON DEV $bt" | tee -a suntest.log
            exit 1
        else
            echo "PASSED"
        fi

    done
done

# success
exit 0
