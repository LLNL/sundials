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
# SUNDIALS pull-request testing driver script
#
# Usage: ./test_pr.sh <num build threads>
#
# Inputs (all optional):
#   <num build threads> = number of threads in parallel build
# ------------------------------------------------------------------------------

# Number of build threads for parallel build (default 1)
bt=1
if [ "$#" -ge 1 ]; then
    bt=$1
fi

# ==============================================================================
# Test SUNDIALS tarball with TPLs and DEVTESTS ON
# ==============================================================================

# run tests
./suntest_tarscript.sh sundials all both both ON DEV $bt

# check return flag
if [ $? -ne 0 ]; then
    echo "FAILED: ./suntest_tarscript.sh sundials all both both ON DEV $bt" | tee -a suntest.log
    exit 1
fi

# success
exit 0
