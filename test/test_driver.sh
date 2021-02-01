#!/bin/bash
# ------------------------------------------------------------------------------
# Programmer(s): David J. Gardner @ LLNL
# ------------------------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2021, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# ------------------------------------------------------------------------------
# SUNDIALS testing driver script
#
# Usage: ./test_driver.sh <num build threads> <test type> <branch name>
#
# Inputs (all optional):
#   <num build threads> = number of threads in parallel build (default 1)
#   <test type>         = test type to run (defalut BRANCH)
#   <branch name>       = the name of the branch (default none)
#
# Test Types:
#   BRANCH  = Quick tests that cover a few configurations
#   PR      = Create tarballs and test more configurations
#   RELEASE = Create tarballs and test many more configurations
# ------------------------------------------------------------------------------

# print input args
echo "./test_driver $@" | tee suntest.log

# Number of threads for parallel build (default 1)
buildthreads=1
if [ "$#" -ge 1 ]; then
    buildthreads=$1
fi

# Set the test type (default BRANCH). When called from Jenkins this input is
# either the branch name or the pull-request ID (PR-#). Based on this input the
# test type is set, but it may be overridden by the next input check for release
# branches. Alternatively, this input can be set to the desired test type when
# running outside of the Jenkins CI system.
tmptype=$2
case "$tmptype" in
    BRANCH|branch|Branch)
        testtype="BRANCH"
        ;;
    PR|pr)
        testtype="PR"
        ;;
    RELEASE|release|Release)
        testtype="RELEASE"
        ;;
    *)
        if [[ "${tmptype:0:2}" == "PR" ]]; then
            # pull-request test
            testtype="PR"
        elif [[ "${tmptype:0:7}" == "RELEASE" ||
                "${tmptype:0:7}" == "Release" ||
                "${tmptype:0:7}" == "release" ]]; then
            # release test
            testtype="RELEASE"
        else
            testtype="BRANCH"
        fi
        ;;
esac

# Set the branch name (no default). When called from Jenkins this input is the
# git branch name. If this is a relase branch, the test type is set to RELEASE
# to ensure a PR for a release branch runs the release tests rather than the PR
# tests. When using this script outside of the Jenkins CI system this input
# does not need to be set.
if [ "$#" -ge 3 ]; then
    branchname=$3
    if [[ "${branchname:0:7}" == "RELEASE" ||
          "${branchname:0:7}" == "Release" ||
          "${branchname:0:7}" == "release" ]]; then
        # release test
        testtype="RELEASE"
    fi
fi

# initialize failure counter (0 = success)
passfail=0

# remove old test directories and logs
\rm -rf build*/ install*/ *.log

# ------------------------------------------------------------------------------
# Run regression tests
# ------------------------------------------------------------------------------

# Print test header
echo "--------------------------------------------------" | tee -a suntest.log
echo "SUNDIALS $testtype test: $branchname " | tee -a suntest.log
date | tee -a suntest.log
echo "--------------------------------------------------" | tee -a suntest.log
git log -1 --pretty='Commit: %H%nAuthor: %an <%ae>%nDate:   %ad%n%n%s' | tee -a suntest.log
echo "--------------------------------------------------" | tee -a suntest.log

# call a test driver based on the test type
case "$testtype" in
    BRANCH)
        ./test_branch.sh $buildthreads
        passfail=$?
        ;;
    PR)
        ./test_pr.sh $buildthreads
        passfail=$?
        ;;
    RELEASE)
        ./test_release.sh $buildthreads
        passfail=$?
        ;;
    *)
        echo "ERROR: Unknown test type"
        passfail=1
        ;;
esac

# Print test footer
if [ $passfail -eq 0 ]; then
    echo "All tests PASSED" | tee -a suntest.log
fi
echo "--------------------------------------------------" | tee -a suntest.log
echo "SUNDIALS $testtype test: $branchname " | tee -a suntest.log
date | tee -a suntest.log
echo "--------------------------------------------------" | tee -a suntest.log

# exit
if [ $passfail -eq 0 ]; then
    # pass
    exit 0
else
    # fail
    exit 1
fi
