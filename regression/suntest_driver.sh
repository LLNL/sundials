#!/bin/bash
# -------------------------------------------------------------------------------
# Programmer(s): David J. Gardner @ LLNL 
# -------------------------------------------------------------------------------
# LLNS Copyright Start
# Copyright (c) 2014, Lawrence Livermore National Security
# This work was performed under the auspices of the U.S. Department 
# of Energy by Lawrence Livermore National Laboratory in part under 
# Contract W-7405-Eng-48 and in part under Contract DE-AC52-07NA27344.
# Produced at the Lawrence Livermore National Laboratory.
# All rights reserved.
# For details, see the LICENSE file.
# LLNS Copyright End
# -------------------------------------------------------------------------------
# SUNDIALS regression testing driver script
# -------------------------------------------------------------------------------

# number of threads for parallel build (defaults to 1)
if [ "$#" -ge 1 ]; then
    buildthreads=$1
else
    buildthreads=1
fi

# set test name (branch name or pull-request) if provided, otherwise leave empty
if [ "$#" -ge 2 ]; then
    testname=$2
fi

# initialize failure counter (0 = success)
nfail=0

# real and index types to test
# NOTE: may need to create answer files for different realtypes
# NOTE: master branch will ignore indextype
realtype=('double')
indextype=('int32_t' 'int64_t')

# real and index types to test
# NOTE: remove above realtype and indextype arrays and uncomment the 
# following arrays this after this file is merged from master to develop
# realtype=('single' 'double' 'extended')
# indextype=('int64_t' 'int32_t')

# remove old test directories and logs
\rm -rf build*/ install*/ *.log

# ------------------------------------------------------------------------------
# Run regression tests
echo "--------------------------------------------------" | tee -a suntest.log
echo "SUNDIALS regression tests: $testname " | tee -a suntest.log
date | tee -a suntest.log
echo "--------------------------------------------------" | tee -a suntest.log
git log -1 | tee -a suntest.log
echo "--------------------------------------------------" | tee -a suntest.log

# loop over build options
for ((i=0; i<${#realtype[@]}; i++)); do
    for ((j=0; j<${#indextype[@]}; j++)); do

        # print test label for Jenkins section collapsing
        echo -e "TEST: ./suntest.sh ${realtype[i]} ${indextype[j]} $buildthreads \n"

        # run tests
        ./suntest.sh ${realtype[i]} ${indextype[j]} $buildthreads

        # check return flag
        if [ $? -ne 0 ]; then
            let nfail+=1
            echo "FAILED: ${realtype[i]} ${indextype[j]}" | tee -a suntest.log
        else
            echo "PASSED: ${realtype[i]} ${indextype[j]}" | tee -a suntest.log
        fi

        # print test label for Jenkins section collapsing
        echo -e "TEST: ./suntest_xsdk.sh ${realtype[i]} ${indextype[j]} $buildthreads \n"

        # run tests using xSDK CMake options
        ./suntest_xsdk.sh ${realtype[i]} ${indextype[j]} $buildthreads

        # check return flag
        if [ $? -ne 0 ]; then
            let nfail+=1
            echo "FAILED: xSDK ${realtype[i]} ${indextype[j]}" | tee -a suntest.log
        else
            echo "PASSED: xSDK ${realtype[i]} ${indextype[j]}" | tee -a suntest.log
        fi               

    done
done

# ------------------------------------------------------------------------------
# Report test results
echo "--------------------------------------------------" | tee -a suntest.log
echo "SUNDIALS regression tests: $testname " | tee -a suntest.log
if [ $nfail -ne 0 ]; then
    echo "FAILED: $nfail failures." | tee -a suntest.log
else
    echo "PASSED" | tee -a suntest.log
fi
date | tee -a suntest.log
echo "--------------------------------------------------" | tee -a suntest.log

# exit
if [ $nfail -ne 0 ]; then
    # failure
    exit 1
else
    # success
    exit 0
fi
