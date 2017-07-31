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

# check number of inputs
if [ "$#" -lt 2 ]; then
    echo "ERROR: Illegal number of parameters, branch name and test URL required"
    exit 1
fi
BRANCHNAME=$1 
TESTURL=$2

# add newer python install to path
export PATH=/usr/apps/python/latest/bin:$PATH

# use newer version of Git (same as Jenkins uses) to path
source /usr/apps/git/2.9.4/setup.sh

# number of threads for parallel builds (optional, if empty will use all threads)
buildthreads=4

# initialize failure counter (0 = success)
nfail=0

# real and index types to test
# NOTE: may need to create answer files for different realtypes
# NOTE: master branch will ignore indextype
realtype=('double')
indextype=('signed_64bit')

# real and index types to test
# NOTE: remove above realtype and indextype arrays and uncomment the 
# following arrays this after this file is merged from master to develop
# realtype=('single' 'double' 'extended')
# indextype=('signed_32bit' 'signed_64bit' 'unsigned_32bit' 'unsigned_64bit')

# remove old test directories and logs
\rm -rf suntest*/ *.log

# ------------------------------------------------------------------------------
# Run regression tests
echo "--------------------------------------------------" | tee -a suntest.log
echo "SUNDIALS regression tests on $BRANCHNAME branch   " | tee -a suntest.log
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

    done
done

# ------------------------------------------------------------------------------
# Return overall pass/fail
echo "--------------------------------------------------" | tee -a suntest.log
echo "SUNDIALS regression tests on $BRANCHNAME branch   " | tee -a suntest.log
if [ $nfail -ne 0 ]; then
    echo "FAILED: $nfail failures." | tee -a suntest.log
else
    echo "PASSED" | tee -a suntest.log
fi
date | tee -a suntest.log
echo "--------------------------------------------------" | tee -a suntest.log

# set to fail to test email notification
nfail=1

# ------------------------------------------------------------------------------
# Email notification
./suntest_notify.py $nfail $BRANCHNAME $TESTURL

# return pass/fail
exit $nfail
