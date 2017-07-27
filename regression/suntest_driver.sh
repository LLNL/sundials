#!/bin/bash
# ------------------------------------------------------------------------------
# SUNDIALS Regression Tests Driver Script
# ------------------------------------------------------------------------------

# number of threads for parallel builds (optional)
nbt=4

# initialize failure counter
nfail=0

# real and index types to test
# NOTE: need to create answers for different realtypes
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

# output to stdout and log file
echo "--------------------------------------------------" | tee -a suntest.log
echo "SUNDIALS Regression Tests" | tee -a suntest.log
date | tee -a suntest.log
echo "--------------------------------------------------" | tee -a suntest.log

# ------------------------------------------------------------------------------
# Minimal setup: no external solvers enabled
for ((i=0; i<${#realtype[@]}; i++)); do
    for ((j=0; j<${#indextype[@]}; j++)); do

        # run tests
        ./suntest_minimal.sh ${realtype[i]} ${indextype[j]} $nbt

        # check number of errors
        nerr=$?
        if [ $nerr -ne 0 ]; then
            let nfail+=nerr
            echo "FAILED: $nerr failures, suntest_minimal.sh ${realtype[i]} ${indextype[j]}" \
                | tee -a suntest.log
        else
            echo "PASSED: $nerr failures, suntest_minimal.sh ${realtype[i]} ${indextype[j]}" \
                | tee -a suntest.log
        fi
    done
done
echo "--------------------------------------------------" | tee -a suntest.log

# ------------------------------------------------------------------------------
# Maximal setup: all external solvers enabled
for ((i=0; i<${#realtype[@]}; i++)); do
    for ((j=0; j<${#indextype[@]}; j++)); do

        # run tests
        ./suntest_maximal.sh ${realtype[i]} ${indextype[j]} $nbt

        # check number of errors
        nerr=$?
        if [ $nerr -ne 0 ]; then
            let nfail+=nerr
            echo "FAILED: $nerr failures, suntest_maximal.sh ${realtype[i]} ${indextype[j]}" \
                | tee -a suntest.log
        else
            echo "PASSED: $nerr failures, suntest_maximal.sh ${realtype[i]} ${indextype[j]}" \
                | tee -a suntest.log
        fi
    done
done
echo "--------------------------------------------------" | tee -a suntest.log

# ------------------------------------------------------------------------------
# Return pass/fail
echo "--------------------------------------------------" | tee -a suntest.log
if [ $nfail -ne 0 ]; then
    echo "FAILED: $nfail failures." | tee -a suntest.log
else
    echo "PASSED: $nfail failures." | tee -a suntest.log
fi
echo "--------------------------------------------------" | tee -a suntest.log

# notify users
# NEED TO ADD

exit $nfail
