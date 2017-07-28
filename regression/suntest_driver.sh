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

# ------------------------------------------------------------------------------
# Run regression tests
echo "--------------------------------------------------" | tee -a suntest.log
echo "SUNDIALS Regression Tests" | tee -a suntest.log
date | tee -a suntest.log
echo "--------------------------------------------------" | tee -a suntest.log

# loop over build options
for ((i=0; i<${#realtype[@]}; i++)); do
    for ((j=0; j<${#indextype[@]}; j++)); do

        # run tests
        ./suntest.sh ${realtype[i]} ${indextype[j]} $nbt

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
echo "Regression Tests Completed" | tee -a suntest.log
date | tee -a suntest.log
if [ $nfail -ne 0 ]; then
    echo "FAILED: $nfail failures." | tee -a suntest.log
else
    echo "PASSED" | tee -a suntest.log
fi
echo "--------------------------------------------------" | tee -a suntest.log

# notify users
# NEED TO ADD

exit $nfail
