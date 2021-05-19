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
# SUNDIALS tarball regression testing script
#
# Usage: ./suntest_tarscript.sh <package> <lib type> <real type> <index size>
#                               <TPL status> <test type> <build threads>
#
# Required Inputs:
#   <package>    = Which tarball to make and test:
#                    arkode   : create ARKode tarball only
#                    cvode    : create CVODE tarball only
#                    cvodes   : create CVODES tarball only
#                    ida      : create IDA tarball only
#                    idas     : create IDAS tarball only
#                    kinsol   : create KINSOL tarball only
#                    sundials : create sundials tarball containing all packages
#                    all      : all of the above options
#   <real type>  = SUNDIALS real type to build/test with:
#                    single   : single precision only
#                    double   : double precision only
#                    extended : extended (quad) precision only
#                    all      : all of the above options
#   <index size> = SUNDIALS index size to build/test with:
#                    32       : 32-bit indices only
#                    64       : 64-bit indices only
#                    both     : both of the above options
#   <lib type>   = Which library type to test:
#                    static   : only build static libraries
#                    shared   : only build shared libraries
#                    each     : build static and shared separately
#                    both     : build static and shared simultaneously
#   <TPL status> = Enable/disable third party libraries:
#                    ON       : All possible TPLs enabled
#                    OFF      : No TPLs enabled
#   <test type>  = Test type to run:
#                    STD      : standard tests
#                    DEV      : development tests
#                    NONE     : no test, configure and compile only
#
# Optional Inputs:
#   <build threads> = number of threads to use in parallel build (default 1)
# ------------------------------------------------------------------------------

# exit the script if a command fails
set -e

# check number of inputs
if [ "$#" -lt 6 ]; then
    echo "ERROR: SIX (6) inputs required"
    echo "package      : [arkode|cvode|cvodes|ida|idas|kinsol|sundials|all]"
    echo "real type    : [single|double|extended|all]"
    echo "index size   : [32|64|both]"
    echo "library type : [static|shared|each|both]"
    echo "TPLs         : [ON|OFF]"
    echo "test type    : [STD|DEV|NONE]"
    exit 1
fi

package=$1      # sundials package to test
realtype=$2     # precision for realtypes
indexsize=$3    # integer size for indices
libtype=$4      # library type to build
tplstatus=$5    # enable/disable third party libraries
testtype=$6     # run standard tests, dev tests, or no tests (compile only)
buildthreads=1  # default number threads for parallel builds

# check if the number of build threads was set
if [ "$#" -gt 6 ]; then
    buildthreads=$7
fi

# ------------------------------------------------------------------------------
# Check inputs
# ------------------------------------------------------------------------------

# check package option
case $package in
    arkode|cvode|cvodes|ida|idas|kinsol|sundials|all) ;;
    *)
        echo "ERROR: Unknown package option: $package"
        exit 1
        ;;
esac

# set real types to test
case "$realtype" in
    SINGLE|Single|single)       realtype=( "single" );;
    DOUBLE|Double|double)       realtype=( "double" );;
    EXTENDED|Extended|extended) realtype=( "extended" );;
    ALL|All|all)                realtype=( "single" "double" "extended" );;
    *)
        echo "ERROR: Unknown real type option: $realtype"
        exit 1
        ;;
esac

# set index sizes to test
case "$indexsize" in
    32)   indexsize=( "32" );;
    64)   indexsize=( "64" );;
    both) indexsize=( "32" "64" );;
    *)
        echo "ERROR: Unknown index size option: $indexsize"
        exit 1
        ;;
esac

# set library types to test
case "$libtype" in
    STATIC|Static|static) libtype=( "static" );;
    SHARED|Shared|shared) libtype=( "shared" );;
    EACH|Each|each)       libtype=( "static" "shared") ;;
    BOTH|Both|both)       libtype=( "both" ) ;;
    *)
        echo "ERROR: Unknown library type option: $libtype"
        exit 1
        ;;
esac

# set TPL status
case "$tplstatus" in
    ON|On|on)
        TPLs=ON
        ;;
    OFF|Off|off)
        TPLs=OFF
        ;;
    *)
        echo "ERROR: Unknown third party library status: $tplstatus"
        exit 1
        ;;
esac

# which tests to run (if any)
case "$testtype" in
    STD|std|Std)
        # only run standard tests
        testtype=STD
        ;;
    DEV|dev|Dev)
        # only run development tests
        testtype=DEV
        ;;
    NONE|none|None)
        # only compile sundials, do not test or install
        testtype=NONE
        ;;
    *)
        echo "ERROR: Unknown test option: $testtype"
        exit 1
        ;;
esac

# ------------------------------------------------------------------------------
# Create tarballs
# ------------------------------------------------------------------------------

# location of testing directory
testdir=`pwd`

# create directory for tarballs
\rm -rf $testdir/tarballs
mkdir $testdir/tarballs

# run tarscript to create tarballs
cd ../scripts

echo "START TARSCRIPT"
./tarscript $package | tee -a tar.log

# check tarscript return code
rc=${PIPESTATUS[0]}
echo -e "\ntarscript returned $rc\n" | tee -a tar.log
if [ $rc -ne 0 ]; then exit 1; fi

# relocate log and tarballs
mv tar.log $testdir/tarballs/.
mv ../tarballs/* $testdir/tarballs/.

# ------------------------------------------------------------------------------
# Test tarballs
# ------------------------------------------------------------------------------

# move to tarball directory
cd $testdir/tarballs

# loop over tarballs and test each one
for tarball in *.tar.gz; do

    # get package name
    package=${tarball%.tar.gz}

    echo "START UNTAR"
    tar -xvzf $tarball 2>&1 | tee -a tar.log

    # check tar return code
    rc=${PIPESTATUS[0]}
    echo -e "\ntar -xzvf returned $rc\n" | tee -a tar.log
    if [ $rc -ne 0 ]; then exit 1; fi

    # move log to package directory
    mv tar.log $package/.

    # move to the extracted package's test directory
    cd $package/test

    # copy environment and testing scripts from original test directory
    if [ -f "$testdir/env.sh" ]; then
        cp $testdir/env.sh .
    fi
    cp -r $testdir/env .
    cp $testdir/suntest.sh .

    # loop over build options
    for rt in "${realtype[@]}"; do
        for is in "${indexsize[@]}"; do
            for lt in "${libtype[@]}"; do

                # print test label for Jenkins section collapsing
                echo "TEST: $rt $is $lt $TPLs $testtype $buildthreads"

                ./suntest.sh $rt $is $lt $TPLs $testtype $buildthreads

                # check return flag
                if [ $? -ne 0 ]; then
                    echo "FAILED: $rt $is $lt $TPLs $testtype"
                    cd $testdir
                    exit 1
                else
                    echo "PASSED"
                fi

            done
        done
    done

    # return to tarball directory
    cd $testdir/tarballs

done

# ------------------------------------------------------------------------------
# Return
# ------------------------------------------------------------------------------

# if we make it here all tests have passed
cd $testdir
exit 0
