#!/usr/bin/env bash

# variables

# Build and test first

source "$BUILD_ROOT/.gitlab/build_and_test.sh" --no-clean

# Embed this into `make benchmark`
nresg=$((BENCHMARK_NNODES * 4)) # Lassen has 4 GPUs per node
nresc=$((BENCHMARK_NNODES * 40)) # Lassen has 40 cores per node

echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
echo "~~~~~ Benchmarking SUNDIALS"
echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"

# We use synchronous kernel launches for benchmarking
# because we are interested in individual region timings.
export CUDA_LAUNCH_BLOCKING=1

if [[ -d ${build_dir} ]]
then 

    # configure for multinode
    $cmake_exe \
        -C "${hostconfig_path}" \
        -DSUNDIALS_BENCHMARK_NUM_CPUS=${nresc} \
        -DSUNDIALS_BENCHMARK_NUM_GPUS=${nresg} \
        "${project_dir}"
    
    cd ${build_dir}
    make benchmark

fi

echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
echo "~~~~~ CLEAN UP"
echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
make clean

