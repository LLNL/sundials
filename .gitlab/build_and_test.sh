#!/usr/bin/bash

# make sure lmod is loaded
if test -e /usr/share/lmod/lmod/init/bash
then
  source /usr/share/lmod/lmod/init/bash
fi

set -o errexit

option=${1:-""}
hostname="$(hostname)"
project_dir="$(pwd)"

build_root=${BUILD_ROOT:-""}
hostconfig=${HOST_CONFIG:-""}
spec=${SPEC:-""}
job_unique_id=${CI_JOB_ID:-""}

sys_type=${SYS_TYPE:-""}
py_env_path=${PYTHON_ENVIRONMENT_PATH:-""}

shared_spack=${SHARED_SPACK:-"UPSTREAM"}

# Dependencies
date

echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
echo "~~~~~ Inputs"
echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
echo "options       = ${option}"
echo "hostname      = ${hostname}"
echo "project_dir   = ${project_dir}"
echo "build_root    = ${build_root}"
echo "hostconfig    = ${hostconfig}"
echo "spec          = ${spec}"
echo "job_unique_id = ${job_unique_id}"
echo "sys_type      = ${sys_type}"
echo "py_env_path   = ${py_env_path}"
echo "shared_spack  = ${shared_spack}"

# remove tailing number from hostname
hostname=${hostname%%[0-9]*}

# number of parallel build jobs
BUILD_JOBS=${BUILD_JOBS:-"1"}

# load newer python to try the clingo concretizer
echo "module load python/3.8.2"
module load python/3.8.2

if [[ "${option}" != "--build-only" && "${option}" != "--test-only" ]]
then
    echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    echo "~~~~~ Building Dependencies"
    echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"

    if [[ -z ${spec} ]]
    then
        echo "SPEC is undefined, aborting..."
        exit 1
    fi

    prefix_opt=""

    if [[ -d /dev/shm ]]
    then
        prefix="/dev/shm/${hostname}"
        if [[ -z ${job_unique_id} ]]; then
          job_unique_id=manual_job_$(date +%s)
          while [[ -d ${prefix}/${job_unique_id} ]] ; do
              sleep 1
              job_unique_id=manual_job_$(date +%s)
          done
        fi

        prefix="${prefix}/${job_unique_id}"
        mkdir -p "${prefix}"
        prefix_opt="--prefix=${prefix}"
    fi

    if [[ -d /usr/workspace/sundials ]]
    then
        upstream="/usr/workspace/sundials/spack_installs/${hostname}"
        mkdir -p "${upstream}"
        upstream_opt="--upstream=${upstream}"
    fi

    if [[ "${shared_spack}" == "UPSTREAM" ]]
    then
        python3 scripts/uberenv/uberenv.py --spec="${spec}" "${prefix_opt}" "${upstream_opt}"
    elif [[ "${shared_spack}" == "ON" ]]
    then
        python3 scripts/uberenv/uberenv.py --spec="${spec}" --prefix="${upstream}"
    else
        python3 scripts/uberenv/uberenv.py --spec="${spec}" "${prefix_opt}"
    fi

    # Ensure correct CUDA module is loaded, only works for module naming
    # convention on Lassen. Only on needed for CUDA 11 (unclear why).
    if [[ -n "${CUDA_SPEC}" ]]; then
        cuda_version="${CUDA_SPEC##*@}"
        echo "module load cuda/${cuda_version}"
        module load cuda/"${cuda_version}"
    fi
fi
date

# Host config file
if [[ -z ${hostconfig} ]]
then
    # If no host config file was provided, we assume it was generated.
    # This means we are looking of a unique one in project dir.
    hostconfigs=( $( ls "${project_dir}/"hc-*.cmake ) )
    if [[ ${#hostconfigs[@]} == 1 ]]
    then
        hostconfig_path=${hostconfigs[0]}
        echo "Found host config file: ${hostconfig_path}"
    elif [[ ${#hostconfigs[@]} == 0 ]]
    then
        echo "No result for: ${project_dir}/hc-*.cmake"
        echo "Spack generated host-config not found."
        exit 1
    else
        echo "More than one result for: ${project_dir}/hc-*.cmake"
        echo "${hostconfigs[@]}"
        echo "Please specify one with HOST_CONFIG variable"
        exit 1
    fi
else
    # Using provided host-config file.
    hostconfig_path="${project_dir}/host-configs/${hostconfig}"
fi

# Build Directory
if [[ -z ${build_root} ]]
then
    build_root=$(pwd)
fi

build_dir="${build_root}/build_${job_unique_id}"

# Build
if [[ "${option}" != "--deps-only" && "${option}" != "--test-only" ]]
then
    date
    echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    echo "~ Host-config: ${hostconfig_path}"
    echo "~ Build Dir:   ${build_dir}"
    echo "~ Project Dir: ${project_dir}"
    echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"

    echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    echo "~~~~~ Building SUNDIALS"
    echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"

    # If building, then delete everything first
    rm -rf "${build_dir}" 2>/dev/null
    mkdir -p "${build_dir}" && cd "${build_dir}"

    date
    cmake --version

    # configure
    cmake -C "${hostconfig_path}" "${project_dir}"

    # build
    VERBOSE_BUILD=${VERBOSE_BUILD:-"OFF"}
    if [[ "${VERBOSE_BUILD}" == "ON" ]]; then
        verbose_build='--verbose'
    fi

    cmake --build . -j ${BUILD_JOBS} ${verbose_build}
    date
fi

# Test
if [[ "${option}" != "--build-only" ]]
then
    echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    echo "~~~~~ Testing SUNDIALS"
    echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"

    if [[ ! -d ${build_dir} ]]
    then
        echo "ERROR: Build directory not found : ${build_dir}" && exit 1
    fi

    cd "${build_dir}"

    VERBOSE_TEST=${VERBOSE_TEST:-"ON"}
    if [[ "${VERBOSE_TEST}" == "ON" ]]; then
        verbose_test='--verbose'
    fi

    date
    ctest --output-on-failure ${verbose_test} -T test 2>&1 | tee tests_output.txt
    date

    no_test_str="No tests were found!!!"
    if [[ "$(tail -n 1 tests_output.txt)" == "${no_test_str}" ]]
    then
        echo "ERROR: No tests were found" && exit 1
    fi

    if grep -q "Errors while running CTest" ./tests_output.txt
    then
        echo "ERROR: failure(s) while running CTest" && exit 1
    fi

fi

# Clean
if [[ "${option}" != "--no-clean" ]]
then
    echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    echo "~~~~~ CLEAN UP"
    echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    make clean
fi
