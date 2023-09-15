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

spack_prefix=${SHARED_SPACK_PREFIX:-"v0.19.1"}
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
# machine specific loads
if [[ "${hostname}" == "lassen" ]]; then
    echo "module load python/3.8.2"
    module load python/3.8.2
elif [[ "${hostname}" == "corona" ]]; then
    echo "module load python/3.9.12"
    module load python/3.9.12
else
    echo "module load python"
    module load python
fi

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
        mkdir -p ${prefix}
        prefix_opt="--prefix=${prefix}"

        # We force Spack to put all generated files (cache and configuration of
        # all sorts) in a unique location so that there can be no collision
        # with existing or concurrent Spack.
        spack_user_cache="${prefix}/spack-user-cache"
        export SPACK_DISABLE_LOCAL_CONFIG=""
        export SPACK_USER_CACHE_PATH="${spack_user_cache}"
        mkdir -p ${spack_user_cache}
    fi

    if [[ -d /usr/workspace/sundials ]]
    then
        upstream="/usr/workspace/sundials/spack_installs/${spack_prefix}/${hostname}"
        mkdir -p "${upstream}"
        upstream_opt="--upstream=${upstream}"
    fi

    if [[ "${shared_spack}" == "UPSTREAM" ]]
    then
        python3 .gitlab/uberenv/uberenv.py --spec="${spec}" "${prefix_opt}" "${upstream_opt}"
    elif [[ "${shared_spack}" == "ON" ]]
    then
        python3 .gitlab/uberenv/uberenv.py --spec="${spec}" --prefix="${upstream}"
    else
        python3 .gitlab/uberenv/uberenv.py --spec="${spec}" "${prefix_opt}"
    fi

    # Ensure correct CUDA module is loaded, only works for module naming
    # convention on Lassen. Only needed for CUDA 11 (unclear why).
    if [[ -n "${CUDA_SPEC}" ]]; then
        cuda_version="${CUDA_SPEC##*@}"
        echo "module load cuda/${cuda_version}"
        module load cuda/"${cuda_version}"
    fi

    module load cmake/3.23
fi
date

# Host config file
if [[ -z ${hostconfig} ]]
then
    # If no host config file was provided, we assume it was generated.
    # This means we are looking of a unique one in project dir.
    hostconfigs=( $( ls "${project_dir}/"*.cmake ) )
    if [[ ${#hostconfigs[@]} == 1 ]]
    then
        hostconfig_path=${hostconfigs[0]}
        echo "Found host config file: ${hostconfig_path}"
    elif [[ ${#hostconfigs[@]} == 0 ]]
    then
        echo "No result for: ${project_dir}/*.cmake"
        echo "Spack generated host-config not found."
        exit 1
    else
        echo "More than one result for: ${project_dir}/*.cmake"
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
    build_root="$(pwd)"
else
    build_root="${build_root}"
fi

build_dir="${build_root}/build_${job_unique_id}_${hostconfig//.cmake/}"
install_dir="${build_root}/install_${job_unique_id}_${hostconfig//.cmake/}"

cmake_exe=`grep 'CMake executable' ${hostconfig_path} | cut -d ':' -f 2 | xargs`

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

    $cmake_exe --version

    # configure
    if [[ "${CI_COMMIT_BRANCH}" == "main" ]]
    then
        # redirect caliper files to release directory
        sundials_version=$(cd ${project_dir}; git describe --abbrev=0)
        $cmake_exe \
            -C "${hostconfig_path}" \
            -DCMAKE_INSTALL_PREFIX=${install_dir} \
            -DSUNDIALS_CALIPER_OUTPUT_DIR="${CALIPER_DIR}/Release/${hostname}/${sundials_version}" \
            "${project_dir}"

    else
        $cmake_exe \
            -C "${hostconfig_path}" \
            -DCMAKE_INSTALL_PREFIX=${install_dir} \
            "${project_dir}"
    fi

    # build
    VERBOSE_BUILD=${VERBOSE_BUILD:-"OFF"}
    if [[ "${VERBOSE_BUILD}" == "ON" ]]; then
        verbose_build='--verbose'
    fi
    $cmake_exe --build . -j ${BUILD_JOBS} ${verbose_build}

    # install
    make install

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
