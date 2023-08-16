..
   -----------------------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2002-2023, Lawrence Livermore National Security
   and Southern Methodist University.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   -----------------------------------------------------------------------------

GitHub CI Testing 
=================

There are two types of CI testing that we run on GitHub via `GitHub actions <https://github.com/LLNL/sundials/actions>`_:

1. Comprehensive (excluding GPUs)
2. Minimal/Short

The comprehensive testing is run only on one platform (Ubuntu) and utilizes Docker + Spack + the
latest E4S release to build SUNDIALS in many different configurations with and without third-party
libraries enabled.

The minimal/short testing runs on more platforms: Windows (MinGW and MSVC), and MacOS but it runs in
only one configuration (using 64-bit indices and double precision) and without third-party
libraries.


Building the Docker containers for CI
-------------------------------------

Original Setup
^^^^^^^^^^^^^^

These are the steps that were originally performed by Cody Balos
to build the Docker container(s) used for the comprehensive CI testing.

1. Create a `spack.yaml` for each configuration of SUNDIALS. E.g., int64 and double:

.. code-block:: yaml

  # This is a Spack Environment file.
  #
  # It describes a set of packages to be installed, along with
  # configuration settings.
  spack:
    packages:
      all:
        providers:
          blas: [openblas]
          mpi: [openmpi]
    # add package specs to the `specs` list
    specs:
    - hypre+int64~internal-superlu
    - petsc+double+int64
    - openmpi
    - openblas+ilp64
    - suite-sparse
    - superlu-dist+int64 ^parmetis+int64
    - trilinos+tpetra gotype=long_long
    config: {}
    modules:
      enable: []
    repos: []
    upstreams: {}
    container:
      images:
        os: ubuntu:20.04

2. Run `spack containerize > Dockerfile` in the directory of the `spack.yaml`

3. The Dockerfile produced in step 2 was then manually modified to leverage
   Spack's buildcache feature and Docker's 'cache' bind-mount. The gist is that
   if the Spack build fails while building the Docker image, the buildcache
   makes it possible to reuse the binaries for the packages that were already installed
   before the build failed. Without the buildcache, the spack build failing would
   result in all packages needing to be built again when re-attempting to build the Docker image.

4. Run `DOCKER_BUILDKIT docker build -t sundials-ci/<index-size>-<precision>:<tag>`

5. Push

Automated building of new containers
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We currently have six different containers, one for each combination of {int32, int64} and {single,
double, extended} precision. These containers are pinned to an E4S release. When E4S does a release,
we can rebuild these containers to use the packages from it. We add E4S as a mirror in the Spack
environment so that its buildcache can be leveraged.

We also maintain two containers for the {int32, double} pair that are built automatically (in a
GitHub action) every week against the latest Spack develop commit. This allows us to test against
the latest versions of dependencies reguarly and detect interface breakages.


GitLab CI Testing
=================

This section provides an overview of the GitLab continuous integration (CI)
features and describes the SUNDIALS GitLab CI setup. For more information see
the official `GitLab CI documentation <https://docs.gitlab.com/ee/ci/>`_.

For information specific to the LLNL GitLab CI see:

* `LLNL GitLab CI <https://lc.llnl.gov/confluence/display/GITLAB/GitLab+CI>`_

* `LLNL GitLab Runner Tags <https://lc.llnl.gov/gitlab/public-info/gitlab-ci/-/wikis/Gitlab-CI-Basic-Information>`_


CI Pipelines and Jobs
---------------------

Pipelines are the main component of the GitLab CI system and are comprised of
two components, stages and jobs. Stages define the logical steps in a pipeline
i.e., *when* to do something. Jobs are the basic configuration component and
define *what* to do. Each job is associated with a stage and multiple jobs in
the same stage are executed in parallel (if there are enough resources). If all
the jobs in a stage succeed, the pipeline moves on to the next stage. If any job
in a stage fails, the next stage is usually (see below) not executed and the
pipeline ends early.


CI Pipeline
^^^^^^^^^^^

The YAML file ``.gitlab-ci.yml`` at the top-level of the repository defines the
GitLab CI pipeline. The SUNDIALS CI configuration file is organized as follows:

* The ``variables:`` keyword defines environment variables shared by all stages
  and jobs. Variables can be used to alter the CI pipeline behavior when
  launching a job from the GitLab "Run pipeline" UI or from the GitLab pipeline
  scheduling UI. See the ``.gitlab-ci.yml`` file for details about all of the
  available variables and what they do.

  .. code-block:: YAML

     variables:
       GIT_SUBMODULE_STRATEGY: recursive
       ALLOC_NAME: ${CI_PROJECT_NAME}_ci_${CI_PIPELINE_ID}
       BUILD_ROOT: ${CI_PROJECT_DIR}
       # ...

* The ``stages:`` keyword defines independent CI stages targeting a specific
  test machine following the prefix naming convention:

  * ``q_`` jobs run on Quartz
  * ``l_`` jobs run on Lassen

  .. code-block:: YAML

     stages:
       - q_build_and_test
       - l_build_and_test
       # ...

* Several hidden job templates (job names start with ``.``) are defined for
  specific architectures and batch queue systems. These jobs provide the batch
  system command to run the ``build_and_test.sh`` script that configures,
  builds, and tests SUNDIALS.

  .. code-block:: YAML

     .build_toss_3_x86_64_ib_script:
       script:
       - echo ${ALLOC_NAME}
       - srun -p pdebug -N 1 -n ${NCPUS} --interactive -t ${DEFAULT_TIME}
         --job-name=${ALLOC_NAME} .gitlab/build_and_test.sh

     # ...

* The ``include:`` keyword loads YAML files defining the jobs for specific
  machines.

  .. code-block:: YAML

     include:
       - local: .gitlab/quartz-jobs.yml
       - local: .gitlab/lassen-jobs.yml
       # ...


CI Jobs
^^^^^^^

As noted above, each stage in the CI pipeline corresponds to testing on a
specific machine. For example, jobs run on Lassen are associated with the
``l_build_and_test`` stage. The actual jobs to run are defined in the YAML
file ``.gitlab/lassen-jobs.yml``.

The Lassen build and test jobs inherit from three job templates:

* ``.build_blueos_3_ppc64le_ib_script`` executes the LSF command to run the
  testing script.

* ``.on_lassen`` defines the tags (``tags:`` keyword) to select a shell runner
  on Lassen and the rules (``rules:`` keyword) for when a job should run.

* ``.lassen_build_and_test`` inherits from the prior two job templates using the
  ``extends:`` keyword and acts as the base jobs that all other Lassen jobs
  inherit from. The base template includes:

  * The ``stage:`` keyword defines which stage the jobs run in.

  * The ``needs:`` keyword lists the job dependencies. Normally, GitLab stages
    are blocking however, by providing the dependencies we can break the
    ordering of stages, in favor of using a DAG. This allows jobs to be run
    out-of-order rather than waiting on the jobs in other stages to complete.

  * The ``artifacts:`` keyword defines ``files:`` and directories (``paths:``)
    created by the job that should be retained and ``when:`` they should be
    attached to the job.

The Lassen tests are defined by jobs that extend the ``.lassen_build_and_test``
template and use the naming convention ``lassen_<compiler>_<test identifiers>``.
For example, tests using GCC, CUDA, and third-party libraries enabled are
defined by the job:

.. code-block:: YAML

   lassen_gcc_cuda_tpls:
     parallel:
       matrix:
         - COMPILER_SPEC: gcc@7.3.1
           CUDA_SPEC: [cuda@10.1.243, cuda@11.2.0]
     variables:
       SPEC: "%${COMPILER_SPEC} precision=double ~int64 +openmp +cuda +raja cuda_arch=70 \
              ^raja+cuda~examples~exercises cuda_arch=70 ^${CUDA_SPEC}"
     extends: .lassen_build_and_test

The ``parallel:`` and ``matrix:`` keywords enable creating multiple jobs with
different variable values for each instance of the job i.e., one job using
``cuda@10.1.243`` and another using ``cuda@11.2.0``. These variables values are
then used to create an environment variable ``SPEC`` with a Spack spec used by
``build_and_test.sh`` when configuring SUNDIALS.

Benchmark Jobs
^^^^^^^^^^^^^^

See :ref:`SUNDIALS Continuous Performance Testing (CPT)<CPT>` for more details.

GitLab CI Test Script
---------------------

The GitLab CI uses the script ``.gitlab/build_and_test.sh``, and when
benchmarking ``.gitlab/build_and_bench.sh``, to configure,
build, and test SUNDIALS. This script leverages two Git submodules:

* `uberenv <https://github.com/LLNL/uberenv>`_ -- automates using a package
  manager (e.g., Spack) to configure and build software. The top-level file
  ``.uberenv_config.json`` defines information need by uberenv including the
  the Spack commit to utilize and the location of Spack config and package
  files.

* `radiuss-spack-configs <https://github.com/sundials-codes/radiuss-spack-configs.git>`_
  -- is the SUNDIALS fork of the `LLNL radiuss-spack-configs <https://github.com/LLNL/radiuss-spack-configs>`_
  repository that provides spack configuration files for various LLNL platfornms
  i.e., ``compilers.yaml`` and ``packages.yaml`` files for Quartz, Lassen, etc.

These submodules work in conjunction with ``spack_packages/sundials/package.py``
to configure and build any third-party libraries needed by the SUNDIALS
configuration and generates an initial CMake cache file for building SUNDIALS.
Other packages can be added to ``spack_packages/<package name>/package.py``
if the default Spack package needs to be overriden. We do this currently for
Caliper, as we need a newer version than in the Spack commit currently used.

Updating Spack
--------------

To update the spack commit used for the CI:

1. The first thing to do is update the spack commit in the
``.uberenv_config.json`` file.
2. Then, a pipeline should be manually launched from the GitLab UI with the
``SHARED_SPACK`` CI variable set to ``ON`` and the ``SPACK_PREFIX`` variable to
the version of spack being set in the uberenv_config.json.

This will create a new spack installation and rebuild all of the specs. 

