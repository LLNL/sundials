..
   -----------------------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2002-2024, Lawrence Livermore National Security
   and Southern Methodist University.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   -----------------------------------------------------------------------------

.. _CI:

GitHub CI Testing
=================

There are two categories of CI testing that we run on GitHub via `GitHub actions <https://github.com/LLNL/sundials/actions>`_:

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
the latest versions of dependencies regularly and detect interface breakages.


GitLab CI Testing
=================

This section provides an overview of the GitLab continuous integration (CI)
features and describes the SUNDIALS GitLab CI setup. For more information see
the official `GitLab CI documentation <https://docs.gitlab.com/ee/ci/>`_.

For information specific to the LLNL GitLab CI see:

* `LLNL GitLab CI <https://lc.llnl.gov/confluence/display/GITLAB/GitLab+CI>`_

* `LLNL GitLab Runner Tags <https://lc.llnl.gov/gitlab/public-info/gitlab-ci/-/wikis/Gitlab-CI-Basic-Information>`_


SUNDIALS utilizes the GitLab CI pipeline code repository shared by LLNL RADIUSS
projects. The `docs <https://radiuss-shared-ci.readthedocs.io/en/latest/>`__ for
the shared project should be reviewed before reading the docs below.


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

Some pipelines are run automatically on new commits (after they are mirrored
from GitHub to LC GitLab). Other pipelines, such as the benchmarking pipeline,
are run on a schedule that is configured through the
`GitLab UI <https://lc.llnl.gov/gitlab/sundials/sundials/-/pipeline_schedules>`__.

Structure
^^^^^^^^^

As previously stated, most of our code for the LC GitLab CI pipelines is sourced from
the templates provided in the
`radiuss-shared-ci <https://radiuss-shared-ci.readthedocs.io/en/latest/>`__ repo.
Here we briefly outline the relevant files:

The ``.gitlab-ci.yml`` file in the root of the repository is the starting point for
defining the SUNDIALS GitLab CI pipelines. The only thing that is typically changed
in this file is the ``SPACK_REF`` variable in the ``variables`` section (this
is done when we update the version of Spack we use for installing dependencies).

The ``.gitlab/subscribed-pipelines.yml`` defines which machines we will test on.
This file may be modified if you need to add a new machine to test on.

The ``.gitlab/custom-jobs-and-variables.yml`` defines variables available in all
pipelines and jobs. This file may be modified if you need to add a new variable
that needs to be accessible to all pipelines and jobs.

The ``.gitlab/jobs/<machine>.yaml`` files define the jobs for a specific machine.
A "hidden" job of the form `.sundials_job_on_<machine>` is defined first in these
files and typically defines variables specific to that machine, such as what command
to use for executing MPI programs. The rest of the jobs in the file extend the
`.sundials_job_on_<machine>` and define the Spack spec that we will build and test.
Take for example, this Tioga job:

.. code-block:: YAML

  tioga_rocmcc_620_tpls:
    parallel:
      matrix:
        - COMPILER_SPEC: rocmcc@6.2.0
          AMDGPU_TARGET: [gfx90a]
    variables:
      # have to use ginkgo@master because our spack version does not have ginkgo@1.8.0: yet (which seems to be needed)
      # similarly, we need a newer magma than available to compile with 'rocm@6:' so we turn it off
      SPEC: "%${COMPILER_SPEC} cstd=99 cxxstd=14 precision=double amdgpu_target=${AMDGPU_TARGET} +rocm+mpi~magma+ginkgo+kokkos ^ginkgo@master+rocm amdgpu_target=${AMDGPU_TARGET} ^kokkos+rocm amdgpu_target=${AMDGPU_TARGET}"
    before_script:
      - module load rocmcc/6.2.0-magic
    extends: [.sundials_job_on_tioga]

The ``parallel:`` and ``matrix:`` keywords could be used to enable creating multiple jobs
with different variable values for each instance of the job, e.g., one job using
``rocmcc@6.2.0`` and another using ``rocmcc@6.2.1``. However, right now they only create
a single job (hence why ``COMPILER_SPEC`` and ``AMDGPU_TARGET`` only have one value). These
variables values are then used to create an environment variable ``SPEC`` which is the Spack spec
used by ``build_and_test.sh`` (discussed below) to configure and build SUNDIALS and the
necessary dependencies.

Disabling a Job
^^^^^^^^^^^^^^^

A job can be disabled by adding the variable ``.ON_<machine>: "OFF"`` to the ``variables:``
section of the job, e.g.,


.. code-block:: YAML

  tioga_rocmcc_620_tpls:
    parallel:
      matrix:
        - COMPILER_SPEC: rocmcc@6.2.0
          AMDGPU_TARGET: [gfx90a]
    variables:
      ON_TIOGA: "OFF" # disable this job
      SPEC: "%${COMPILER_SPEC} cstd=99 cxxstd=14 precision=double amdgpu_target=${AMDGPU_TARGET} +rocm+mpi~magma+ginkgo+kokkos ^ginkgo@master+rocm amdgpu_target=${AMDGPU_TARGET} ^kokkos+rocm amdgpu_target=${AMDGPU_TARGET}"
    before_script:
      - module load rocmcc/6.2.0-magic
    extends: [.sundials_job_on_tioga]

These variables can also be set when manually or scheduling a pipeline in the GitLab UI.

Updating Spack
^^^^^^^^^^^^^^

To update the spack commit used for the CI simply replace the commit hash in the
``SPACK_REF`` variable inside the ``.gitlab-ci.yml`` file with the new commit hash.
This will create a new spack build cache in ``/usr/workspace/sundials/ci/spack_stuff``
and rebuild all of the specs.


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
  repository that provides spack environment files for various LLNL platfornms
  i.e., ``spack.yaml`` for Dane, Tioga, etc.

These submodules work in conjunction with ``scripts/sundials/package.py``
to configure and build any third-party libraries needed by the SUNDIALS
configuration and generates an initial CMake cache file for building SUNDIALS.
Other packages can be added to ``spack/packages`` if the default Spack package
needs to be overridden.

Running Locally
^^^^^^^^^^^^^^^

It is possible to run these scripts locally on an LC machine. First set a ``SPACK_REF``
environment variable to a spack commit that you want to use, and then set a ``SPEC``
environment variable with a SUNDIALS spack spec that you want to test.
