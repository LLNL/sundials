..
   -----------------------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2002-2022, Lawrence Livermore National Security
   and Southern Methodist University.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   -----------------------------------------------------------------------------

.. _GitHub:

GitHub
======

[GitHub actions](https://github.com/LLNL/sundials/actions) is also used for SUNDIALS CI. There are
two types of CI testing that we run on GitHub:

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

```yaml
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
  - superlu-dist+int64 ^metis+int64
  - trilinos+tpetra gotype=long_long
  config: {}
  modules:
    enable: []
  repos: []
  upstreams: {}
  container:
    images:
      os: ubuntu:20.04
```

2. Run `spack containerize > Dockerfile` in the directory of the `spack.yaml`

3. The Dockerfile produced in step 2 was then manually modified to leverage
   Spack's buildcache feature and Docker's 'cache' bind-mount. The gist is that
   if the Spack build fails while building the Docker image, the buildcache
   makes it possible to reuse the binaries for the packages that were already installed
   before the build failed. Without the buildcache, the spack build failing would
   result in all packages needing to be built again when re-attempting to build the Docker image.

4. Run `DOCKER_BUILDKIT docker build -t sundials-ci/<index-size>-<precision>-<tag>`

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

