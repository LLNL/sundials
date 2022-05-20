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

Testing Locally
===============

Using sundials-ci Docker Containers
-----------------------------------

It is possible to use the sundials-ci Docker containers available at
`https://hub.docker.com/r/balos1/sundials-ci <https://hub.docker.com/r/balos1/sundials-ci>_`
for testing locally. This will allow you to get as close as possible to running tests
in the CI environment on your own machine.

Developers with Linux machines can use Docker to pull the image(s) and then
run the test suite.

.. code-block:: shell

   docker pull balos1/sundials-ci:int32-double-latest
   docker run -t -d --name sundialsci-int32-double-latest --cap-add SYS_PTRACE -v "/path/to/your/sundials/development/repo":/sundials balos1/sundials-ci:int32-single-latest
   docker exec -it sundialsci-int32-double-latest bash # drops you into a docker container
   cd /sundials/test
   ./test_driver.sh --testtype CUSTOM --env env/docker.sh --tpls --realtype double --indexsize 32


On Macs, it is recommended to use `Podman <https://podman.io/>`_ (and then the same steps above apply. using ``podman`` instead of ``docker``).