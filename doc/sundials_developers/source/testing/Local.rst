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

Testing Locally
===============

Using SUNDIALS CI Docker Containers
-----------------------------------

It is possible to use the SUNDIALS CI containers available at
`https://github.com/orgs/LLNL/packages?ecosystem=container <https://github.com/orgs/LLNL/packages?ecosystem=container>`_
for testing locally. This will allow you to get as close as possible to running tests
in the CI environment on your own machine.

If you have Docker or `Podman <https://podman.io/>`_ installed on your machine already,
than the easiest way to use the containers is via the CMake targets ``setup_local_ci``
and ``test_local_ci``, for example:

.. code-block:: shell

   $ cd builddir
   $ make setup_local_ci 
   ...
   $ make test_local_ci 
   ...

The ``setup_local_ci`` target will pull the container for the configured ``SUNDIALS_PRECISION``
and ``SUNDIALS_INDEX_SIZE``. Then the ``test_local_ci`` target will run the full test suite
in the container.

Alternatively, if you want to work with Docker directly, you can pull the image(s) and then
run the test suite manually. The ``run`` command will pull the image and start the container:

.. code-block:: shell

   docker run -t -d --name sundialsci-int32-double-latest --cap-add SYS_PTRACE -v "/path/to/your/sundials/development/repo":/sundials ghcr.io/llnl/sundials-ci-int32-double:latest

The ``exec`` command can then be used to execute the test script:

.. code-block:: shell

   docker exec -w /sundials/test sundialsci-int32-double-latest ./test_driver.sh --testtype CUSTOM --env env/docker.sh --tpls --realtype double --indexsize 32

Alternatively, you can drop into a bash shell inside the container to run specific examples:

.. code-block:: shell

   docker exec -it sundialsci-int32-double-latest bash

On Macs, it is recommended to use `Podman <https://podman.io/>`_ (and then the
same steps above apply using ``podman`` instead of ``docker``). Podman is
useful on Linux too, as it can run rootless easily.

.. note::

   Its possible that the difference in your local machine architecture, and the
   one used to build the docker container(s), results in different answers and
   failing tests. You can provide the path to your own directory with answer
   files by setting the environment variable ``SUNDIALS_TEST_ANSWER_DIR`` with
   the path, and adding the argument ``-e SUNDIALS_TEST_ANSWER_DIR`` to the
   ``docker run`` command above.
