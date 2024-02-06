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


.. _CPT:

Continuous Performance Testing (CPT)
====================================

In order to protect against performance regression of SUNDIALS at all scales,
we leverage the GitLab CI infrastructure setup for CI testing to perform
continuous performance testing.

The CPT suite consists of the :ref:`SUNDIALS benchmark programs<Benchmarks>`
programs, which can scale up to full supercomputers, and the normal SUNDIALS
examples program suite (i.e., the programs in the ``examples`` directory of
the SUNDIALS repo).

The CI suite can run the regular SUNDIALS CI tests, or it can run the SUNDIALS
CPT suite by setting the ``BENCHMARK`` variable to ``ON`` when running a
pipeline from the GitLab CI UI.
The benchmark problems are run with Caliper and a report for Spot and a
human-readable runtime-report are generated.
The runtime-report is printed to the stdout and can be viewed in the GitLab
CI job output. The Spot output files are made available as job artifacts.
We maintain a pipeline that runs the CPT suite weekly on a schedule,
see `<https://lc.llnl.gov/gitlab/sundials/sundials/-/pipeline_schedules>`_.

Performance over time is tracked with `Caliper <https://lc.llnl.gov/confluence/display/CALI/Spot+DB>`_
and `SPOT web framework <https://lc.llnl.gov/confluence/display/SpotDoc/Spot+Documentation>`_.

Locally Building and Running the CPT
------------------------------------

The SUNDIALS example suite can be run with Caliper profiling enabled and
Adiak enabled for the CPT suite by setting the CMake options

.. code-block:: bash

  $ cmake \
  > -DSUNDIALS_BUILD_WITH_PROFILING=ON \
  > -DENABLE_CALIPER=ON \
  > -DCaliper_DIR=/path/to/caliper \
  > -DENABLE_ADIAK=ON \
  > -Dadiak_DIR=/path/to/adiak/lib/cmake/adiak \
  > -DSUNDIALS_TEST_DEVTESTS=ON \
  > -DSUNDIALS_TEST_PROFILE=ON \

This command will result in ``--profile`` option being passed to the SUNDIALS
test runner Python script, ``test/testRunner``, which will in turn set the
``CALI_CONFIG`` environment variable before running every test so that when
you run ``make test`` the examples will produce `.cali` output files
documenting the performance.

Note: Caliper prints to the `.out` files by default. Ensure all Caliper configs
requested have the `output` option defined to ensure output data is saved in a
separate location from the test output. Otherwise, the `.out` files for each
test will contain the output and the tests will fail.

Turning on the ``BUILD_BENCHMARKS`` option will build benchmarks. Running
``make benchmark`` will execute all the available benchmarks and produce
`.cali` output files for each one. To change what parameters benchmarks are run
with, edit the respective `CMakeLists.txt`. The ``BENCHMARK_VARS`` variable
determines how many tests to run with different parameters. Arguments passed
into the ``sundials_add_benchmark`` macro change how the benchmark is run.

To specify where `.cali` output files are placed, define the CMake option
SUNDIALS_CALIPER_OUTPUT_DIR with the directory path. By default `.cali` output
files are placed in the build directory (:ref:`Installation`) under
``Benchmarking/output`` and ``Testing/output``.

Refer to section :ref:`Benchmarks` for details on instructions on building
and running the ``benchmarks/`` programs locally. Refer to section <<examples>>
for instructions on build and running examples locally.
