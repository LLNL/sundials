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

Continuous Performance Testing (CPT)
====================================

In order to protect against performance regression of SUNDIALS at all scales, we leverage the GitLab CI infrastructure setup for CI testing to perform continuous performance testing. 

The CPT suite consists of the :ref:`SUNDIALS benchmark programs<Benchmarks>` programs, which can scale up to full supercomputers, and the normal SUNDIALS examples program suite (i.e., the programs in the ``examples`` directory of the SUNDIALS repo).

The CI suite can run the regular SUNDIALS CI tests, or it can run the SUNDIALS CPT suite by setting the ``BENCHMARK`` variable to ``ON`` when running a pipeline from the GitLab CI UI.
The benchmark problems are run with Caliper and a report for Spot and a human-readable runtime-report are generated. 
The runtime-report is printed to the stdout and can be viewed in the GitLab CI job output. The Spot output files are made available as job artifacts.
We maintain a pipeline that runs the CPT suite weekly on a schedule, see `<https://lc.llnl.gov/gitlab/sundials/sundials/-/pipeline_schedules>`_.

FUTURE: WE WILL TRACK PERFORMANCE OVER TIME WITH `CALIPER <https://lc.llnl.gov/confluence/display/CALI/Spot+DB>`_ AND THE `SPOT WEB FRAMEWORK <https://lc.llnl.gov/confluence/display/SpotDoc/Spot+Documentation>`_.


Locally Building and Running the CPT
------------------------------------

The SUNDIALS example suite can be run with Caliper profiling enabled for the CPT suite by setting the CMake options

.. code-block:: cmake

  cmake -DSUNDIALS_BUILD_WITH_PROFILING=ON -DENABLE_CALIPER=ON -DCaliper_DIR=/path/to/caliper -DSUNDIALS_TEST_DEVTESTS=ON -DSUNDIALS_TEST_PROFILE=ON 

This command will result in ``--profile`` option being passed to the SUNDIALS test runner Python script, ``test/testRunner``, which will in turn set the ``CALI_CONFIG`` environment variable before running every test so that when you run ``make test`` the examples will produce `.cali` output files documenting the performance.

Refer to documentation section :ref:`Benchmarks` for instructions on building and running the ``benchmarks/`` programs locally.
