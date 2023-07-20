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

SPOT Performance Visualization Tool
====================================

In order to protect against performance regression of SUNDIALS at all scales,
we leverage the GitLab CI infrastructure setup for CI testing to perform
continuous performance testing. 

To track and visualize performance over time, we collect data with
`Caliper <https://lc.llnl.gov/confluence/display/CALI/Spot+DB>`_ and
`Adiak <https://github.com/LLNL/Adiak>`_. Steps to setting up continuous
performance testing(CPT) can be found in
:ref:`Continuous Performance Testing<CPT>`.

Documentation to using the SPOT performance visualization tool can be found
`here <https://lc.llnl.gov/confluence/display/SpotDoc/Spot+Documentation>`_.

Setting Up Your Own SPOT Visualizations
---------------------------------------

To display data from non-GitLab jobs or create a local collection of runs,
input into the SPOT search bar the absolute path to the directory containing the `.cali` files and refresh the page.
By default `.cali`` files will be generated in the build directory under
``Benchmarking/output`` and ``Testing/output``. To specify where `.cali` output
files are placed, define the CMake option SUNDIALS_TEST_OUTPUT_DIR with the
desired directory path.

To retain the same filters as a given SPOT visualization link in
:ref:`Bookmarks`, swap out the ``sf`` value in the URL with the 
directory path containing the `.cali` files.

.. _Bookmarks:

Bookmarks to Notable SPOT Visualizations
----------------------------------------------

Links will be added as informative SPOT visualizations are found.

.. future links to add: aggregations of examples
.. 