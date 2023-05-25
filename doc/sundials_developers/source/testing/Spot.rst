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

In order to protect against performance regression of SUNDIALS at all scales, we leverage the GitLab CI infrastructure setup for CI testing to perform continuous performance testing. 

To track and visualize performance over time, we collect data with `Caliper <https://lc.llnl.gov/confluence/display/CALI/Spot+DB>`_ and `Adiak <https://github.com/LLNL/Adiak>`_. Steps to setting up continuous performance testing(CPT) can be found at :ref:`Continuous Performance Testing<CPT>`.

Documentation to using the SPOT performance viisualization tool can be found at `here<https://lc.llnl.gov/confluence/display/SpotDoc/Spot+Documentation>`_.