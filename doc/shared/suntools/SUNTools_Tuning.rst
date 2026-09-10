..
   ----------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2025-2026, Lawrence Livermore National Security,
   University of Maryland Baltimore County, and the SUNDIALS contributors.
   Copyright (c) 2013-2025, Lawrence Livermore National Security
   and Southern Methodist University.
   Copyright (c) 2002-2013, Lawrence Livermore National Security.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   ----------------------------------------------------------------

.. _SUNTOOLS.Tuning:

**********
Autotuning
**********

The ``suntools tune`` command searches over SUNDIALS ``SetOptions`` parameters
by appending parameter names and values to an executable command. Parameters
can be supplied on the command line:

.. code-block:: bash

   suntools tune \
      --params arkode.table_names choice:TABLE_A,TABLE_B \
      -- ./arkstep_app

For repeatable or larger searches, put the executable, parameters, search
settings, and objective in a YAML configuration file:

.. code-block:: yaml

   executable:
     command: ./arkstep_app
   parameters:
     - name: arkode.table_names
       type: choice
       values: [TABLE_A, TABLE_B]
   search:
     max_evals: 40
     workers: 1
     output_dir: tune-results
   objective:
     metric: wall_time
     direction: minimize

Run the search with:

.. code-block:: bash

   suntools tune --config tune.yaml

The results directory contains the baseline, best, and worst trial records,
along with CSV and JSON-lines summaries of all trials. Objective and constraint
metrics can be extracted from an executable's output using regular expressions.
