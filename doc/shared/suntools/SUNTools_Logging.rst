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

.. _SUNTOOLS.Logging:

****************
Log Manipulation
****************

The ``suntools.logs`` module provides Python helpers for reading the structured
messages emitted by :c:type:`SUNLogger`. Parsed records retain the logging
level, MPI rank, scope, label, and payload, making them convenient for scripts
that inspect solver behavior or generate plots.

The command-line interface can filter a log by solver region. For example, the
following command keeps integrator, nonlinear-solver, and linear-solver output
(the default filter):

.. code-block:: bash

   suntools parse_logs path/to/sundials.log

Select a comma-separated subset with ``--filter`` and invert the selection with
``--invert``. A path of ``-`` reads from standard input:

.. code-block:: bash

   suntools parse_logs --filter nonlinear,linear path/to/sundials.log
   cat path/to/sundials.log | suntools parse_logs --invert --filter linear -

The Python parser can be used directly when a structured representation is
needed:

.. code-block:: python

   from suntools import logs

   records = logs.log_file_to_list("path/to/sundials.log")

The :mod:`suntools.table` and :mod:`suntools.csv` modules provide helpers for
reading table-formatted and CSV-formatted SUNDIALS statistics, respectively.
