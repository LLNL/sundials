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

.. _SUNTools.API:

suntools API
============

The API documentation below is extracted from the RST-formatted docstrings in
the ``suntools`` source.  The parser functions return ordinary Python
dictionaries and lists, so their output can be passed directly to analysis or
plotting packages.

Log parsing
-----------

.. autofunction:: suntools.logs.log_file_to_list

.. autofunction:: suntools.logs.print_log

.. autofunction:: suntools.logs.get_history

.. autoclass:: suntools.logs.StepData
   :members:

CSV and table parsing
---------------------

.. autofunction:: suntools.csv.keys

.. autofunction:: suntools.csv.read

.. autofunction:: suntools.csv.write

.. autofunction:: suntools.table.parse_table

.. autofunction:: suntools.utils.str2num

Command-line API
----------------

.. autofunction:: suntools.cli.build_parser

.. autofunction:: suntools.cli.main

.. autofunction:: suntools.tune.cli.run_from_args
