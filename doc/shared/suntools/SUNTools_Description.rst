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

###########################
SUNDIALS Python tools
###########################

``suntools`` is a Python package containing utilities for analyzing SUNDIALS
output and tuning SUNDIALS applications.  It is distributed with SUNDIALS and
can be installed from the ``suntools`` directory with

.. code-block:: bash

   python -m pip install ./suntools

The package currently provides the following modules:

* :mod:`suntools.logs` parses and filters files produced by
  :c:type:`SUNLogger`.
* :mod:`suntools.csv` reads SUNDIALS statistics in CSV format.
* :mod:`suntools.table` reads SUNDIALS statistics in table format.
* :mod:`suntools.tune` configures and runs parameter searches for SUNDIALS
  executables.

The ``suntools`` command-line program provides the ``parse_logs`` and
``tune`` subcommands.  The complete Python API is documented in
:ref:`SUNTools.API`.

.. note::

   The ``tune`` command requires the tuning dependencies listed by the
   ``suntools`` project.  The Ytopt backend is optional and can be installed
   with ``python -m pip install "suntools[ytopt]"``.
