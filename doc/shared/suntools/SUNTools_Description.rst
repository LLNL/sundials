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

.. _SUNTOOLS:

*********************************
suntools - utilities for SUNDIALS
*********************************

``suntools`` is a Python package containing utilities for analyzing SUNDIALS
output and tuning SUNDIALS applications.  It is distributed with SUNDIALS and
can be installed from the ``suntools`` directory with

.. code-block:: bash

   python -m pip install ./suntools

The optional ``ytopt`` tuning backend can be installed with

.. code-block:: bash

   python -m pip install "suntools[ytopt]"

The optional GPTune tuning backend can be installed with

.. code-block:: bash

   python -m pip install "suntools[gptune]"

Both optional backends can be installed together with

.. code-block:: bash

   python -m pip install "suntools[gptune,ytopt]"

The GPTune backend requires a complete GPTune runtime installation.  The
``gptune`` PyPI wheel does not include all of GPTune's runtime dependencies,
including the ``autotune`` package, so installing this extra may not provide a
complete runtime.  Install GPTune using its documented source or Spack
installation, including its Python dependencies.  A quick check is

.. code-block:: bash

   python -c "from autotune.problem import TuningProblem; from GPTune.gptune import GPTune"

The package currently provides the following modules:

* ``suntools.logs`` parses and filters files produced by
  :c:type:`SUNLogger`.
* ``suntools.csv`` reads SUNDIALS statistics in CSV format.
* ``suntools.table`` reads SUNDIALS statistics in table format.
* ``suntools.tune`` configures and runs parameter searches for SUNDIALS
  executables.

The ``suntools`` command-line program provides the ``parse_logs`` and
``tune`` subcommands.  The complete Python API is documented in
:ref:`SUNTools.API`.
