..
   -----------------------------------------------------------------------------
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
   -----------------------------------------------------------------------------

.. _SUNTOOLS:

********
SUNTOOLS
********

``suntools`` is a Python package containing utilities for working with
SUNDIALS applications. It can parse SUNDIALS logging and statistics output and
can tune ``SetOptions`` parameters for SUNDIALS executables.

The package source and installation metadata are in the ``suntools`` directory
of the SUNDIALS source tree. Install it in an environment with:

.. code-block:: bash

   python -m pip install -e suntools

The optional ``ytopt`` tuning backend can be installed with:

.. code-block:: bash

   python -m pip install -e "suntools[ytopt]"

After installation, the Python modules can be imported as ``suntools`` and the
command-line interface is available as ``suntools``.

.. toctree::
   :maxdepth: 1

   logging
   tuning
   SUNTools_links.rst
