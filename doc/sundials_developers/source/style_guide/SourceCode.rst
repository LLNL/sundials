..
   Author(s): David J. Gardner @ LLNL
   -----------------------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2002-2021, Lawrence Livermore National Security
   and Southern Methodist University.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   -----------------------------------------------------------------------------

.. _Style.Code:

Coding Style
============

Formatting
----------

Indentation
^^^^^^^^^^^

Spaces not tabs

Line Length
^^^^^^^^^^^

Loops
^^^^^

Comments
--------

Function Comments
^^^^^^^^^^^^^^^^^

Implementation Comments
^^^^^^^^^^^^^^^^^^^^^^^

TODO Comments
^^^^^^^^^^^^^

Following the Google Style Guide [GoogleStyle]_, TODO comments are used to note
code that is "temporary, a short-term solution, or good-enough but not perfect."

A consistent TODO comment format provides an easy to search for keyword with
details on how to get more information. TODO comments should start with ``TODO``
followed by a unique identifier, enclosed in parentheses, for the person most
knowledgeable about the issue and a brief description of the TODO item.
Generally, these comments should be used sparingly and are not a substitute for
creating an issue or bug report. When applicable, the comment should include the
relevant issue or bug report number.

Examples:

.. code-block:: c

   /* TODO(DJG): Update to new API in the next major release (Issue #256) */
