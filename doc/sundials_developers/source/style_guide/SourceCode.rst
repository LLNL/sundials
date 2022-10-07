..
   Author(s): David J. Gardner @ LLNL
   -----------------------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2002-2022, Lawrence Livermore National Security
   and Southern Methodist University.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   -----------------------------------------------------------------------------

.. _Style.Code:

Coding Style
============

Naming conventions
------------------

All exported symbols that will be publically available must be namespaced
appropriately!

- ``SUN_`` or ``SUNDIALS_`` for macros
- ``sun`` for typedef's to native types (e.g., ``sunindextype``)
- ``SUN`` for public functions that are not class functions (see
  :numref:`Style.Classes` for class/struct naming conventions)
- ``sun`` for private functions that are non-native.
- ``sundials::`` for C++ code (nesting under `sundials::` is OK)
- ``sundials::<somename>::impl`` for C++ code that is private (implementation
  only)

Themes
^^^^^^

Generally Pascal case (e.g. ``DoSomething``) is used for public names and
camelcase for private names (e.g. ``doSomething``).

Macros/Constants
^^^^^^^^^^^^^^^^

Upper case should be used for macros and constants.

Variable names
^^^^^^^^^^^^^^

Snake case is preferred for local variable names e.g. ``foo_bar``.

C function names
^^^^^^^^^^^^^^^^

Functions should have a descriptive name that lets a reader know what it does.
For functions that return a boolean value, prefer the convention
``Is<statement>``, e.g. ``IsOutputRank``. Pascal case (with the appropriate
namespace prefix) should be used for all public function names that are not
class functions (see :ref:`Style.Classes` for class naming conventions).

C++ function names
^^^^^^^^^^^^^^^^^^

All functions should be in the ``sundials::`` namespace. Pascal case should be
used for public function names. Camelcase should be used for private function
names.


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
