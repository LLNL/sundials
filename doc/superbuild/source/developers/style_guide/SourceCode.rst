..
   Author(s): David J. Gardner @ LLNL
   -----------------------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2002-2023, Lawrence Livermore National Security
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

All new code added to SUNDIALS should be linted with  `clang-tidy
<https://clang.llvm.org/extra/clang-tidy/>`_ and formatted with `clang-format
<https://clang.llvm.org/docs/ClangFormat.html>`_. The ``.clang-tidy`` and
``.clang-format`` files in the root of the project define our configurations
for the tools respectively.

It may be necessary to override clang-tidy at times. This can be done with the
``NOLINT`` magic comments e.g.,

.. code-block:: cpp

  template<class GkoSolverType, class GkoMatrixType>
  int SUNLinSolFree_Ginkgo(SUNLinearSolver S)
  {
    auto solver{static_cast<LinearSolver<GkoSolverType, GkoMatrixType>*>(S->content)};
    delete solver; // NOLINT
    return SUNLS_SUCCESS;
  }

  class BaseObject {
  protected:
    // NOLINTNEXTLINE(cppcoreguidelines-non-private-member-variables-in-classes)
    SUNContext sunctx_{};
  };

See the clang-tidy documentation for more details.

Indentation
^^^^^^^^^^^

Spaces not tabs

Comments
--------

TODO Comments
^^^^^^^^^^^^^

Following the `Google Style Guide <https://google.github.io/styleguide/>`_ , TODO comments are used
to note code that is "temporary, a short-term solution, or good-enough but not perfect."

A consistent TODO comment format provides an easy to search for keyword with details on how to get
more information. TODO comments should start with ``TODO`` followed by a unique identifier, enclosed
in parentheses, for the person most knowledgeable about the issue and a brief description of the
TODO item. Generally, these comments should be used sparingly and are not a substitute for creating
an issue or bug report. When applicable, the comment should include the relevant issue or bug report
number.

Examples:

.. code-block:: c

   /* TODO(DJG): Update to new API in the next major release (Issue #256) */
