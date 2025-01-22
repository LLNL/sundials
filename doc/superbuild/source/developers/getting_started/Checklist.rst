..
   Author(s): David J. Gardner @ LLNL
   -----------------------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2002-2025, Lawrence Livermore National Security
   and Southern Methodist University.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   -----------------------------------------------------------------------------

.. _DevelopmentChecklist:

Development Checklist
=====================

When adding new functionality to SUNDIALS, in the form of bugfixes, new solvers,
new solver options, new test problems, modifications to the SUNDIALS build
system, etc. developers should adhere to the following checklist.

#. When adding new source code files, ensure that the corresponding
   ``CMakeLists.txt`` file(s) is/are updated to include your file(s). If your
   files depend on external libraries, emulate the CMake options used for other
   third party libraries to ensure that these files are only included when
   SUNDIALS is configured to use that library.

#. It can be helpful to configure SUNDIALS using the C flags ``-Wall -Werror``.
   When building, modify your file to remove any error/warning messages output
   during compilation of your code. This can be done automatically with the
   advanced CMake options ``ENABLE_ALL_WARNINGS`` and
   ``ENABLE_WARNINGS_AS_ERRORS``.

#. Configure your build with a minimal set of configuration options enabled
   (serial). Run ``make``, ``make test``, ``make install``, and
   ``make test_install`` to ensure that everything functions smoothly, when
   no external libraries are supplied.

#. Configure your build with a maximal set of configuration options enabled
   (OpenMP, Pthreads, MPI, Lapack, KLU, etc.) and run ``make``, ``make test``,
   ``make install``, and ``make test_install`` to ensure that everything
   functions smoothly when external libraries are supplied.

#. When implementing a bug-fix to an existing package/solver:

   * If this change affects the user interface, update the documentation to
     describe the updated interface.
   * For any problems that now fail when running "make test", verify that the
     bugfix either made only cosmetic changes or improved the results, and
     update test output file in the ``examples/`` directory.

#. When adding new solvers or new solver options:

   * Update the documentation to include descriptions of your work.
   * Add a new example problem (or multiple problems) to the ``examples/``
     directory to demonstrate how to use your solver/option, and to include in
     SUNDIALS' automated nightly tests.
   * For any problems that now fail when running ``make test``, verify that your
     updates either made only cosmetic changes or improved the results, and
     update test output in the ``examples/`` directory.

#. When adding new test or example problems to the ``examples/`` directory:

   * Ensure that your new file has a unique name.
   * After running ``make install``, change to the installed ``examples/``
     directory and ensure that ``make`` succeeds, since the CMake-generated
     Makefile system differs from how the examples are built within SUNDIALS.
   * Ensure that the reference output is included e.g., if a file ``foo.c`` is
     added, also add ``foo.out``. 
   * Update the example problem documentation for to include a description of
     the new problem.

#. When adding any new files, update the corresponding package script in the
   ``scripts/`` directory to include your file(s) within the distribution.

#. If answer files changed, and it is expected/desired, then update the `.out` files
   that are embedded in the `examples/` directory AND the
   `"answers" repository <https://github.com/sundials-codes/answers>`_. 

#. If you changed any header files, re-run SWIG to generate updated fortran interfaces.
   This is done by navigating to the `swig/` directory and running `make all32 all64`.
   If you do not have `swig` installed on your system, you can obtain a git patch file
   from the SWIG GitHub action that we run on all pull requests. The patch can be found
   under the job artifacts (if there were in fact changes that required updates
   to the Fortran).
