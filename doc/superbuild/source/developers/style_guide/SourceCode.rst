..
   Author(s): David J. Gardner, Cody J. Balos @ LLNL
   -----------------------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2002-2023, Lawrence Livermore National Security
   and Southern Methodist University.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   -----------------------------------------------------------------------------

.. _Style.Naming:

Naming
======

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

Generally Pascal case (e.g. ``DoSomething``) is used for public names and
camelcase for private names (e.g. ``doSomething``).

Macros/Constants
----------------

Upper case should be used for macros and constants.

Variable names
--------------

Snake case is preferred for local variable names e.g. ``foo_bar``.

C function names
----------------

Functions should have a descriptive name that lets a reader know what it does.
For functions that return a boolean value, prefer the convention
``Is<statement>``, e.g. ``IsOutputRank``. Pascal case (with the appropriate
namespace prefix) should be used for all public function names that are not
class functions (see :ref:`Style.Classes` for class naming conventions).

C++ function names
------------------

All functions should be in the ``sundials::`` namespace. Pascal case should be
used for public function names. Camelcase should be used for private function
names.

Names for Vectors, Matrices, and Solvers
----------------------------------------

The SUNDIALS vector, matrix, linear solver, and nonlinear solver classes use the
naming convention ``<short class name><method>`` for base class methods where
each component of the name uses Pascal case. See
:numref:`Style.Table.OldBaseClassMethodNaming` for examples.

.. note::

   This naming convention *only* applies to the vector, matrix, and solver
   classes. All other classes should follow the naming convention described in
   :ref:`Style.Naming.NewClasses`.

.. _Style.Table.OldBaseClassMethodNaming:

.. Table:: SUNDIALS base class naming convention examples for vectors, matrices,
           linear solvers and nonlinear solvers.

   +-----------------------+------------------+------------+-----------------------+
   | Base Class            | Short Name       | Operation  | Method                |
   +-----------------------+------------------+------------+-----------------------+
   | ``N_Vector``          | ``N_V``          | Linear Sum | ``N_VLinearSum``      |
   +-----------------------+------------------+------------+-----------------------+
   | ``SUNMatrix``         | ``SUNMat``       | Zero       | ``SUNMatZero``        |
   +-----------------------+------------------+------------+-----------------------+
   | ``SUNLinearSolver``   | ``SUNLinSol``    | Setup      | ``SUNLinSolSetup``    |
   +-----------------------+------------------+------------+-----------------------+
   | ``SUNNonlinarSolver`` | ``SUNNonlinSol`` | Solve      | ``SUNNonlinSolSolve`` |
   +-----------------------+------------------+------------+-----------------------+

Derived class implementations of the base class methods should follow the naming
convention ``<short class name><method>_<implementation>``. See
:numref:`Style.Table.OldDerivedClassMethodNaming` for examples.

.. _Style.Table.OldDerivedClassMethodNaming:

.. Table:: SUNDIALS derived class naming convention examples for vectors,
           matrices, linear solvers and nonlinear solvers.

   +---------------+-----------------------+------------------------------+
   | Derived Class | Base Class Method     | Method Implementation        |
   +---------------+-----------------------+------------------------------+
   | Serial        | ``N_VLinearSum``      | ``N_VLinearSum_Serial``      |
   +---------------+-----------------------+------------------------------+
   | Dense         | ``SUNMatZero``        | ``SUNMatZero_Dense``         |
   +---------------+-----------------------+------------------------------+
   | SPGMR         | ``SUNLinSolSetup``    | ``SUNLinSolSetup_SPGMR``     |
   +---------------+-----------------------+------------------------------+
   | Newton        | ``SUNNonlinSolSolve`` | ``SUNNonlinSolSolve_Newton`` |
   +---------------+-----------------------+------------------------------+

Implementation specific methods do not currently have a consistent naming
convention across the different derived classes. When adding new methods to an
existing class, follow the naming style used within that class. When adding a
new derived class, use the same style as above for implementations of the base
class method i.e., ``<short class name><method>_<implementation>``.

.. _Style.Naming.NewClasses:

Names for New Classes
---------------------

All new base classes should use the naming convention ``<class name>_<method>``
for the base class methods. See
:numref:`Style.Table.NewBaseClassMethodNaming` for examples.

.. _Style.Table.NewBaseClassMethodNaming:

.. Table:: SUNDIALS naming conventions for methods in new base classes.

   +-----------------------+------------+---------------------------+
   | Base Class            | Operation  | Method                    |
   +-----------------------+------------+---------------------------+
   | ``SUNMemoryHelper``   | Alloc      | ``SUNMemoryHelper_Alloc`` |
   +-----------------------+------------+---------------------------+

Derived class implementations of the base class methods should follow the naming
convention  ``<class name>_<method>_<implementation>``. See
:numref:`Style.Table.NewDerivedClassMethodNaming` for examples.

.. _Style.Table.NewDerivedClassMethodNaming:

.. Table:: SUNDIALS naming conventions for derived class implementations of
           methods in new base classes.

   +---------------+---------------------------+--------------------------------+
   | Derived Class | Base Class Method         | Method Implementation          |
   +---------------+---------------------------+--------------------------------+
   | CUDA          | ``SUNMemoryHelper_Alloc`` | ``SUNMemoryHelper_Alloc_Cuda`` |
   +---------------+---------------------------+--------------------------------+

For destructor functions, use ``Destroy`` rather than ``Free`` or some other alternative.


.. _Style.Classes.Cpp:

Naming Convention for C++ Classes
---------------------------------

C++ classes should have a descriptive name. The class name should not be
prefixed with ``SUN``, but it should reside in the ``sundials::`` namespace.
Public C++ class functions should use Pascal case (e.g. ``DoSomething``).
Private C++ class functions should use camelcase (e.g. ``doSomething``).

C++ private class members should use snake case with a trailing underscore
(e.g. ``some_var_``).

.. _Style.Code:

Coding Conventions and Rules 
============================

#. Do not use language features that are not compatible with C99, C++14,
   and MSVC v1900+ (Visual Studio 2015). Examples of such features include
   variable-length arrays. Exceptions are allowed when interfacing with a
   library which requires a newer standard.

#. All new code added to SUNDIALS should be linted with  `clang-tidy
   <https://clang.llvm.org/extra/clang-tidy/>`_ and formatted with `clang-format
   <https://clang.llvm.org/docs/ClangFormat.html>`_. See :ref:`Style.Formatting`
   for details.

#. If clang-format breaks lines in a way that is unreadable, use ``//`` to break the line. For example,
   sometimes (mostly in C++ code) you may have code like this:

   .. code-block:: cpp

      MyClass::callAFunctionOfSorts::doSomething().doAnotherThing().doSomethingElse();

   That you would like to format as (for readability):

   .. code-block:: cpp

      MyObject::callAFunctionOfSorts()
            .doSomething()
            .doAnotherThing()
            .doSomethingElse();

   Clang-format might produce something like:
   
   .. code-block:: cpp

      MyObject::callAFunctionOfSorts().doSomething().doAnotherThing()
            .doSomethingElse();
   ```

   unless you add the `//`.
  
   .. code-block:: cpp

      MyObject::callAFunctionOfSorts()
            .doSomething()       //
            .doAnotherThing()    //
            .doSomethingElse();  //
  
   There are other scenarios (e.g., a function call with a lot of parameters) where doing this type of line break is useful for readability too.

#. It may be necessary to override clang-tidy at times. This can be done with
   the ``NOLINT`` magic comments e.g.,

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

#. Spaces not tabs.

#. Comments should use proper spelling and grammar.

#. Following the Google Style Guide [GoogleStyle]_, TODO comments are used to note
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

#. All SUNDIALS data structures should hold onto a ``SUNContext`` object. Exceptions 
   are the ``SUNLogger`` and ``SUNProfiler`` classes.

#. All SUNDIALS functions should return a ``SUNErrCode``. Many older functions
   do not do this and are exceptions to the rule for backwards compatiblilty. 
   In addition, internal helper functions may or may-not return a ``SUNErrCode``.

#. All SUNDIALS functions, with the exception of some functions
   that do not have access to a ``SUNContext``, should begin with a call to 
   ``SUNFunctionBegin()``. The argument to ``SUNFunctionBegin()`` is a ``SUNContext``
   which should come from the first object in the function parameter list that has a
   ``SUNContext``.  This macro is used for error handling and declares a 
   local variable access via the macro ``SUNCTX_``.

   .. code-block:: c

      SUNErrCode N_VLinearCombination_Serial(int nvec, realtype* c, N_Vector* X, N_Vector z)
      {
         SUNFunctionBegin(X[0]->sunctx); // Correct
         
         int          i;
         sunindextype j, N;
         realtype*    zd=NULL;
         realtype*    xd=NULL;

         /* invalid number of vectors */
         SUNAssert(nvec >= 1, SUN_ERR_ARG_OUTOFRANGE);

         // ...
      }

      SUNErrCode N_VLinearCombination_Serial(int nvec, realtype* c, N_Vector* X, N_Vector z)
      {
         int          i;
         sunindextype j, N;
         realtype*    zd=NULL;
         realtype*    xd=NULL;

         SUNFunctionBegin(X[0]->sunctx); // Incorrect, SUNFunctionBegin should occur as early as possible  

         /* invalid number of vectors */
         SUNAssert(nvec >= 1, SUN_ERR_ARG_OUTOFRANGE);

         // ...
      }


      int CVodeGetEstLocalErrors(void *cvode_mem, N_Vector ele)
      {
         CVodeMem cv_mem;

         if (cvode_mem==NULL) {
            cvProcessError(NULL, CV_MEM_NULL, __LINE__, __func__, __FILE__, MSGCV_NO_MEM);
            return(CV_MEM_NULL);
         }

         cv_mem = (CVodeMem) cvode_mem;

         SUNFunctionBegin(cv_mem->sunctx); // Correct - this is as early as possible to call SUNFunctionBegin
         SUNFunctionBegin(ele->sunctx); // Incorrect - cvode_mem is first in the function parameter list

         // ...
      }


#. All references to ``SUNContext`` objects should be done via the ``SUNCTX_``
   macro. The only exceptions are functions in the ``SUNContext`` class. 

#. All calls to SUNDIALS functions that return a ``SUNErrCode`` should have
   their return value checked with a macro from the ``SUNCheckCall`` family.
   These macros are documented in the header file ``sundials/priv/sundials_errors_impl.h``
   and ``sundials/priv/sundials_mpi_errors_impl.h``.

   .. code-block:: c

    SUNCheckCall(N_VLinearCombination(...)); // Correct

   Avoid storing the return value and then checking the stored value except when absolutely necessary.

   .. code-block:: c

    SUNErrCode err;
    err = N_VLinearCombination(...);
    SUNCheckCall(err); // Avoid except when absolutely necessary.

#. All calls to SUNDIALS functions that *do not* return a ``SUNErrCode`` should
   be followed by checking the last error stored in the ``SUNContext``. 
   The exception to this rule is for internal helper functions. 
   These should not be checked unless they return a ``SUNErrCode``. 
   These checks are done with the ``SUNCheckLastErr`` macros. 

   .. code-block:: c

    // Correct 
    (void) N_VLinearSum(...); SUNCheckLastErr(); 

    // Incorrect - we must check for errors in N_VDotProd before calling a second function
    sunrealtype norm2 = SUNRsqrt(N_VDotProd(...)); SUNCheckLastErr(); 

    // Correct
    sunrealtype norm2 = N_VDotProd(...); SUNCheckLastErr(); 
    norm2 = SUNRsqrt(norm2);

#. Programmer errors should be checked with the ``SUNAssert`` macro, that verifies whether its
   argument evaluates to "true", and specifies an error flag otherwise. By programmer errors we
   mean, for example, illegal inputs such as mismatching dimensions or a ``NULL`` value for
   something that should not be. 

   .. code-block:: c
      
      SUNLinearSolver SUNLinSol_Band(N_Vector y, SUNMatrix A, SUNContext sunctx)
      {
         SUNFunctionBegin(sunctx);
         SUNLinearSolver S;
         SUNLinearSolverContent_Band content;
         sunindextype MatrixRows;

         // Correct - check these with SUNAssert
         SUNAssert(SUNMatGetID(A) == SUNMATRIX_BAND, SUN_ERR_ARG_WRONGTYPE);
         SUNAssert(SUNBandMatrix_Rows(A) == SUNBandMatrix_Columns(A), SUN_ERR_ARG_DIMSMISMATCH);
         SUNAssert(y->ops->nvgetarraypointer, SUN_ERR_ARG_INCOMPATIBLE);

         // ...
      }

#. If statements and loops should always have braces even if they are one line.

#. Return statements should not unecessarily use parentheses. Prefer ``return
   x;`` to ``return(x);``. Note, however, lots of older SUNDIALS source code
   uses ``return(x);``. 


.. _Style.Formatting:

Formatting
----------

All new code added to SUNDIALS should be linted with  `clang-tidy
<https://clang.llvm.org/extra/clang-tidy/>`_ and formatted with `clang-format
<https://clang.llvm.org/docs/ClangFormat.html>`_. The ``.clang-tidy`` and
``.clang-format`` files in the root of the project define our configurations
for the tools respectively. To run both tools on a directory, you can run

.. code-block::

   cd scripts
   ./tidy_and_format.sh <path to directory to tidy and format> <path to builddir>

For clang-tidy to work, the build directory must have a ``compile_commands.json``, which 
can be generated by setting the CMake option ``CMAKE_EXPORT_COMPILE_COMMANDS``
to ``ON``. If you just want to run the formatter, then you can run

.. code-block::

   cd scripts
   ./format.sh <path to directory to format>

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
