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

.. _Style.Classes:

Classes
========

.. _Style.Classes.Old:

Naming Convention for Vectors, Matrices, and Solvers
----------------------------------------------------

The SUNDIALS vector, matrix, linear solver, and nonlinear solver classes use the
naming convention ``<short class name><method>`` for base class methods where
each component of the name uses camel case. See
:numref:`Style.Table.OldBaseClassMethodNaming` for examples.

.. note::

   This naming convention *only* applies to the vector, matrix, and solver
   classes. All other classes should follow the naming convention described in
   :ref:`Style.Classes.New`.

.. _Style.Table.OldBaseClassMethodNaming:

.. Table:: SUNDIALS base class naming convention examples for vectors, matrices,
           linear solvers and nonlinear solvers.

   +-----------------------+------------+------------+-----------------------+
   | Base Class            | Short Name | Operation  | Method                |
   +-----------------------+------------+------------+-----------------------+
   | ``N_Vector``          | ``N_V``    | Linear Sum | ``N_VLinearSum``      |
   +-----------------------+------------+------------+-----------------------+
   | ``SUNMatrix``         | ``SUNMat`` | Zero       | ``SUNMatZero``        |
   +-----------------------+------------+------------+-----------------------+
   | ``SUNLinearSolver``   | ``Setup``  | Setup      | ``SUNLinSolSetup``    |
   +-----------------------+------------+------------+-----------------------+
   | ``SUNNonlinarSolver`` | ``Solve``  | Solve      | ``SUNNonlinSolSolve`` |
   +-----------------------+------------+------------+-----------------------+

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

.. _Style.Classes.New:

Naming Convention for New Classes
---------------------------------

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
