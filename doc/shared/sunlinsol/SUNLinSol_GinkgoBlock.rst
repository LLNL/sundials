..
   ----------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2002-2025, Lawrence Livermore National Security
   and Southern Methodist University.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   ----------------------------------------------------------------

.. _SUNLinSol.GinkgoBlock:

The SUNLINEARSOLVER_GINKGOBLOCK Module
=======================================

.. versionadded:: X.Y.Z

The SUNLINEARSOLVER_GINKGOBLOCK implementation of the ``SUNLinearSolver`` API provides an
interface to the batched linear solvers from the Ginkgo linear algebra library :cite:p:`ginkgo-toms-2022`.
Like SUNLINEARSOLVER_GINKGO, this module is written in C++14 and is distributed as a header file.
To use the SUNLINEARSOLVER_GINKGOBLOCK ``SUNLinearSolver``, users will need to include 
``sunlinsol/sunlinsol_ginkgoblock.hpp``. The module is meant to be used with the SUNMATRIX_GINKGOBLOCK 
module described in :numref:`SUNMatrix.GinkgoBlock`. Instructions on building SUNDIALS  with Ginkgo enabled are given
in :numref:`Installation.Options.Ginkgo`.  For instructions on building and using Ginkgo itself, refer to the
`Ginkgo website and documentation <https://ginkgo-project.github.io/>`_.

.. note::

  It is assumed that users of this module are aware of how to use Ginkgo. This module does not
  try to encapsulate Ginkgo linear solvers, rather it provides a lightweight iteroperability layer
  between Ginkgo and SUNDIALS. Most, if not all, of the Ginkgo linear solver should work with this
  interface.

.. _SUNLinSol.GinkgoBlock.Usage:

Using SUNLINEARSOLVER_GINKGOBLOCK
---------------------------------

After choosing a compatible ``N_Vector`` (see :numref:`SUNMatrix.GinkgoBlock.CompatibleNVectors`) and creating a Ginkgo-enabled ``SUNMatrix`` (see
:numref:`SUNMatrix.GinkgoBlock`) to use the SUNLINEARSOLVER_GINKGOBLOCK module, we create the linear solver object:

.. code-block:: cpp

   using GkoBatchMatrixType = gko::batch::matrix::Csr<sunrealtype, sunindextype>;
   using GkoBatchSolverType = gko::batch::solver::Bicgstab<sunrealtype>;
   using SUNGkoMatrixType   = sundials::ginkgo::BlockMatrix<GkoBatchMatrixType>;
   using SUNGkoLinearSolverType =
      sundials::ginkgo::BlockLinearSolver<GkoBatchSolverType, GkoBatchMatrixType>;

   SUNGkoLinearSolverType LS{gko_exec, gko::batch::stop::tolerance_type::absolute,
                             precond_factory, num_batches, sunctx};

Next, we can pass the instance of ``sundials::ginkgo::BlockLinearSolver`` to any function
expecting a ``SUNLinearSolver`` object through the implicit conversion operator or explicit conversion function.

.. code-block:: cpp

   // Attach linear solver and matrix to CVODE.
   //
   // Implicit conversion from sundials::ginkgo::BlockLinearSolver<GkoBatchSolverType, GkoBatchMatrixType>
   // to a SUNLinearSolver object is done.
   //
   // For details about creating A see the SUNMATRIX_GINKGOBLOCK module.
   CVodeSetLinearSolver(cvode_mem, LS, A);

   // Alternatively with explicit conversion of LS to a SUNLinearSolver
   // and A to a SUNMatrix:
   CVodeSetLinearSolver(cvode_mem, LS.Convert(), A.Convert());

After attaching the linear solver to the SUNDIALS integrator, one must change the norm factor the integrator uses
since the Ginkgo linear solver will take norms over individual batches, not the entire system.

.. code-block:: cpp

   // When using ARKODE:
   ARKodeSetLSNormFactor(arkode_mem, std::sqrt(batch_size));

   // When using CVODE:
   CVodeSetLSNormFactor(cvode_mem, std::sqrt(batch_size));

   // When using IDA:
   IDASetLSNormFactor(ida_mem, std::sqrt(batch_size));

.. warning:: 

   Setting the linear solver norm factor is essential. If this is not set, you will likely see lots of linear solver convergence failures.

.. warning::

  :c:func:`SUNLinSolFree` should never be called on a ``SUNLinearSolver`` that was created via conversion
  from a :cpp:type:`sundials::ginkgo::BlockLinearSolver`. Doing so may result in a double free.


.. _SUNLinSol.GinkgoBlock.API:

SUNLINEARSOLVER_GINKGOBLOCK API
-------------------------------

In this section we list the public API of the :cpp:type:`sundials::ginkgo::BlockLinearSolver` class.

.. cpp:class:: template<class GkoBatchSolverType, class GkoBatchMatType> \
               sundials::ginkgo::BlockLinearSolver : public sundials::ConvertibleTo<SUNLinearSolver>

   .. cpp:function:: BlockLinearSolver(std::shared_ptr<const gko::Executor> gko_exec, sunindextype num_blocks, SUNContext sunctx)

      Constructs a new BlockLinearSolver with default tolerance type and max iterations.

      :param gko_exec: The `gko::Executor` to use
      :param num_blocks: Number of batches (block systems)
      :param sunctx: The SUNDIALS simulation context (:c:type:`SUNContext`)

   .. cpp:function:: BlockLinearSolver(std::shared_ptr<const gko::Executor> gko_exec, gko::batch::stop::tolerance_type tolerance_type, sunindextype num_blocks, SUNContext sunctx)

      Constructs a new BlockLinearSolver with specified tolerance type.

      :param gko_exec: The `gko::Executor` to use
      :param tolerance_type: Ginkgo batch solver tolerance type
      :param num_blocks: Number of batches (block systems)
      :param sunctx: The SUNDIALS simulation context (:c:type:`SUNContext`)

   .. cpp:function:: BlockLinearSolver(std::shared_ptr<const gko::Executor> gko_exec, std::shared_ptr<gko::batch::BatchLinOpFactory> precon_factory, sunindextype num_blocks, SUNContext sunctx)

      Constructs a new BlockLinearSolver with a preconditioner factory.

      :param gko_exec: The `gko::Executor` to use
      :param precon_factory: Ginkgo batch preconditioner factory
      :param num_blocks: Number of batches (block systems)
      :param sunctx: The SUNDIALS simulation context (:c:type:`SUNContext`)

   .. cpp:function:: BlockLinearSolver(std::shared_ptr<const gko::Executor> gko_exec, int max_iters, sunindextype num_blocks, SUNContext sunctx)

      Constructs a new BlockLinearSolver with a maximum number of iterations.

      :param gko_exec: The `gko::Executor` to use
      :param max_iters: Maximum number of iterations
      :param num_blocks: Number of batches (block systems)
      :param sunctx: The SUNDIALS simulation context (:c:type:`SUNContext`)

   .. cpp:function:: BlockLinearSolver(std::shared_ptr<const gko::Executor> gko_exec, gko::batch::stop::tolerance_type tolerance_type, int max_iters, sunindextype num_blocks, SUNContext sunctx)

      Constructs a new BlockLinearSolver with specified tolerance type and maximum iterations.

      :param gko_exec: The `gko::Executor` to use
      :param tolerance_type: Ginkgo batch solver tolerance type
      :param max_iters: Maximum number of iterations
      :param num_blocks: Number of batches (block systems)
      :param sunctx: The SUNDIALS simulation context (:c:type:`SUNContext`)

   .. cpp:function:: BlockLinearSolver(std::shared_ptr<const gko::Executor> gko_exec, std::shared_ptr<gko::batch::BatchLinOpFactory> precon_factory, int max_iters, sunindextype num_blocks, SUNContext sunctx)

      Constructs a new BlockLinearSolver with a preconditioner factory and maximum iterations.

      :param gko_exec: The `gko::Executor` to use
      :param precon_factory: Ginkgo batch preconditioner factory
      :param max_iters: Maximum number of iterations
      :param num_blocks: Number of batches (block systems)
      :param sunctx: The SUNDIALS simulation context (:c:type:`SUNContext`)

   .. cpp:function:: BlockLinearSolver(std::shared_ptr<const gko::Executor> gko_exec, gko::batch::stop::tolerance_type tolerance_type, std::shared_ptr<gko::batch::BatchLinOpFactory> precon_factory, sunindextype num_blocks, SUNContext sunctx)

      Constructs a new BlockLinearSolver with specified tolerance type and preconditioner factory.

      :param gko_exec: The `gko::Executor` to use
      :param tolerance_type: Ginkgo batch solver tolerance type
      :param precon_factory: Ginkgo batch preconditioner factory
      :param num_blocks: Number of batches (block systems)
      :param sunctx: The SUNDIALS simulation context (:c:type:`SUNContext`)

   .. cpp:function:: BlockLinearSolver(std::shared_ptr<const gko::Executor> gko_exec, gko::batch::stop::tolerance_type tolerance_type, std::shared_ptr<gko::batch::BatchLinOpFactory> precon_factory, int max_iters, sunindextype num_blocks, SUNContext sunctx)

      Constructs a new BlockLinearSolver with all options specified.

      :param gko_exec: The `gko::Executor` to use
      :param tolerance_type: Ginkgo batch solver tolerance type
      :param precon_factory: Ginkgo batch preconditioner factory
      :param max_iters: Maximum number of iterations
      :param num_blocks: Number of batches (block systems)
      :param sunctx: The SUNDIALS simulation context (:c:type:`SUNContext`)

   .. cpp:function:: BlockLinearSolver(BlockLinearSolver&& that_solver) noexcept

      Move constructor.

   .. cpp:function:: BlockLinearSolver& operator=(BlockLinearSolver&& rhs)

      Move assignment.

   .. cpp:function:: ~BlockLinearSolver() override = default

      Default destructor.

   .. cpp:function:: operator SUNLinearSolver() override

      Implicit conversion to a :c:type:`SUNLinearSolver`.

   .. cpp:function:: operator SUNLinearSolver() const override

      Implicit conversion to a :c:type:`SUNLinearSolver`.

   .. cpp:function:: SUNLinearSolver Convert() override

      Explicit conversion to a :c:type:`SUNLinearSolver`.

   .. cpp:function:: SUNLinearSolver Convert() const override

      Explicit conversion to a :c:type:`SUNLinearSolver`.

   .. cpp:function:: std::shared_ptr<const gko::Executor> GkoExec() const

      Get the ``gko::Executor`` associated with the Ginkgo solver.

   .. cpp:function:: std::shared_ptr<typename GkoBatchSolverType::Factory> GkoFactory()

      Get the underlying Ginkgo solver factory.

   .. cpp:function:: GkoBatchSolverType* GkoSolver()

      Get the underlying Ginkgo solver.

      .. note::

         This will be `nullptr` until the linear solver setup phase.

   .. cpp:function:: int AvgNumIters() const

      Get the number of linear solver iterations in the most recent solve.

   .. cpp:function:: int StddevNumIters() const

      Get the standard deviation of the number of iterations across the batches during the last solve.

   .. cpp:function:: int SumAvgNumIters() const

      Get the running sum of the average number of iterations in this solvers lifetime.

   .. cpp:function:: void SetScalingVectors(N_Vector s1, N_Vector s2)

      Sets the left and right scaling vector to be used.

   .. cpp:function:: GkoBatchSolverType* Setup(BlockMatrix<GkoBatchMatType>* A)

      Setup the linear system.

      :param A: the linear system matrix

      :returns: Pointer to the Ginkgo solver generated from the factory

   .. cpp:function:: gko::LinOp* Solve(N_Vector b, N_Vector x, sunrealtype tol)

      Solve the linear system Ax = b to the specified tolerance.

      :param b: the right-hand side vector
      :param x: the solution vector
      :param tol: the tolerance to solve the system to

      :returns: ``gko::LinOp*`` the solution
