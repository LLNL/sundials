..
   Programmer(s): Daniel R. Reynolds @ SMU
   ----------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2002-2021, Lawrence Livermore National Security
   and Southern Methodist University.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   ----------------------------------------------------------------

.. _Usage.ARKStep.Skeleton:

A skeleton of the user's main program
============================================

The following is a skeleton of the user's main program (or calling
program) for the integration of an ODE IVP using the ARKStep module.
Most of the steps are independent of the NVECTOR, SUNMATRIX, SUNLINSOL
and SUNNONLINSOL implementations used.  For the steps that are not,
refer to sections :numref:`NVectors`, :numref:`SUNMatrix`,
:numref:`SUNLinSol`, and  :numref:`SUNNonlinSol` for the specific name of
the function to be called or macro to be referenced.

.. index:: User main program

#. Initialize parallel or multi-threaded environment, if appropriate.

   For example, call ``MPI_Init`` to initialize MPI if used, or set
   ``num_threads``, the number of threads to use within the threaded
   vector functions, if used.

#. Create the SUNDIALS simulation context object.

   Call :c:func:`SUNContext_Create` to create a ``SUNContext`` object,
   providing an MPI communicator if using MPI. See the section
   :ref:`SUNDIALS.SUNContext`` for more information.

   This step may look like

   .. code-block:: c

      SUNContext ctx;
      SUNContext_Create(NULL, &ctx);

   when not using MPI, or like

   .. code-block:: c

      MPI_Comm comm = MPI_COMM_WORLD;
      SUNContext ctx;
      SUNContext_Create((void*) &comm, &ctx);

   when using MPI.

#. Set problem dimensions, etc.

   This generally includes the problem size, ``N``, and may include
   the local vector length ``Nlocal``.

   .. note::

      The variables ``N`` and ``Nlocal`` should be of type
      ``sunindextype``.

#. Set vector of initial values

   To set the vector ``y0`` of initial values, use the appropriate
   functions defined by the particular NVECTOR implementation.

   For native SUNDIALS vector implementations (except the CUDA and
   RAJA based ones), use a call of the form

   .. code-block:: c

      y0 = N_VMake_***(..., ydata);

   if the ``realtype`` array ``ydata`` containing the initial values of
   :math:`y` already exists.  Otherwise, create a new vector by making
   a call of the form

   .. code-block:: c

      y0 = N_VNew_***(...);

   and then set its elements by accessing the underlying data where it
   is located with a call of the form

   .. code-block:: c

      ydata = N_VGetArrayPointer_***(y0);

   See sections :numref:`NVectors.NVSerial` through
   :numref:`NVectors.Pthreads` for details.

   For the HYPRE and PETSc vector wrappers, first create and initialize
   the underlying vector, and then create the NVECTOR wrapper with a call
   of the form

   .. code-block:: c

      y0 = N_VMake_***(yvec);

   where ``yvec`` is a HYPRE or PETSc vector.  Note that calls like
   ``N_VNew_***(...)`` and ``N_VGetArrayPointer_***(...)`` are not
   available for these vector wrappers.  See sections
   :numref:`NVectors.ParHyp` and :numref:`NVectors.NVPETSc` for details.

   If using either the CUDA- or RAJA-based vector implementations use
   calls to the module-specific routines

   .. code-block:: c

      y0 = N_VMake_***(...);

   as applicable.  See sections
   :numref:`NVectors.CUDA` and :numref:`NVectors.RAJA` for details.

#. Create ARKStep object

   Call ``arkode_mem = ARKStepCreate(...)`` to create the ARKStep memory
   block. :c:func:`ARKStepCreate()` returns a ``void*`` pointer to
   this memory structure. See section
   :numref:`Usage.ARKStep.Initialization` for details.

#. Specify integration tolerances

   Call :c:func:`ARKStepSStolerances()` or
   :c:func:`ARKStepSVtolerances()` to specify either a scalar relative
   tolerance and scalar absolute tolerance, or a scalar relative
   tolerance and a vector of absolute tolerances,
   respectively.  Alternatively, call :c:func:`ARKStepWFtolerances()`
   to specify a function which sets directly the weights used in
   evaluating WRMS vector norms. See section
   :numref:`Usage.ARKStep.Tolerances` for details.

   If a problem with non-identity mass matrix is used, and the
   solution units differ considerably from the equation units,
   absolute tolerances for the equation residuals (nonlinear and
   linear) may be specified separately through calls to
   :c:func:`ARKStepResStolerance()`, :c:func:`ARKStepResVtolerance()`, or
   :c:func:`ARKStepResFtolerance()`.

#. Create matrix object

   If a nonlinear solver requiring a linear solver will be used (e.g.,
   a Newton iteration) and the linear solver will be a matrix-based linear
   solver, then a template Jacobian matrix must be created by using the
   appropriate functions defined by the particular SUNMATRIX
   implementation.

   For the SUNDIALS-supplied SUNMATRIX implementations, the
   matrix object may be created using a call of the form

   .. code-block:: c

      SUNMatrix A = SUNBandMatrix(...);

   or similar for the other matrix modules (see section :numref:`SUNMatrix` for
   further information).

   Similarly, if the problem involves a non-identity mass matrix, and
   the mass-matrix linear systems will be solved using a direct linear
   solver, then a template mass matrix must be created by using the
   appropriate functions defined by the particular SUNMATRIX
   implementation.

   .. note::

      The dense, banded, and sparse matrix objects are usable only in a
      serial or threaded environment..


#. Create linear solver object

   If a nonlinear solver requiring a linear solver will be used (e.g.,
   a Newton iteration), or if the problem involves a non-identity mass
   matrix, then the desired linear solver object(s) must be created by
   using the appropriate functions defined by the particular SUNLINSOL
   implementation.

   For any of the SUNDIALS-supplied SUNLINSOL implementations, the
   linear solver object may be created using a call of the form

   .. code-block:: c

      SUNLinearSolver LS = SUNLinSol_*(...);

   where ``*`` can be replaced with "Dense", "SPGMR", or other
   options, as discussed in section :numref:`SUNLinSol`.

#. Set linear solver optional inputs

   Call ``*Set*`` functions from the selected linear solver module
   to change optional inputs specific to that linear solver.  See the
   documentation for each SUNLINSOL module in section
   :numref:`SUNLinSol` for details.

#. Attach linear solver module

   If a linear solver was created above for implicit stage solves,
   initialize the ARKLS linear solver interface by attaching the
   linear solver object (and Jacobian matrix object, if applicable)
   with the call (for details see section :numref:`Usage.ARKStep.LinearSolvers`):

   .. code-block:: c

      ier = ARKStepSetLinearSolver(...);

   Similarly, if the problem involves a non-identity mass matrix,
   initialize the ARKLS mass matrix linear solver interface by
   attaching the mass linear solver object (and mass matrix object,
   if applicable) with the call (for details see section
   :numref:`Usage.ARKStep.LinearSolvers`):

   .. code-block:: c

      ier = ARKStepSetMassLinearSolver(...);

#. Create nonlinear solver object

   If the problem involves an implicit component, and if a non-default
   nonlinear solver object will be used for implicit stage solves
   (see section :numref:`Usage.ARKStep.NonlinearSolvers`),
   then the desired nonlinear solver object must be created by using
   the appropriate functions defined by the particular SUNNONLINSOL
   implementation (e.g., ``NLS = SUNNonlinSol_***(...);`` where
   ``***`` is the name of the nonlinear solver (see section
   :numref:`SUNNonlinSol` for details).

   For the SUNDIALS-supplied SUNNONLINSOL implementations, the
   nonlinear solver object may be created using a call of the form

   .. code-block:: c

      SUNNonlinearSolver NLS = SUNNonlinSol_*(...);

   where ``*`` can be replaced with "Newton", "FixedPoint", or other
   options, as discussed in section :numref:`SUNNonlinSol`.

#. Attach nonlinear solver module

   If a nonlinear solver object was created above, then it must be
   attached to ARKStep using the call (for details see section
   :numref:`Usage.ARKStep.NonlinearSolvers`):

   .. code-block:: c

      ier = ARKStepSetNonlinearSolver(...);

#. Set nonlinear solver optional inputs

   Call the appropriate set functions for the selected nonlinear
   solver module to change optional inputs specific to that nonlinear
   solver.  These *must* be called after attaching the nonlinear
   solver to ARKStep, otherwise the optional inputs will be
   overridden by ARKStep defaults.  See section
   :numref:`SUNNonlinSol` for more information on optional inputs.

#. Set optional inputs

   Call ``ARKStepSet*`` functions to change any optional inputs that
   control the behavior of ARKStep from their default values. See
   section :numref:`Usage.ARKStep.OptionalInputs` for details.

#. Specify rootfinding problem

   Optionally, call :c:func:`ARKStepRootInit()` to initialize a rootfinding
   problem to be solved during the integration of the ODE system. See
   section :numref:`Usage.ARKStep.RootFinding` for general details, and
   section :numref:`Usage.ARKStep.OptionalInputs` for relevant optional
   input calls.

#. Advance solution in time

   For each point at which output is desired, call

   .. code-block:: c

      ier = ARKStepEvolve(arkode_mem, tout, yout, &tret, itask);

   Here, ``itask`` specifies the return mode. The vector ``yout``
   (which can be the same as the vector ``y0`` above) will contain
   :math:`y(t_\text{out})`. See section
   :numref:`Usage.ARKStep.Integration` for details.

#. Get optional outputs

   Call ``ARKStepGet*`` functions to obtain optional output. See
   section :numref:`Usage.ARKStep.OptionalOutputs` for details.

#. Deallocate memory for solution vector

   Upon completion of the integration, deallocate memory for the
   vector ``y`` (or ``yout``) by calling the destructor function:

   .. code-block:: c

      N_VDestroy(y);

#. Free solver memory

   Call :c:func:`ARKStepFree()` to free the memory allocated for
   the ARKStep module (and any nonlinear solver module).

#. Free linear solver and matrix memory

   Call :c:func:`SUNLinSolFree()` and (possibly)
   :c:func:`SUNMatDestroy()` to free any memory allocated for the
   linear solver and matrix objects created above.

#. Free nonlinear solver memory

   If a user-supplied ``SUNNonlinearSolver`` was provided to ARKStep,
   then call :c:func:`SUNNonlinSolFree()` to free any memory allocated
   for the nonlinear solver object created above.

#. Finalize MPI, if used

   Call ``MPI_Finalize`` to terminate MPI.
