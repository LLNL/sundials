..
   Programmer(s): David J. Gardner @ LLNL
                  Daniel R. Reynolds @ SMU
   ----------------------------------------------------------------
   Based on ERKStep by Daniel R. Reynolds @ SMU
   ----------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2002-2020, Lawrence Livermore National Security
   and Southern Methodist University.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   ----------------------------------------------------------------

:tocdepth: 3


.. _MRIStep_CInterface.Skeleton:

A skeleton of the user's main program
============================================

The following is a skeleton of the user's main program (or calling
program) for the integration of an ODE IVP using the MRIStep module.
Most of the steps are independent of the NVECTOR, SUNMATRIX, SUNLINSOL
and SUNNONLINSOL implementations used.
For the steps that are not, refer to the sections :ref:`NVectors`, :ref:`SUNMatrix`,
:ref:`SUNLinSol`, and  :ref:`SUNNonlinSol` for
the specific name of the function to be called or macro to be
referenced.

.. index:: User main program

#. Initialize parallel or multi-threaded environment, if appropriate.

   For example, call ``MPI_Init`` to initialize MPI if used, or set
   ``num_threads``, the number of threads to use within the threaded
   vector functions, if used.

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

   See the sections :ref:`NVectors.NVSerial` through
   :ref:`NVectors.Pthreads` for details.

   For the HYPRE and PETSc vector wrappers, first create and initialize
   the underlying vector, and then create the NVECTOR wrapper with a call
   of the form

   .. code-block:: c

      y0 = N_VMake_***(yvec);

   where ``yvec`` is a HYPRE or PETSc vector.  Note that calls like
   ``N_VNew_***(...)`` and ``N_VGetArrayPointer_***(...)`` are not
   available for these vector wrappers.  See the sections
   :ref:`NVectors.ParHyp` and :ref:`NVectors.NVPETSc` for details.

   If using either the CUDA- or RAJA-based vector implementations use
   calls to the module-specific routines

   .. code-block:: c

      y0 = N_VMake_***(...);

   as applicable.  See the sections
   :ref:`NVectors.CUDA` and :ref:`NVectors.RAJA` for details.

#. Create an ARKStep object for the fast (inner) integration

   Call ``inner_arkode_mem = ARKStepCreate(...)`` to create the ARKStep memory
   block. :c:func:`ARKStepCreate()` returns a ``void*`` pointer to
   this memory structure. See the section
   :ref:`ARKStep_CInterface.Initialization` for details.

#. Configure the fast (inner) integrator

   Specify tolerances, create and attach matrix and/or solver objects,
   or call ``ARKStepSet*`` functions to configure the fast integrator
   as desired. See sections :ref:`ARKStep_CInterface.Skeleton` and
   :ref:`ARKStep_CInterface.OptionalInputs` for details on configuring
   ARKStep.

   *Notes on using ARKStep as a fast integrator:*

   It is the user's responsibility to create, configure, and attach the
   ``inner_arkode_mem`` to the MRIStep module.  User-specified options
   regarding how this fast integration should be performed (e.g., adaptive
   versus fixed time step, explicit/implicit/ImEx partitioning, algebraic
   solvers, etc.) will be respected during integration of the fast time scales
   during MRIStep integration.

   If a *user_data* pointer needs to be passed to user functions called by the
   fast (inner) integrator then it should be attached here by calling
   :c:func:`ARKStepSetUserData()`. This *user_data* pointer will only be passed
   to user-supplied functions that are attached to the fast (inner) integrator.
   To supply a *user_data* pointer to user-supplied functions called by the slow
   (outer) integrator the desired pointer should be attached by calling
   :c:func:`MRIStepSetUserData()` after creating the MRIStep memory below. Note
   the *user_data* pointers attached to the inner and outer integrators may
   be the same or different depending on what is required by the user code.

   Specifying a rootfinding problem for the fast integration is not
   supported. Rootfinding problems should be created and initialized with
   the slow integrator. See the steps below and :c:func:`MRIStepRootInit()`
   for more details.

   We note that due to the algorithms supported in MRIStep, the
   ARKStep module used for the fast time scale must be configured
   with an identity mass matrix.

#. Create an MRIStep object for the slow (outer) integration

   Call ``arkode_mem = MRIStepCreate(..., inner_arkode_mem)`` to create the 
   MRIStep memory block. :c:func:`MRIStepCreate()` returns a ``void*`` pointer 
   to this memory structure. See the section
   :ref:`MRIStep_CInterface.Initialization` for details.

#. Set the slow step size

   Call :c:func:`MRIStepSetFixedStep()` to specify the slow time step
   size.

   Specifically, if MRIStep is configured to use an implicit solver for the
   slow time scale, then the following steps are recommended:

#. Create and configure implicit solvers

   If MRIStep is configured to use an implicit solver for the slow time scale, then:
   
   #. Specify integration tolerances

      Call :c:func:`MRIStepSStolerances()` or :c:func:`MRIStepSVtolerances()` to
      specify either a scalar relative tolerance and scalar absolute tolerance,
      or a scalar relative tolerance and a vector of absolute tolerances,
      respectively.  Alternatively, call :c:func:`MRIStepWFtolerances()`
      to specify a function which sets directly the weights used in
      evaluating WRMS vector norms. See the section
      :ref:`MRIStep_CInterface.Tolerances` for details.

   #. Create nonlinear solver object

      If a non-default nonlinear solver object is desired for implicit
      MRI stage solves (see the section :ref:`MRIStep_CInterface.NonlinearSolvers`),
      then that nonlinear solver object must be created by using
      the appropriate functions defined by the particular SUNNONLINSOL
      implementation (e.g., ``NLS = SUNNonlinSol_***(...);`` where
      ``***`` is the name of the nonlinear solver (see the section
      :ref:`SUNNonlinSol` for details).

      For the SUNDIALS-supplied SUNNONLINSOL implementations, the
      nonlinear solver object may be created using a call of the form

      .. code-block:: c

         SUNNonlinearSolver NLS = SUNNonlinSol_*(...);

      where ``*`` can be replaced with "Newton", "FixedPoint", or other
      options, as discussed in the sections
      :ref:`ARKStep_CInterface.NonlinearSolvers` and :ref:`SUNNonlinSol`.

      Note: by default, MRIStep will use the Newton nonlinear solver
      (see section :ref:`SUNNonlinSol_Newton`), so a custom nonlinear solver
      object is only needed when using a *different* solver, or for the user
      to exercise additional controls over the Newton solver.
           
   #. Attach nonlinear solver module

      If a nonlinear solver object was created above, then it must be
      attached to MRIStep using the call (for details see the
      section :ref:`MRIStep_CInterface.NonlinearSolvers`):

      .. code-block:: c
  
         ier = MRIStepSetNonlinearSolver(...);

   #. Set nonlinear solver optional inputs

      Call the appropriate set functions for the selected nonlinear
      solver module to change optional inputs specific to that nonlinear
      solver.  These *must* be called after attaching the nonlinear
      solver to MRIStep, otherwise the optional inputs will be
      overridden by MRIStep defaults.  See the section
      :ref:`SUNNonlinSol` for more information on optional inputs.

   #. Create matrix object

      If a nonlinear solver requiring a linear solver will be used (e.g.,
      a Newton iteration) and if that linear solver will be matrix-based,
      then a template Jacobian matrix must be created by using the
      appropriate functions defined by the particular SUNMATRIX
      implementation.

      For the SUNDIALS-supplied SUNMATRIX implementations, the
      matrix object may be created using a call of the form

      .. code-block:: c

         SUNMatrix A = SUNBandMatrix(...);

      or

      .. code-block:: c
                      
         SUNMatrix A = SUNDenseMatrix(...);

      or

      .. code-block:: c

         SUNMatrix A = SUNSparseMatrix(...);

      or similarly for the CUDA and SuperLU_DIST matrix modules (see the
      sections :ref:`SUNMatrix_cuSparse` or :ref:`SUNMatrix_SLUNRloc` for
      further information).

      NOTE: The dense, banded, and sparse matrix objects are usable only in a
      serial or threaded environment.

   #. Create linear solver object

      If a nonlinear solver requiring a linear solver will be used (e.g.,
      a Newton iteration), then the desired linear solver object(s) must be
      created by using the appropriate functions defined by the particular
      SUNLINSOL implementation.

      For any of the SUNDIALS-supplied SUNLINSOL implementations, the
      linear solver object may be created using a call of the form

      .. code-block:: c
   
         SUNLinearSolver LS = SUNLinSol_*(...);

      where ``*`` can be replaced with "Dense", "SPGMR", or other
      options, as discussed in the sections
      :ref:`MRIStep_CInterface.LinearSolvers` and :ref:`SUNLinSol`.

   #. Set linear solver optional inputs

      Call ``*Set*`` functions from the selected linear solver module
      to change optional inputs specific to that linear solver.  See the
      documentation for each SUNLINSOL module in the section
      :ref:`SUNLinSol` for details.

   #. Attach linear solver module

      If a linear solver was created above for implicit MRI stage solves,
      initialize the ARKLS linear solver interface by attaching the
      linear solver object (and Jacobian matrix object, if applicable)
      with the call (for details see the section
      :ref:`MRIStep_CInterface.LinearSolvers`):

      .. code-block:: c
   
         ier = MRIStepSetLinearSolver(...);

#. Set optional inputs

   Call ``MRIStepSet*`` functions to change any optional inputs that
   control the behavior of MRIStep from their default values. See the
   section :ref:`MRIStep_CInterface.OptionalInputs` for details.

#. Specify rootfinding problem

   Optionally, call :c:func:`MRIStepRootInit()` to initialize a rootfinding
   problem to be solved during the integration of the ODE system. See
   the section :ref:`MRIStep_CInterface.RootFinding` for general details, and
   the section :ref:`MRIStep_CInterface.OptionalInputs` for relevant optional
   input calls.

#. Advance solution in time

   For each point at which output is desired, call

   .. code-block:: c

      ier = MRIStepEvolve(arkode_mem, tout, yout, &tret, itask);

   Here, ``itask`` specifies the return mode. The vector ``yout``
   (which can be the same as the vector ``y0`` above) will contain
   :math:`y(t_\text{out})`. See the section
   :ref:`MRIStep_CInterface.Integration` for details.

#. Get optional outputs

   Call ``MRIStepGet*`` and/or ``ARKStepGet*`` functions to obtain optional
   output from the slow or fast integrators respectively. See
   the section :ref:`MRIStep_CInterface.OptionalOutputs` and
   :ref:`ARKStep_CInterface.OptionalOutputs` for details.

#. Deallocate memory for solution vector

   Upon completion of the integration, deallocate memory for the
   vector ``y`` (or ``yout``) by calling the NVECTOR destructor
   function:

   .. code-block:: c

      N_VDestroy(y);

#. Free solver memory

    Call ``ARKStepFree(&inner_arkode_mem)`` and ``MRIStepFree(&arkode_mem)`` to
    free the memory allocated for fast and slow integration modules respectively.

#. Free linear solver and matrix memory

    Call :c:func:`SUNLinSolFree()` and (possibly)
    :c:func:`SUNMatDestroy()` to free any memory allocated for any
    linear solver and/or matrix objects created above for either the fast or
    slow integrators.

#. Free nonlinear solver memory

   If a user-supplied ``SUNNonlinearSolver`` was provided to MRIStep,
   then call :c:func:`SUNNonlinSolFree()` to free any memory allocated
   for the nonlinear solver object created above.

#. Finalize MPI, if used

    Call ``MPI_Finalize`` to terminate MPI.


    
SUNDIALS provides some linear solvers only as a means for users to get
problems running and not as highly efficient solvers.  For example, if
solving a dense system, we suggest using the LAPACK solvers if the
size of the linear system is :math:`> 50,000` (thanks to A. Nicolai
for his testing and recommendation).  See the table
:ref:`ARKStep_CInterface.solver-vector` for a listing of the 
linear solver interfaces available as ``SUNLinearSolver`` modules and
the vector implementations required for use.  
