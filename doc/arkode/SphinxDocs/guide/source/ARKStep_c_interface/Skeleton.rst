..
   Programmer(s): Daniel R. Reynolds @ SMU
   ----------------------------------------------------------------
   Copyright (c) 2013, Southern Methodist University.
   All rights reserved.
   For details, see the LICENSE file.
   ----------------------------------------------------------------

:tocdepth: 3


.. _ARKStep_CInterface.Skeleton:

A skeleton of the user's main program
============================================

The following is a skeleton of the user's main program (or calling
program) for the integration of an ODE IVP using the ARKStep module.
Most of the steps are independent of the NVECTOR, SUNMATRIX, and
SUNLINSOL implementations used.  For the steps that are not, refer to
the sections :ref:`NVectors`, :ref:`SUNMatrix`  and :ref:`SUNLinSol`
for the specific name of the function to be called or macro to be
referenced.

.. index:: User main program

1. Initialize parallel or multi-threaded environment, if appropriate.

   For example, call ``MPI_Init`` to initialize MPI if used, or set
   ``num_threads``, the number of threads to use within the threaded
   vector functions, if used.

2. Set problem dimensions, etc.

   This generally includes the problem size, ``N``, and may include
   the local vector length ``Nlocal``.

   .. note::

      The variables ``N`` and ``Nlocal`` should be of type
      ``sunindextype``.

3. Set vector of initial values

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
   a call of the form

   .. code-block:: c

      y0 = N_VMake_***(..., c);

   where ``c`` is a pointer to a ``suncudavec`` or ``sunrajavec``
   vector class if this class already exists.  Otherwise, create a new
   vector by making a call of the form

   .. code-block:: c

      N_VGetDeviceArrayPointer_***

   or

   .. code-block:: c

      N_VGetHostArrayPointer_***

   Note that the vector class will allocate memory on both the host
   and device when instantiated.  See the sections
   :ref:`NVectors.CUDA` and :ref:`NVectors.RAJA` for details.

4. Create ARKStep object

   Call ``arkode_mem = ARKStepCreate(...)`` to create the ARKStep memory
   block. :c:func:`ARKStepCreate()` returns a ``void*`` pointer to
   this memory structure. See the section
   :ref:`ARKStep_CInterface.Initialization` for details.

5. Specify integration tolerances

   Call :c:func:`ARKStepSStolerances()` or
   :c:func:`ARKStepSVtolerances()` to specify either a scalar relative
   tolerance and scalar absolute tolerance, or a scalar relative
   tolerance and a vector of absolute tolerances,
   respectively.  Alternatively, call :c:func:`ARKStepWFtolerances()`
   to specify a function which sets directly the weights used in
   evaluating WRMS vector norms. See the section
   :ref:`ARKStep_CInterface.Tolerances` for details.

   If a problem with non-identity mass matrix is used, and the
   solution units differ considerably from the equation units,
   absolute tolerances for the equation residuals (nonlinear and
   linear) may be specified separately through calls to
   :c:func:`ARKStepResStolerance()`, :c:func:`ARKStepResVtolerance()`, or
   :c:func:`ARKStepResFtolerance()`.

6. Set optional inputs

   Call ``ARKStepSet*`` functions to change any optional inputs that
   control the behavior of ARKStep from their default values. See the
   section :ref:`ARKStep_CInterface.OptionalInputs` for details.

7. Create matrix object

   If a matrix-based linear solver is to be used within a Newton
   iteration or for solving non-identity mass matrix systems, then a
   template Jacobian and/or mass matrix must be created by using the
   appropriate functions defined by the particular SUNMATRIX
   implementation.

   NOTE: The dense, banded, and sparse matrix objects are usable only in a
   serial or threaded environment.

8. Create linear solver object

   If a Newton iteration is chosen, or if the problem involves a
   non-identity mass matrix, then the desired linear solver object(s)
   must be created by using the appropriate functions defined by the
   particular SUNLINSOL implementation.

9. Set linear solver optional inputs

   Call ``*Set*`` functions from the selected linear solver module
   to change optional inputs specific to that linear solver.  See the
   documentation for each SUNLINSOL module in the section
   :ref:`SUNLinSol` for details.

10. Attach linear solver module

    If a Newton iteration is chosen for implicit or ImEx methods,
    initialize the ARKDLS or ARKSPILS linear solver interface by
    attaching the linear solver object (and Jacobian matrix object, if
    applicable) with one of the following calls (for details see the
    section :ref:`ARKStep_CInterface.LinearSolvers`):

    .. code-block:: c

       ier = ARKDlsSetLinearSolver(...);

       ier = ARKSpilsSetLinearSolver(...);

    Similarly, if the problem involves a non-identity mass matrix,
    initialize the ARKDLS or ARKSPILS mass matrix linear solver
    interface by attaching the mass linear solver object (and mass
    matrix object, if applicable) with one of the following calls (for
    details see the section :ref:`ARKStep_CInterface.LinearSolvers`):

    .. code-block:: c

       ier = ARKDlsSetMassLinearSolver(...);

       ier = ARKSpilsSetMassLinearSolver(...);

11. Set linear solver interface optional inputs

    Call ``ARKDlsSet*`` or ``ARKSpilsSet*`` functions to change
    optional inputs specific to that linear solver interface. See the
    section :ref:`ARKStep_CInterface.OptionalInputs` for details.

12. Create nonlinear solver object

    If the problem involves an implicit component, and if a non-default
    nonlinear solver object will be used for implicit stage solves,
    then the desired nonlinear solver object must be created by using
    the appropriate functions defined by the particular SUNNONLINSOL 
    implementation. 

13. Attach nonlinear solver module

    If a nonlinear solver object was created above, then it must be
    attached to ARKStep using the call (for details see the
    section :ref:`ARKStep_CInterface.NonlinearSolvers`):

    .. code-block:: c

       ier = ARKStepSetNonlinearSolver(...);

14. Set nonlinear solver optional inputs

    Call the appropriate set functions for the selected nonlinear
    solver module to change optional inputs specific to that nonlinear
    solver.  These *must* be called after attaching the nonlinear
    solver to ARKStep, otherwise the optional inputs will be
    overridden by ARKStep defaults.  See the section
    :ref:`SUNNonlinSol` for more information on optional inputs.

12. Specify rootfinding problem

    Optionally, call :c:func:`ARKStepRootInit()` to initialize a rootfinding
    problem to be solved during the integration of the ODE system. See
    the section :ref:`ARKStep_CInterface.RootFinding` for general details, and
    the section :ref:`ARKStep_CInterface.OptionalInputs` for relevant optional
    input calls.

13. Advance solution in time

    For each point at which output is desired, call

    .. code-block:: c

       ier = ARKStepEvolve(arkode_mem, tout, yout, &tret, itask);

    Here, ``itask`` specifies the return mode. The vector ``yout``
    (which can be the same as the vector ``y0`` above) will contain
    :math:`y(t_\text{out})`. See the section
    :ref:`ARKStep_CInterface.Integration` for details.

15. Get optional outputs

    Call ``ARKStepGet*`` functions to obtain optional output. See
    the section :ref:`ARKStep_CInterface.OptionalOutputs` for details.

16. Deallocate memory for solution vector

    Upon completion of the integration, deallocate memory for the
    vector ``y`` (or ``yout``) by calling the destructor function:

    .. code-block:: c

       N_VDestroy(y);

17. Free solver memory

    Call ``ARKStepFree(&arkode_mem)`` to free the memory allocated for
    the ARKStep module (and any nonlinear solver module).

18. Free linear solver and matrix memory

    Call :c:func:`SUNLinSolFree()` and (possibly)
    :c:func:`SUNMatDestroy()` to free any memory allocated for the
    linear solver and matrix objects created above.

19. Finalize MPI, if used

    Call ``MPI_Finalize`` to terminate MPI.



SUNDIALS provides some linear solvers only as a means for users to get
problems running and not as highly efficient solvers.  For example, if
solving a dense system, we suggest using the LAPACK solvers if the
size of the linear system is :math:`> 50,000` (thanks to A. Nicolai
for his testing and recommendation).  The table below shows the
linear solver interfaces available as ``SUNLinearSolver`` modules and
the vector implementations required for use.  As an example, one
cannot use the dense direct solver interfaces with the MPI-based
vector implementation.  However, as discussed in section
:ref:`SUNLinSol` the SUNDIALS packages operate on generic
``SUNLinearSolver`` objects, allowing a user to develop their own
solvers should they so desire.



.. _ARKStep_CInterface.solver-vector:

SUNDIALS linear solver interfaces and vector implementations that can be used for each
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. cssclass:: table-bordered


+---------------+--------+----------+--------+----------+---------+--------+------+------+----------+
| Linear Solver | Serial | Parallel | OpenMP | pThreads | *hypre* | PETSc  | CUDA | RAJA | User     |
| Interface     |        | (MPI)    |        |          | Vec.    | Vec.   |      |      | Suppl.   |
+===============+========+==========+========+==========+=========+========+======+======+==========+
| Dense         | X      |          | X      | X        |         |        |      |      | X        |
+---------------+--------+----------+--------+----------+---------+--------+------+------+----------+
| Band          | X      |          | X      | X        |         |        |      |      | X        |
+---------------+--------+----------+--------+----------+---------+--------+------+------+----------+
| LapackDense   | X      |          | X      | X        |         |        |      |      | X        |
+---------------+--------+----------+--------+----------+---------+--------+------+------+----------+
| LapackBand    | X      |          | X      | X        |         |        |      |      | X        |
+---------------+--------+----------+--------+----------+---------+--------+------+------+----------+
| KLU           | X      |          | X      | X        |         |        |      |      | X        |
+---------------+--------+----------+--------+----------+---------+--------+------+------+----------+
| SuperLU_MT    | X      |          | X      | X        |         |        |      |      | X        |
+---------------+--------+----------+--------+----------+---------+--------+------+------+----------+
| SPGMR         | X      | X        | X      | X        | X       | X      | X    | X    | X        |
+---------------+--------+----------+--------+----------+---------+--------+------+------+----------+
| SPFGMR        | X      | X        | X      | X        | X       | X      | X    | X    | X        |
+---------------+--------+----------+--------+----------+---------+--------+------+------+----------+
| SPBCGS        | X      | X        | X      | X        | X       | X      | X    | X    | X        |
+---------------+--------+----------+--------+----------+---------+--------+------+------+----------+
| SPTFQMR       | X      | X        | X      | X        | X       | X      | X    | X    | X        |
+---------------+--------+----------+--------+----------+---------+--------+------+------+----------+
| PCG           | X      | X        | X      | X        | X       | X      | X    | X    | X        |
+---------------+--------+----------+--------+----------+---------+--------+------+------+----------+
| User supplied | X      | X        | X      | X        | X       | X      | X    | X    | X        |
+---------------+--------+----------+--------+----------+---------+--------+------+------+----------+
