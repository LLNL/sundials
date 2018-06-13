..
   Programmer(s): Daniel R. Reynolds @ SMU
   ----------------------------------------------------------------
   Copyright (c) 2013, Southern Methodist University.
   All rights reserved.
   For details, see the LICENSE file.
   ----------------------------------------------------------------

:tocdepth: 3


.. _CInterface.Skeleton:

A skeleton of the user's main program
============================================

The following is a skeleton of the user's main program (or calling
program) for the integration of an ODE IVP.  Most of the steps are
independent of the NVECTOR, SUNMATRIX, and SUNLINSOL implementations
used.  For the steps that are not, refer to the sections
:ref:`NVectors`, :ref:`SUNMatrix`  and :ref:`SUNLinSol` for the
specific name of the function to be called or macro to be referenced.

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
       
4. Create ARKode object

   Call ``arkode_mem = ARKodeCreate()`` to create the ARKode memory
   block. :c:func:`ARKodeCreate()` returns a pointer to the ARKode memory
   structure. See the section :ref:`CInterface.Initialization` for
   details.  

5. Create ARKode time-stepping module

   Call either :c:func:`ARKStepCreate()` or :c:func:`ERKStepCreate()`
   to provide required problem specifications, allocate internal
   memory for ARKode, and initialize ARKode.  These routines return a
   flag, the value of which indicates either success or an illegal
   argument value. See the section :ref:`CInterface.Initialization`
   for details.

6. Specify integration tolerances

   Call :c:func:`ARKodeSStolerances()` or
   :c:func:`ARKodeSVtolerances()` to specify either a scalar relative
   tolerance and scalar absolute tolerance, or a scalar relative
   tolerance and a vector of absolute tolerances,
   respectively.  Alternatively, call :c:func:`ARKodeWFtolerances()` 
   to specify a function which sets directly the weights used in
   evaluating WRMS vector norms. See the section
   :ref:`CInterface.Tolerances` for details. 

   If a problem with non-identity mass matrix is used, and the
   solution units differ considerably from the equation units,
   absolute tolerances for the equation residuals (nonlinear and
   linear) may be specified separately through calls to 
   :c:func:`ARKodeResStolerance()`, :c:func:`ARKodeResVtolerance()` or
   :c:func:`ARKodeResFtolerance()`.

7. Set optional inputs 

   Call ``ARKodeSet*``, ``ARKStepSet*`` or ``ERKStepSet*`` functions
   to change any optional inputs that control the behavior of ARKode
   from their default values. See the section
   :ref:`CInterface.OptionalInputs` for details.

8. Create matrix object
 
   If a direct linear solver is to be used within a Newton iteration
   or for solving non-identity mass matrix systems, then a template
   Jacobian and/or mass matrix must be created by using the
   appropriate functions defined by the particular SUNMATRIX
   implementation.
  
   NOTE: The dense, banded, and sparse matrix objects are usable only in a
   serial or threaded environment.

9. Create linear solver object
 
   If a Newton iteration is chosen, or if the problem involves a
   non-identity mass matrix, then the desired linear solver object(s)
   must be created by using the appropriate functions defined by the
   particular SUNLINSOL implementation.
  
10. Set linear solver optional inputs

    Call ``*Set*`` functions from the selected linear solver module
    to change optional inputs specific to that linear solver.  See the
    documentation for each SUNLINSOL module in the section
    :ref:`SUNLinSol` for details. 

11. Attach linear solver module

    If a Newton iteration is chosen for implicit or ImEx methods,
    initialize the ARKDLS or ARKSPILS linear solver interface by
    attaching the linear solver object (and Jacobian matrix object, if
    applicable) with one of the following calls (for details see the
    section :ref:`CInterface.LinearSolvers`): 

    .. code-block:: c
   
       ier = ARKDlsSetLinearSolver(...);

       ier = ARKSpilsSetLinearSolver(...);

    Similarly, if the problem involves a non-identity mass matrix,
    initialize the ARKDLS or ARKSPILS mass matrix linear solver
    interface by attaching the mass linear solver object (and mass
    matrix object, if applicable) with one of the following calls (for
    details see the section :ref:`CInterface.LinearSolvers`): 

    .. code-block:: c
   
       ier = ARKDlsSetMassLinearSolver(...);

       ier = ARKSpilsSetMassLinearSolver(...);

12. Set linear solver interface optional inputs 

    Call ``ARKDlsSet*`` or ``ARKSpilsSet*`` functions to change
    optional inputs specific to that linear solver interface. See the
    section :ref:`CInterface.OptionalInputs` for details. 

13. Specify rootfinding problem

    Optionally, call :c:func:`ARKodeRootInit()` to initialize a rootfinding
    problem to be solved during the integration of the ODE system. See
    the section :ref:`CInterface.RootFinding` for general details, and
    the section :ref:`CInterface.OptionalInputs` for relevant optional
    input calls. 

14. Advance solution in time

    For each point at which output is desired, call 

    .. code-block:: c
   
       ier = ARKode(arkode_mem, tout, yout, &tret, itask);

    Here, ``itask`` specifies the return mode. The vector ``yout``
    (which can be the same as the vector ``y0`` above) will contain
    :math:`y(t_\text{out})`. See the section
    :ref:`CInterface.Integration` for details. 

15. Get optional outputs 

    Call ``ARK*Get*`` functions to obtain optional output. See
    the section :ref:`CInterface.OptionalOutputs` for details.  

16. Deallocate memory for solution vector 

    Upon completion of the integration, deallocate memory for the
    vector ``y`` (or ``yout``) by calling the destructor function
    defined by the NVECTOR implementation:

    .. code-block:: c
   
       N_VDestroy_***(y);

17. Free solver memory 

    Call ``ARKodeFree(&arkode_mem)`` to free the memory allocated for ARKode. 

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



.. _CInterface.solver-vector:

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
