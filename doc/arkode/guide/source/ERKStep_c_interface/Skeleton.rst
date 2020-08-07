..
   Programmer(s): Daniel R. Reynolds @ SMU
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


.. _ERKStep_CInterface.Skeleton:

A skeleton of the user's main program
============================================

The following is a skeleton of the user's main program (or calling
program) for the integration of an ODE IVP using the ERKStep module.
Most of the steps are independent of the NVECTOR implementation used.
For the steps that are not, refer to the section :ref:`NVectors` for
the specific name of the function to be called or macro to be
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

4. Create ERKStep object

   Call ``arkode_mem = ERKStepCreate(...)`` to create the ERKStep memory
   block. :c:func:`ERKStepCreate()` returns a ``void*`` pointer to
   this memory structure. See the section
   :ref:`ERKStep_CInterface.Initialization` for details.

5. Specify integration tolerances

   Call :c:func:`ERKStepSStolerances()` or
   :c:func:`ERKStepSVtolerances()` to specify either a scalar relative
   tolerance and scalar absolute tolerance, or a scalar relative
   tolerance and a vector of absolute tolerances,
   respectively.  Alternatively, call :c:func:`ERKStepWFtolerances()`
   to specify a function which sets directly the weights used in
   evaluating WRMS vector norms. See the section
   :ref:`ERKStep_CInterface.Tolerances` for details.

6. Set optional inputs

   Call ``ERKStepSet*`` functions to change any optional inputs that
   control the behavior of ERKStep from their default values. See the
   section :ref:`ERKStep_CInterface.OptionalInputs` for details.

7. Specify rootfinding problem

   Optionally, call :c:func:`ERKStepRootInit()` to initialize a rootfinding
   problem to be solved during the integration of the ODE system. See
   the section :ref:`ERKStep_CInterface.RootFinding` for general details, and
   the section :ref:`ERKStep_CInterface.OptionalInputs` for relevant optional
   input calls.

8. Advance solution in time

   For each point at which output is desired, call

   .. code-block:: c

      ier = ERKStepEvolve(arkode_mem, tout, yout, &tret, itask);

   Here, ``itask`` specifies the return mode. The vector ``yout``
   (which can be the same as the vector ``y0`` above) will contain
   :math:`y(t_\text{out})`. See the section
   :ref:`ERKStep_CInterface.Integration` for details.

9. Get optional outputs

   Call ``ERKStepGet*`` functions to obtain optional output. See
   the section :ref:`ERKStep_CInterface.OptionalOutputs` for details.

10. Deallocate memory for solution vector

    Upon completion of the integration, deallocate memory for the
    vector ``y`` (or ``yout``) by calling the NVECTOR destructor
    function:

    .. code-block:: c

       N_VDestroy(y);

11. Free solver memory

    Call ``ERKStepFree(&arkode_mem)`` to free the memory allocated for
    the ERKStep module.

12. Finalize MPI, if used

    Call ``MPI_Finalize`` to terminate MPI.
