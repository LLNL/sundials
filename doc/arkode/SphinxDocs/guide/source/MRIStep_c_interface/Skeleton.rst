..
   Programmer(s): David J. Gardner @ LLNL
   ----------------------------------------------------------------
   Based on ERKStep by Daniel R. Reynolds @ SMU
   ----------------------------------------------------------------
   LLNS/SMU Copyright Start
   Copyright (c) 2018, Southern Methodist University and
   Lawrence Livermore National Security

   This work was performed under the auspices of the U.S. Department
   of Energy by Southern Methodist University and Lawrence Livermore
   National Laboratory under Contract DE-AC52-07NA27344.
   Produced at Southern Methodist University and the Lawrence
   Livermore National Laboratory.

   All rights reserved.
   For details, see the LICENSE file.
   LLNS/SMU Copyright End
   ----------------------------------------------------------------

:tocdepth: 3


.. _MRIStep_CInterface.Skeleton:

A skeleton of the user's main program
============================================

The following is a skeleton of the user's main program (or calling
program) for the integration of an ODE IVP using the MRIStep module.
Most of the steps are independent of the NVECTOR implementation used.
For the steps that are not, refer to the section :ref:`NVectors` for
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

#. Create MRIStep object

   Call ``arkode_mem = MRIStepCreate(...)`` to create the MRIStep memory
   block. :c:func:`MRIStepCreate()` returns a ``void*`` pointer to
   this memory structure. See the section
   :ref:`MRIStep_CInterface.Initialization` for details.

#. Set the slow and fast step sizes

   Call :c:func:`MRIStepSetFixedStep()` to specify the slow and fast time step
   sizes. 
           
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

   Call ``MRIStepGet*`` functions to obtain optional output. See
   the section :ref:`MRIStep_CInterface.OptionalOutputs` for details.

#. Deallocate memory for solution vector

   Upon completion of the integration, deallocate memory for the
   vector ``y`` (or ``yout``) by calling the NVECTOR destructor
   function:

   .. code-block:: c

      N_VDestroy(y);

#. Free solver memory

    Call ``MRIStepFree(&arkode_mem)`` to free the memory allocated for
    the MRIStep module.

#. Finalize MPI, if used

    Call ``MPI_Finalize`` to terminate MPI.
