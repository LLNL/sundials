.. ----------------------------------------------------------------
   Programmer(s): Steven B. Roberts @ LLNL
   ----------------------------------------------------------------
   Based on MRIStep by David J. Gardner @ LLNL
   Daniel R. Reynolds @ SMU
   ----------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2002-2024, Lawrence Livermore National Security
   and Southern Methodist University.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   ----------------------------------------------------------------

.. _ARKODE.Usage.SplittingStep.Skeleton:

A skeleton of the user's main program
============================================

While SplittingStep usage generally follows the same pattern as the rest of
ARKODE, since it is the composition of other steppers, we summarize the
differences in using SplittingStep here.  Steps that are unchanged from the
skeleton program presented in :numref:`ARKODE.Usage.Skeleton` are *italicized*.

.. index:: SplittingStep user main program

#. *Initialize parallel or multi-threaded environment, if appropriate.*

#. *Create the SUNDIALS simulation context object*

#. *Set problem dimensions, etc.*

#. *Set vector of initial values*

#. Create a stepper object for each problem partition

   * If using ARKStep as an inner integrator, create the ARKStep object with
     :c:func:`ARKStepCreate` and configure the integrator as desired for
     evolving the partition. See sections :numref:`ARKODE.Usage.Skeleton`,
     :numref:`ARKODE.Usage.OptionalInputs`, and
     :numref:`ARKODE.Usage.ARKStep.OptionalInputs` for details on configuring
     ARKStep.

     Once the ARKStep object is setup, create a ``SUNStepper`` object with
     :c:func:`ARKStepCreateSUNStepper`.

   * If supplying a user-defined inner integrator, create the ``SUNStepper``
     object as described in section TODO(SBR).

   .. note::

      When using ARKStep as an inner integrator it is the user's responsibility
      to create and configure the integrator. User-specified options regarding
      how the integration should be performed (e.g., adaptive vs. fixed time
      step, explicit/implicit/ImEx partitioning, algebraic solvers, etc.) will
      be respected during evolution of a partition during SplittingStep
      integration.

      If a *user_data* pointer needs to be passed to user functions called by
      an inner integrator then it should be attached here by calling
      :c:func:`ARKodeSetUserData()`. This *user_data* pointer will only be
      passed to user-supplied functions that are attached to an inner
      integrator. To supply a *user_data* pointer to user-supplied functions
      called by the outer integrator the desired pointer should be attached by
      calling :c:func:`ARKodeSetUserData()` after creating the SplittingStep
      memory below. The *user_data* pointers attached to the inner and outer
      integrators may be the same or different depending on what is required by
      the user code.

      Specifying a rootfinding problem for an inner integrator is not supported.
      Rootfinding problems should be created and initialized with the outer
      integrator. See the steps below and :c:func:`ARKodeRootInit()` for more
      details.

#. Create a SplittingStep object for the outer integration

   Create the SplittingStep object by calling :c:func:`SplittingStepCreate`. One
   of the inputs to :c:func:`SplittingStepCreate` is an array of ``SUNStepper``
   objects with one to evolve each partition.

#. Set the SplittingStep step size

   Call :c:func:`ARKodeSetFixedStep()` on the SplittingStep object to specify
   the outer time step size.

#. *Set optional inputs*

#. *Specify rootfinding problem*

#. *Advance solution in time*

#. *Get optional outputs*

#. *Deallocate memory for solution vector*

#. Free solver memory

   * If ARKStep was used as an inner IVP integrator, call
     :c:func:`SUNStepper_Free` and :c:func:`ARKodeFree` to free the memory
     allocated for that inner integrator.

   * If a user-defined inner integrator was supplied, free the integrator
     content and call :c:func:`SUNStepper_Free` to free the ``SUNStepper``
     object.

   * Call :c:func:`ARKodeFree` to free the memory allocated for the
     SplittingStep outer integration object.

#. *Free the SUNContext object*

#. *Finalize MPI, if used*
