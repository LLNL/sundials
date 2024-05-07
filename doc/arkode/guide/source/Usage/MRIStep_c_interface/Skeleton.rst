.. ----------------------------------------------------------------
   Programmer(s): David J. Gardner @ LLNL
                  Daniel R. Reynolds @ SMU
   ----------------------------------------------------------------
   Based on ERKStep by Daniel R. Reynolds @ SMU
   ----------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2002-2024, Lawrence Livermore National Security
   and Southern Methodist University.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   ----------------------------------------------------------------

.. _ARKODE.Usage.MRIStep.Skeleton:

A skeleton of the user's main program
============================================

While MRIStep usage generally follows the same pattern as the rest of
ARKODE, since it involves the solution of both MRIStep for the slow
time scale and another time integrator for the fast time scale, we
summarize the differences in using MRIStep here.  Steps that are
unchanged from the skeleton program presented in
:numref:`ARKODE.Usage.Skeleton` are *italicized*.

.. index:: MRIStep user main program

#. *Initialize parallel or multi-threaded environment, if appropriate.*

#. *Create the SUNDIALS simulation context object*

#. *Set problem dimensions, etc.*

#. *Set vector of initial values*

#. Create an inner stepper object to solve the fast (inner) IVP

   * If using ARKStep as the fast (inner) integrator, create the ARKStep object
     with :c:func:`ARKStepCreate` and configure the integrator as desired for
     evolving the fast time scale. See sections :numref:`ARKODE.Usage.Skeleton`,
     :numref:`ARKODE.Usage.OptionalInputs`, and
     :numref:`ARKODE.Usage.ARKStep.OptionalInputs` for details on configuring
     ARKStep.

     Once the ARKStep object is setup, create an ``MRIStepInnerStepper`` object
     with :c:func:`ARKStepCreateMRIStepInnerStepper`.

   * If supplying a user-defined fast (inner) integrator, create the
     ``MRIStepInnerStepper`` object as described in section
     :numref:`ARKODE.Usage.MRIStep.CustomInnerStepper`.

   .. note::

      When using ARKStep as a fast (inner) integrator it is the user's
      responsibility to create, configure, and attach the integrator to the
      MRIStep module. User-specified options regarding how this fast integration
      should be performed (e.g., adaptive vs. fixed time step,
      explicit/implicit/ImEx partitioning, algebraic solvers, etc.) will be
      respected during evolution of the fast time scale during MRIStep
      integration.

      Due to the algorithms supported in MRIStep, the ARKStep module used for
      the fast time scale must be configured with an identity mass matrix.

      If a *user_data* pointer needs to be passed to user functions called by
      the fast (inner) integrator then it should be attached here by calling
      :c:func:`ARKodeSetUserData()`. This *user_data* pointer will only be
      passed to user-supplied functions that are attached to the fast (inner)
      integrator. To supply a *user_data* pointer to user-supplied functions
      called by the slow (outer) integrator the desired pointer should be
      attached by calling :c:func:`ARKodeSetUserData()` after creating the
      MRIStep memory below. The *user_data* pointers attached to the inner and
      outer integrators may be the same or different depending on what is
      required by the user code.

      Specifying a rootfinding problem for the fast integration is not
      supported. Rootfinding problems should be created and initialized with
      the slow integrator. See the steps below and :c:func:`ARKodeRootInit()`
      for more details.

#. Create an MRIStep object for the slow (outer) integration

   Create the MRIStep object by calling  :c:func:`MRIStepCreate`. One of the
   inputs to :c:func:`MRIStepCreate` is the ``MRIStepInnerStepper`` object for
   solving the fast (inner) IVP created in the previous step.

#. Set the slow step size

   Call :c:func:`ARKodeSetFixedStep()` on the MRIStep object to specify the
   slow time step size.

#. Create and configure implicit solvers (*as appropriate*)

   Specifically, if MRIStep is configured with an implicit slow right-hand side
   function in the prior step, then the following steps are recommended:

   #. *Specify integration tolerances*

   #. *Create matrix object*

   #. *Create linear solver object*

   #. *Set linear solver optional inputs*

   #. *Attach linear solver module*

   #. *Create nonlinear solver object*

   #. *Attach nonlinear solver module*

   #. *Set nonlinear solver optional inputs*

#. *Set optional inputs*

#. *Specify rootfinding problem*

#. *Advance solution in time*

#. *Get optional outputs*

#. *Deallocate memory for solution vector*

#. Free solver memory

   * If ARKStep was used as the fast (inner) IVP integrator, call
     :c:func:`MRIStepInnerStepper_Free` and :c:func:`ARKodeFree` to free the
     memory allocated for the fast (inner) integrator.

   * If a user-defined fast (inner) integrator was supplied, free the integrator
     content and call :c:func:`MRIStepInnerStepper_Free` to free the
     ``MRIStepInnerStepper`` object.

   * Call :c:func:`ARKodeFree` to free the memory allocated for the MRIStep
     slow integration object.

#. *Free linear solver and matrix memory (as appropriate)*

#. *Free nonlinear solver memory (as appropriate)*

#. *Free the SUNContext object*

#. *Finalize MPI, if used*
