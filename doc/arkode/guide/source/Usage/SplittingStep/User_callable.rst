.. ----------------------------------------------------------------
   Programmer(s): Steven B. Roberts @ LLNL
   ----------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2002-2024, Lawrence Livermore National Security
   and Southern Methodist University.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   ----------------------------------------------------------------

.. _ARKODE.Usage.SplittingStep.UserCallable:

SplittingStep User-callable functions
=====================================

This section describes the SplittingStep-specific functions that may be called
by the user to setup and then solve an IVP using the SplittingStep time-stepping
module.

As discussed in the main :ref:`ARKODE user-callable function introduction
<ARKODE.Usage.UserCallable>`, each of ARKODE's time-stepping modules
clarifies the categories of user-callable functions that it supports.
SplittingStep does not support any of the categories beyond the functions that
apply for all time-stepping modules.


.. _ARKODE.Usage.SplittingStep.Initialization:

SplittingStep initialization functions
--------------------------------------

.. c:function:: void* SplittingStepCreate(SUNStepper* steppers, int partitions, sunrealtype t0, N_Vector y0, SUNContext sunctx)

   This function allocates and initializes memory for a problem to be solved
   using the SplittingStep time-stepping module in ARKODE.

   **Arguments:**
      * *steppers* -- ?
      * *partitions* -- ?
      * *t0* -- the initial value of :math:`t`.
      * *y0* -- the initial condition vector :math:`y(t_0)`.
      * *sunctx* -- the :c:type:`SUNContext` object (see :numref:`SUNDIALS.SUNContext`)

   **Return value:**
      If successful, a pointer to initialized problem memory of type ``void*``,
      to be passed to all user-facing SplittingStep routines listed below. If
      unsuccessful, a ``NULL`` pointer will be returned, and an error message
      will be printed to ``stderr``.

   **Example usage:**

      .. code-block:: C

         /* fast (inner) and slow (outer) ARKODE objects */
         void *inner_arkode_mem = NULL;
         void *outer_arkode_mem = NULL;

         /* MRIStepInnerStepper to wrap the inner (fast) ARKStep object */
         MRIStepInnerStepper stepper = NULL;

         /* create an ARKStep object, setting fast (inner) right-hand side
            functions and the initial condition */
         inner_arkode_mem = ARKStepCreate(ffe, ffi, t0, y0, sunctx);

         /* setup ARKStep */
         . . .

         /* create MRIStepInnerStepper wrapper for the ARKStep memory block */
         flag = ARKStepCreateMRIStepInnerStepper(inner_arkode_mem, &stepper);

         /* create an MRIStep object, setting the slow (outer) right-hand side
            functions and the initial condition */
         outer_arkode_mem = MRIStepCreate(fse, fsi, t0, y0, stepper, sunctx)

   **Example codes:**
      * ``examples/arkode/C_serial/ark_advection_diffusion_reaction_splitting.c``
      * ``examples/arkode/C_serial/ark_analytic_splitting.c``
   
   .. versionadded:: x.y.z


Optional inputs for IVP method selection
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. c:function:: int SplittingStep_SetCoefficients(void* arkode_mem, SplittingStepCoefficients coefficients)

   Specifies a customized set of coefficients for the operator splitting method.

   **Arguments:**

   * *arkode_mem* -- pointer to the SplittingStep memory block.

   * *coefficients* -- the table of coupling coefficients for the MRI method.

   **Return value:**

   * *ARK_SUCCESS* if successful

   * *ARK_MEM_NULL* if the SplittingStep memory is ``NULL``

   * *ARK_ILL_INPUT* if an argument has an illegal value

   **Notes:**

   For a description of the :c:type:`SplittingStepCoefficients` type and related
   functions for creating splitting coefficients see
   :numref:`ARKODE.Usage.SplittingStep.SplittingStepCoefficients`.

   **Warning:**

   This should not be used with :c:func:`ARKodeSetOrder`.
