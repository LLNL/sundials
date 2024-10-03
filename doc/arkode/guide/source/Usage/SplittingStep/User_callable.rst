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

   :param steppers: an array of :c:type:`SUNStepper` with one for each
      partition of the IVP.
   :param partitions: the number :math:`P > 1` of partitions in the IVP.
   :param t0: the initial value of :math:`t`.
   :param y0: the initial condition vector :math:`y(t_0)`.
   :param sunctx: the :c:type:`SUNContext` object (see :numref:`SUNDIALS.SUNContext`)

   :return: If successful, a pointer to initialized problem memory of type
      ``void*``, to be passed to all user-facing SplittingStep routines listed
      below. If unsuccessful, a ``NULL`` pointer will be returned, and an error
      message will be printed to ``stderr``.

   **Example usage:**

      .. code-block:: C

         /* inner ARKODE objects for integrating individual partitions */
         void *partition_mem[] = {NULL, NULL};

         /* SUNSteppers to wrap the inner ARKStep objects */
         SUNStepper steppers[] = {NULL, NULL};

         /* create ARKStep objects, setting right-hand side functions and the
            initial condition */
         partition_mem[0] = ARKStepCreate(fe1, fi1, t0, y0, sunctx);
         partition_mem[1] = ARKStepCreate(fe2, fi2, t0, y0, sunctx);

         /* setup ARKStep */
         . . .

         /* create SUNStepper wrappers for the ARKStep memory blocks */
         flag = ARKodeCreateSUNStepper(partition_mem[0], &stepper[0]);
         flag = ARKodeCreateSUNStepper(partition_mem[1], &stepper[1]);

         /* create a SplittingStep object with two partitions */
         arkode_mem = SplittingStepCreate(steppers, 2, t0, y0, sunctx);

   **Example codes:**
      * ``examples/arkode/C_serial/ark_advection_diffusion_reaction_splitting.c``
      * ``examples/arkode/C_serial/ark_analytic_partitioned.c``
   
   .. versionadded:: x.y.z


Optional inputs for IVP method selection
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. c:function:: int SplittingStep_SetCoefficients(void* arkode_mem, SplittingStepCoefficients coefficients)

   Specifies a customized set of coefficients for the operator splitting method.

   :param arkode_mem: pointer to the SplittingStep memory block.
   :param coefficients: the splitting coefficients for the method.

   :retval ARK_SUCCESS: if successful
   :retval ARK_MEM_NULL: if the SplittingStep memory is ``NULL``
   :retval ARK_ILL_INPUT: if an argument has an illegal value

   **Notes:**

   For a description of the :c:type:`SplittingStepCoefficients` type and related
   functions for creating splitting coefficients see
   :numref:`ARKODE.Usage.SplittingStep.SplittingStepCoefficients`.

   **Warning:**

   This should not be used with :c:func:`ARKodeSetOrder`.
   
   .. versionadded:: x.y.z


.. _ARKODE.Usage.SplittingStep.OptionalOutputs:


Optional output functions
------------------------------

.. c:function:: int SplittingStep_GetNumEvolves(void* arkode_mem, int partition, long int *evolves)

   Returns the number of times the :c:type:`SUNStepper` for the given partition
   index has been evolved (so far).

   :param arkode_mem: pointer to the SplittingStep memory block.
   :param partition: index of the partition between 0 and :math:`P - 1` or a
      negative number to indicate the total number across all
      partitions.
   :param evolves: number of :c:type:`SUNStepper` evolves.

   :retval ARK_SUCCESS: if successful
   :retval ARK_MEM_NULL: if the SplittingStep memory was ``NULL``
   :retval ARK_ILL_INPUT: if *partition* was out of bounds
   
   .. versionadded:: x.y.z
