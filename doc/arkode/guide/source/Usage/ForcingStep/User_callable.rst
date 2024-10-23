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

.. _ARKODE.Usage.ForcingStep.UserCallable:

ForcingStep User-callable functions
===================================

This section describes the ForcingStep-specific functions that may be called
by the user to setup and then solve an IVP using the ForcingStep time-stepping
module.

As discussed in the main :ref:`ARKODE user-callable function introduction
<ARKODE.Usage.UserCallable>`, each of ARKODE's time-stepping modules
clarifies the categories of user-callable functions that it supports.
ForcingStep does not support any of the categories beyond the functions that
apply for all time-stepping modules.


.. _ARKODE.Usage.ForcingStep.Initialization:

ForcingStep initialization functions
------------------------------------

.. c:function:: void* ForcingStepCreate(SUNStepper stepper1, SUNStepper stepper2, sunrealtype t0, N_Vector y0, SUNContext sunctx)

   This function allocates and initializes memory for a problem to be solved
   using the ForcingStep time-stepping module in ARKODE.

   :param stepper1: A :c:type:`SUNStepper` to integrate partition one. All
      :c:type:`SUNStepper` operations are required to be implemented except
      :c:func:`SUNStepper_SetForcing`.
   :param stepper2: A :c:type:`SUNStepper` to integrate partition two
      including the forcing from partition one. All :c:type:`SUNStepper`
      operations are required to be implemented.
   :param t0: The initial value of :math:`t`.
   :param y0: The initial condition vector :math:`y(t_0)`.
   :param sunctx: The :c:type:`SUNContext` object (see
      :numref:`SUNDIALS.SUNContext`)

   :return: If successful, a pointer to initialized problem memory of type
      ``void*``, to be passed to all user-facing ForcingStep routines listed
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

         /* create a ForcingStep object */
         arkode_mem = ForcingStepCreate(steppers[0], steppers[1], t0, y0, sunctx);

   **Example codes:**
      * ``examples/arkode/C_serial/ark_analytic_partitioned.c``
   
   .. versionadded:: x.y.z


.. _ARKODE.Usage.ForcingStep.OptionalOutputs:


Optional output functions
------------------------------

.. c:function:: int ForcingStep_GetNumEvolves(void* arkode_mem, int partition, long int *evolves)

   Returns the number of times the :c:type:`SUNStepper` for the given partition
   index has been evolved (so far).

   :param arkode_mem: pointer to the ForcingStep memory block.
   :param partition: index of the partition (0 or 1) or a negative number to
     indicate the total number across both partitions.
   :param evolves: number of :c:type:`SUNStepper` evolves.

   :retval ARK_SUCCESS: if successful
   :retval ARK_MEM_NULL: if the ForcingStep memory was ``NULL``
   :retval ARK_ILL_INPUT: if *partition* was out of bounds
   
   .. versionadded:: x.y.z