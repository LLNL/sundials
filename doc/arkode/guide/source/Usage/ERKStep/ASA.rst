.. _ARKODE.Usage.ERKStep.ASA:

Adjoint Sensitivity Analysis
============================

The previous sections discuss using ARKStep for the integration of forward ODE models.
This section discusses how to use ARKStep for adjoint sensitivity analysis as introduced
in :numref:`ARKODE.Mathematics.ASA`. To use ARKStep for adjoint sensitivity analysis (ASA), users simply setup the forward
integration as usual (following :numref:`ARKODE.Usage.Skeleton`) with one exception:
a :c:type:`SUNAdjointCheckpointScheme` object must be created and passed to
:c:func:`ARKodeSetAdjointCheckpointScheme` before the call to the :c:func:`ARKodeEvolve`
function. After the forward model integration code, a :c:type:`SUNAdjointStepper` object
can be created for the adjoint model integration by calling :c:func:`ERKStepCreateAdjointStepper`.
The code snippet below demonstrates these steps in brief and the example code
``examples/arkode/C_serial/ark_lotka_volterra_asa.c`` demonstrates these steps in detail.

.. code-block:: C

   // 1. Create a SUNAdjointCheckpointScheme object

   // 2. Setup ERKStep for forward integration

   // 3. Attach the SUNAdjointCheckpointScheme

   // 4. Evolve the forward model

   // 5. Create the SUNAdjointStepper

   // 6. Setup the adjoint model

   // 7. Evolve the adjoint model

   // 8. Cleanup



User Callable Functions
-----------------------

This section describes ERKStep-specific user-callable functions for performing
adjoint sensitivity analysis with methods with ERKStep.

.. c:function:: int ERKStepCreateAdjointStepper(void* arkode_mem, N_Vector sf, SUNAdjointStepper* adj_stepper_ptr)

   Creates a :c:type:`SUNAdjointStepper` object compatible with the provided ARKStep instance for
   integrating the adjoint sensitivity system :eq:`ARKODE_DISCRETE_ADJOINT`.

   :param arkode_mem: a pointer to the ARKStep memory block.
   :param sf: the sensitivity vector holding the adjoint system terminal condition.
      This must be an instance of the ManyVector ``N_Vector`` implementation with at
      least one subvector (depending on if sensitivities to parameters should be computed).
      The first subvector must be :math:`\partial g_y(y(t_f)) \in \mathbb{R}^N`. If sensitivities to parameters should be computed, then the second subvector must be :math:`g_p(y(t_f), p) \in \mathbb{R}^{N_s}`.
   :param adj_stepper_ptr: the newly created :c:type:`SUNAdjointStepper` object.

   :return:
      * ``ARK_SUCCESS`` if successful
      * ``ARK_MEM_FAIL`` if a memory allocation failed
      * ``ARK_ILL_INPUT`` if an argument has an illegal value.
