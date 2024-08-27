.. _ARKODE.Usage.ARKStep.ASA:

Adjoint Sensitivity Analysis
============================

This section describes ARKStep-specific user-callable functions for performing
adjoint sensitivity analysis with methods with ARKStep.


.. c:function:: int ARKStepCreateAdjointStepper(void* arkode_mem, N_Vector sf, SUNAdjointStepper* adj_stepper_ptr)

   Creates a :c:type:`SUNAdjointStepper` object compatible with the provided ARKStep instance for 
   integrating the adjoint sensitivity system :eq:`ARKODE_ADJOINT_ODE`.

   :param arkode_mem: a pointer to the ARKStep memory block.
   :param sf: the sensitivity vector holding the adjoint system terminal condition.
      This must be an instance of the ManyVector ``N_Vector`` implementation with two subvectors. 
      The first subvector must hold :math:`\partial g/\partial y |_{t=t_f} \in \mathbb{R}^N` and the second subvector must
      hold  :math:`\partial g / \partial p |_{t=t_f} \in \mathbb{R}^d`.
   :param adj_stepper_ptr: the newly created :c:type:`SUNAdjointStepper` object.

   :return:
      * ``ARK_SUCCESS`` if successful
      * ``ARK_MEM_FAIL`` if a memory allocation failed
      * ``ARK_ILL_INPUT`` if an argument has an illegal value.
