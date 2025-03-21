.. _ARKODE.Usage.ASA:

Adjoint Sensitivity Analysis
============================

The previous sections discuss using ARKODE for the integration of forward ODE models. This section
discusses how to use ARKODE for adjoint sensitivity analysis as introduced in
:numref:`ARKODE.Mathematics.ASA`. To use ARKStep for adjoint sensitivity analysis (ASA), users
simply setup the forward integration as usual (following :numref:`ARKODE.Usage.Skeleton`) with a few
differences. Below we provide an updated version of the ARKode usage in section
:numref:`ARKODE.Usage.Skeleton` where steps that are unchanged from are *italicized*. The example
code ``examples/arkode/C_serial/ark_lotka_volterra_asa.c`` demonstrates these steps in detail.

.. index:: Adjoint Sensitivity Analysis user main program

#. *Initialize parallel or multi-threaded environment, if appropriate.*

#. *Create the SUNDIALS simulation context object*

#. *Set problem dimensions, etc.*

#. *Set vector of initial values*

#. *Create a stepper object*

#. *Set optional inputs*

#. *Specify rootfinding problem*

#. Create a :c:type:`SUNAdjointCheckpointScheme` object

   Create the :c:type:`SUNAdjointCheckpointScheme` object by calling ``SUNAdjointCheckpointScheme_Create_*``.
   Available :c:type:`SUNAdjointCheckpointScheme` implementations are found in
   section :numref:`SUNAdjoint.CheckpointScheme`.
   
#. Attach the checkpoint scheme object to ARKODE 
   
   Call :c:func:`ARKodeSetAdjointCheckpointScheme()`.

#. *Advance solution in time*

#. *Get optional outputs*

#. Create the sensitivities vector with the terminal condition

   The sensitivities vector must be an instance of the :ref:`ManyVector N_Vector implementation <NVectors.ManyVector>`.
   You will have one subvector for the initial condition sensitivities and an additional subvector if you
   want sensitivities with respect to parameters.
   The vector should contain the terminal condition from the forward problem.

#. Create the :c:type:`SUNAdjointStepper` object

   Call :c:func:`ERKStepCreateAdjointStepper` or :c:func:`ARKStepCreateAdjointStepper`.

#. Set either the Jacobian-vector product, vector-Jacobian product, or Jacobian evaluation functions

   Users must supply one of:
   
   * :math:`(\partial f/\partial y)^*v`,
   * :math:`v^*(\partial f/\partial y)`,
   * :math:`(\partial f/\partial y)`,

   and, if sensitivities with respect to the parameters is desired, one of

   * :math:`(\partial f/\partial p)^*v`,
   * :math:`v^*(\partial f/\partial p)`,
   * :math:`(\partial f/\partial p)`.

   These user-supplied routines can be set with :c:func:`SUNAdjointStepper_SetJacTimesVecFn`, or
   :c:func:`SUNAdjointStepper_SetJacFn`.

#. Set optional ASA input

   Refer to :numref:`SUNAdjoint.Stepper` for options.

#. Advance the adjoint sensitivity analysis ODE

   Call :c:func:`SUNAdjointStepper_Evolve`.

#. Get optional ASA outputs

   Refer to :numref:`SUNAdjoint.Stepper` for options.

#. Deallocate memory for ASA objects

   Deallocate the sensitivities vector, :c:type:`SUNAdjointStepper`, 
   and :c:type:`SUNAdjointCheckpointScheme` objects.

#. *Deallocate memory for solution vector*

#. Free solver memory

   * If an ARKODE stepper module was used as a partition IVP integrator, call
     :c:func:`SUNStepper_Destroy` and :c:func:`ARKodeFree` to free the memory
     allocated for that integrator.

   * If a user-defined partition integrator was supplied, free the integrator
     content and call :c:func:`SUNStepper_Destroy` to free the :c:type:`SUNStepper`
     object.

   * Call :c:func:`ARKodeFree` to free the memory allocated for the
     SplittingStep integration object.

#. *Free the SUNContext object*

#. *Finalize MPI, if used*



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
      The first subvector must be :math:`\partial g_y(y(t_f)) \in \mathbb{R}^N`. If sensitivities to parameters should be computed, then the second subvector must be :math:`g_p(y(t_f), p) \in \mathbb{R}^{N_s}`,
      otherwise only one subvector should be provided.
   :param adj_stepper_ptr: the newly created :c:type:`SUNAdjointStepper` object.

   :return:
      * ``ARK_SUCCESS`` if successful
      * ``ARK_MEM_FAIL`` if a memory allocation failed
      * ``ARK_ILL_INPUT`` if an argument has an illegal value.

   .. versionadded:: x.y.z

.. c:function:: int ARKStepCreateAdjointStepper(void* arkode_mem, N_Vector sf, SUNAdjointStepper* adj_stepper_ptr)

   Creates a :c:type:`SUNAdjointStepper` object compatible with the provided ARKStep instance for
   integrating the adjoint sensitivity system :eq:`ARKODE_DISCRETE_ADJOINT`.

   :param arkode_mem: a pointer to the ARKStep memory block.
   :param sf: the sensitivity vector holding the adjoint system terminal condition.
      This must be an instance of the ManyVector ``N_Vector`` implementation with at
      least one subvector (depending on if sensitivities to parameters should be computed).
      The first subvector must be :math:`\partial g_y(y(t_f)) \in \mathbb{R}^N`. If sensitivities to parameters should be computed, then the second subvector must be :math:`g_p(y(t_f), p) \in \mathbb{R}^{N_s}`,
      otherwise only one subvector should be provided.
   :param adj_stepper_ptr: the newly created :c:type:`SUNAdjointStepper` object.

   :return:
      * ``ARK_SUCCESS`` if successful
      * ``ARK_MEM_FAIL`` if a memory allocation failed
      * ``ARK_ILL_INPUT`` if an argument has an illegal value.

   .. versionadded:: x.y.z

   .. note:: 

      Currently only explicit methods are supported for ASA.
      