.. ----------------------------------------------------------------
   Programmer(s): Daniel R. Reynolds @ SMU
   ----------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2002-2024, Lawrence Livermore National Security
   and Southern Methodist University.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   ----------------------------------------------------------------

.. _ARKODE.Constants:

================
ARKODE Constants
================

Below we list all input and output constants used by the main solver,
timestepper, and linear solver modules, together with a short
description of their meaning.  :numref:`ARKODE.Constants.in_constants`
contains the ARKODE input constants, and :numref:`ARKODE.Constants.out_constants`
contains the ARKODE output constants.

.. _ARKODE.Constants.in_constants:
.. table:: ARKODE input constants
   :widths: 38 52

   +-----------------------------------------------+------------------------------------------------------------+
   | **Shared input constants**                    |                                                            |
   +-----------------------------------------------+------------------------------------------------------------+
   | :index:`ARK_NORMAL`                           | Solver should return at a specified output time.           |
   +-----------------------------------------------+------------------------------------------------------------+
   | :index:`ARK_ONE_STEP`                         | Solver should return after each successful step.           |
   +-----------------------------------------------+------------------------------------------------------------+
   |                                               |                                                            |
   +-----------------------------------------------+------------------------------------------------------------+
   | **Full right-hand side evaluation constants** |                                                            |
   +-----------------------------------------------+------------------------------------------------------------+
   | :index:`ARK_FULLRHS_START`                    | Calling the full right-hand side function at the           |
   |                                               | start of the integration.                                  |
   +-----------------------------------------------+------------------------------------------------------------+
   | :index:`ARK_FULLRHS_END`                      | Calling the full right-hand side function at the end of    |
   |                                               | a step.                                                    |
   +-----------------------------------------------+------------------------------------------------------------+
   | :index:`ARK_FULLRHS_OTHER`                    | Calling the full right-hand side function at the some      |
   |                                               | other point e.g., for dense output.                        |
   +-----------------------------------------------+------------------------------------------------------------+
   |                                               |                                                            |
   +-----------------------------------------------+------------------------------------------------------------+
   | **Interpolation module input constants**      |                                                            |
   +-----------------------------------------------+------------------------------------------------------------+
   | :index:`ARK_INTERP_NONE`                      | Disables polynomial interpolation for dense output.        |
   +-----------------------------------------------+------------------------------------------------------------+
   | :index:`ARK_INTERP_HERMITE`                   | Specifies use of the Hermite polynomial interpolation      |
   |                                               | module (for non-stiff problems).                           |
   +-----------------------------------------------+------------------------------------------------------------+
   | :index:`ARK_INTERP_LAGRANGE`                  | Specifies use of the Lagrange polynomial interpolation     |
   |                                               | module (for stiff problems).                               |
   +-----------------------------------------------+------------------------------------------------------------+
   | :index:`ARK_INTERP_MAX_DEGREE`                | Maximum possible interpolating polynomial degree.          |
   +-----------------------------------------------+------------------------------------------------------------+
   |                                               |                                                            |
   +-----------------------------------------------+------------------------------------------------------------+
   | **Relaxtion module input constants**          |                                                            |
   +-----------------------------------------------+------------------------------------------------------------+
   | :index:`ARK_RELAX_BRENT`                      | Specifies Brent's method as the relaxation nonlinear       |
   |                                               | solver.                                                    |
   +-----------------------------------------------+------------------------------------------------------------+
   | :index:`ARK_RELAX_NEWTON`                     | Specifies Newton's method as the relaxation nonlinear      |
   |                                               | solver.                                                    |
   +-----------------------------------------------+------------------------------------------------------------+
   |                                               |                                                            |
   +-----------------------------------------------+------------------------------------------------------------+
   | **MRI method types**                          |                                                            |
   +-----------------------------------------------+------------------------------------------------------------+
   | :index:`MRISTEP_EXPLICIT`                     | Use an explicit (at the slow time scale) MRI method.       |
   +-----------------------------------------------+------------------------------------------------------------+
   | :index:`MRISTEP_IMPLICIT`                     | Use an implicit (at the slow time scale) MRI method.       |
   +-----------------------------------------------+------------------------------------------------------------+
   | :index:`MRISTEP_IMEX`                         | Use an ImEx (at the slow time scale) MRI method.           |
   +-----------------------------------------------+------------------------------------------------------------+



.. _ARKODE.Constants.out_constants:
.. table:: ARKODE output constants
   :widths: 25 5 60

   +-------------------------------------+------+------------------------------------------------------------+
   | **Shared output constants**                                                                             |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_SUCCESS`                | 0    | Successful function return.                                |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_TSTOP_RETURN`           | 1    | ARKODE succeeded by reaching the specified stopping point. |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_ROOT_RETURN`            | 2    | ARKODE succeeded and found one more more roots.            |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_WARNING`                | 99   | ARKODE succeeded but an unusual situation occurred.        |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_TOO_MUCH_WORK`          | -1   | The solver took ``mxstep`` internal steps but could not    |
   |                                     |      | reach ``tout``.                                            |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_TOO_MUCH_ACC`           | -2   | The solver could not satisfy the accuracy                  |
   |                                     |      | demanded by the user for some internal step.               |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_ERR_FAILURE`            | -3   | Error test failures occurred too many times during one     |
   |                                     |      | internal time step, or the minimum step size was reached.  |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_CONV_FAILURE`           | -4   | Convergence test failures occurred too many times during   |
   |                                     |      | one internal time step, or the minimum step size was       |
   |                                     |      | reached.                                                   |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_LINIT_FAIL`             | -5   | The linear solver's initialization function failed.        |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_LSETUP_FAIL`            | -6   | The linear solver's setup function failed in an            |
   |                                     |      | unrecoverable manner.                                      |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_LSOLVE_FAIL`            | -7   | The linear solver's solve function failed in an            |
   |                                     |      | unrecoverable manner.                                      |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_RHSFUNC_FAIL`           | -8   | The right-hand side function failed in an                  |
   |                                     |      | unrecoverable manner.                                      |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_FIRST_RHSFUNC_ERR`      | -9   | The right-hand side function failed at the first call.     |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_REPTD_RHSFUNC_ERR`      | -10  | The right-hand side function had repeated recoverable      |
   |                                     |      | errors.                                                    |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_UNREC_RHSFUNC_ERR`      | -11  | The right-hand side function had a recoverable error, but  |
   |                                     |      | no recovery is possible.                                   |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_RTFUNC_FAIL`            | -12  | The rootfinding function failed in an unrecoverable        |
   |                                     |      | manner.                                                    |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_LFREE_FAIL`             | -13  | The linear solver's memory deallocation function failed.   |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_MASSINIT_FAIL`          | -14  | The mass matrix linear solver's initialization function    |
   |                                     |      | failed.                                                    |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_MASSSETUP_FAIL`         | -15  | The mass matrix linear solver's setup function failed in   |
   |                                     |      | an unrecoverable manner.                                   |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_MASSSOLVE_FAIL`         | -16  | The mass matrix linear solver's solve function failed in   |
   |                                     |      | an unrecoverable manner.                                   |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_MASSFREE_FAIL`          | -17  | The mass matrix linear solver's memory deallocation        |
   |                                     |      | function failed.                                           |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_MASSMULT_FAIL`          | -18  | The mass matrix-vector product function failed.            |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_CONSTR_FAIL`            | -19  | The inequality constraint test failed repeatedly or        |
   |                                     |      | failed with the minimum step size.                         |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_MEM_FAIL`               | -20  | A memory allocation failed.                                |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_MEM_NULL`               | -21  | The ``arkode_mem`` argument was ``NULL``.                  |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_ILL_INPUT`              | -22  | One of the function inputs is illegal.                     |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_NO_MALLOC`              | -23  | The ARKODE memory block was not allocated by               |
   |                                     |      | a call to :c:func:`ARKStepCreate`,                         |
   |                                     |      | :c:func:`ERKStepCreate`, or :c:func:`MRIStepCreate`.       |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_BAD_K`                  | -24  | The derivative order :math:`k` is larger than allowed.     |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_BAD_T`                  | -25  | The time :math:`t` is outside the last step taken.         |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_BAD_DKY`                | -26  | The output derivative vector is ``NULL``.                  |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_TOO_CLOSE`              | -27  | The output and initial times are too close to each other.  |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_VECTOROP_ERR`           | -28  | An error occurred when calling an :c:type:`N_Vector`       |
   |                                     |      | routine.                                                   |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_NLS_INIT_FAIL`          | -29  | An error occurred when initializing a SUNNonlinSol module. |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_NLS_SETUP_FAIL`         | -30  | A non-recoverable error occurred when setting up a         |
   |                                     |      | SUNNonlinSol module.                                       |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_NLS_SETUP_RECVR`        | -31  | A recoverable error occurred when setting up a             |
   |                                     |      | SUNNonlinSol module.                                       |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_NLS_OP_ERR`             | -32  | An error occurred when calling a set/get routine in a      |
   |                                     |      | SUNNonlinSol module.                                       |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_INNERSTEP_ATTACH_ERR`   | -33  | An error occurred when attaching the inner stepper module. |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_INNERSTEP_FAIL`         | -34  | An error occurred in the inner stepper module.             |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_PREINNERFN_FAIL`        | -35  | An error occurred in the MRIStep pre inner integrator      |
   |                                     |      | function.                                                  |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_POSTINNERFN_FAIL`       | -36  | An error occurred in the MRIStep post inner integrator     |
   |                                     |      | function.                                                  |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_INTERP_FAIL`            | -40  | An error occurred in the ARKODE polynomial interpolation   |
   |                                     |      | module.                                                    |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_INVALID_TABLE`          | -41  | An invalid Butcher or MRI table was encountered.           |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_CONTEXT_ERR`            | -42  | An error occurred with the SUNDIALS context object         |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_RELAX_FAIL`             | -43  | An error occurred in computing the relaxation parameter    |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_RELAX_MEM_FAIL`         | -44  | The relaxation memory structure is ``NULL``                |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_RELAX_FUNC_FAIL`        | -45  | The relaxation function returned an unrecoverable error    |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_RELAX_JAC_FAIL`         | -46  | The relaxation Jacobian function returned an unrecoverable |
   |                                     |      | error                                                      |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_CONTROLLER_ERR`         | -47  | An error with a SUNAdaptController object was encountered. |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_STEPPER_UNSUPPORTED`    | -48  | An operation was not supported by the current              |
   |                                     |      | time-stepping module.                                      |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_UNRECOGNIZED_ERROR`     | -99  | An unknown error was encountered.                          |
   +-------------------------------------+------+------------------------------------------------------------+
   |                                                                                                         |
   +-------------------------------------+------+------------------------------------------------------------+
   | **ARKLS linear solver module output constants**                                                         |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARKLS_SUCCESS`              | 0    | Successful function return.                                |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARKLS_MEM_NULL`             | -1   | The ``arkode_mem`` argument was ``NULL``.                  |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARKLS_LMEM_NULL`            | -2   | The ARKLS linear solver interface has not been             |
   |                                     |      | initialized.                                               |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARKLS_ILL_INPUT`            | -3   | The ARKLS solver interface is not compatible with          |
   |                                     |      | the current :c:type:`N_Vector` module, or an input value   |
   |                                     |      | was illegal.                                               |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARKLS_MEM_FAIL`             | -4   | A memory allocation request failed.                        |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARKLS_PMEM_NULL`            | -5   | The preconditioner module has not been initialized.        |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARKLS_MASSMEM_NULL`         | -6   | The ARKLS mass-matrix linear solver interface has not been |
   |                                     |      | initialized.                                               |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARKLS_JACFUNC_UNRECVR`      | -7   | The Jacobian function failed in an unrecoverable manner.   |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARKLS_JACFUNC_RECVR`        | -8   | The Jacobian function had a recoverable error.             |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARKLS_MASSFUNC_UNRECVR`     | -9   | The mass matrix function failed in an unrecoverable        |
   |                                     |      | manner.                                                    |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARKLS_MASSFUNC_RECVR`       | -10  | The mass matrix function had a recoverable error.          |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARKLS_SUNMAT_FAIL`          | -11  | An error occurred with the current :c:type:`SUNMatrix`     |
   |                                     |      | module.                                                    |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARKLS_SUNLS_FAIL`           | -12  | An error occurred with the current                         |
   |                                     |      | :c:type:`SUNLinearSolver` module.                          |
   +-------------------------------------+------+------------------------------------------------------------+

.. c:enum:: ARKRelaxSolver

   Nonlinear solver identifiers used to specify the method for solving
   :eq:`ARKODE_RELAX_NLS` when relaxation is enabled.

   .. c:enumerator:: ARK_RELAX_NEWTON

      Newton's method

   .. c:enumerator:: ARK_RELAX_BRENT

      Brent's method
..
   Commented-out table rows:

      +-------------------------------------+------+------------------------------------------------------------+
      | :index:`ARK_POSTPROCESS_STEP_FAIL`  | -37  | An error occurred when calling the user-provided           |
      |                                     |      | step-based :c:func:`ARKPostProcessFn` routine.             |
      +-------------------------------------+------+------------------------------------------------------------+
      | :index:`ARK_POSTPROCESS_STAGE_FAIL` | -38  | An error occurred when calling the user-provided           |
      |                                     |      | stage-based :c:func:`ARKPostProcessFn` routine.            |
      +-------------------------------------+------+------------------------------------------------------------+
