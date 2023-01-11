.. -----------------------------------------------------------------------------
   Programmer(s): David J. Gardner @ LLNL
   -----------------------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2002-2022, Lawrence Livermore National Security
   and Southern Methodist University.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   -----------------------------------------------------------------------------

.. _ARKODE.Usage.ARKStep.Relaxation:

Relaxation Methods
==================

This section describes user-callable functions for utilizing relaxation methods
with ARKStep. For more information on relaxation Runge-Kutta methods see
:numref:`ARKODE.Mathematics.Relaxation`.

Enabling or Disabling Relaxation
--------------------------------

.. c:function:: int ARKStepSetRelaxFn(void* arkode_mem, int nrfn, ARKRelaxFn rfn, ARKRelaxJacFn rjac)

   Attaches the user supplied functions for evaluating the relaxation function
   (``rfn``) and its Jacobian (``rjac``) and specifies the number of relaxation
   functions (``nrfn``).

   :param arkode_mem: the ARKStep memory structure
   :param nrfn: the number of relaxation functions
   :param rfn: the user-defined function to compute the relaxation function
               :math:`\xi_i(y)`
   :param rjac: the user-defined function to compute the relaxation Jacobian
                :math:`\xi'_i(y)`

   :retval ARK_SUCCESS: the function exited successfully
   :retval ARK_MEM_NULL: ``arkode_mem`` was ``NULL``
   :retval ARK_ILL_INPUT: an invalid input combination was provided see the
                          output error message for more details
   :retval ARK_MEM_FAIL: a memory allocation failed

   .. note::

      If ``nrfn = 0`` and both ``rfn = rjac = NULL`` relaxation is disabled.

Optional Input Functions
------------------------

This section describes optional input functions used to control the relaxation
method.

.. c:function:: int ARKStepSetRelaxEtaFail(void* arkode_mem, sunrealtype eta_rf)

   Sets the step size reduction factor applied after a failed relaxation solve.
   The default value is 0.25. Input values :math:`\geq 1` will result in the
   default value being used.

   :param arkode_mem: the ARKStep memory structure
   :param eta_rf: the step size reduction factor

   :retval ARK_SUCCESS: the value was successfully set
   :retval ARK_MEM_NULL: ``arkode_mem`` was ``NULL``
   :retval ARK_RELAX_MEM_NULL: the internal relaxation memory structure was
                               ``NULL``

.. c:function:: int ARKStepSetRelaxLowerBound(void* arkode_mem, sunrealtype lower)

   Sets the smallest acceptable value for the relaxation parameter. Values
   smaller than the lower bound will result in a relaxation solve failure and
   the step will be repeated with a smaller step size (determined by
   :c:func:`ARKStepSetRelaxEtaFail`). The default value is 0.8. Input values
   :math:`\geq 1` will result in the default value being used.

   :param arkode_mem: the ARKStep memory structure
   :param lower: the relaxation parameter lower bound

   :retval ARK_SUCCESS: the value was successfully set
   :retval ARK_MEM_NULL: ``arkode_mem`` was ``NULL``
   :retval ARK_RELAX_MEM_NULL: the internal relaxation memory structure was
                               ``NULL``

.. c:function:: int ARKStepSetRelaxMaxFails(void* arkode_mem, int max_fails)

   Sets the maximum number relaxation failures allowed in a single step attempt
   before the integration is halted with an error. The default value is 10.
   Input values :math:`\leq 0` will result in the default value being used.

   :param arkode_mem: the ARKStep memory structure
   :param max_iters: the maximum number relaxtion failures allowed in a step

   :retval ARK_SUCCESS: the value was successfully set
   :retval ARK_MEM_NULL: ``arkode_mem`` was ``NULL``
   :retval ARK_RELAX_MEM_NULL: the internal relaxation memory structure was
                               ``NULL``

.. c:function:: int ARKStepSetRelaxMaxIters(void* arkode_mem, int max_iters)

   Sets the maximum number of nonlinear iterations allowed when solving for the
   relaxation parameter. If the maximum number of iterations is reached before
   meeting the solve tolerance (determined by :c:func:`ARKStepSetRelaxTol`), the
   step will be repeated with a smaller step size (determined by
   :c:func:`ARKStepSetRelaxEtaFail`). The default value is 5. Input values
   :math:`\leq 0` will result in the default value being used.

   :param arkode_mem: the ARKStep memory structure
   :param max_iters: the maximum number of solver iterations allowed

   :retval ARK_SUCCESS: the value was successfully set
   :retval ARK_MEM_NULL: ``arkode_mem`` was ``NULL``
   :retval ARK_RELAX_MEM_NULL: the internal relaxation memory structure was
                               ``NULL``

.. c:function:: int ARKSteSetRelaxSolver(void* arkode_mem, ARKRelaxSolver solver)

   Sets the nonlinear solver method used to compute the relaxation parameter.
   The default value is ``ARK_RELAX_NEWTON``.

   :param arkode_mem: the ARKStep memory structure
   :param solver: the nonlinear solver to use

   :retval ARK_SUCCESS: the value was successfully set
   :retval ARK_MEM_NULL: ``arkode_mem`` was ``NULL``
   :retval ARK_RELAX_MEM_NULL: the internal relaxation memory structure was
                               ``NULL``

.. c:function:: int ARKStepSetRelaxTol(void* arkode_mem, sunrealtype tol)

   Sets the nonlinear solver tolerance to use when computing the relaxation
   parameter. If the tolerance is not reached within the maximum number of
   iterations (determined by :c:func:`ARKStepSetRelaxMaxIters`), the step will
   be repeated with a smaller step size (determined by
   :c:func:`ARKStepSetRelaxEtaFail`). The default value is 1.0e-14. Input values
   :math:`\leq 0.0` will result in the default value being used.

   :param arkode_mem: the ARKStep memory structure
   :param tol: the nonlinear solver tolerance to use

   :retval ARK_SUCCESS: the value was successfully set
   :retval ARK_MEM_NULL: ``arkode_mem`` was ``NULL``
   :retval ARK_RELAX_MEM_NULL: the internal relaxation memory structure was
                               ``NULL``

.. c:function:: int ARKStepSetRelaxUpperBound(void* arkode_mem, sunrealtype upper)

   Sets the largest acceptable value for the relaxation parameter. Values
   larger than the upper bound will result in a relaxation solve failure and
   the step will be repeated with a smaller step size (determined by
   :c:func:`ARKStepSetRelaxEtaFail`). The default value is 1.2. Input values
   :math:`\leq 1` will result in the default value being used.

   :param arkode_mem: the ARKStep memory structure
   :param eta_rf: the relaxation parameter upper bound

   :retval ARK_SUCCESS: the value was successfully set
   :retval ARK_MEM_NULL: ``arkode_mem`` was ``NULL``
   :retval ARK_RELAX_MEM_NULL: the internal relaxation memory structure was
                               ``NULL``

Optional Output Functions
-------------------------

This section describes optional output functions used to retrieve information
about the performance of the relaxation method.

.. c:function:: int ARKStepGetNumRelaxFnEvals(void* arkode_mem, long int* r_evals)

   Get the number of times the user's relaxation function was evaluated.

   :param arkode_mem: the ARKStep memory structure
   :param r_evals: the number of relaxation function evaluations

   :retval ARK_SUCCESS: the value was successfully set
   :retval ARK_MEM_NULL: ``arkode_mem`` was ``NULL``
   :retval ARK_RELAX_MEM_NULL: the internal relaxation memory structure was
                               ``NULL``

.. c:function:: int ARKStepGetNumRelaxJacEvals(void* arkode_mem, long int* J_evals)

   Get the number of times the user's relaxation Jacobian was evaluated.

   :param arkode_mem: the ARKStep memory structure
   :param J_evals: the number of relaxation Jacobian evaluations

   :retval ARK_SUCCESS: the value was successfully set
   :retval ARK_MEM_NULL: ``arkode_mem`` was ``NULL``
   :retval ARK_RELAX_MEM_NULL: the internal relaxation memory structure was
                               ``NULL``

.. c:function:: int ARKStepGetNumRelaxFails(void* arkode_mem, long int* fails)

   Get the total number of relaxation failures.

   :param arkode_mem: the ARKStep memory structure
   :param fails: the total number of relaxation failures

   :retval ARK_SUCCESS: the value was successfully set
   :retval ARK_MEM_NULL: ``arkode_mem`` was ``NULL``
   :retval ARK_RELAX_MEM_NULL: the internal relaxation memory structure was
                               ``NULL``

.. c:function:: int ARKStepGetNumRelaxSolveFails(void* arkode_mem, long int* fails)

   Get the number of times the relaxation parameter nonlinear solver failed.

   :param arkode_mem: the ARKStep memory structure
   :param fails: the number of relaxation nonlinear solver failures

   :retval ARK_SUCCESS: the value was successfully set
   :retval ARK_MEM_NULL: ``arkode_mem`` was ``NULL``
   :retval ARK_RELAX_MEM_NULL: the internal relaxation memory structure was
                               ``NULL``

.. c:function:: int ARKStepGetNumRelaxSolveIters(void* arkode_mem, long int* iters)

   Get the number of relaxation parameter nonlinear solver iterations.

   :param arkode_mem: the ARKStep memory structure
   :param iters: the number of relaxation nonlinear solver iterations

   :retval ARK_SUCCESS: the value was successfully set
   :retval ARK_MEM_NULL: ``arkode_mem`` was ``NULL``
   :retval ARK_RELAX_MEM_NULL: the internal relaxation memory structure was
