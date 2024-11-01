
.. _SUNAdjoint:

############################
Adjoint Sensitivity Analysis
############################

This section presents the ``SUNAdjointStepper`` and ``SUNAdjointCheckpointScheme`` classes.
The ``SUNAdjointStepper`` represents a generic adjoint sensitivity analysis (ASA) procedure
to obtain the adjoint sensitivities of an IVP of the form

.. math::
   \dot{y}(t) = f(t, y, p), \qquad y(t_0) = y_0, \qquad y \in \mathbb{R}^N,
   :label: SUNADJOINT_IVP

where :math:`p` is some set of :math:`N_s` problem parameters.

.. note::
  The API itself does not implement ASA, but it provides a common
  interface for ASA capabilities implemented in the SUNDIALS packages. Right now it supports :ref:`the
  ASA capabilities in ARKODE <ARKODE.Mathematics.ASA>`, while the ASA capabilities in :ref:`CVODES
  <CVODES.Mathematics.ASA>` and :ref:`IDAS <IDAS.Mathematics.ASA>` must be used directly.

Suppose we have a functional :math:`g(y(t_f),p)` for which we would like to compute the gradients
:math:`\partial g(y(t_f),p)/\partial y(t_0)` and/or :math:`\partial g(y(t_f),p)/\partial p`.  This
most often arises in the form of an optimization problem such as

.. math::
   \min_{\xi} \bar{\Psi}(\xi) = g(y(t_f), p)
   :label: SUNADJOINT_OPTIMIZATION_PROBLEM

where :math:`\xi \subset \{y(t_0), p\}`. The adjoint method is one approach to obtaining the
gradients that is particularly efficient when there are relatively few functionals and a large
number of parameters. While :ref:`CVODES <CVODES.Mathematics.ASA>` and
:ref:`IDAS <IDAS.Mathematics.ASA>` *continuous* adjoint methods
(differentiate-then-discretize), ARKODE provides *discrete* adjoint methods
(discretize-then-differentiate). For the continuous approach, we derive and solve the adjoint IVP
backwards in time

.. math::
   \lambda'(t) &= -f_y^T(t, y, p) \lambda,\quad \lambda(t_F) = g_y^T(y(t_f), p), \\
   \mu'(t) &= -f_p^T(t, y, p) \mu,\quad \mu(t_F) = g_p^T(y(t_f), p), \quad t_f \geq t \geq t_0, \\
   :label: SUNADJOINT_CONTINUOUS_ADJOINT_IVP

where :math:`\lambda(t) \in \mathbb{R}^N`, :math:`\mu(t) \in \mathbb{R}^{N_s}`
:math:`f_y \equiv \partial f/\partial y \in \mathbb{R}^{N \times N}` is the Jacobian with respect to the dependent variable,
and :math:`f_p \equiv \partial f/\partial p \in \mathbb{R}^{N \times N_s}` is the Jacobian with respect to the parameters
(:math:`N` is the size of the original IVP, :math:`N_s` is the number of parameters).
When solved with a numerical time integration scheme, the solution to the continuous adjoint IVP
are numerical approximations of the continuous adjoint sensitivities

.. math::
   \lambda(t_0) \approx  g_y^T(y(t_0), p),\quad \mu(t_0) \approx g_p^T(y(t_0), p)
   :label: SUNADJOINT_CONTINUOUS_ADJOINT_SOLUTION

For the discrete adjoint approach, we first numerically discretize the original IVP :eq:`SUNADJOINT_IVP`
using either a time integration scheme :math:`\varphi` so that

.. math::
   y_0 = y(t_0),\quad y_n = \varphi(y_{n-k}, \cdots, y_{n-1}, p), \quad k = n, \cdots, 1.
   :label: SUNADJOINT_DISCRETE_IVP

For linear multistep methods :math:`k \geq 1` and for one step methods :math:`k = 1`.
Reformulating the optimization problem for the discrete case, we have

.. math::
   \min_{\xi} \Psi(\xi) = g(y_n, p)
   :label: SUNADJOINT_DISCRETE_OPTIMIZATION_PROBLEM

The gradients of :eq:`SUNADJOINT_DISCRETE_OPTIMIZATION_PROBLEM` can be computed using the transposed chain
rule backwards in time to obtain the discete adjoint variables :math:`\lambda_n, \lambda_{n-1}, \cdots, \lambda_0`
and :math:`\mu_n, \mu_{n-1}, \cdots, \mu_0`,

.. math::
   \lambda_n &= g_y^T(y_n, p), \quad \lambda_k = 0, \quad \mu_n = g_y^T(y_n, p), \quad \mu_k = 0, \quad k = n - 1, \cdots, 0, \\
   \lambda_{\ell} &= \lambda_{\ell} + \left(\frac{\partial \varphi}{\partial y_{\ell}}(y_0, \cdots, y_{k-1}, p)\right)^T \lambda_{k},
   \quad \mu_{\ell} = \mu_{\ell} + \left(\frac{\partial \varphi}{\partial p}(y_0, \cdots, y_{k-1}, p)\right)^T \lambda_{k}, \\
   \quad & \quad \ell = k - 1, \cdots, 0, \quad k = n, \cdots, 0.
   :label: SUNADJOINT_DISCRETE_ADJOINT

The solution of the discrete adjoint equations :eq:`SUNADJOINT_DISCRETE_ADJOINT` is the sensitivities of the discrete cost function
:eq:`SUNADJOINT_DISCRETE_OPTIMIZATION_PROBLEM` with respect to changes in the discretized IVP :eq:`SUNADJOINT_DISCRETE_IVP`.

.. math::
   \lambda_0 = g_y^T(y_0, p), \quad \mu_0 = g_p^T(y_0, p).
   :label: SUNADJOINT_DISCRETE_ADJOINT_SOLUTION
