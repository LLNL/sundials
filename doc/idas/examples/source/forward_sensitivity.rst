..
   -----------------------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2025-2026, Lawrence Livermore National Security,
   University of Maryland Baltimore County, and the SUNDIALS contributors.
   Copyright (c) 2013-2025, Lawrence Livermore National Security
   and Southern Methodist University.
   Copyright (c) 2002-2013, Lawrence Livermore National Security.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   -----------------------------------------------------------------------------

.. _IDAS.Examples.FSA:

Forward sensitivity analysis examples
=====================================

All IDAS examples support the ``IDA_SIMULTANEOUS`` and ``IDA_STAGGERED``
sensitivity methods. Sensitivities may be included in or excluded from the
error test through :c:func:`IDASetSensErrCon`. This section describes one
serial and one MPI example; the source files document the remaining examples.

.. _IDAS.Examples.FSA.SliderCrank:

A serial dense example: idasSlCrank_FSA_dns
-------------------------------------------

This example uses forward sensitivity analysis for a multibody dynamics
problem. The system consists of a crank and connecting rod with a
translational spring-damper (TSD) and a constant force acting on the rod. It
has one degree of freedom and is represented with three generalized
coordinates: the crank angle, horizontal position of the translational joint,
and connecting-rod angle. The crank has length :math:`a`, mass :math:`m_1`,
and moment of inertia :math:`J_1`; the connecting rod has length 2, mass
:math:`m_2`, and moment of inertia :math:`J_2`.

.. figure:: figs/idas/slider_crank.png
   :align: center
   :width: 90%

   Slider-crank mechanism modeled with three generalized coordinates.

The equations of motion are

.. math::

   M(y) \ddot y &= Q(y,\dot y) - \Phi_y^T(y) \lambda, \\
   \Phi(y) &= 0,

where :math:`y \in \mathbb{R}^3`, :math:`M` is the generalized mass matrix,
:math:`Q` is the generalized applied force, :math:`\Phi \in \mathbb{R}^2`
contains the position constraints, and :math:`\lambda \in \mathbb{R}^2`
contains their Lagrange multipliers.

For IDAS, this index-three DAE is converted to the stabilized index-two
Gear--Gupta--Leimkuhler formulation :cite:p:`GGL:85` by introducing two
additional multipliers :math:`\mu` and appending the velocity constraints:

.. math::
   :label: IDAS.Examples.GGL

   \dot y &= v - \Phi_y^T(y)\mu, \\
   M(y)\dot v &= Q(y,v) - \Phi_y^T(y)\lambda, \\
   \Phi(y) &= 0, \\
   \Phi_y(y)v &= 0.

The position constraints are

.. math::

   \Phi(y) = \begin{bmatrix}
      y_2-a\cos(y_1)-a\cos(y_3) \\
      a\sin(y_1)+\sin(y_3)
   \end{bmatrix},

and the generalized force is

.. math::

   Q(y,v) = \begin{bmatrix}
   -(f/\ell)a[\sin(y_3-y_1)/2+y_2\sin(y_1)]/2 \\
   (f/\ell)[\cos(y_3)/2-y_2+a\cos(y_1)/2]+F \\
   -(f/\ell)[y_2\sin(y_3)-a\sin(y_3-y_1)/2]/2-F\sin(y_3)
   \end{bmatrix},

where

.. math::

   f &= k(\ell-\ell_0)+c\ell', \\
   \ell^2 &= y_2^2-y_2[\cos(y_3)+a\cos(y_1)]+(1+a^2)/4
              +a\cos(y_3-y_1)/2, \\
   2\ell\ell' &= 2y_2v_2-v_2[\cos(y_3)+a\cos(y_1)]
      +y_2[\sin(y_3)v_3+a\sin(y_1)v_1]
      -a\sin(y_3-y_1)(v_3-v_1)/2.

The mass matrix is :math:`M=\operatorname{diag}(J_1,m_2,J_2)`. The example
uses :math:`a=0.5`, :math:`J_1=1`, :math:`J_2=2`, :math:`m_1=m_2=1`,
:math:`F=k=c=\ell_0=1`, and final time :math:`t_f=10`.

The state is :math:`Y=[y,v,\lambda,\mu]\in\mathbb{R}^{10}`. At :math:`t=0`,
consistent initial conditions are

.. math::

   y_1 &= \pi/2, & y_3 &= \arcsin(-a), & y_2 &= \cos(y_3), \\
   v_1 &= v_2=v_3=0, & \lambda_1 &= \lambda_2=\mu_1=\mu_2=0, \\
   \dot y_1 &= \dot y_2=\dot y_3=0, \\
   \dot v_1 &= Q_1(0)/J_1, & \dot v_2 &= Q_2(0)/m_2,
      & \dot v_3 &= Q_3(0)/J_2.

The remaining multiplier derivatives are zero. The relative and scalar
absolute tolerances are :math:`10^{-6}` and :math:`10^{-7}`. The algebraic
variables :math:`\lambda` and :math:`\mu` are excluded from the error test
using :c:func:`IDASetId` and :c:func:`IDASetSuppressAlg`.

Sensitivities with respect to the TSD parameters :math:`k` and :math:`c` are
used to estimate the gradient of the integrated kinetic energy

.. math::

   G=\int_{t_0}^{t_f}\left(\frac12J_1\dot y_1^2+
   \frac12m_2\dot y_2^2+\frac12J_2\dot y_3^2\right)\,dt.

The result is compared with backward, forward, and centered finite-difference
approximations. IDAS internally approximates the sensitivity residuals, and
quadrature sensitivities provide the gradient of :math:`G`.

.. figure:: figs/idas/x2sensi.png
   :align: center
   :width: 90%

   Sensitivities of :math:`y_2` with respect to the TSD parameters.

The following output uses simultaneous sensitivities and full error control:

.. literalinclude:: ../../../../examples/idas/serial/idasSlCrank_FSA_dns.out
   :language: none

.. _IDAS.Examples.FSA.Brusselator:

An MPI example using IDABBDPRE: idasBruss_FSA_kry_bbd_p
-------------------------------------------------------

This program solves the two-species time-dependent Brusselator PDE with GMRES
and the :ref:`IDABBDPRE <IDAS.Usage.precond.idabbdpre>` preconditioner:

.. math::

   \frac{\partial u}{\partial t} &= \epsilon_1(u_{xx}+u_{yy})
      +u^2v-(B+1)u+A, \\
   \frac{\partial v}{\partial t} &= \epsilon_2(v_{xx}+v_{yy})-u^2v+Bu.

The domain is the unit square, :math:`0\leq t\leq t_f=1`, with
:math:`\epsilon_1=\epsilon_2=0.002`, :math:`A=1`, and :math:`B=3.4`.
Homogeneous Neumann boundary conditions are used, with initial conditions

.. math::

   u=1-0.5\cos(\pi y), \qquad v=3.5-2.5\cos(\pi x).

The PDE is centrally differenced on a uniform 2D mesh. Each process owns a
submesh, and boundary conditions are implemented by copying the first interior
line into ghost values.

IDAS computes sensitivities with respect to :math:`\epsilon_i`. Spatially
integrating these gives sensitivities of the final spatial average

.. math::

   g=\iint u(x,y,t_f)\,dx\,dy, \qquad
   \frac{dg}{d\epsilon_i}=\iint
      \frac{\partial u(x,y,t_f)}{\partial\epsilon_i}\,dx\,dy.

A four-process run using simultaneous sensitivities and full error control is
``mpirun -np 4 idasBruss_FSA_kry_bbd_p -sensi sim t``. Its output is:

.. literalinclude:: ../../../../examples/idas/parallel/idasBruss_FSA_kry_bbd_p.out
   :language: none
