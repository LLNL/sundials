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

.. _IDAS.Examples.ASA:

Adjoint sensitivity analysis examples
=====================================

This section describes one serial and one MPI adjoint sensitivity example.
The source files document the remaining examples.

.. _IDAS.Examples.ASA.AkzoNobel:

A serial dense example: idasAkzoNob_ASAi_dns
--------------------------------------------

This program solves the index-one Akzo-Nobel chemical kinetics problem. The
model describes a process in which two species are mixed while carbon dioxide
is continuously added. It has the form

.. math::
   :label: IDAS.Examples.AkzoNobelDAE

   y' = f(y,z), \qquad 0=g(y,z),

with :math:`y\in\mathbb{R}^5` and :math:`z\in\mathbb{R}`. The differential
right-hand side is

.. math::

   f(y,z)=\begin{bmatrix}
   -2r_1+r_2-r_3-r_4 \\
   -\frac12r_1-r_4-\frac12r_5+F_{in} \\
   r_1-r_2+r_3 \\
   -r_2+r_3-2r_4 \\
   r_2-r_3+r_5
   \end{bmatrix},

where

.. math::

   r_1 &= k_1y_1^4y_2^{1/2}, & r_2 &= k_2y_3y_4, \\
   r_3 &= (k_2/K)y_1y_5, & r_4 &= k_3y_1y_4^{1/2}, \\
   r_5 &= k_4z^2y_2^{1/2}, &
   F_{in} &= klA\left(\frac{p(CO_2)}{H}-y_2\right).

The algebraic equation is

.. math::

   g(y,z)=K_sy_1y_4-z.

Since :math:`\partial g/\partial z` is nonsingular, the DAE has differentiation
index one. The dense linear solver is used with the default difference-quotient
Jacobian approximation.

The adjoint capability computes gradients, with respect to the initial values
of :math:`y`, of

.. math::

   G=\int_0^{t_f}y_1\,dt.

The initial value of :math:`z` is not a free parameter because it is determined
by :math:`y`. The first five components of the adjoint solution at :math:`t=0`
give the sensitivity of :math:`G` to the initial values of :math:`y`.

.. literalinclude:: ../../../../examples/idas/serial/idasAkzoNob_ASAi_dns.out
   :language: none

.. _IDAS.Examples.ASA.Brusselator:

An MPI example using IDABBDPRE: idasBruss_ASAp_kry_bbd_p
---------------------------------------------------------

This program solves the same Brusselator PDE as
:ref:`IDAS.Examples.FSA.Brusselator` and uses adjoint sensitivity analysis to
compute gradients of the model output

.. math::

   g(t)=\iint u(t,x,y)\,dx\,dy.

For perturbations :math:`\delta u_0` and :math:`\delta v_0` in the initial
profiles, the final output perturbation is

.. math::

   \delta g(t_f)=\iint[\lambda(0,x,y)\delta u_0+
      \mu(0,x,y)\delta v_0],dx\,dy,

where :math:`\lambda` and :math:`\mu` solve the adjoint PDEs

.. math::

   \frac{\partial\lambda}{\partial t} &=
      -\epsilon_1(\lambda_{xx}+\lambda_{yy})
      -(2uv-B-1)\lambda+(2uv-B)\mu, \\
   \frac{\partial\mu}{\partial t} &=
      -\epsilon_2(\mu_{xx}+\mu_{yy})-u^2\lambda+u^2\mu,

with homogeneous Neumann boundary conditions and final-time conditions

.. math::

   \lambda(t_f,x,y)=1, \qquad \mu(t_f,x,y)=0.

The adjoint PDEs are discretized and solved in the same manner as the forward
Brusselator equations.

.. literalinclude:: ../../../../examples/idas/parallel/idasBruss_ASAp_kry_bbd_p.out
   :language: none
