..
   Copyright (c) 2002-2026, Lawrence Livermore National Security and
   the SUNDIALS contributors. SPDX-License-Identifier: BSD-3-Clause

.. _CVODES.Examples.ASA:

Adjoint sensitivity analysis examples
=====================================

A serial dense example: cvsRoberts_ASAi_dns
-------------------------------------------

For Robertson kinetics with parameters :math:`p`, this example computes the
gradient of

.. math::

   G(p)=\int_{t_0}^{T}y_3(t,p)\,dt.

The adjoint satisfies

.. math::

   \dot\lambda=-f_y^T\lambda-g_y^T,qquad \lambda(T)=0,

and backward quadratures satisfy

.. math::

   \dot\xi=g_p^T+f_p^T\lambda,qquad \xi(T)=0,qquad
   \frac{dG}{dp}=-\xi^T(t_0).

The forward phase uses BDF, a dense solver, a user Jacobian, and a quadrature
for :math:`G`. ``CVodeAdjInit`` enables checkpointing with Hermite
interpolation. After ``CVodeF``, the program creates a backward problem with
``CVodeCreateB`` and ``CVodeInitB``, attaches its dense solver and Jacobian,
initializes backward quadratures, and integrates with ``CVodeB``. It then
demonstrates reinitializing the backward problem for a second final time.

.. literalinclude:: ../../../../examples/cvodes/serial/cvsRoberts_ASAi_dns.out
   :language: none

An MPI nonstiff example: cvsAdvDiff_ASAp_non_p
----------------------------------------------

For :math:`u_t=p_1u_{xx}+p_2u_x` on :math:`[0,2]`, this example computes
gradients of :math:`g(t_f)=\int u(t_f,x)\,dx`. The adjoint PDE is

.. math::

   \mu_t+p_1\mu_{xx}-p_2\mu_x=0,qquad
   \mu(t_f,x)=1,qquad\mu(t,0)=\mu(t,2)=0.

Backward quadratures integrate :math:`\mu u_{xx}` and :math:`\mu u_x` to
obtain derivatives with respect to :math:`p_1` and :math:`p_2`.
:math:`\mu(t_0,x)` also gives the response to perturbations in the initial
profile. Forward and backward systems use Adams and fixed-point iteration.

.. figure:: figs/cvodes/cvsadjnonx_p.png
   :width: 75%
   :align: center

   Forward solution and initial-time adjoint sensitivity.

.. literalinclude:: ../../../../examples/cvodes/parallel/cvsAdvDiff_ASAp_non_p.out
   :language: none

An MPI CVBBDPRE example: cvsAtmDisp_ASAi_kry_bbd_p
--------------------------------------------------

This example models atmospheric transport in two or three dimensions:

.. math::

   c_t-k\nabla^2c+v\cdot\nabla c+S=0,

with homogeneous flux boundaries and zero initial concentration. It computes
the gradient, with respect to distributed-source parameters, of

.. math::

   G(p)=\frac12\int_0^T\int_\Omega \lVert c(t,x)\rVert^2\,d\Omega\,dt.

The forward and adjoint PDEs are centrally differenced with ghost cells and
solved using BDF, GMRES, and
:ref:`CVBBDPRE <CVODES.Usage.SIM.precond.cvbbdpre>`. Forward and backward
quadratures compute :math:`G` and its gradient. Processes exchange boundary
lines in 2D or surfaces in 3D.

.. figure:: figs/cvodes/cvsadjkryx_p2D.png
   :width: 90%
   :align: center

   Two-dimensional source gradient and nominal source distribution.

.. figure:: figs/cvodes/cvsadjkryx_p3Dcf.png
   :width: 70%
   :align: center

   Nominal three-dimensional source parameters.

.. figure:: figs/cvodes/cvsadjkryx_p3Dgrad.png
   :width: 70%
   :align: center

   Isosurfaces of the three-dimensional source gradient.

.. literalinclude:: ../../../../examples/cvodes/parallel/cvsAtmDisp_ASAi_kry_bbd_p.out
   :language: none
