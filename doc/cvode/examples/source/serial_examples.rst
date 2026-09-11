..
   Copyright (c) 2002-2026, Lawrence Livermore National Security and
   the SUNDIALS contributors. SPDX-License-Identifier: BSD-3-Clause

.. _CVODE.Examples.Serial:

Serial example problems
=======================

.. _CVODE.Examples.Serial.Roberts:

A dense example: cvRoberts_dns
------------------------------

This Robertson kinetics example :cite:p:`Rob:66` solves

.. math::

   y'_1 &= -0.04y_1+10^4y_2y_3,\\
   y'_2 &= 0.04y_1-10^4y_2y_3-3\cdot10^7y_2^2,\\
   y'_3 &= 3\cdot10^7y_2^2,

with :math:`y(0)=(1,0,0)` from :math:`t=0` to :math:`4\cdot10^{10}`.
Rootfinding locates :math:`y_1=10^{-4}` and :math:`y_3=0.01`.

The program creates serial vectors and initializes BDF integration with
:c:func:`CVodeCreate`, :c:func:`CVodeInit`, and :c:func:`CVodeSVtolerances`.
:c:func:`CVodeRootInit` registers two root functions. A dense matrix and linear
solver are attached with :c:func:`CVodeSetLinearSolver`, and
:c:func:`CVodeSetJacFn` registers the analytic Jacobian. The output loop calls
:c:func:`CVode` in ``CV_NORMAL`` mode and handles ``CV_ROOT_RETURN`` with
:c:func:`CVodeGetRootInfo`.

.. literalinclude:: ../../../../examples/cvode/serial/cvRoberts_dns.out
   :language: none

.. _CVODE.Examples.Serial.AdvDiff:

A banded example: cvAdvDiff_bnd
-------------------------------

This program solves a two-dimensional advection-diffusion equation on
:math:`0\le x\le2`, :math:`0\le y\le1`:

.. math::

   u_t=0.5u_{xx}+0.75u_{yy}+u_x,

with homogeneous Dirichlet boundaries and initial condition
:math:`u(0,x,y)=x(2-x)y(1-y)e^{5xy}`. Central differencing on a
:math:`10\times5` interior mesh produces a 50-equation ODE system and a banded
Jacobian. The program uses BDF, Newton iteration, a band matrix and linear
solver, and a user-supplied Jacobian.

.. literalinclude:: ../../../../examples/cvode/serial/cvAdvDiff_bnd.out
   :language: none

.. _CVODE.Examples.Serial.Diurnal:

A Krylov example: cvDiurnal_kry
-------------------------------

This example models two reacting species transported by horizontal advection
and vertical diffusion in a 2D domain:

.. math::

   \frac{\partial c_i}{\partial t}=K_h\frac{\partial^2c_i}{\partial x^2}
   +V\frac{\partial c_i}{\partial x}
   +\frac{\partial}{\partial y}\left(K_v(y)\frac{\partial c_i}{\partial y}\right)
   +R_i(c_1,c_2,t),\quad i=1,2.

The diurnal reaction rates depend on sunlight through time-varying rate
coefficients. The PDE is centrally differenced on a rectangular mesh with
homogeneous Neumann conditions. CVODE uses BDF and GMRES with a left
preconditioner formed from the block-diagonal part of the Newton matrix. The
setup routine saves and conditionally reuses the block-diagonal Jacobian.

.. literalinclude:: ../../../../examples/cvode/serial/cvDiurnal_kry.out
   :language: none
