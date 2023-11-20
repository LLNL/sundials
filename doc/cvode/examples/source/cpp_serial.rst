..
   Programmer(s): Daniel M. Margolis @ SMU
   ----------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2002-2023, Lawrence Livermore National Security
   and Southern Methodist University.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   ----------------------------------------------------------------

:tocdepth: 3


.. _serial_cpp:

====================================
Serial C++ example problems
====================================



.. _cv_heat2D:

cv_heat2D
===============================================

Description
------------

This test example simulates a simple anisotropic 2D heat equation,

.. math::

   u_t = k_x u_{xx} + k_y u_{yy} + b,

for :math:`t` in :math:`[0, 1]` and :math:`(x, y)` in :math:`[0, 1]^2`,
with initial conditions

.. math::

   u(0, x, y) = \sin^2 (\pi x) \sin^2 (\pi y) + 1,

stationary boundary conditions

.. math::

   u_t (t, 0, y) = u_t (t, 1, y) = u_t (t, x, 0) = u_t (t, x, 1) = 0,

and the heat source

.. math::

   b(t, x, y) = -2 \pi \sin^2(\pi x) \sin^2(\pi y) \sin(\pi t) \cos(\pi t)
                -k_x 2 \pi^2 (\cos^2(\pi x) - \sin^2(\pi x)) \sin^2(\pi y) \cos^2(\pi t)
                -k_y 2 \pi^2 (\cos^2(\pi y) - \sin^2(\pi y)) \sin^2(\pi x) \cos^2(\pi t).

Under this setup, the problem has the analytical solution

.. math::

   u(t, x, y) =\sin^2(\pi x) \sin^2(\pi y) \cos^2(\pi t) + 1.

The spatial derivatives are computed using second-order centered differences,
with the data distributed over :math:`nx * ny` points on a uniform spatial grid. The
problem is advanced in time with BDF methods using an inexact Newton method
paired with the PCG or SPGMR linear solver. Several command line options are
available to change the problem parameters and CVODE settings. Use the flag
``--help`` for more information.


Problem output
---------------

.. literalinclude:: ../../../../examples/cvode/CXX_serial/cv_heat2D.out
   :language: text


Numerical method
----------------

The example routine solves this problem using a Backwards Differentiation
Formula in fixed-leading coefficient form.  Each stage is solved using the
built-in modified Newton iteration.  Internally, Newton will use the
SUNLINSOL_PCG linear solver or the SUNLINSOL_SPGMR linear solver.  The
example file contains functions to evaluate both
:math:`f(t,x,y)` and :math:`Pre\_Jac(t,x,y)`.

We specify the relative and absolute tolerances, :math:`rtol=0`
and :math:`atol=10^{-8}`, respectively.  Aside from these choices,
this problem uses only the default CVODE solver parameters.

12 outputs are printed, and run statistics are printed at the end as well.



.. _cv_kpr:

cv_kpr
================

Description
------------

This example problem is otherwise known as the Kvaerno-Prothero-Robinson
ODE test problem, structured in the following way

.. math::

   \begin{bmatrix} u' \\ v' \end{bmatrix} =
   \begin{bmatrix} a & b \\ c & d \end{bmatrix}
   \begin{bmatrix} \frac{-1 + u^2 - r(t)}{2u} \\
                   \frac{-2 + v^2 - s(t)}{2v} \end{bmatrix} +
   \begin{bmatrix} \frac{r'(t)}{2u} \\ \frac{s'(t)}{2v} \end{bmatrix}.

This problem has an analytical solution given by

.. math::

   u(t) = \sqrt{1 + r(t)}
   v(t) = \sqrt{2 + s(t)}

where, in this test, we use the functions

.. math::

   r(t) = 0.5 \cos(t)
   s(t) = \cos(2t)



Problem output
---------------

.. literalinclude:: ../../../../examples/cvode/CXX_serial/cv_kpr.out
   :language: text


Numerical method
-----------------

The example routine solves this problem using a Backwards Differentiation
Formula in fixed-leading coefficient form.  Each stage is solved using the
built-in modified Newton iteration.  Internally, Newton will use the
SUNLINSOL_DENSE linear solver.  The example file contains functions to
evaluate both :math:`f(t,u,v)` and :math:`J(t,u,v)`.

We specify the relative and absolute tolerances, :math:`rtol=10^{-6}`
and :math:`atol=10^{-10}`, respectively.  Aside from these choices,
this problem uses only the default CVODE solver parameters.

10 outputs are printed, and run statistics are printed at the end as well.
