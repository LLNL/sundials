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

.. _IDA.Examples.Serial:

Serial example problems
=======================

.. _IDA.Examples.Serial.Roberts:

A dense example: idaRoberts_dns
-------------------------------

This example, due to Robertson :cite:p:`Rob:66`, models a three-species
chemical kinetics system in DAE form. Differential equations determine
:math:`y_1` and :math:`y_2`, while an algebraic equation determines
:math:`y_3`:

.. math::
   :label: IDA.Examples.RobertsDAE

   y'_1 &= -0.04y_1+10^4y_2y_3, \\
   y'_2 &= 0.04y_1-10^4y_2y_3-3\cdot10^7y_2^2, \\
   0 &= y_1+y_2+y_3-1.

The initial values are :math:`(y_1,y_2,y_3)=(1,0,0)`, and the system is
integrated from :math:`t=0` to :math:`4\cdot10^{10}`. Rootfinding locates the
points where :math:`y_1=10^{-4}` and :math:`y_3=0.01`.

After creating serial vectors, the program initializes IDA with
:c:func:`IDACreate`, :c:func:`IDAInit`, and :c:func:`IDASVtolerances`.
:c:func:`IDARootInit` registers the two root functions. A dense matrix and
linear solver are attached with :c:func:`IDASetLinearSolver`, and
:c:func:`IDASetJacFn` registers the analytic Jacobian.

The output loop calls :c:func:`IDASolve` in ``IDA_NORMAL`` mode. A return of
``IDA_ROOT_RETURN`` triggers :c:func:`IDAGetRootInfo`; successful integration
advances the output time by a factor of ten. Finally, the program reports
integration, residual, Jacobian, nonlinear iteration, error-test, convergence,
and rootfinding statistics.

.. literalinclude:: ../../../../examples/ida/serial/idaRoberts_dns.out
   :language: none

.. _IDA.Examples.Serial.FoodWeb:

A banded example: idaFoodWeb_bnd
--------------------------------

This example models a multispecies food web :cite:p:`Bro:86` with
predator-prey interactions and diffusion on the unit square. For :math:`s=2p`
species, prey species satisfy differential equations and predator species
satisfy algebraic equations:

.. math::
   :label: IDA.Examples.FoodWebPDE

   \frac{\partial c^i}{\partial t} &= R_i(x,y,c)+d_i(c^i_{xx}+c^i_{yy}),
      &&i=1,\ldots,p, \\
   0 &= R_i(x,y,c)+d_i(c^i_{xx}+c^i_{yy}), &&i=p+1,\ldots,s,

where

.. math::

   R_i(x,y,c)=c^i\left(b_i+\sum_{j=1}^s a_{ij}c^j\right),

.. math::

   a_{ij}=\begin{cases}
   -1&i=j,\\ -0.5\cdot10^{-6}&i\le p,\ j>p,\\
   10^4&i>p,\ j\le p,\\0&\text{otherwise},
   \end{cases}

.. math::

   b_i=\begin{cases}
   1+\alpha xy+\beta\sin(4\pi x)\sin(4\pi y)&i\le p,\\
   -[1+\alpha xy+\beta\sin(4\pi x)\sin(4\pi y)]&i>p,
   \end{cases}
   \quad
   d_i=\begin{cases}1&i\le p,\\0.5&i>p.\end{cases}

Homogeneous Neumann boundary conditions are used on :math:`0\le x,y\le1`,
with :math:`0\le t\le1`, :math:`\alpha=50`, and :math:`\beta=1000`. Initial
prey profiles are :math:`10+i[16x(1-x)y(1-y)]^2`; predator values are
:math:`10^5`.

Central differences on a :math:`20\times20` mesh with :math:`p=1` produce a
system of size 800 and a banded Jacobian with half-bandwidth 40. The program
uses a band matrix and band linear solver. It identifies differential and
algebraic components with :c:func:`IDASetId`, then calls :c:func:`IDACalcIC`
with ``IDA_YA_YDP_INIT`` to correct the algebraic states and differential
derivatives before integration.

.. literalinclude:: ../../../../examples/ida/serial/idaFoodWeb_bnd.out
   :language: none

.. _IDA.Examples.Serial.Heat2D:

A Krylov example: idaHeat2D_kry
-------------------------------

This example solves the heat equation on the unit square,

.. math::
   :label: IDA.Examples.Heat2DPDE

   \frac{\partial u}{\partial t}=u_{xx}+u_{yy}\quad (x,y)\in\Omega,
   \qquad u=0\quad (x,y)\in\partial\Omega,

for :math:`0\le t\le10.24`, with :math:`u(x,y,0)=16x(1-x)y(1-y)`. Central
differences on a :math:`10\times10` mesh yield 100 equations; the discrete
boundary conditions are algebraic equations.

The program imposes nonnegativity constraints, creates a GMRES solver with
modified Gram--Schmidt orthogonalization and five restarts, and supplies a
diagonal left preconditioner through :c:func:`IDASetPreconditioner`. It solves
the problem again after reinitialization using classical Gram--Schmidt. The
residual function applies boundary conditions and central differences; the
preconditioner retains the inverse diagonal of
:math:`J=\partial F/\partial u+c_j\partial F/\partial u'`.

.. literalinclude:: ../../../../examples/ida/serial/idaHeat2D_kry.out
   :language: none

.. _IDA.Examples.Serial.Heat2DKLU:

A sparse direct example: idaHeat2D_klu
--------------------------------------

This example solves the same heat problem with the KLU sparse direct solver and
a compressed sparse column matrix. The ``jacHeat3`` function handles the
special three-point grid case; ``jacHeat`` constructs the general pattern for
grids of at least four points. Both fill the column pointers, row indices, and
matrix data for :math:`J=\partial F/\partial u+c_j\partial F/\partial u'`.
The related ``idaHeat2D_sps`` example demonstrates SuperLU_MT.

.. literalinclude:: ../../../../examples/ida/serial/idaHeat2D_klu.out
   :language: none
