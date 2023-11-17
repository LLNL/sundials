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

.. _parallel_c:

====================================
Parallel C example problems
====================================



.. _cvAdvDiff_diag_p:

cvAdvDiff_diag_p
======================

Description
------------

This is a simple example problem, with the program for its solution 
provided by CVODE. The problem is the semi-discrete form
of the advection-diffusion equation in 1-D:

.. math::

   \frac{du}{dt} = \frac{d^2 u}{dx^2} + 0.5 \frac{du}{dx}

on the interval :math:`0 \leq x \leq 2`, and the time interval
:math:`0 \leq t \leq 5`.  Homogeneous Dirichlet boundary conditions
are posed, and the initial condition is the following:

.. math::

   u(x,t=0) = x \cdot (2 - x) \cdot e^{2x}.

The PDE is discretized on a uniform grid of size :math:`Mx+2` with
central differencing, and with boundary values eliminated,
leaving an ODE system of size :math:`NEQ = Mx`.

This program solves the problem with the ADAMS integration method,
and with Newton iteration using diagonal approximate Jacobians.
It uses scalar relative and absolute tolerances.

Output is printed at :math:`t = 0.5, 1.0, \ldots, 5`.  Run
statistics (optional outputs) are printed at the end.

This version uses MPI for user routines.

Execute with Number of Processors :math:`= N`,  with 
:math:`1 \leq N \leq Mx`.


Problem output
---------------

.. include:: ../../../../examples/cvode/parallel/cvAdvDiff_diag_p.out
   :literal:


Numerical method
-----------------

The numerical method is similar to the previous implementation,
except that we now use SUNDIALS' MPI-enabled vector kernel module,
NVECTOR_PARALLEL, and have similarly parallelized the supplied right-hand
side residual.

The example routine solves this problem using the Adams-Moulton module.
Each stage is solved using the built-in modified Newton iteration.
Internally, Newton will use the CVDIAG linear solver via the CVode interface.
The example file contains functions to evaluate :math:`f(t,u)`.

We specify the relative and absolute tolerances, :math:`reltol=0`
and :math:`abstol=10^{-5}`, respectively.  Aside from these choices,
this problem uses only the default MPI and CVode solver parameters.


.. _cvAdvDiff_non_p:

cvAdvDiff_non_p
=====================

Description
------------

This problem is exactly the same as ``cvAdvDiff_diag_p`` above except
that this program solves the problem with the option for nonstiff
systems: ADAMS method along with fixed-point iteration.


Problem output
---------------

.. include:: ../../../../examples/cvode/parallel/cvAdvDiff_non_p.out
   :literal:


Numerical method
-----------------

Similar to ``cvAdvDiff_diag_p`` above solve except that CVDIAG is
replaced by SUNNONLINSOL_FIXEDPOINT.


.. _cvDiurnal_kry_bbd_p:

cvDiurnal_kry_bbd_p
===================================================


As described in [HSR2017]_, this problem
models a two-species diurnal kinetics advection-diffusion PDE system
in two spatial dimensions,

.. math::

   \frac{\partial c_i}{\partial t} = 
     K_h \frac{\partial^2 c_i}{\partial x^2} + 
     V \frac{\partial     c_i}{\partial x} + 
     \frac{\partial}{\partial y}\left( K_v(y) 
     \frac{\partial c_i}{\partial y}\right) + 
     R_i(c_1,c_2,t),\qquad i=1,2 

where

.. math::

   R_1(c_1,c_2,t) &= -q_1*c_1*c_3 - q_2*c_1*c_2 + 2*q_3(t)*c_3 + q_4(t)*c_2, \\
   R_2(c_1,c_2,t) &=  q_1*c_1*c_3 - q_2*c_1*c_2 - q_4(t)*c_2, \\
   K_v(y) &= K_{v0} e^{y/5}.

Here :math:`K_h`, :math:`V`, :math:`K_{v0}`, :math:`q_1`, :math:`q_2`,
and :math:`c_3` are constants, and :math:`q_3(t)` and :math:`q_4(t)`
vary diurnally.  The problem is posed on the square spatial domain
:math:`(x,y) \in [0,20]\times[30,50]`, with homogeneous Neumann
boundary conditions, and for time interval :math:`t\in [0,86400]` sec
(1 day).

We enforce the initial conditions 

.. math::

   c_1(x,y) &=  10^6 \chi(x)\eta(y) \\
   c_2(x,y) &=  10^{12} \chi(x)\eta(y) \\
   \chi(x) &= 1 - \sqrt{\frac{x - 10}{10}} + \frac12 \sqrt[4]{\frac{x - 10}{10}} \\
   \eta(y) &= 1 - \sqrt{\frac{y - 40}{10}} + \frac12 \sqrt[4]{\frac{x - 10}{10}}.

The PDE system is treated by central differences on a uniform
mesh, with simple polynomial initial profiles.


Problem output
---------------

.. include:: ../../../../examples/cvode/parallel/cvDiurnal_kry_bbd_p.out
   :literal:


Numerical method
----------------

We employ a method of lines approach, wherein we first
semi-discretize in space to convert the system of 2 PDEs into a larger
system of ODEs.  To this end, the spatial derivatives are computed
using second-order centered differences, with the data distributed
over :math:`Mx*My` points on a uniform spatial grid.  As a result, CVode
approaches the problem as one involving :math:`2*Mx*My` coupled ODEs.

The problem is decomposed in parallel into uniformly-sized subdomains,
with two subdomains in each direction (four in total), and where each
subdomain has five points in each direction (i.e. :math:`Mx=My=10`).

This program solves the problem with a DIRK method, using a Newton
iteration with the preconditioned SUNLINSOL_SPGMR iterative linear
solver through the CVode interface.

The preconditioner matrix used is block-diagonal, with banded blocks,
constructed using the CVBBDPRE module.  Each block is generated using
difference quotients, with half-bandwidths ``mudq = mldq = 10``, but
the retained banded blocks have half-bandwidths ``mukeep = mlkeep = 2``.
A copy of the approximate Jacobian is saved and conditionally reused
within the preconditioner routine. 

Two runs are made for this problem, first with left and then with
right preconditioning.

Performance data and sampled solution values are printed at
selected output times, and all performance counters are printed
on completion.




.. _cvDiurnal_kry_p:

cvDiurnal_kry_p
===================================================

As described in [HSR2017]_, this test problem
models a two-species diurnal kinetics advection-diffusion PDE system
in two spatial dimensions,

.. math::

   \frac{\partial c_i}{\partial t} = 
     K_h \frac{\partial^2 c_i}{\partial x^2} + 
     V \frac{\partial     c_i}{\partial x} + 
     \frac{\partial}{\partial y}\left( K_v(y) 
     \frac{\partial c_i}{\partial y}\right) + 
     R_i(c_1,c_2,t),\qquad i=1,2 

where

.. math::

   R_1(c_1,c_2,t) &= -q_1*c_1*c_3 - q_2*c_1*c_2 + 2*q_3(t)*c_3 + q_4(t)*c_2, \\
   R_2(c_1,c_2,t) &=  q_1*c_1*c_3 - q_2*c_1*c_2 - q_4(t)*c_2, \\
   K_v(y) &= K_{v0} e^{y/5}.

Here :math:`K_h`, :math:`V`, :math:`K_{v0}`, :math:`q_1`, :math:`q_2`,
and :math:`c_3` are constants, and :math:`q_3(t)` and :math:`q_4(t)`
vary diurnally.  The problem is posed on the square spatial domain
:math:`(x,y) \in [0,20]\times[30,50]`, with homogeneous Neumann
boundary conditions, and for time interval :math:`t\in [0,86400]` sec
(1 day).

We enforce the initial conditions 

.. math::

   c^1(x,y) &=  10^6 \chi(x)\eta(y) \\
   c^2(x,y) &=  10^{12} \chi(x)\eta(y) \\
   \chi(x) &= 1 - \sqrt{\frac{x - 10}{10}} + \frac12 \sqrt[4]{\frac{x - 10}{10}} \\
   \eta(y) &= 1 - \sqrt{\frac{y - 40}{10}} + \frac12 \sqrt[4]{\frac{x - 10}{10}}.

The PDE system is treated by central differences on a uniform
mesh, with simple polynomial initial profiles.


Problem output
---------------

.. include:: ../../../../examples/cvode/parallel/cvDiurnal_kry_p.out
   :literal:


Numerical method
----------------

We employ a method of lines approach, wherein we first semi-discretize
in space to convert the system of 2 PDEs into a larger system of ODEs.
To this end, the spatial derivatives are computed using second-order
centered differences, with the data distributed over :math:`Mx*My`
points on a uniform spatial grid.  As a result, CVode approaches the
problem as one involving :math:`2*Mx*My` coupled ODEs. 

The problem is decomposed in parallel into uniformly-sized subdomains,
with two subdomains in each direction (four in total), and where each
subdomain has five points in each direction (i.e. :math:`Mx=My=10`).

This program solves the problem with a DIRK method, using a Newton
iteration with the preconditioned SUNLINSOL_SPGMR iterative linear
solver, through the CVode interface.

The preconditioner matrix used is block-diagonal, with block-diagonal
portion of the Newton matrix used as a left preconditioner.  A copy of
the block-diagonal portion of the Jacobian is saved and conditionally
reused within the preconditioner routine. 

Performance data and sampled solution values are printed at
selected output times, and all performance counters are printed
on completion.
