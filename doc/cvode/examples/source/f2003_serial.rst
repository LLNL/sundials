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


.. _serial_f2003:

=======================================
Serial Fortran 2003 example problems
=======================================



.. _cv_advdiff_bnd:

cv_advdiff_bnd
===================================================

Description
------------

This simple example problem is a serial Fortran 2003 implementation of
previous serial C example ``cvAdvDiff_bnd``.


Problem output
---------------

.. literalinclude:: ../../../../examples/cvode/F2003_serial/cv_advdiff_bnd_f2003.out
   :language: text


Numerical method
----------------

The numerical method is identical to the previous implementation,
except that we now use Fortran 2003.



.. _cv_analytic_fp:

cv_analytic_fp
===================================================

Description
------------

This is a simple example problem with an analytical solution.
Our problem is then:

.. math::

   \frac{dy}{dt} = \lambda \cdot y + \frac{1}{1 + t^2} - \lambda \cdot \text{atan} (t)

for :math:`t` in the interval :math:`[0.0, 10.0]`, with initial
condition :math:`y = 0`.

The stiffness of the problem is directly proportional to the value
of :math:`\lambda`. The value of :math:`\lambda` should be negative
to result in a well-posed ODE; for values with magnitude larger than
100, the problem becomes quite stiff.


Problem output
---------------

.. literalinclude:: ../../../../examples/cvode/F2003_serial/cv_analytic_fp_f2003.out
   :language: text


Numerical method
----------------

The example routine solves this problem using a Adams-Moulton
methods.  Each stage is solved using the SUNDIALS Fixed-point iteration.
The example file contains a function to evaluate :math:`f(t,y)`.

We specify the relative and absolute tolerances, :math:`rtol = 10^{-6}`
and :math:`atol = 10^{-10}`, respectively.  Aside from these choices,
this problem uses only the default CVODE solver parameters.

Output is printed every 1.0 units of time (10 total).
Run statistics (optional outputs) are printed at the end.



.. _cv_analytic_sys_dns:

cv_analytic_sys_dns
===================================================

Description
------------

This is a simple example problem with an analytical solution.
Our problem is then:

.. math::

   \frac{dy}{dt} = A * y

where :math:`A = V * D * V^{-1}` and

.. math::

   A &=\begin{pmatrix} \frac{\lambda}{4} - \frac{23}{40} &
                       \frac{\lambda}{4} - \frac{3}{40} &
                       \frac{\lambda}{4} - \frac{13}{40} \\
                       \frac{\lambda}{4} + \frac{21}{40} &
                       \frac{\lambda}{4} + \frac{1}{40} &
                       \frac{\lambda}{4} - \frac{11}{40} \\
                       \frac{\lambda}{2} + \frac{1}{20} &
                       \frac{\lambda}{2} + \frac{1}{20} &
                       \frac{\lambda}{2} - \frac{1}{20} &
       \end{pmatrix} \\
   V &=\begin{pmatrix} 1 & -1 & 1 \\
                       -1 & 2 & 1 \\
                       0 & -1 & 2
       \end{pmatrix} \quad
   D = \begin{pmatrix} -\frac{1}{2} & 0 & 0 \\
                       0 & -\frac{1}{10} & 0 \\
                       0 & 0 & \lambda
       \end{pmatrix} \quad
   V^{-1} = \begin{pmatrix} \frac{5}{4} & \frac{1}{4} & -\frac{3}{4} \\
                            \frac{1}{2} & \frac{1}{2} & -\frac{1}{2} \\
                            \frac{1}{4} & \frac{1}{4} & \frac{1}{4}
            \end{pmatrix}

and :math:`\lambda` is a large negative number. The analytical solution
to this problem is:

.. math::

   y(t) = V * e^{D \cdot t} * V^{-1} * y_0

for :math:`t` in the interval :math:`[0.0, 0.05]`, with initial condition:
:math:`y(0) = \begin{bmatrix} 1 \\ 1 \\ 1 \end{bmatrix}`.

The stiffness of the problem is directly proportional to the value of :math:`\lambda`.
The value of :math:`\lambda` should be negative to result in a well-posed ODE;
for values with magnitude larger than 100 the problem becomes quite stiff.

In this example, we choose :math:`\lambda = -100`.

Output is printed every 1.0 units of time (10 total).
Run statistics (optional outputs) are printed at the end.



Problem output
---------------

.. literalinclude:: ../../../../examples/cvode/F2003_serial/cv_analytic_sys_dns_f2003.out
   :language: text


Numerical method
----------------

The example routine solves this problem using a Backwards Differentiation
Formula in fixed-leading coefficient form.  Each stage is solved using the
built-in modified Newton iteration.  Internally, Newton will use the
SUNLINSOL_DENSE linear solver.  The example file
contains a function to evaluate :math:`f(t,y)`.

We specify the relative and absolute tolerances, :math:`rtol=10^{-6}`
and :math:`atol=10^{-10}`, respectively.  Aside from these choices,
this problem uses only the default CVODE solver parameters.


.. _cv_analytic_sys_dns_jac:

cv_analytic_sys_dns_jac
===================================================

Description
------------

This problem is exactly the same as :ref:`cv_analytic_sys_dns` above except
that here, as the title suggests, we provide a user-supplied Jacobian function.


Problem output
---------------

.. literalinclude:: ../../../../examples/cvode/F2003_serial/cv_analytic_sys_dns_jac_f2003.out
   :language: text


Numerical method
----------------

As previously stated, this problem is the same as ``cv_analytic_sys_dns``
apart from the example file containing another function to evaluate our
Jacobian denoted :math:`J(t,y)`.



.. _cv_analytic_sys_klu:

cv_analytic_sys_klu
===================================================

Description
------------

This problem is exactly the same as :ref:`cv_analytic_sys_dns` above except
that here, as the title suggests, we use sparse matrices to operate the
KLU linear solver module.


Problem output
---------------

.. literalinclude:: ../../../../examples/cvode/F2003_serial/cv_analytic_sys_klu_f2003.out
   :language: text


Numerical method
----------------

Again, this problem is the same as ``cv_analytic_sys_dns`` apart from
using SUNLINSOL_KLU and consequentially, SUNMATRIX_SPARSE as well.


.. _cv_brusselator_dns:

cv_brusselator_dns
===================================================

Description
------------

This problem has 3 components, :math:`Y = \begin{bmatrix} u \\ v \\ w \end{bmatrix}`,
that depend on the independent variable :math:`t` satisfying the
equations,

.. math::

   \frac{du}{dt} &= a - (w + 1) u + v u^2, \\
   \frac{dv}{dt} &= w u - v u^2, \\
   \frac{dw}{dt} &= \frac{b - w}{\varepsilon} - w u.

We integrate over the interval :math:`0 \leq t \leq 10`, with the
initial conditions :math:`Y_0 = \begin{bmatrix} u(0) = 3.9 \\ v(0) = 1.1 \\ w(0) = 2.8 \end{bmatrix}`,
and parameters :math:`a = 1.2`, :math:`b = 2.5` and :math:`\varepsilon = 10^{-5}`.

Here, all three components exhibit a rapid transient change during
the first :math:`0.2` time units, followed by a slow and smooth evolution.
After each unit time interval, the solution is output to the screen.


Problem output
---------------

.. literalinclude:: ../../../../examples/cvode/F2003_serial/cv_brusselator_dns_f2003.out
   :language: text


Numerical method
----------------

Since this driver and utility functions are written in Fortran 2003,
this example demonstrates the use of the FCVODE interface for the
CVODE solver.  For time integration, the solver uses a BDF method.
The implicit systems are solved using the built-in modified Newton
iteration, with the SUNMATRIX_DENSE matrix module and accompanying
SUNLINSOL_DENSE linear solver module.
Both the Jacobian routine and right-hand side functions are supplied
by functions provided in the example file.

The only non-default solver options are the tolerances
:math:`atol=10^{-10}` and :math:`rtol=10^{-6}`, adaptivity method 2 (I
controller), a maximum of 8 Newton iterations per step, a nonlinear
solver convergence coefficient :math:`nlscoef=10^{-8}`, and a maximum
of :math:`1000` internal time steps.



.. _cv_diurnal_kry_bp:

cv_diurnal_kry_bp
===================================================

This problem models a two-species diurnal
kinetics advection-diffusion PDE system in two spatial dimensions,

.. math::

   \frac{\partial c_i}{\partial t} =
     K_h \frac{\partial^2 c_i}{\partial x^2} +
     V \frac{\partial     c_i}{\partial x} +
     \frac{\partial}{\partial y}\left( K_v(y)
     \frac{\partial c_i}{\partial y}\right) +
     R_i(c_1,c_2,t),\quad i=1,2

where

.. math::

   R_1(c_1,c_2,t) &= -q_1*c_1*c_3 - q_2*c_1*c_2 + 2*q_3(t)*c_3 + q_4(t)*c_2, \\
   R_2(c_1,c_2,t) &=  q_1*c_1*c_3 - q_2*c_1*c_2 - q_4(t)*c_2, \\
   K_v(y) &= K_{v0} e^{y/5}.

Here :math:`K_h`, :math:`V`, :math:`K_{v0}`, :math:`q_1`, :math:`q_2`,
and :math:`c_3` are constants, and :math:`q_3(t)` and :math:`q_4(t)`
vary diurnally.  The problem is posed on the square spatial domain
:math:`(x,y) \in [0,20] \times [30,50]`, with homogeneous Neumann
boundary conditions, and for time interval :math:`t \in [0,86400]` sec
(1 day).

We enforce the initial conditions

.. math::

   c^1(x,y) &=  10^6 \chi(x) \eta(y) \\
   c^2(x,y) &=  10^{12} \chi(x) \eta(y) \\
   \chi(x) &= 1 - \sqrt{\frac{x - 10}{10}} + \frac{1}{2} \sqrt[4]{\frac{x - 10}{10}} \\
   \eta(y) &= 1 - \sqrt{\frac{y - 40}{10}} + \frac{1}{2} \sqrt[4]{\frac{x - 10}{10}}.



Problem output
---------------

.. literalinclude:: ../../../../examples/cvode/F2003_serial/cv_diurnal_kry_bp_f2003.out
   :language: text


Numerical method
----------------

We employ a method of lines approach, wherein we first semi-discretize
in space to convert the system of 2 PDEs into a larger system of ODEs.
To this end, the spatial derivatives are computed using second-order
centered differences, with the data distributed over :math:`Mx*My`
points on a uniform spatial grid.  As a result, CVODE approaches the
problem as one involving :math:`2*Mx*My` coupled ODEs. In this
problem, we use a relatively coarse uniform mesh with
:math:`Mx = My = 10`.

This program solves the problem with a BDF method, using a Newton
iteration with the preconditioned SUNLINSOL_SPGMR iterative linear
solver module.

The left preconditioner used is a banded matrix, constructed using
the CVBANDPRE module.  The banded preconditioner matrix is generated using
difference quotients, with half-bandwidths ``mu = ml = 2``.

Performance data and sampled solution values are printed at
selected output times, and all performance counters are printed
on completion.



.. _cv_diurnal_kry:

cv_diurnal_kry
===================================================

Description
------------

This problem is almost exactly similar to :ref:`cv_diurnal_kry_bp` above,
except here, we don't use the built in Banded Preconditioner, instead
supplying one. Interestingly, in doing so, we estimate a banded
preconditioner with ``mu = ml = 2`` as before using only a :math:`2 \times 2`-
block-diagonal matrix, treating each :math:`2 \times 2` block separately,
and finding the inverse of each block as they're added to the identity.
This method is employed to setup our user-supplied Preconditioner.


Problem output
---------------

.. literalinclude:: ../../../../examples/cvode/F2003_serial/cv_diurnal_kry_f2003.out
   :language: text


Numerical method
----------------

Instead of relying on CVBANDPRE as previously, we use CVSETPRECONDITIONER instead.
Interestingly, in doing so, we estimate a banded preconditioner with ``mu = ml = 2``
as before using only a :math:`2 \times 2`-block-diagonal matrix, treating each
:math:`2 \times 2` block separately, and finding the inverse of each block as they're
added to the identity.  This method is employed to setup our user-supplied Preconditioner.
For our Preconditioner's solver, we employ :math:`2 \times 2` matrix-vector multiplication.



.. _cv_roberts_dns:

cv_roberts_dns
===================================================

Description
------------

This problem is exactly the same problem as the C serial problem
``cvRoberts_dns`` aside from using FCVODE as a consequence of running
in Fortran 2003.  The details are significant enough to mention again.

This example is a simple problem, with the coding needed for its
solution completed via CVODE.  The problem is from chemical kinetics,
and consists of the following three rate equations:

.. math::

   \frac{dy_1}{dt} &= -0.04 \cdot y_1 + 10^4 \cdot y_2 \cdot y_3 \\
   \frac{dy_2}{dt} &= 0.04 \cdot y_1 - 10^4 \cdot y_2 \cdot y_3 - 3 \times 10^7 \cdot y_2^2 \\
   \frac{dy_3}{dt} &= 3 \times 10^7 \cdot y_2^2

on the interval from :math:`t = 0.0` to :math:`t = 4 \times 10^{10}`,
with initial conditions: :math:`y_1 = 1.0`, :math:`y_2 = y_3 = 0`.
The problem is stiff.

While integrating the system, we also use the rootfinding feature to
find the points at which :math:`y_1 = 10^{-4}` or at which :math:`y_3 = 0.01`.
This program solves the problem with the BDF method, Newton iteration
with the dense linear solver, and a user-supplied Jacobian routine.


Problem output
---------------

.. literalinclude:: ../../../../examples/cvode/F2003_serial/cv_roberts_dns_f2003.out
   :language: text


Numerical method
----------------

The example routine solves this problem using a Backwards Differentiation
Formula in fixed-leading coefficient form.  Each stage is solved using the
built-in modified Newton iteration.  Internally, Newton will use the
SUNLINSOL_DENSE linear solver.  The example file
contains functions to evaluate both :math:`f(t, y_1, y_2, y_3)` and
:math:`J(t, y_1, y_2, y_3)`.  Additionally, a root-finding function will,
as previously mentioned, find roots at :math:`y_3 = 0.01` and
:math:`y_1 = 10^{-4}` using the built-in root-finding mechanism in CVODE.

We specify the scalar relative and vector-valued absolute tolerances, :math:`rtol=10^{-4}`
and :math:`atol=\begin{Bmatrix} 10^{-8} \\ 10^{-14} \\ 10^{-6} \end{Bmatrix}`
, respectively.  Aside from these choices, this problem uses only the default
CVODE solver parameters.

11 normal + 2 root output times are printed at multiplicatively equally-spaced
points as well as at the two roots, and run statistics are printed at the end.


.. _cv_roberts_dns_constraints:

cv_roberts_dns_constraints
===================================================

Description
------------

This example problem is the same as :ref:`cv_roberts_dns` above except that the
constraint :math:`y_i \geq 0` is posed for all components :math:`i = 1, 2, 3`.


Problem output
---------------

.. literalinclude:: ../../../../examples/cvode/F2003_serial/cv_roberts_dns_constraints_f2003.out
   :language: text


Numerical method
----------------

Here, we specify the scalar relative and vector-valued absolute tolerances,
:math:`rtol=10^{-4}` and :math:`atol=\begin{Bmatrix} 10^{-6} \\ 10^{-11} \\ 10^{-5} \end{Bmatrix}`
, respectively.

Aside from this, see ``cv_roberts_dns`` above.


.. _cv_roberts_dnsL:

cv_roberts_dnsL
===================================================

This problem models
the kinetics of a three-species autocatalytic reaction.  This is an
ODE system with 3 components, :math:`Y = \begin{bmatrix} y_1 \\ y_2 \\ y_3 \end{bmatrix}`,
satisfying the equations,

.. math::

   \frac{d y_1}{dt} &= -0.04 y_1 + 10^4 y_2 y_3, \\
   \frac{d y_2}{dt} &= 0.04 y_1 - 10^4 y_2 y_3 - 3 \cdot 10^7 y_2^2, \\
   \frac{d y_3}{dt} &= 3 \cdot 10^7 y_2^2.

We integrate over the interval :math:`0 \leq t \leq 4 \cdot 10^{10}`, with initial
conditions  :math:`Y(0) = \begin{bmatrix} 1 \\ 0 \\ 0 \end{bmatrix}`.

Additionally, we supply the following two root-finding equations:

.. math::

   g_1(u) = u - 10^{-4}, \\
   g_2(w) = w - 10^{-2}.

While these are not inherently difficult nonlinear equations, they
easily serve the purpose of determining the times at which our
solutions attain desired target values.


Problem output
---------------

.. literalinclude:: ../../../../examples/cvode/F2003_serial/cv_roberts_dnsL_f2003.out
   :language: text


Numerical method
----------------

This program solves the problem with a BDF method, using a Newton
iteration with the SUNLINSOL_LAPACKDENSE linear solver module.

As with the :ref:`cv_roberts_dns` problem, we enable CVODE's
rootfinding module to find the times at which either :math:`u=10^{-4}`
or :math:`w=10^{-2}`.

Performance data and solution values are printed at
selected output times, along with additional output at root-finding
events.  All performance counters are printed on completion.





.. _cv_roberts_klu:

cv_roberts_klu
===================================================

Description
------------

This example problem is the same as :ref:`cv_roberts_dns` above except that
here we use the KLU sparse linear solver instead of the dense linear
solver.


Problem output
---------------

.. literalinclude:: ../../../../examples/cvode/F2003_serial/cv_roberts_klu_f2003.out
   :language: text


Numerical method
----------------

Aside from using SUNLINSOL_KLU, see ``cv_roberts_dns`` above.
