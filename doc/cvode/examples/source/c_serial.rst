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


.. _serial_c:

====================================
Serial C example problems
====================================



.. _cvAdvDiff_bnd:

cvAdvDiff_bnd
====================================

Description
------------

This is a very simple C example showing how to use the CVode solver
interface with a banded Jacobian.

The problem is the semi-discrete form of the advection-diffusion
equation in 2-D:

.. math::

   \frac{du}{dt} = \frac{d^2 u}{dx^2} + 0.5 \frac{du}{dx} + \frac{d^2 u}{dy^2}

on the rectangle :math:`0 \leq x \leq 2`, :math:`0 \leq y \leq 1`, and
the time interval :math:`0 \leq t \leq 1`.  Homogeneous Dirichlet
boundary conditions are posed, and the initial condition is

.. math::

   u(x, y, t=0) = x (2 - x) y (1 - y) e^{5xy}.

The PDE is discretized on a uniform :math:`mx + 2` by :math:`my + 2`
grid with central differencing, and with boundary values eliminated,
leaving an ODE system of size :math:`neq = mx \cdot my`.

This example solves the problem with the BDF method, Newton iteration
with the SUNBAND linear solver, and a user-supplied Jacobian routine.

It uses scalar relative and absolute tolerances.  Output is printed at
:math:`t = 0.1, 0.2, \ldots, 1.0`.  Run statistics (optional outputs)
are printed at the end.

Following the initial comment block, this program has a number
of ``#include`` lines, which allow access to useful items in CVODE
header files.  This is true for all examples.

Problem output
---------------

.. include:: ../../../../examples/cvode/serial/cvAdvDiff_bnd.out
   :literal:

Numerical method
----------------

The example routine solves this problem using a Backwards Differentiation
Formula in fixed-leading coefficient form.  Each stage is solved using the 
built-in modified Newton iteration.  Internally, Newton will use the
SUNLINSOL_DENSE linear solver via the CVode interface, which in the case
of this scalar-valued problem is just division.  The example file contains
functions to evaluate both :math:`f(t,u)` and :math:`J(t,u)`.  

We specify the relative and absolute tolerances, :math:`reltol=0`
and :math:`abstol=10^{-5}`, respectively.  Aside from these choices,
this problem uses only the default CVode solver parameters.



.. _cvAdvDiff_bndL:

cvAdvDiff_bndL
====================================

Description
------------

This example problem is the same as the previous one aside from the
use of the LAPACK_BAND linear solver. 

Problem output
---------------

.. include:: ../../../../examples/cvode/serial/cvAdvDiff_bndL.out
   :literal:



.. _cvAnalytic_mels:

cvAnalytic_mels
==============================================

Description
------------

This example problem is a simple example problem with an analytical
solution, represented as follows,

.. math::

   \frac{dy}{dt} = \lambda \cdot y + \frac{1}{1 + t^2} - \lambda \cdot \text{atan} (t)

for :math:`t` in the interval :math:`[0.0, 10.0]` with initial
condition: :math:`y = 0`.

The stiffness of the problem is directly proportional to the value
of ':math:`\lambda`'.  The value of :math:`\lambda` should be negative
in order to result in a well-posed ODE; for values with magnitudes
larger than 100, the problem becomes quite stiff.

Problem output
---------------

.. include:: ../../../../examples/cvode/serial/cvAnalytic_mels.out
   :literal:

Numerical method
----------------

This program solves the problem with the BDF method, Newton iteration,
and a custom 'Matrix-embedded' SUNLinearSolver.  Output is printed
every 1.0 units of time (10 total).  Run statistics (optional outputs)
are printed at the end.



.. _cvDirectDemo_ls:

cvDirectDemo_ls
================================================

Description
------------

Demonstration program for CVODE -- direct linear solvers.  Two
separate problems are solved using both the CV_ADAMS and CV_BDF
linear multistep methods in combination with the SUNNONLINSOL_FIXEDPOINT
and SUNNONLINSOL_NEWTON nonlinear solver modules:

*  Problem 1: Van der Pol oscillator

   *  :math:`\ddot x - 3 \cdot (1 - x^2) \cdot \dot x + x = 0` where
      :math:`x(0) = 2` and :math:`\dot x(0) = 0`.
   
   *  This second-order ODE is converted to a first-order system by
      defining :math:`y_0 = x` and :math:`y_1 = \dot x`.
      
   *  The NEWTON iteration cases use the following types of Jacobian
      approximation: (1) dense, user-supplied, (2) dense, difference
      quotient approximation, (3) diagonal approximation

*  Problem 2: :math:`\dot y = A \cdot y`, where :math:`A` is a banded
   lower triangular matrix derived from a 2-D advection PDE.

   *  The NEWTON iteration cases use the following types of Jacobian
      approximation: (1) band, user-supplied, (2) band, difference
      quotient approximation, (3) diagonal approximation.
   
   *  For each problem, in the series of eight runs, CVodeInit is
      called only once, for the first run, whereas CVodeReInit is
      called for each of the remaining seven runs.

*  Notes: This program demonstrates the usage of the sequential
   macros NV_Ith_S, SM_ELEMENT_D, SM_COLUMN_B, and SM_COLUMN_ELEMENT_B.

   *  The NV_Ith_S macro is used to reference the components of an N_Vector.

      *  It works for any size :math:`N = NEQ`, but due to efficiency concerns
         it should only by used when the problem size is small.

      *  The Problem 1 right hand side and Jacobian functions :math:`f_1` and
         :math:`Jac_1` both use NV_Ith_S.

      *  The N_VGetArrayPointer function gives the user access to the
         memory used for the component storage of an N_Vector. In the
         sequential case, the user may assume that this is one contiguous
         array of reals. The N_VGetArrayPointer function gives a more
         efficient means (than the NV_Ith_S macro) to access the components
         of an N_Vector and should be used when the problem size is large.

      *  The Problem 2 right hand side function :math:`f_2` uses the
         N_VGetArrayPointer function.

   *  The SM_ELEMENT_D macro used in :math:`Jac_1` gives access to an element
      of a dense SUNMatrix. It should be used only when the problem size is
      small (the size of a Dense SUNMatrix is :math:`NEQ \times NEQ`) due to
      efficiency concerns.

   *  For larger problem sizes, the macro SM_COLUMN_D can be used in order
      to work directly with a column of a Dense SUNMatrix.

   *  The SM_COLUMN_B and SM_COLUMN_ELEMENT_B allow efficient columnwise access
      to the elements of a Banded SUNMatix. These macros are used in the 
      :math:`Jac_2` function.


Problem output
---------------

.. include:: ../../../../examples/cvode/serial/cvDirectDemo_ls.out
   :literal:


Numerical method
----------------

This program solves the two separate problems first with ADAMS and 
then with BDF method, using a Fixed-point iteration and a Newton
iteration for their nonlinear solvers, interchangeably.  The first problem
uses the SUNLINSOL_DENSE linear solver module via the CVode interface.  
Similarly, the second problem uses the SUNLINSOL_BAND linear solver
module via the CVode interface.  Additionally, this example occasionally
provides a routine to CVode to compute the dense and banded Jacobian,
respectively. 

The problems are run using scalar relative and absolute tolerances of
:math:`reltol=0` and :math:`abstol=10^{-6}`, respectively.

4 or 5 outputs are printed at intermittant or multiplicatively equal
intervals, and run statistics are printed at the end.



.. _cvDisc_dns:

cvDisc_dns
=============

Description
------------

This test problem is a sample of two simple 1D examples to illustrate
integrating over discontinuities. For instance,

*  Discontinuity in solution

   .. math::

      y' = -y \qquad ; \enspace y(0) = 1 \qquad ; \enspace t = [0, 1] \\
      y' = -y \qquad ; \enspace y(1) = 1 \qquad ; \enspace t = [1, 2]

*  Discontinuity in RHS (:math:`y'`)

   .. math::

      y' = -y \qquad ; \enspace y(0) = 1 \qquad ; \enspace t = [0, 1] \\
      z' = -5 \cdot z \quad ; \ z(1) = y(1) \quad ; \enspace t = [1, 2]

   *  This case is solved twice, first by explicitly treating the
      discontinuity point and secondly by letting the integrator
      deal with the discontinuity.


Problem output
---------------

.. include:: ../../../../examples/cvode/serial/cvDisc_dns.out
   :literal:


Numerical method
----------------

This program solves the two separate sample problems both with BDF methods,
using the built-in Newton iteration for their nonlinear solvers.  Both
use the SUNLINSOL_DENSE linear solver module via the CVode interface.
Additionally, this example only advances one step in time using CV_ONE_STEP,
for both problems. The second problem is run twice, once explicitly, and once
via CVODE autonomously.

The problems are run using scalar relative and absolute tolerances of
:math:`reltol=10^{-3}` and :math:`abstol=10^{-4}`, respectively.

22 or 38 outputs are printed, and run statistics are printed at the end.



.. _cvDiurnal_kry:

cvDiurnal_kry
================================================

Description
------------

We now investigate a time-dependent system of discretized partial
differential equations.  This example problem is an ODE system generated
from the following 2-species diurnal kinetics advection-diffusion PDE system
in 2 space dimensions:

.. math::

   \frac{dc(i)}{dt} &= K_h\left(\frac{d}{dx}\right)^2 c(i) + V \frac{dc(i)}{dx} + \frac{d}{dy}\left(K_v(y) \frac{dc(i)}{dy}\right) + R_i(c_1,c_2,t) \quad \text{ for } i = 1,2, \quad \text{ where} \\
   R_1(c_1,c_2,t) &= -q_1 c_1 c_3 - q_2 c_1 c_2 + 2 q_3(t) c_3 + q_4(t) c_2 , \\
   R_2(c_1,c_2,t) &=  q_1 c_1 c_3 - q_2 c_1 c_2 - q_4(t) c_2 , \\
   K_v(y) &= K_{v_0} exp\left(\frac{y}{5}\right) ,

:math:`K_h`, :math:`V`, :math:`K_{v_0}`, :math:`q_1`, :math:`q_2`, and
:math:`c_3` are constants, and :math:`q_3(t)` and :math:`q_4(t)` vary
diurnally.  The problem is posed on the square :math:`0 \leq x \leq 20`,
:math:`30 \leq y \leq 50` (all in km), with homogeneous Neumann boundary
conditions, and for time :math:`0 \leq t \leq 86400` seconds (1 day).

The PDE system is treated by central differences on a uniform
:math:`10 \times 10` mesh, with simple polynomial initial profiles.

The problem is solved with CVODE, using the BDF/GMRES method (i.e.
using the SUNLinSol_SPGMR linear solver) and the block-diagonal part
of the Newton matrix as a left preconditioner. A copy of the block-
diagonal part of the Jacobian is saved and conditionally reused within
the Preconditioner routine.


Problem output
---------------

.. include:: ../../../../examples/cvode/serial/cvDiurnal_kry.out
   :literal:


Numerical method
----------------

The example routine solves this problem using a Backwards Differentiation
Formula in fixed-leading coefficient form.  Each stage is solved using the 
built-in modified Newton iteration.  Internally, Newton will use the
SUNLINSOL_SPGMR linear solver via the CVode interface, which in the case
of this scalar-valued problem is just division.  The example file contains
functions to evaluate both :math:`f(t,u)` and :math:`J(t,u)`.  

We specify the relative and absolute tolerances, :math:`reltol=10^{-5}`
and :math:`abstol=10^{-3}`, respectively.  Aside from these choices,
this problem uses only the default CVode solver parameters.

12 outputs are printed, and run statistics are printed at the end.


.. _cvDiurnal_kry_bp:

cvDiurnal_kry_bp
=====================

This example problem is a duplicate of the :ref:`cvDiurnal_kry` problem
above, except here we use the module CVBANDPRE and solve with left
and right preconditioning.


Problem output
---------------

.. include:: ../../../../examples/cvode/serial/cvDiurnal_kry_bp.out
   :literal:


Numerical method
----------------

See ``cvDiurnal_kry`` for details.
   


.. _cvHeat2D_klu:

cvHeat2D_klu
===================

Description
------------

This (future) example problem is a 2D heat equation in a serial 
environment using sparse matrices, and is based on ``idaHeat2D_klu.c`` and
``cvRoberts_klu.c``.

This example solves a discretized 2D heat equation problem, such
that this version uses the KLU sparse solver.

The PDE system solved is a spatial discretization of the PDE

.. math::

   \frac{du}{dt} = \frac{d^2 u}{dx^2} + \frac{d^2 u}{dy^2}

on the unit square. The boundary condition is :math:`u = 0` on
all edges.  Initial conditions are given by 

.. math::
   
   u = 16 \cdot x \cdot (1 - x) \cdot y \cdot (1 - y).

The PDE is treated with central differences on a uniform 
:math:`MGRID \times MGRID` grid.  The values of u at the interior
points satisfy ODEs, and equations :math:`u = 0` at the boundaries
are appended, to form a DAE system of size :math:`N = MGRID^2`. 
Here :math:`MGRID = 10`.

The system is solved with CVODE using the direct sparse linear system
solver, half-bandwidths equal to :math:`M`, and default difference-quotient
Jacobian.

Output is taken at :math:`t = 0, 0.01, 0.02, 0.04, \ldots, 10.24`.


Numerical method
----------------

The example routine solves this problem using a Backwards Differentiation
Formula in fixed-leading coefficient form.  Each stage is solved using the 
built-in modified Newton iteration.  Internally, Newton will use the
SUNLINSOL_KLU linear solver via the CVode interface.  The example file
contains functions to evaluate both :math:`f(t,x,y)` and :math:`J(t,x,y)`.  

We specify the relative and absolute tolerances, :math:`rtol=0`
and :math:`atol=10^{-8}`, respectively.  Aside from these choices,
this problem uses only the default CVode solver parameters.

12 outputs (will be) printed, and run statistics are printed at the end.



.. _cvKrylovDemo_ls:

cvKrylovDemo_ls
========================

Description
------------

Demonstration program for CVODE -- This example loops through the 
available iterative linear solvers:

*  SPGMR -- Scaled, Preconditioned, Generalized Minimum Residual
*  SPFGMR -- Scaled, Preconditioned, Flexible, Generalized Minimum Residual
*  SPBCGS -- Scaled, Preconditioned, Bi-Conjugate Gradient, Stabilized
*  SPTFQMR -- Scaled, Preconditioned, Transpose-Free Quasi-Minimum Residual

The example problem is the same problem as is in :ref:`cvDiurnal_kry_bp`
and :ref:`cvDiurnal_kry` as seen earlier; However, in this case, the
problem is solved with CVODE, with the BDF/GMRES, BDF/FGMRES
BDF/Bi-CGStab, and BDF/TFQMR methods (i.e. using the SUNLinSol_SPGMR,
SUNLinSol_SPFGMR, SUNLinSol_SPBCGS, and SUNLinSol_SPTFQMR linear solvers)
and the block-diagonal part of the Newton matrix as a left preconditioner.

A copy of the block-diagonal part of the Jacobian is saved and
conditionally reused within the Preconditioner routine.


Problem output
---------------

.. include:: ../../../../examples/cvode/serial/cvKrylovDemo_ls.out
   :literal:


Numerical method
----------------

This program solves the same problem four times via the Backwards 
Differentiation Formula method, using a Newton iteration for their
nonlinear solvers.  The first portion of this example uses the
SUNLINSOL_SPGMR linear solver module via the CVode GMRES interface,
while the second portion uses the SUNLINSOL_SPFGMR linear solver module
via the CVode FGMRES interface, the third portion uses the SUNLinSol_SPBCGS
linear solver module via the CVode Bi-CGStab interface, and the fourth
portion uses the SUNLINSOL_SPTFQMR linear solver module via the CVode
TFQMR interface.  Additionally, this example provides a routine to CVode
to compute the densely-referenced, block-diagonal Jacobian. 

The problems are run using scalar relative and absolute tolerances of
:math:`reltol=10^{-5}` and :math:`abstol=10^{-3}`, respectively.

12 outputs are printed and run statistics are printed at the end.



.. _cvKrylovDemo_prec:

cvKrylovDemo_prec
========================

Description
------------

Demonstration program for CVODE - Krylov linear solver.
ODE system from :math:`ns`-species interaction PDE in 2 dimensions.

This program solves a stiff ODE system that arises from a system
of partial differential equations. The PDE system is a food web
population model, with predator-prey interaction and diffusion on
the unit square in two dimensions. The dependent variable vector is:

.. math::

   c = (c^1, c^2, \ldots, c^{ns})

and the PDEs are as follows:

.. math::

   \frac{dc^i}{dt} = d(i) \cdot \left( c_{xx}^i + c_{yy}^i \right) +
      f_i (x, y, c) \qquad (i = 1, \ldots, ns)

where

.. math::
   
   f_i (x, y, c) = c^i \cdot \left( b(i) + \sum_{j = 1}^{ns} a(i, j) c^j \right)

The number of species is :math:`ns = 2 \cdot np`, with the first
:math:`np` being predators. The coefficients :math:`a(i, j)`, :math:`b(i)`,
and :math:`d(i)` are:

.. math::

   \begin{cases}
   a(i, j) = -g \qquad &i \leq np, \enspace j > np \\
   a(i, j) = e \qquad &i > np, \enspace j \leq np \\
   a(i, j) = -a \qquad &\text{otherwise} \\
   \end{cases} \\
   \begin{cases}
   b(i) = b \cdot (1 + \alpha \cdot x \cdot y) &i \leq np \\
   b(i) = -b \cdot (1 + \alpha \cdot x \cdot y) &i > np \\
   \end{cases} \\
   \begin{cases}
   d(i) = D_{\text{prey}} &i \leq np \\
   d(i) = D_{\text{predator}} &i > np \\
   \end{cases}

The spatial domain is the unit square. The final time is 10.  The
boundary conditions are: normal derivative = 0.  A polynomial in
:math:`x` and :math:`y` is used to set the initial conditions.

The PDEs are discretized by central differencing on an 
:math:`MX \text{ by } MY` mesh.

The resulting ODE system is stiff.

The ODE system is solved using Newton iteration and the
SUNLinSol_SPGMR linear solver (scaled preconditioned GMRES).

The preconditioner matrix used is the product of two matrices:
(1) A matrix, only defined implicitly, based on a fixed number of
Gauss-Seidel iterations using the diffusion terms only.  (2) A
block-diagonal matrix based on the partial derivatives of the
interaction terms :math:`f` only, using block-grouping (computing
only a subset of the :math:`ns` by :math:`ns` blocks).

Four different runs are made for this problem.  The product
preconditoner is applied on the left and on the right.  In each
case, both the modified and classical Gram-Schmidt options are
tested. In the series of runs, CVodeInit, SUNLinSol_SPGMR, and
CVSetLinearSolver are called only for the first run, whereas
CVodeReInit, SUNLinSol_SPGMRSetPrecType, and
SUNLinSol_SPGMRSetGSType are called for each of the remaining
three runs.

A problem description, performance statistics at selected output
times, and final statistics are written to standard output.  On
the first run, solution values are also printed at output times.
Error and warning messages are written to standard error, but
there should be no such messages.

Note: This program requires the dense linear solver functions
SUNDlsMat_newDenseMat, SUNDlsMat_newIndexArray,
SUNDlsMat_denseAddIdentity, SUNDlsMat_denseGETRF,
SUNDlsMat_denseGETRS, SUNDlsMat_destroyMat and
SUNDlsMat_destroyArray.

Note: This program assumes the sequential implementation for the
type N_Vector and uses the N_VGetArrayPointer function to gain
access to the contiguous array of components of an N_Vector.


Problem output
---------------

.. include:: ../../../../examples/cvode/serial/cvKrylovDemo_prec.out
   :literal:


Numerical method
-----------------

The example routine solves this problem using a Backwards Differentiation
Formula in fixed-leading coefficient form.  Each stage is solved using the 
built-in modified Newton iteration.  Internally, Newton will use the
SUNLINSOL_SPGMR linear solver via the CVode interface.  The example file
contains functions to evaluate both :math:`f(t,x,y)` and :math:`J(t,x,y)`.  

We specify the relative and absolute tolerances, :math:`rtol=10^{-5}`
and :math:`atol=10^{-5}`, respectively.  The Gram-Schmidt type is set
manually for the SPGMR linear solver, and the Preconditioner is applied
on the left and the right.

6 species matrix-outputs are printed, and run statistics are printed at
the end.  A total of 18 time-steps are taken throughout the problem's
evolution.


References
-----------

Peter N. Brown and Alan C. Hindmarsh, Reduced Storage Matrix Methods
in Stiff ODE Systems, J. Appl. Math. & Comp., 31 (1989), pp. 40-91.
Also available as Lawrence Livermore National Laboratory Report
UCRL-95088, Rev. 1, June 1987.



.. _cvParticle_dns:

cvParticle_dns
=====================

Description
------------

This example problem solves the equation for a particle moving 
counterclockwise with velocity :math:`\alpha` on the unit circle in
the :math:`xy`-plane. The ODE system is given by:

.. math::

   x' &= -\alpha \cdot y \\
   y' &= \alpha \cdot x

where :math:`x` and :math:`y` are subject to the constraint:

.. math::

   x^2 + y^2 - 1 = 0

with initial condition :math:`x = 1` and :math:`y = 0` at :math:`t = 0`.
The system has the analytic solution:

.. math::

   x(t) &= \cos (\alpha \cdot t) \\
   y(t) &= \sin (\alpha \cdot t)

For a description of the command line options for this example, run
the program with the ``--help`` flag.


Problem output
---------------

.. include:: ../../../../examples/cvode/serial/cvParticle_dns.out
   :literal:


Numerical method
----------------

The example routine solves this problem using a Backwards Differentiation
Formula in fixed-leading coefficient form.  Each stage is solved using the 
built-in modified Newton iteration.  Internally, Newton will use the
SUNLINSOL_DENSE linear solver via the CVode interface.  The example file
contains functions to evaluate both :math:`f(t,x,y)` and :math:`J(t,x,y)`.  

We specify the relative and absolute tolerances, :math:`rtol=10^{-4}`
and :math:`atol=10^{-9}`, respectively.  Aside from these choices,
this problem uses only the default CVode solver parameters.

2 output times are printed, and run statistics are printed at the end.




.. _cvPendulum_dns:

cvPendulum_dns
======================

Description
------------

This example problem solves a simple pendulum equation in Cartesian
coordinates where the pendulum bob has mass 1 and is suspended from the
origin with a rod of length 1.  The governing equations are:

.. math::

   x' &= v_x \\
   y' &= v_y \\
   v_x' &= -x \cdot T \\
   v_y' &= -y \cdot T - g

with the constraints:

.. math::

   &x^2 + y^2 - 1 = 0 \\
   &x \cdot v_x + y \cdot v_y = 0

where :math:`x` and :math:`y` are the pendulum bob position, :math:`v_x` and
:math:`v_y` are the bob velocity in the :math:`x` and :math:`y` directions
respectively, T is the tension in the rod, and :math:`g` is acceleration due
to gravity chosen such that the pendulum has period 2.  The initial condition
at :math:`t = 0` is :math:`x = 1`, :math:`y = 0`, :math:`v_x = 0`, and
:math:`v_y = 0`.

A reference solution is computed using the pendulum equation in terms of the
angle between the :math:`x`-axis and the pendulum rod i.e., :math:`\theta in [0, -\pi]`.
The governing equations are:

.. math::

   \theta' &= v_\theta \\
   v_\theta' &= -g \cdot \cos (\theta)

where :math:`\theta` is the angle from the :math:`x`-axis, :math:`v_\theta` is
the angular velocity, and :math:`g` the same acceleration due to gravity from above.
The initial condition at :math:`t = 0` is :math:`\theta = 0` and :math:`v_theta = 0`.

The Cartesian formulation is run to a final time :math:`t_f` (default 30) with
and without projection for various integration tolerances. The error in the
position and velocity at :math:`t_f` compared to the reference solution, the
error in the position constraint equation, and various integrator statistics are
printed to the screen for each run.

When projection is enabled a user-supplied function is used to project the
position, velocity, and error to the constraint manifold.

Optional command line inputs may be used to change the final simulation time
(default 30), the initial tolerance (default :math:`10^{-5}`), the number of outputs
(default 1), or disable error projection. Use the option ``--help`` for a list
of the command line flags.

 
Problem output
---------------

.. include:: ../../../../examples/cvode/serial/cvPendulum_dns.out
   :literal:


Numerical method
----------------

The example routine solves this problem using a Backwards Differentiation
Formula in fixed-leading coefficient form.  Each stage is solved using the 
built-in modified Newton iteration.  Internally, Newton will use the
SUNLINSOL_DENSE linear solver via the CVode interface.  The example file
contains functions to evaluate both :math:`f(t,x,y)` and :math:`J(t,x,y)`.  

We specify the relative and absolute tolerances for several runs to be

*  :math:`rtol=10^{-5}` and :math:`atol=10^{-5}`,
*  :math:`rtol=10^{-6}` and :math:`atol=10^{-6}`,
*  :math:`rtol=10^{-7}` and :math:`atol=10^{-7}`,
*  :math:`rtol=10^{-8}` and :math:`atol=10^{-8}`,
*  :math:`rtol=10^{-9}` and :math:`atol=10^{-9}`,

respectively.  Aside from these choices, this problem uses only the
default CVode solver parameters.

2 output times are printed for each solve, one with a user-supplied
Projection function, and one without.



.. _cvRoberts_dns:

cvRoberts_dns
============================================

Description
------------

This example is a simple problem, with the coding needed for its
solution completed via CVODE. The problem is from chemical kinetics,
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

.. include:: ../../../../examples/cvode/serial/cvRoberts_dns.out
   :literal:


Numerical method
----------------

The example routine solves this problem using a Backwards Differentiation
Formula in fixed-leading coefficient form.  Each stage is solved using the 
built-in modified Newton iteration.  Internally, Newton will use the
SUNLINSOL_DENSE linear solver via the CVode interface.  The example file
contains functions to evaluate both :math:`f(t, y_1, y_2, y_3)` and 
:math:`J(t, y_1, y_2, y_3)`.  Additionally, a root-finding function will,
as previously mentioned, find roots at :math:`y_3 = 0.01` and 
:math:`y_1 = 10^{-4}` using the built-in root-finding mechanism in CVode.

We specify the scalar relative and vector-valued absolute tolerances, :math:`rtol=10^{-4}`
and :math:`atol=\begin{pmatrix} 10^{-8} \\ 10^{-14} \\ 10^{-6} \end{pmatrix}`
, respectively.  Aside from these choices, this problem uses only the default
CVode solver parameters.

11 normal + 2 root output times are printed at multiplicatively equally-spaced
points as well as at the two roots, and run statistics are printed at the end.



.. _cvRoberts_dns_constraints:

cvRoberts_dns_constraints
==========================================

Description
------------

This example problem is the same as :ref:`cvRoberts_dns` above except that the
constraint :math:`y_i \geq 0` is posed for all components :math:`i = 1, 2, 3`.


Problem output
---------------

.. include:: ../../../../examples/cvode/serial/cvRoberts_dns_constraints.out
   :literal:


Numerical method
----------------

Here, we specify the scalar relative and vector-valued absolute tolerances,
:math:`rtol=10^{-4}` and :math:`atol=\begin{pmatrix} 10^{-6} \\ 10^{-11} \\ 10^{-5} \end{pmatrix}`
, respectively.

Aside from this, see ``cvRoberts_dns`` above.
   


.. _cvRoberts_dns_negsol:

cvRoberts_dns_negsol
================================

Description
------------

This example problem is the same as :ref:`cvRoberts_dns` above except that here
we allow negative solutions to take form in the first run, and intercept them
via CVode on the second run.


Problem output
---------------

.. include:: ../../../../examples/cvode/serial/cvRoberts_dns_negsol.out
   :literal:


Numerical method
----------------

Here, we specify the scalar relative and vector-valued absolute tolerances,
:math:`rtol=10^{-4}` and :math:`atol=\begin{pmatrix} 10^{-7} \\ 10^{-13} \\ 10^{-5} \end{pmatrix}`
, respectively.

Aside from this, see ``cvRoberts_dns`` above.



.. _cvRoberts_dns_uw:

cvRoberts_dns_uw
======================

Description
------------

This example problem is the same as :ref:`cvRoberts_dns` above except that
here it uses a user-supplied function to compute the error weights
required for the WRMS norm calculations.


Problem output
---------------

.. include:: ../../../../examples/cvode/serial/cvRoberts_dns_uw.out
   :literal:


Numerical method
-----------------

See ``cvRoberts_dns`` above.



.. _cvRoberts_dnsL:

cvRoberts_dnsL
====================

Description
------------

This example problem is the same as :ref:`cvRoberts_dns` above except that
here we use the LAPACK dense linear solver instead of the dense linear
solver.


Problem output
---------------

.. include:: ../../../../examples/cvode/serial/cvRoberts_dnsL.out
   :literal:


Numerical method
-----------------

Aside from using SUNLINSOL_LAPACKDENSE, see ``cvRoberts_dns`` above.



.. _cvRoberts_klu:

cvRoberts_klu
===================

Description
------------

This example problem is the same as :ref:`cvRoberts_dns` above except that
here we use the KLU sparse linear solver instead of the dense linear
solver.


Problem output
---------------

.. include:: ../../../../examples/cvode/serial/cvRoberts_klu.out
   :literal:


Numerical method
-----------------

Aside from using SUNLINSOL_KLU, see ``cvRoberts_dns`` above.



.. _cvRoberts_block_klu:

cvRoberts_block_klu
==============================

Description
------------

This example problem is the same as :ref:`cvRoberts_klu` above except
that in this example problem, we simulate a scenario where a set of
independent ODEs are grouped together to form a larger system.

This program takes one optional argument, the number of groups
of independent ODE systems:

.. math::

   \qquad \texttt{./cvRoberts_block_klu [number of groups]}

The problem is comparable to the CUDA version -
``cvRoberts_block_cusolversp_batchqr.cu``.  As previously suggested,
it was based off of the ``cvRoberts_klu.c`` example.
 

Problem output
---------------

.. include:: ../../../../examples/cvode/serial/cvRoberts_block_klu.out
   :literal:


Numerical method
----------------

See ``cvRoberts_klu`` and ``cvRoberts_dns`` above.



.. _cvRoberts_sps:

cvRoberts_sps
===================

Description
------------

This example problem is the same as :ref:`cvRoberts_dns` above except
that in this example problem, we solve the problem with the SuperLUMT
sparse direct linear solver.


Problem output
---------------

.. include:: ../../../../examples/cvode/serial/cvRoberts_sps.out
   :literal:


Numerical method
-----------------

Aside from using SUNLINSOL_SUPERLUMT along with the necessary sparse
matrices, see ``cvRoberts_dns`` above.



.. _cvRocket_dns:

cvRocket_dns
=================

Description
------------

The following is a simple example problem, with the coding needed
for its solution provided by CVODE. The problem is a simpliflied model
of a rocket, ascending vertically, with mass decreasing over time.  The
system (of size 2) is given by:

.. math::

   y_1 &= \text{rocket height} \enspace H, &&\text{where } y_1(0) = 0, \\
   y_2 &= \text{rocket velocity} \enspace v, &&\text{where } y_2(0) = 0, \\
   \frac{dH}{dt} &= v, \\
   \frac{dv}{dt} &= a(t,v).

The upward acceleration :math:`a(t,v)` is given by:

.. math::

   a(t,v) = \frac{F}{M_r + M_f} - D \cdot v - g,

where :math:`F =` engine thrust force (constant), :math:`M_r =` rocket mass
without fuel, :math:`M_f =` fuel mass :math:`= M_{f_0} - r \cdot t`,
:math:`r =` fuel burn rate, :math:`D =` drag coefficient, and :math:`g =`
gravitational acceleration.

The engine force is reset to 0 when the fuel mass reaches 0, or when
:math:`H` reaches a preset height :math:`H_c`, whichever happens first.
Root-finding is used to locate the time at which :math:`M_f = 0` or
:math:`H = H_c`, and also the time at which the rocket reaches its maximum
height, given by the condition :math:`v = 0` and :math:`t > 0`.

The problem is solved with the BDF method and Dense linear solver.

Run statistics (optional outputs) are printed at the end.


Problem output
---------------

.. include:: ../../../../examples/cvode/serial/cvRocket_dns.out
   :literal:


Numerical method
-----------------

The example routine solves this problem using a Backwards Differentiation
Formula in fixed-leading coefficient form.  Each stage is solved using the 
built-in modified Newton iteration.  Internally, Newton will use the
SUNLINSOL_DENSE linear solver via the CVode interface.  The example file
contains functions to evaluate both :math:`f(t, y_1, y_2)` and 
:math:`J(t, y_1, y_2)`.  Additionally, a root-finding function will,
as previously mentioned, find roots at either :math:`M_f = 0` or :math:`H = H_c` 
as well as at :math:`v = 0` for :math:`t > 0` using the built-in root-finding
mechanism in CVode.

We specify the scalar relative and vector-valued absolute tolerances, :math:`rtol=10^{-5}`
and :math:`atol=\begin{pmatrix} 0.01 \\ 0.1 \end{pmatrix}`
, respectively.  Aside from these choices, this problem uses only the default
CVode solver parameters.

68 normal + 2 root output times are printed at normatively equally-spaced
points as well as at the two roots, and run statistics are printed at the end.



