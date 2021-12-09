..
   Daniel R. Reynolds @ SMU
   ----------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2002-2021, Lawrence Livermore National Security
   and Southern Methodist University.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   ----------------------------------------------------------------

.. _SUNLinSol:

######################
Solving Linear Systems
######################

For problems that require the solution of linear systems of equations,
the SUNDIALS packages operate using generic linear solver modules
defined through the ``SUNLinearSolver``, or ``SUNLinSol``, API.  This allows SUNDIALS
packages to utilize any valid SUNLinSol implementation that provides
a set of required functions.  These functions can be divided into
three categories.  The first are the core linear solver functions.  The
second group consists of "set" routines to supply the linear solver object
with functions provided by the SUNDIALS package, or for modification
of solver parameters.  The last group consists of "get" routines for
retrieving artifacts (statistics, residual vectors, etc.) from the
linear solver.  All of these functions are defined in the header file
``sundials/sundials_linearsolver.h``.

The implementations provided with SUNDIALS work in coordination
with the SUNDIALS generic ``N_Vector`` and ``SUNMatrix`` modules to
provide a set of compatible data structures and solvers for the
solution of linear systems using direct or iterative (matrix-based or matrix-free)
methods. Moreover, advanced users can provide a customized
``SUNLinearSolver`` implementation to any SUNDIALS package,
particularly in cases where they provide their own ``N_Vector`` and/or
``SUNMatrix`` modules.

Historically, the SUNDIALS packages have been designed to specifically
leverage the use of either *direct linear solvers* or matrix-free,
*scaled, preconditioned, iterative linear solvers*.  However, user-supplied
implementations for matrix-based iterative linear solvers and linear solvers
with 'embedded' matrices are also supported.

The iterative linear solvers packaged with SUNDIALS leverage scaling
and preconditioning, as applicable, to balance error between solution
components and to accelerate convergence of the linear solver.  To
this end, instead of solving the linear system :math:`Ax = b`
directly, these apply the underlying iterative algorithm to the
transformed system

.. math::
   :label: eq:transformed_linear_system

   \tilde{A} \tilde{x} = \tilde{b}

where

.. math::
  :label: eq:transformed_linear_system_components

  \tilde{A} &= S_1 P_1^{-1} A P_2^{-1} S_2^{-1},\\
  \tilde{b} &= S_1 P_1^{-1} b,\\
  \tilde{x} &= S_2 P_2 x,

and where

* :math:`P_1` is the left preconditioner,

* :math:`P_2` is the right preconditioner,

* :math:`S_1` is a diagonal matrix of scale factors for
  :math:`P_1^{-1} b`,

* :math:`S_2` is a diagonal matrix of scale factors for :math:`P_2 x`.

SUNDIALS solvers request that iterative linear solvers stop
based on the 2-norm of the scaled preconditioned residual meeting a
prescribed tolerance

.. math::

   \left\| \tilde{b} - \tilde{A} \tilde{x} \right\|_2  <  \text{tol}.


When provided an iterative SUNLinSol implementation that does not
support the scaling matrices :math:`S_1` and :math:`S_2`, SUNDIALS'
packages will adjust the value of :math:`\text{tol}` accordingly
(see the section :numref:`SUNLinSol.Iterative_Tolerance` for more details).  In
this case, they instead request that iterative linear solvers stop
based on the criteria

.. math::

   \left\| P_1^{-1} b - P_1^{-1} A x \right\|_2  <  \text{tol}.

We note that the corresponding adjustments to :math:`\text{tol}` in
this case are non-optimal, in that they cannot balance error between
specific entries of the solution :math:`x`, only the aggregate error
in the overall solution vector.

We further note that not all of the SUNDIALS-provided iterative linear
solvers support the full range of the above options (e.g., separate
left/right preconditioning), and that some of the SUNDIALS packages
only utilize a subset of these options.  Further details on these
exceptions are described in the documentation for each
``SUNLinearSolver`` implementation, or for each SUNDIALS package.


For users interested in providing their own SUNLinSol module, the
following section presents the SUNLinSol API and its implementation
beginning with the definition of SUNLinSol functions in sections
:numref:`SUNLinSol.CoreFn` -- :numref:`SUNLinSol.GetFn`. This is followed by
the definition of functions supplied to a linear solver implementation in
section :numref:`SUNLinSol.SUNSuppliedFn`. The linear solver return
codes are described in section :numref:`SUNLinSol.ErrorCodes`. The
``SUNLinearSolver`` type and the generic SUNLinSol module are defined
in section :numref:`SUNLininSol.Generic`.  The section
:numref:`SUNLinSol.Compatibility` discusses compatibility between the
SUNDIALS-provided SUNLinSol modules and SUNMATRIX modules.  Section
:numref:`SUNLinSol.API.Custom` lists the requirements for supplying a
custom SUNLinSol module and discusses some intended use cases. Users wishing to
supply their own SUNLinSol module are encouraged to use the SUNLinSol
implementations provided with SUNDIALS as a template for supplying custom
linear solver modules. The SUNLinSol functions required by this SUNDIALS
package as well as other package specific details are given in
section :numref:`SUNLinSol.ARKODE`. The remaining sections of this chapter
present the SUNLinSol modules provided with SUNDIALS.


.. toctree::
   :maxdepth: 1

   SUNLinSol_API_link.rst
   ARKode_requirements
   SUNLinSol_links.rst
