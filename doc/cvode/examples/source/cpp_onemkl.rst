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


.. _onemkl_cpp:

====================================
OneMKL C++ example problems
====================================

.. _cvRoberts_blockdiag_onemkl:

cvRoberts_blockdiag_onemkl
==================================

Description
------------

This example is a simple problem, with the coding needed for its
solution completed via CVODE. The problem is from chemical kinetics,
and consists of the following three rate equations:

.. math::

   \frac{dy_1}{dt} &= -0.04 \cdot y_1 + 10^4 \cdot y_2 \cdot y_3 \\
   \frac{dy_2}{dt} &= -\frac{dy_1}{dt} - \frac{dy_3}{dt} \\
   \frac{dy_3}{dt} &= 3 \times 10^7 \cdot y_2^2

on the interval from :math:`t = 0.0` to :math:`t = 4 \times 10^{10}`,
with initial conditions: :math:`y_1 = 1.0`, :math:`y_2 = y_3 = 0`.
The problem is stiff.

This program solves the problem with the BDF method, Newton iteration
with the dense linear solver, and a user-supplied Jacobian routine.
since the grouping of the independent systems results in a block diagonal linear
system, with the oneMKL SUNLinearSolver.  Alternatively, the SPGMR linear
solver may be used with a Jacobi preconditioner.

The problem uses a scalar relative tolerance and a vector absolute tolerance.
Output is printed in decades from :math:`t = 0.1` to :math:`10^6`. Run statistics
(optional outputs) are printed at the end.

The program takes three optional argument: the number of independent ODE
systems (default 100), a flag to use a direct (1, default) or iterative (0)
linear solver, and a flag to enable (1, default) or disable (0) solution
output:

.. math::

   /texttt{./cvRoberts\_blockdiag\_onemkl [number of groups] [solver type] [output]}

This problem is comparable to the ``cvRoberts_block_klu.c`` example.


Problem output
---------------

.. include:: ../../../../examples/cvode/CXX_onemkl/cvRoberts_blockdiag_onemkl.out
   :literal:


Numerical method
-----------------

The example routine solves this problem using a Backwards Differentiation
Formula in fixed-leading coefficient form.  Each stage is solved using the
built-in modified Newton iteration.  Internally, Newton will use the
SUNLINSOL_ONEMKLDENSE linear solver via the CVode interface.  Alternatively,
Newton can use SUNLINSOL_SPGMR linear solver via the CVode interface.  The
example file contains functions to evaluate both :math:`f(t, y_1, y_2, y_3)`
and :math:`J(t, y_1, y_2, y_3)`.

We specify the scalar relative and vector-valued absolute tolerances, :math:`rtol=10^{-4}`
and :math:`atol=\begin{pmatrix} 10^{-8} \\ 10^{-14} \\ 10^{-6} \end{pmatrix}`
, respectively.  Aside from these choices, this problem uses only the default
CVode solver parameters.

8 times are printed with 10 groups at each multiplicatively equally-spaced time-step,
and run statistics are printed at the end.
