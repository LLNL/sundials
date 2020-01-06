..
   Programmer(s): Daniel R. Reynolds @ SMU
   ----------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2002-2020, Lawrence Livermore National Security
   and Southern Methodist University.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   ----------------------------------------------------------------

:tocdepth: 3



.. _SUNLinSol.ARKode:

ARKode SUNLinearSolver interface
==============================================

In the table below, we list the SUNLinSol module linear solver
functions used within the ARKLS interface.  As with the SUNMATRIX module, we
emphasize that the ARKode user does not need to know detailed usage of linear
solver functions by the ARKode code modules in order to use ARKode. The
information is presented as an implementation detail for the interested reader.

The linear solver functions listed below are marked with "X" to
indicate that they are required, or with "O" to indicate that they are
only called if they are non-``NULL`` in the ``SUNLinearSolver``
implementation that is being used.  Note:

1. :c:func:`SUNLinSolNumIters()` is only used to accumulate overall
   iterative linear solver statistics.  If it is not implemented by
   the ``SUNLinearSolver`` module, then ARKLS will consider all
   solves as requiring zero iterations.

2. Although :c:func:`SUNLinSolResNorm()` is optional, if it is not
   implemented by the ``SUNLinearSolver`` then ARKLS will consider all
   solves a being *exact*.

3. Although ARKLS does not call :c:func:`SUNLinSolLastFlag()`
   directly, this routine is available for users to query linear
   solver failure modes directly.

4. Although ARKLS does not call :c:func:`SUNLinSolFree()`
   directly, this routine should be available for users to call when
   cleaning up from a simulation.



.. cssclass:: table-bordered

===========================  ======  =========  ================
Routine                      DIRECT  ITERATIVE  MATRIX_ITERATIVE
===========================  ======  =========  ================
SUNLinSolGetType             X       X          X
SUNLinSolSetATimes           O       X          O
SUNLinSolSetPreconditioner   O       O          O
SUNLinSolSetScalingVectors   O       O          O
SUNLinSolInitialize          X       X          X
SUNLinSolSetup               X       X          X
SUNLinSolSolve               X       X          X
SUNLinSolNumIters\ :sup:`1`          O          O
SUNLinSolResNorm\ :sup:`2`           O          O
SUNLinSolLastFlag\ :sup:`3`
SUNLinSolFree\ :sup:`4`
SUNLinSolSpace               O       O          O
===========================  ======  =========  ================


Since there are a wide range of potential SUNLinSol use cases, the following
subsections describe some details of the ARKLS interface, in the case that
interested users wish to develop custom SUNLinSol modules.


.. _SUNLinSol.Lagged_matrix:

Lagged matrix information
---------------------------------------------------

If the SUNLinSol identifies as having type
``SUNLINEARSOLVER_DIRECT`` or ``SUNLINEARSOLVER_MATRIX_ITERATIVE``,
then the SUNLinSol object solves a
linear system *defined* by a SUNMATRIX object. ARKLS will update the
matrix information infrequently according to the strategies outlined in
the section :ref:`Mathematics.Linear.Setup`.  To this end, we
differentiate between the *desired* linear system
:math:`\mathcal A x = b` with :math:`\mathcal A = (M-\gamma J)` 
and the *actual* linear system

.. math::
   \tilde{\mathcal A} \tilde{x} = b \quad\Leftrightarrow\quad (M-\tilde{\gamma} J)\tilde{x} = b.

Since ARKLS updates the SUNMATRIX object infrequently, it is likely
that :math:`\gamma\ne\tilde{\gamma}`, and in turn :math:`\mathcal
A\ne\tilde{\mathcal A}`.  Therefore, after calling the
SUNLinSol-provided :c:func:`SUNLinSolSolve()` routine, we test whether
:math:`\gamma / \tilde{\gamma} \ne 1`, and if this is the case we
scale the solution :math:`\tilde{x}` to obtain the desired linear
system solution :math:`x` via

.. math::
   x = \frac{2}{1 + \gamma / \tilde{\gamma}} \tilde{x}.
   :label: eq:rescaling

The motivation for this selection of the scaling factor
:math:`c = 2/(1 + \gamma/\tilde{\gamma})` follows the derivation in
[BBH1989]_ and [H2000]_.  In short, if we consider a stationary
iteration for the linear system as consisting of a solve with
:math:`\tilde{\mathcal A}` followed by scaling by :math:`c`,
then for a linear constant-coefficient problem, the error in the
solution vector will be reduced at each iteration by the error matrix
:math:`E = I - c \tilde{\mathcal A}^{-1} \mathcal A`, with a
convergence rate given by the spectral radius of :math:`E`.  Assuming
that stiff systems have a spectrum spread widely over the left
half-plane, :math:`c` is chosen to minimize the magnitude of the
eigenvalues of :math:`E`.


.. _SUNLinSol.Iterative_Tolerance:

Iterative linear solver tolerance
---------------------------------------------------

If the SUNLinSol object self-identifies as having type
``SUNLINEARSOLVER_ITERATIVE`` or ``SUNLINEARSOLVER_MATRIX_ITERATIVE``,
then ARKLS will set the input tolerance ``delta`` as described in
:ref:`Mathematics.Error.Linear`.  However, if the iterative linear
solver does not support scaling matrices (i.e., the
:c:func:`SUNLinSolSetScalingVectors()` routine is ``NULL``), then
ARKLS will attempt to adjust the linear solver tolerance to account
for this lack of functionality.  To this end, the following
assumptions are made:

* The units of the IVP solution and linear residual are the same
  (i.e., the error and residual weight vectors in section
  :ref:`Mathematics.Error.Norm` are the same); this is automatically
  satisfied with identity mass matrix, :math:`M=I`, or similar.

* All solution components have similar magnitude; hence the error
  weight vector :math:`w` used in the WRMS norm (see the section
  :ref:`Mathematics.Error.Norm`) should satisfy the assumption

  .. math::
     w_i \approx w_{mean},\quad \text{for}\quad i=0,\ldots,n-1.

* The SUNLinSol object uses a standard 2-norm to measure convergence.

Under these assumptions, ARKLS uses identical left and right scaling matrices,
:math:`S_1 = S_2 = S = \operatorname{diag}(w)`, so the linear solver
convergence requirement is converted as follows
(using the notation from the beginning of this chapter):

.. math::
   &\left\| \tilde{b} - \tilde{A} \tilde{x} \right\|_2  <  \text{tol}\\
   \Leftrightarrow \quad & \left\| S P_1^{-1} b - S P_1^{-1} A x \right\|_2  <  \text{tol}\\
   \Leftrightarrow \quad & \sum_{i=0}^{n-1} \left[w_i \left(P_1^{-1} (b - A x)\right)_i\right]^2  <  \text{tol}^2\\
   \Leftrightarrow \quad & w_{mean}^2 \sum_{i=0}^{n-1} \left[\left(P_1^{-1} (b - A x)\right)_i\right]^2  <  \text{tol}^2\\
   \Leftrightarrow \quad & \sum_{i=0}^{n-1} \left[\left(P_1^{-1} (b - A x)\right)_i\right]^2  <  \left(\frac{\text{tol}}{w_{mean}}\right)^2\\
   \Leftrightarrow \quad & \left\| P_1^{-1} (b - A x)\right\|_2  <  \frac{\text{tol}}{w_{mean}}

Therefore the tolerance scaling factor

.. math::
   w_{mean} = \|w\|_2 / \sqrt{n}

is computed and the scaled tolerance ``delta`` :math:`= \text{tol} /
w_{mean}` is supplied to the SUNLinSol object.
