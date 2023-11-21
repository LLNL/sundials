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


.. _parallel_f2003:

==============================
Parallel Fortran 2003 Examples
==============================



.. _cv_diag_kry_bbd_p:

cv_diag_kry_bbd_p
=========================

Description
------------

This problem models a stiff, linear, diagonal
ODE system,

.. math::

   \frac{\partial y_i}{\partial t} = -\alpha i y_i, \quad i=1,\ldots N.


Here :math:`\alpha=10` and :math:`N=10 * N_P`, where :math:`N_P` is the
number of MPI tasks used for the problem.  The problem has initial
conditions :math:`y_i=1` and evolves for the time interval :math:`t\in
[0,1]`.


Problem output
---------------

.. literalinclude:: ../../../../examples/cvode/F2003_parallel/cv_diag_kry_bbd_f2003.out
   :language: text


Numerical method
----------------

This program solves the problem with a BDF method, using a Newton
iteration with the preconditioned SUNLINSOL_SPGMR iterative linear
solver module.

A diagonal preconditioner matrix is used, formed automatically through
difference quotients within the CVBBDPRE module.  Since CVBBDPRE is
developed for use of a block-banded preconditioner, in this solver
each block is set to have half-bandwidths ``mudq = mldq = 0`` to
retain only the diagonal portion.

Two runs are made for this problem, first with left and then with
right preconditioning (``iPre`` is first set to :math:`1` and then to :math:`2`).

Performance data is printed at equally-spaced output times, and maximum
errors and final performance counters are printed on completion.



.. _cv_diag_kry_p:

cv_diag_kry_p
=================

Description
------------

This is basically the same probelm mathematically as is described above
in ``cv_diag_kry_bbd_p``.


Problem output
---------------

.. literalinclude:: ../../../../examples/cvode/F2003_parallel/cv_diag_kry_f2003.out
   :language: text


Numerical method
-----------------

This program solves the problem with a BDF method, using a Newton
iteration with the preconditioned SUNLINSOL_SPGMR iterative linear
solver module.

A user-supplied preconditioner matrix is used, where in the solver
routine uses a diagonal preconditioner :math:`P = I - \gamma \cdot J`,
where :math:`J` is a diagonal approximation to the true Jacobian, given by:

.. math::

   J = diag(0, 0, 0, -4 \alpha, \ldots, -N \alpha).

The vector :math:`r` is copied to :math:`z`, and the inverse of :math:`P`
(restricted to the local vector segment) is applied to the vecdtor :math:`z`.

Performance data is printed at equally-spaced output times, and maximum
errors and final performance counters are printed on completion.



.. _cv_diag_non_p:

cv_diag_non_p
=================

Description
------------

This problem models a nonstiff, linear, diagonal
ODE system,

.. math::

   \frac{\partial y_i}{\partial t} = -\alpha i y_i, \quad i=1,\ldots N.


Here :math:`\alpha=\frac{10}{N}` and :math:`N=10 N_P`, where :math:`N_P` is the
number of MPI tasks used for the problem.  The problem has initial
conditions :math:`y_i=1` and evolves for the time interval :math:`t\in [0,1]`.



Problem output
---------------

.. literalinclude:: ../../../../examples/cvode/F2003_parallel/cv_diag_non_f2003.out
   :language: text


Numerical method
-----------------

This program solves the problem with an ADAMS method, using a Fixed-point
iteration.

Performance data is printed at equally-spaced output times, and maximum
errors and final performance counters are printed on completion.
