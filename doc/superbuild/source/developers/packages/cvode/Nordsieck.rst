..
   Author(s): David J. Gardner @ LLNL
   -----------------------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2002-2024, Lawrence Livermore National Security
   and Southern Methodist University.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   -----------------------------------------------------------------------------

.. _CVODE.Alg.Nordsieck:

The Nordsieck History Array
---------------------------

A step from time :math:`t_{n-1}` to :math:`t_n` with step size :math:`h_n =
t_{n} - t_{n-1}` using a method of order :math:`q` begins with the polynomial
interpolant, :math:`\pi_{n-1}`, satisfying the :math:`q + 1` conditions for
Adams methods,

.. math::

   \pi_{n-1}(t_{n-1}) &= y_{n-1} \\
   \dot{\pi}_{n-1}(t_{n-j}) &= \dot{y}_{n-j}, \quad j = 1,2,\ldots,q

or the conditions for BDF methods

.. math::

   \pi_{n-1}(t_{n-1}) = y_{n-1}

The solution at :math:`t_n` is obtained by constructing the polynomial
interpolant, :math:`\pi_n`, satisfying the :math:`q + 2` conditions for Adams
methods

.. math::

   \pi_{n}(t_n) &= y_n, \\
   \pi_{n-1}(t_{n-1}) &= y_{n-1}, \\
   \dot{\pi}(t_{n-j}) &= \dot{y}_{n-j}

or the conditions for BDF methods

.. math::

   \pi_{n}(t_n) = y_n

The solution history is represented using the Nordsieck array given by

.. math::

   Z_{n-1} = [ y_{n-1},\, h_n \dot{y}_{n-1},\, \frac{1}{2} h_n^2 \ddot{y}_{n-1},\, \ldots,\, \frac{1}{q!} h_n^q y^{(q)}_{n-1}]

where the derivatives, :math:`y^{(j)}_{n-1}` are approximated by
:math:`\pi^{(j)}_{n-1}(t_{n-1})`. With this formulation, computing the time step
is now a problem of constructing :math:`Z_n` from :math:`Z_{n-1}`. This
computation is done by forming the predicted array, :math:`Z_{n(0)}`, with
columns :math:`\frac{1}{j!} h^{(j)}_{n} \pi^{(j)}_{n-1}(t_n)` and computing the
coefficients :math:`l = [l_0,\, l_1,\, \ldots,\, l_q]` such that

.. math::

   Z_n = Z_{n(0)} + c_n l,

where :math:`c_n` is the correction to the predicted solution i.e.,

.. math::

   c_n = y_n - y_{n(0)} = \pi_n(t_n) - \pi_{n-1}(t_n)

From the second column of ?? we obtain the nonlinear system

.. math::

   F(y_n) \equiv y_n - \gamma f(t_n, y_n) - a_n = 0

for computing :math:`y_n` where :math:`\gamma \equiv h_n / l_1` and :math:`a_n
\equiv y_{n(0)} - \gamma \dot{y}_{n(0)}`.

Predicting the History Array
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The predicted Nordsieck array is computed as the matrix product

.. math::

   Z_{n(0)} = Z_{n-1} A_q,

where :math:`A_q` is an order :math:`q + 1` Pascal triangle matrix. For example,
with :math:`q = 5` the matrix is

.. math::

   A_5 =
   \begin{bmatrix}
   1 &   &    &    &   &   \\
   1 & 1 &    &    &   &   \\
   1 & 2 & 1  &    &   &   \\
   1 & 3 & 3  & 1  &   &   \\
   1 & 4 & 6  & 4  & 1 &   \\
   1 & 5 & 10 & 10 & 5 & 1 \\
   \end{bmatrix}.

The product :math:`Z_{n-1} A_q` is computed in-place using repeated additions

.. code::

   for (k = 1; k <= q; k++)
   {
     for (j = q; j <= k; j++)
     {
       Z[:][j-1] += Z[:][j]
     }
   }

Correcting the History Array
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Recall the BDF is
