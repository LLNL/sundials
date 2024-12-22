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

.. _CVODE.Alg.Resizing:

Resizing
--------

A time step with step size :math:`h_n = t_{n-1} - t_{n}` using a method of order
:math:`q` begins with the polynomial interpolant, :math:`\pi_{n-1}`, satisfying
the :math:`q + 1` conditions for Adams methods

.. math::

   \pi_{n-1}(t_{n-1}) = y_{n-1}
   \dot{\pi}_{n-1}(t_{n-j}) = \dot{y}_{n-j}, \quad j = 1,2,\ldots,q
