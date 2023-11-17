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

.. _openmp_c:

====================================
OpenMP C example problems
====================================




.. _cvAdvDiff_bnd_omp:

cvAdvDiff_bnd_omp
============================================

Description
------------

This simple example problem is an OpenMP implementation of previous
serial C example ``cvAdvDiff_bnd`` showing how to use the CVode solver
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

Optionally, we can set the number of threads from environment
variable or command line. To check the current value for number of
threads from environment:

.. code-block:: bash

   echo $OMP_NUM_THREADS

Execution
----------

To use the default value or the number of threads from the environment
value, run without arguments:

.. code-block:: bash

   ./cvAdvDiff_bnd_omp

The environment variable can be over-ridden with a command line
argument specifying the number of threads to use, e.g:

.. code-block:: bash

   ./cvAdvDiff_bnd_omp 5



Problem output
---------------

.. include:: ../../../../examples/cvode/C_openmp/cvAdvDiff_bnd_omp.out
   :literal:


Numerical method
----------------

The numerical method is identical to the previous implementation,
except that we now use SUNDIALS' OpenMP-enabled vector kernel module,
NVECTOR_OPENMP, and have similarly threaded the supplied right-hand
side residual and banded Jacobian construction functions.
