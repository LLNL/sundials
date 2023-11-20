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

.. _parhyp_c:

====================================
Parallel Hypre C example problems
====================================




.. _cvAdvDiff_non_ph:

cvAdvDiff_non_ph
============================================

This problem is mathematically identical to the parallel C example
problem :ref:`cvAdvDiff_non_p`.  As before, the problem is the
semi-discrete form of the advection-diffusion equation in 1-D:

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

.. literalinclude:: ../../../../examples/cvode/parhyp/cvAdvDiff_non_ph.out
   :language: text


Numerical method
----------------

The numerical method is identical to the previous implementation,
except that we now use the HYPRE parallel vector module,
NVECTOR_PARHYP.  The output of these two examples is identical.
Below, we discuss only the main differences between the two
implementations; familiarity with the HYPRE library is helpful.

We use the HYPRE IJ vector interface to allocate the template vector
and create the parallel partitioning:

.. code-block:: c

   HYPRE_IJVectorCreate(comm, my_pe*local_N, (my_pe + 1)*local_N - 1, &Uij);
   HYPRE_IJVectorSetObjectType(Uij, HYPRE_PARCSR);
   HYPRE_IJVectorInitialize(Uij);

The *initialize* call means that vector elements are ready to be set using
the IJ interface. We choose the initial condition vector :math:`x_0 =
x(t_0)` as the template vector, and we set its values in the
``SetInitialProfiles(...)`` function. We complete assembly of the
HYPRE vector with the calls:

.. code-block:: c

   HYPRE_IJVectorAssemble(Uij);
   HYPRE_IJVectorGetObject(Uij, (void**) &Upar);

The *assemble* call is collective and makes the HYPRE vector ready to
use.  The sets the handle ``Upar`` to the actual HYPRE vector.  The
handle is then passed to the ``N_VMake`` function, which creates the
template ``N_Vector``, ``u``, as a wrapper around the HYPRE vector.
All of the other vectors used in the computation are created by
cloning this template vector.

Furthermore, since the template vector does not own the underlying
HYPRE vector (it was created using the ``HYPRE_IJVectorCreate`` call
above), so it is the user's responsibility to destroy it by calling
``HYPRE_IJVectorDestroy(Uij)`` after the template vector ``u`` has
been destroyed.  This function will destroy both the HYPRE vector
and its IJ interface.

To access individual elements of the solution and derivative vectors
``u`` and ``udot`` in the IVP right-hand side function, ``f``, the
user needs to first extract the HYPRE vector by calling
``N_VGetVector_ParHyp``, and then use HYPRE-specific methods to access
the data from that point on.

.. note::

   Currently, interfaces to HYPRE solvers and preconditioners are not
   available directly through the SUNDIALS interfaces, however these
   could be utilized directly as preconditioners for SUNDIALS' Krylov
   solvers.  Direct interfaces to the HYPRE solvers will be provided
   in subsequent SUNDIALS releases.  The current HYPRE vector
   interface is included in this release mainly for testing purposes
   and as a preview of functionality to come.
