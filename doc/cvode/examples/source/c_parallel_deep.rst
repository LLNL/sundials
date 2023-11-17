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


.. _parallel_deep_c:

====================================================
Parallel C example problems -- Deep dive
====================================================



.. _deep_dive.cvAdvDiff_non_p:

A nonstiff example: cvAdvDiff_non_p [DD]
===============================================

Description -- step-by-step
----------------------------

This problem begins with a simple diffusion-advection equation
for :math:`u = u(t,x)`

.. math::
   :label: PDE1

   \frac{\partial u}{\partial t}=\frac{\partial ^{2}u}{\partial x^{2}}
   + 0.5\frac{\partial u}{\partial x}

for :math:`0 \leq t \leq 5, ~~ 0\leq x \leq 2`, and subject to homogeneous
Dirichlet boundary conditions and initial values given by

.. math::
   :label: BCIC1

   u(t,0) &= 0, \quad u(t,2) = 0, \\
   u(0,x) &= x(2-x)e^{2x}.

A system of :math:`M_X` ODEs is obtained by discretizing the :math:`x`-axis with :math:`M_X + 2`
grid points and replacing the first and second order spatial derivatives
with their *central difference approximations*. Since the value of :math:`u` is
constant at the two endpoints, the semi-discrete equations for those points
can be eliminated.  With :math:`u_i` as the approximation to :math:`u(t,x_i)`,
:math:`x_i = i (\Delta x)`, and :math:`\Delta x = \frac{2}{M_X + 1}`, the resulting system of
ODEs, :math:`\dot u = f(t,u)`, can now be written:

.. math::
   :label: cvAdvDiff_pode

   \dot{u}_i=\frac{u_{i+1}-2u_{i}+u_{i-1}}{(\Delta x)^{2}}
             + 0.5 \frac{u_{i+1}-u_{i-1}}{2(\Delta x)} ~.

This equation holds for :math:`i = 1, 2, \ldots, M_X`, with the understanding
that :math:`u_0 = u_{M_X + 1} = 0.`

In the *parallel processing environment*, we may think of the several
processors as being laid out on a straight line with each processor to
compute its contiguous subset of the solution vector.  Consequently
the computation of the *right hand side* of the Equation (:math:numref:`cvAdvDiff_pode`) requires
that each interior processor must pass the first component of its block of
the solution vector to its left-hand neighbor, acquire the last component of
that neighbor's block, pass the last component of its block of the solution
vector to its right-hand neighbor, and acquire the first component of that
neighbor's block. If the processor is the first (:math:`0^\text{th}`) or last processor,
then communication to the *left or right* (respectively) is not required.

This problem uses the *Adams-Moulton* (non-stiff) integration formula and *fixed-point*
nonlinear solver.  It is unrealistically simple, but serves to illustrate
use of the *parallel version* of **CVODE**.

The :literal:`cvAdvDiff_non_p.c` file begins with ``#include`` declarations for
various required header files, including lines for
:literal:`nvector_parallel` to access the parallel :literal:`N_Vector` type and related
*macros*, and for :literal:`mpi.h` to access **MPI** types and constants. Following
that are definitions of problem constants and a data block for communication
with the :math:`f` routine.  That block includes the *number of PEs*, the *index*
of the *local PE*, and the *MPI communicator*.

The :literal:`main` program begins with **MPI** calls to initialize **MPI** and to set
*multi-processor environment parameters* ``npes`` (*number of PEs*) and
``my_pe`` (*local PE index*).  The *local vector length* is set according
to ``npes`` and the *problem size* ``NEQ`` (which may or may not be
multiple of ``npes``).  The value ``my_base`` is the *base value for
computing global indices* (from :math:`1` to :math:`NEQ`) for the local vectors.
The solution vector :math:`u` is created with a call to :literal:`N_VNew_Parallel`
and loaded with a call to :literal:`SetIC`.  The calls to :literal:`CVodeCreate`,
:literal:`CVodeInit`, and :literal:`CVodeSStolerances` specify a **CVODE** solution with
the *nonstiff method* and *scalar tolerances*.  The call to :literal:`CVodeSetUserdata`
insures that the pointer ``data`` is passed to the :math:`f` routine whenever it
is called.

A heading is printed (if on processor :math:`0`).  In a loop over :math:`t_{out}` *values*,
:literal:`CVode` is called, and the return value checked for errors.  The
*max-norm* of the solution and the total number of time steps so far
are printed at each output point.  Finally, some statistical counters are
printed, memory is freed, and **MPI** is finalized.

The :literal:`SetIC` routine uses the last two arguments passed to it to compute
the set of global indices (:math:`my\_base + 1` to :math:`my\_base + my\_length`)
corresponding to the local part of the solution vector :math:`u`, and then to
load the corresponding initial values.  The :literal:`PrintFinalStats` routine
uses :literal:`CVodeGet***` calls to get various counters, and then prints these.
The counters are:

* ``nst`` (*number of steps*),
* ``nfe`` (*number of* :math:`f` *evaluations*),
* ``nni`` (*number of nonlinear iterations*),
* ``netf`` (*number of error test failures*),
* and ``ncfn`` (*number of nonlinear convergence failures*).

This routine is suitable for general use with **CVODE** applications to nonstiff problems.

The :literal:`f` function is an implementation of the Equation (:math:numref:`cvAdvDiff_pode`),
but preceded by communication operations appropriate for the parallel setting.
It copies the local vector :math:`u` into a larger array :math:`z`, shifted by :math:`1`
to allow for the storage of immediate neighbor components.  The first and last
components of :math:`u` are sent to neighboring processors with :literal:`MPI_Send` calls,
and the immediate neighbor solution values are received from the neighbor
processors with :literal:`MPI_Recv` calls, except that zero is loaded into :math:`z[0]`
or :math:`z[my\_length + 1]` instead if at the actual boundary.  Then the central
difference expressions are easily formed from the :math:`z` array, and loaded into
the data array of the :math:`\dot u =` ``udot`` vector.

The :literal:`cvAdvDiff_non_p.c` file includes a routine :literal:`check_flag` that checks the
return values from calls in :literal:`main`.  This routine was written to be used
by any *parallel* **SUNDIALS** application.

Program output
---------------

The output below is for :literal:`cvAdvDiff_non_p` with :math:`M_X = 10` and four processors.
Varying the number of processors will alter the output, only because
of roundoff-level differences in various vector operations.  The fairly
high value of ``ncfn`` indicates that this problem is on the borderline
of being *stiff*.

.. include:: ../../../../examples/cvode/parallel/cvAdvDiff_non_p.out
   :literal:


.. _deep_dive.cvDiurnal_kry_p:

A user preconditioner example: cvDiurnal_kry_p [DD]
===============================================================

Description -- step-by-step
----------------------------

As an example of using **CVODE** with the *Krylov linear solver*
**SUNLinSol_SPGMR**, **CVode** *linear solver* interface, and
the *parallel* **MPI NVECTOR_PARALLEL** module, we describe a test problem based on
the system of PDEs given above for the :literal:`cvDiurnal_kry` example.
As before, we discretize the PDE system with *central differencing*, to
obtain an ODE system :math:`\dot u = f(t,u)` representing (:math:numref:`cvDiurnalpde`).
But in this case, the discrete solution vector is distributed over
many processors.  Specifically, we may think of the processors as
being laid out in a rectangle, and each processor being assigned a
subgrid of size :math:`M_{X_{SUB}} \times M_{Y_{SUB}}` of the :math:`x-y` grid. If
there are :math:`NPE_X` processors in the :math:`x` direction and :math:`NPE_Y`
processors in the :math:`y` direction, then the overall grid size is
:math:`M_X \times M_Y` with :math:`M_X = NPE_X \times M_{X_{SUB}}` and
:math:`M_Y = NPE_Y \times M_{Y_{SUB}}`, and the size of the ODE system is
:math:`2 \cdot M_X \cdot M_Y`.

To compute :math:`f` in this setting, the processors pass and receive
information as follows.  The solution components for the bottom row of
grid points in the current processor are passed to the processor below
it and the solution for the top row of grid points is received from
the processor below the current processor. The solution for the top
row of grid points for the current processor is sent to the processor
above the current processor, while the solution for the bottom row of
grid points is received from that processor by the current
processor. Similarly the solution for the first column of grid points
is sent from the current processor to the processor to its left and
the last column of grid points is received from that processor by the
current processor. The *communication* for the solution at the right
edge of the processor is similar. If this is the last processor in a
particular direction, then message passing and receiving are bypassed
for that direction.

This code is intended to provide a more realistic example than that in
:literal:`cvAdvDiff_non_p`, and to provide a template for a stiff ODE system
arising from a PDE system. The solution method is *BDF* with *Newton*
iteration and **SPGMR**. The *left preconditioner* is the *block-diagonal*
part of the *Newton matrix*, with :math:`2 \times 2` *blocks*, and the
corresponding *diagonal blocks* of the *Jacobian* are saved each time the
*preconditioner* is generated, for re-use later under certain conditions.

The organization of the :literal:`cvDiurnal_kry_p` program deserves some comments. The
*right-hand side* routine :math:`f` calls two other routines:

* :literal:`ucomm`, which carries out inter-processor communication;
* and :literal:`fcalc`, which operates on local data only and contains the actual
  calculation of :math:`f(t,u)`.

The :literal:`ucomm` function in turn calls three routines which do, respectively,
non-blocking receive operations, blocking send operations, and
receive-waiting. All three use **MPI**, and transmit data from the local :math:`u`
vector into a local working array :math:`u_{ext}`, an extended copy of :math:`u`.
The :literal:`fcalc` function copies :math:`u` into :math:`u_{ext}`, so that the
calculation of :math:`f(t,u)` can be done conveniently by operations on
:math:`u_{ext}` only.  Most other features of :literal:`cvDiurnal_kry_p.c` are the same as
in :literal:`cvDiurnal_kry.c`, except for extra logic involved with distributed
vectors.

Program output
---------------

The following is a sample output from :literal:`cvDiurnal_kry_p`, for four processors
(in a :math:`2 \times 2` array) with a :math:`5 \times 5` subgrid on each.
The output will vary slightly if the number of processors is changed.

.. include:: ../../../../examples/cvode/parallel/cvDiurnal_kry_p.out
   :literal:


.. _deep_dive.cvDiurnal_kry_bbd_p:

A CVBBDPRE preconditioner example: cvDiurnal_kry_bbd_p [DD]
========================================================================

Description -- step-by-step
----------------------------

In this example, :literal:`cvDiurnal_kry_bbd_p`, we solve the same
problem as in :literal:`cvDiurnal_kry_p`
above, but instead of supplying the *preconditioner*, we use the **CVBBDPRE** module,
which generates and uses a *band-block-diagonal preconditioner*.  The
*half-bandwidths* of the *Jacobian block* on each processor are both equal to
:math:`2 \cdot M_{X_{SUB}}`, and that is the value supplied as ``mudq`` and ``mldq``
in the call to :literal:`CVBBDPrecInit`.  But in order to reduce storage and computation
costs for *preconditioning*, we supply the values :math:`mu_{keep} = ml_{keep} = 2`
(:math:`= N_{VARS}`) as the *half-bandwidths* of the retained *band matrix blocks*.
This means that the *Jacobian elements* are computed with a *difference quotient*
scheme using the true bandwidth of the block, but only a narrow band matrix
(bandwidth :math:`5`) is kept as the *preconditioner*.

As in :literal:`cvDiurnal_kry_p.c`, the :math:`f` routine in :literal:`cvDiurnal_kry_bbd_p.c` simply calls a
*communication* routine, :literal:`fucomm`, and then a strictly computational routine,
:literal:`flocal`.  However, the call to :literal:`CVBBDPrecInit` specifies the pair of
routines to be called as :literal:`ucomm` and :literal:`flocal`, where :literal:`ucomm` is ``NULL``.
This is because each call by the solver to :literal:`ucomm` is
preceded by a call to :math:`f` with the same :math:`(t,u)` arguments, and therefore the
*communication* needed for :literal:`flocal` in the solver's calls to it have already been
done.

In :literal:`cvDiurnal_kry_bbd_p.c`, the problem is *solved twice* --- first with *preconditioning
on the left*, and then *on the right*.  Thus prior to the second solution, calls
are made to reset the initial values (:literal:`SetInitialProfiles`), the main solver
memory (:literal:`CVodeReInit`), the **CVBBDPRE** memory (:literal:`CVBBDPrecReInit`),
as well as the *preconditioner* type (:literal:`SUNLinSol_SPGMRSetPrecType`).

Program output
---------------

Sample output from :literal:`cvDiurnal_kry_bbd_p` follows, again using :math:`5 \times 5` *subgrids*
on a :math:`2 \times 2` processor grid.  The performance of the *preconditioner*,
as measured by the number of *Krylov* iterations per *Newton* iteration,
``nli/nni``, is very close to that of :literal:`cvDiurnal_kry_p` when *preconditioning* is *on
the left*, but slightly poorer when it is *on the right*.

.. include:: ../../../../examples/cvode/parallel/cvDiurnal_kry_bbd_p.out
   :literal:
