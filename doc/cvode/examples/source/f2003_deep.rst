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


.. _deep_f2003:

===================================================
Fortran 2003 example problems -- Deep dive
===================================================

The **Fortran** example problem programs supplied with the **CVODE**
package are all written in standard **Fortran 2003** and use *double precision
arithmetic*.  Before running any of these examples, the user should
make sure that the **Fortran** *data types* for *real* and *integer variables*
appropriately match the **C** *types*.  See `SUNDIALS Fortran
<https://sundials.readthedocs.io/en/latest/sundials/Fortran_link.html>`_ in the
**SUNDIALS User Document** for details.



.. _deep_dive.cv_diurnal_kry:

A serial example: cv_diurnal_kry [DD]
==============================================

Description -- step-by-step
----------------------------

The :literal:`cv_diurnal_kry` example is a **Fortran** equivalent of the :literal:`cvDiurnal_kry` problem.
(In fact, it was derived from an *earlier* **Fortran** *example program* for **VODPK**.)

The :literal:`main` program begins with several declarations, creating the **SUNDIALS** *context*
pointer via :literal:`FSUNContext_Create`, creating the **SUNDIALS NVector** pointers
via :literal:`FN_VNew_Serial`, initializing the **SUNDIALS NVector** external pointer
(to assist with *data modification*) via :literal:`FN_VGetArrayPointer`, and dual ``do``
loops to initialize our :math:`u_{vec}` via embedded :math:`x`, :math:`y`, :math:`c_x`,
and :math:`c_y`.  With this we have successfully initialized the **NVECTOR_SERIAL** module via **SUNDIALS**.
:literal:`main` then calls :literal:`FCVodeCreate` and :literal:`FCVodeInit` as per usual, thereby
initializing the **FCVODE** interface with ``CV_BDF``.  :literal:`FCVodeSetMaxNumSteps` and
:literal:`FCVodeSStolerances` are then called to set the *maximum number of steps* per solve to be
:math:`1000` and initialize *scalar tolerances*.  :literal:`FSUNLinSol_SPGMR`, :literal:`FSUNLinSol_SPGMRSetGSType`,
and :literal:`FCVodeSetLinearSolver` are then called to initialize the **SUNLinSol_SPGMR** module with specific
``iGStype`` characteristics to the pointer ``sunls``.  :literal:`main` then calls :literal:`FCVodeSetPreconditioner`
to specify *user-supplied preconditioner setup* and *solve* routines.  Then, via a ``do`` loop,
:literal:`main` calls :literal:`FCVode` in a loop over :math:`t_{out}` values, meanwhile printing
selected *solution values* and *performance data*, getting ``lnst`` from :literal:`FCVodeGetNumSteps`
and ``lh`` from :literal:`FCVodeGetCurrentStep`.  At the end, it prints a number of performance counters,
via the user sub-routine :literal:`CVodeStats`, such as:

* ``nsteps`` (total internal steps taken),
* ``nfe`` (total *right-hand side* function calls),
* ``npe`` (total number of *preconditioner* evaluations),
* ``nps`` (total number of *preconditioner* solves),
* ``netfails`` (number of *error test fails*),
* ``nniters`` (number of *nonlinear iterations*),
* ``nliters`` (number of *linear iterations*),
* ``avdim`` :math:`\bigg(` average *Krylov space dimension* :math:`= \frac{\texttt{nliters}}{\texttt{nniters}} \bigg)`,
* ``ncf`` (number of *nonlinear convergence failures*),
* ``ncfl`` (number of *linear convergence failures*),
* ``lenrw`` and ``leniw`` (main solver *real/int workspace size*)
* ``lenrwls`` and ``leniwls`` (linear solver *real/int workspace size*),

frees memory with calls to :literal:`FCVFree`, :literal:`FN_VDestroy`, :literal:`FSUNLinSolFree`,
and :literal:`FSUNContext_Free`.

In :literal:`cv_diurnal_kry.f90`, or more specifically, in the :literal:`diurnal_mod` module,
we have that the :literal:`RhsFn` function is a straghtforward implementation of
the discretized form of Equations (:math:numref:`cvDiurnalpde`).  In :literal:`PreSet`, the *block-diagonal*
part of the *Jacobian*, :math:`J_{bd}`, is computed (and copied to :math:`p`) if :math:`jok = 0`
via the subroutine :literal:`Prec_Jac`, but is simply copied from :math:`bd` to :math:`p` if
:math:`jok = 1`.  In both cases, the *preconditioner matrix* :math:`p` is formed from 
:math:`J_{bd}` and its :math:`2 \times 2` blocks are *'LU-factored'* via the subroutine
:literal:`Prec_LU`, which serves as a sort of *block-by-block inversion* of :math:`p_{2 \times 2} - I_2`.
In :literal:`PreSolve`, the solution of a linear system :math:`Px = z` is solved by doing *backsolve
operations* on the :math:`2 \times 2` blocks.  This is done by, first, copying :math:`r` into :math:`z`.
Then, the subordinate routine :literal:`Prec_Sol` is used to perform *matrix-vector multiplication* to
finally perform our *backsolve*.  The remainder of :literal:`cv_diurnal_kry.f90` previously consisted of, now
defunct per our problem's functional requisites, routines from **LINPACK** and the **BLAS** needed for
*matrix* and *vector* operations.

Program output
---------------

The following is sample output from :literal:`cv_diurnal_kry`, using a :math:`10 \times 10` mesh.
The performance of **FCVODE** here is quite similar to that of **CVODE** on the :literal:`cvDiurnal_kry`
problem, as expected.

.. include:: ../../../../examples/cvode/F2003_serial/cv_diurnal_kry_f2003.out
   :literal:


.. _deep_dive.cv_diag_kry_bbd_p:

A parallel example: cv_diag_kry_bbd_p [DD]
===================================================

Description -- step-by-step
----------------------------

This example, :literal:`cv_diag_kry_bbd_p`, uses a simple diagonal ODE system to illustrate
the use of **FCVODE** in a *parallel* setting.  The system is then

.. math::
   :label: diagode

   \dot y_i = -\alpha ~i~ y_i ~~~ (i = 1, \ldots, N)

on the time interval :math:`0 \leq t \leq 1`.  In this case, we use :math:`\alpha = 10`
and :math:`N = 10 * nprocs`, where ``nprocs`` is the *number of processors* specified at run time.
The *linear solver* to be used is **SPGMR** with the **CVBBDPRE** (*band-block-diagonal*) *preconditioner*.
Since the system *Jacobian* is *diagonal*, the *half-bandwidths* specified are all *zero*.  The problem
is solved twice --- with *preconditioning on the left*, then *on the right*.

The source file for this problem begins with **MPI** calls to *initialize* **MPI** and to get the *number
of processors* and *local processor index*, or ``nprocs`` and ``nlocal`` respectively.  This is done via
the **MPI** subroutines :literal:`MPI_Init`, :literal:`MPI_Comm_size`, and :literal:`MPI_Comm_rank`.
:literal:`FSUNContext_Create`, :literal:`FCVodeCreate`, and :literal:`FCVodeInit` are called as usual.
:literal:`FN_VNew_Parallel` is then called to initialize the **MPI-parallel NVector** module so that we
can use :literal:`FN_VGetArrayPointer` as per usual, while :literal:`FSUNLinSol_SPGMR` :literal:`FCVodeSetLinearSolver`,
and :literal:`FSUNLinSol_SPGMRSetGSType` are called as per usual to initialize the **SPGMR SUNLinSol** module
and attach the *linear solver* to **FCVODE**.
Following the call to :literal:`FCVodeSStolerances`, the *preconditioner* is attached to **FCVODE** with calls to
:literal:`FCVBBDPrecInit`.  Then, via a ``do`` loop over ``iPreType``, we have several ``if`` statements, specific
to ``iPreType``.  In the first one, we consider :code:`iPreType == 2`, calling :literal:`FCVodeReInit`,
:literal:`FCVBBDPrecReInit`, and :literal:`FSUNLinSol_SPGMRSetPrecType` to *re-initialize* with *preconditioning
on the right*, as stated in the *print* statement.  This is all that is needed to prepare for the second run of our
Diagonal PDE.  Meanwhile, in the second one, while :code:`iPreType == 1`, we only *print* out the noted fact that
we now have *preconditioning on the left*.  After this, in a loop over :math:`t_{out}` values, it calls
:literal:`FCVode` and *prints* the *step* and :math:`f` *evaluation counters*. After that, it *computes* and *prints*
the *maximum global error* via :literal:`MPI_Reduce`. Afterwards, it *computes* and *prints* all the relevant performance
counters, such as:

* ``nst`` (internal *solver steps*),
* ``nfe`` (total *right-hand side* evaluations),
* ``npre`` (total *preconditioner* setups),
* ``npsol`` (total *preconditioner* solves),
* ``nni`` (total *nonlinear iterations*),
* ``nli`` (total *linear iterations*),
* ``avdim`` :math:`\bigg(` average *Krylov subspace dimension* :math:`= \frac{\texttt{nliters}}{\texttt{nniters}} \bigg)`,
* ``ncfn`` (total *convergence failures* -- *nonlinear*),
* ``ncfl`` (total *convergence failures* -- *linear*),
* ``netf`` (total number of *error test failures*),
* ``lenrw`` and ``leniw`` (main solver *real/int workspace sizes*),
* ``lenrwls`` and ``leniwls`` (linear solver *real/int workspace sizes*),
* ``lenrwbbd`` and ``leniwbbd`` (*BBD preconditioner real/int workspace sizes*),
* ``ngebbd`` (total number of :math:`g` *evaluations* by *BBD preconditioner*).

Finally, it *frees memory* and terminates **MPI** via the subroutines :literal:`FCVodeFree`, :literal:`FN_VDestroy`,
:literal:`MPI_Barrier`, and :literal:`MPI_Finalize`.  Notice that in the :literal:`firhs` function, outside of
the main :literal:`driver` routine, the local processor index ``myid`` and the local vector size ``nlocal``
are used to form the *module-specific index values* needed to evaluate the *right-hand side* of our Equation (:math:numref:`diagode`). 

Program output
---------------

The following is a sample output from :literal:`cv_diag_kry_bbd_p`, with `nprocs = 4`.
As expected, the performance is identical for *left* vs *right preconditioning*.

.. include:: ../../../../examples/cvode/F2003_parallel/cv_diag_kry_bbd_f2003.out
   :literal:

