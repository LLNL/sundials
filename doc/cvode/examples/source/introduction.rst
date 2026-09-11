..
   Copyright (c) 2002-2026, Lawrence Livermore National Security and
   the SUNDIALS contributors. SPDX-License-Identifier: BSD-3-Clause

.. _CVODE.Examples.Intro:

Introduction
============

CVODE includes examples covering its integration methods, nonlinear and linear
solvers, preconditioners, and vector implementations. Except for demo programs,
names generally have the form ``[slv][PbName]_[ls]_[prec]_[p]``. Here ``cv``
identifies CVODE, ``ls`` identifies the linear solver (``non`` means fixed-point
iteration without one), ``bp`` and ``bbd`` identify CVBANDPRE and CVBBDPRE, and
``p`` denotes an MPI parallel vector.

Serial examples
---------------

The ``examples/cvode/serial`` directory contains:

- Robertson chemical kinetics variants using dense, LAPACK dense, KLU, and
  SuperLU_MT solvers, with examples of constraints, custom error weights,
  recoverable right-hand-side errors, and rootfinding.
- ``cvRocket_dns``, demonstrating stop and restart at a discontinuity defined
  by a root.
- ``cvAdvDiff_bnd`` and ``cvAdvDiff_bndL``, solving a 2D
  advection-diffusion equation with native and LAPACK band solvers.
- ``cvDiurnal_kry`` and ``cvDiurnal_kry_bp``, solving a two-species diurnal
  advection-diffusion problem with GMRES and user or CVBANDPRE preconditioning.
- Direct and Krylov demonstration programs exercising multiple integration,
  nonlinear-solver, linear-solver, Jacobian, preconditioning, and
  orthogonalization choices.
- Discontinuity, analytic-solution, particle-projection, and constrained
  pendulum examples.

Distributed and accelerator examples
------------------------------------

The ``parallel`` directory contains MPI advection-diffusion and diurnal
kinetics examples, including user and CVBBDPRE preconditioning. Related
implementations use OpenMP, MPIManyVector, hypre, CUDA, HIP, RAJA, SYCL,
oneMKL, PETSc, SuperLU_DIST, Ginkgo, Kokkos, and MAGMA. C++ examples include
serial and parallel heat equations and the Kvaerno--Prothero--Robinson test.

Fortran examples
----------------

The ``F2003_serial`` examples cover analytic systems, the Brusselator,
advection-diffusion, diurnal kinetics, and Robertson kinetics with several
linear solvers and options. ``F2003_parallel`` includes nonstiff and stiff
diagonal systems, with and without CVBBDPRE.

.. note::

   Example output may vary slightly between platforms. Solution differences
   should remain within tolerances, while cumulative counters may differ by
   10--20%. Conditional code supporting multiple SUNDIALS configurations may
   not be needed in a typical application.
