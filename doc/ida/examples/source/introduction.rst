..
   -----------------------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2025-2026, Lawrence Livermore National Security,
   University of Maryland Baltimore County, and the SUNDIALS contributors.
   Copyright (c) 2013-2025, Lawrence Livermore National Security
   and Southern Methodist University.
   Copyright (c) 2002-2013, Lawrence Livermore National Security.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   -----------------------------------------------------------------------------

.. _IDA.Examples.Intro:

Introduction
============

IDA includes examples illustrating different integration and linear solver
options and several vector implementations. Except for "demo" programs,
example names generally have the form ``[slv][PbName]_[ls]_[prec]_[p]``:

- ``[slv]`` identifies the solver (``ida``).
- ``[PbName]`` identifies the problem.
- ``[ls]`` identifies the linear solver.
- ``[prec]`` identifies a preconditioner; ``bbd`` denotes
  :ref:`IDABBDPRE <IDA.Usage.CC.precond.idabbdpre>`.
- ``[p]`` denotes an example using the MPI parallel vector.

Serial examples
---------------

The ``examples/ida/serial`` directory contains examples using the
:ref:`serial vector <NVectors.NVSerial>`:

- ``idaRoberts_dns`` solves the Robertson chemical kinetics problem
  :cite:p:`Rob:66`, including rootfinding, with a dense linear solver and
  user-supplied Jacobian. ``idaRoberts_klu`` and ``idaRoberts_sps`` use KLU
  and SuperLU_MT instead.
- ``idaSlCrank_dns`` solves index-two DAEs for a planar slider-crank mechanism
  obtained through a stabilized reduction of an index-three formulation.
- ``idaHeat2D_bnd`` solves a 2D heat equation with a band solver, uses
  :c:func:`IDACalcIC` to correct boundary values, and imposes positivity.
  ``idaHeat2D_kry`` uses GMRES and a diagonal preconditioner, while
  ``idaHeat2D_klu`` uses KLU.
- ``idaFoodWeb_bnd`` solves a 2D predator-prey reaction-diffusion problem with
  a band solver. ``idaFoodWeb_kry`` uses GMRES and a user preconditioner.
- ``idaKrylovDemo_ls`` demonstrates GMRES, BiCGSTAB, and TFQMR on the heat
  problem.
- ``idaAnalytic_mels`` uses a custom matrix-embedded linear solver on a problem
  with a known solution.

MPI and accelerator examples
----------------------------

The ``examples/ida/parallel`` directory contains MPI versions of the heat and
food-web problems, including ``idaHeat2D_kry_p``,
``idaHeat2D_kry_bbd_p``, ``idaFoodWeb_kry_p``, and
``idaFoodWeb_kry_bbd_p``. The ``C_openmp``, ``cuda``, and ``mpicuda``
directories provide OpenMP, CUDA, and MPI+CUDA variants.

The ``petsc`` directory contains ``idaHeat2D_petsc_spgmr`` and
``idaHeat2D_petsc_snes``. The C++ ``raja``, ``mpiraja``, and ``trilinos``
directories provide corresponding RAJA, MPI+RAJA, and Tpetra variants.

Fortran examples
----------------

The Fortran interface examples include ``idaRoberts_dns_f2003`` and
``idaHeat2D_kry_f2003`` in ``F2003_serial``,
``idaHeat2D_kry_bbd_f2003`` in ``F2003_parallel``, and
``idaHeat2D_kry_omp_f2003`` in ``F2003_openmp``.

.. note::

   Example output can differ slightly between machines. Solution values may
   differ within the requested tolerances, while cumulative counters such as
   step and nonlinear iteration counts may differ by 10--20%.

   The examples support multiple SUNDIALS configuration options and therefore
   contain conditional code that a typical application may not require.
