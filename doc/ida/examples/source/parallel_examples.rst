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

.. _IDA.Examples.Parallel:

Parallel example problems
=========================

.. _IDA.Examples.Parallel.Heat2D:

A user-preconditioner example: idaHeat2D_kry_p
----------------------------------------------

This MPI example solves the same heat equation as
:ref:`IDA.Examples.Serial.Heat2D` using the
:ref:`MPI parallel vector <NVectors.NVParallel>`, GMRES, and a user-defined
preconditioner.

Processes form an ``NPEX`` by ``NPEY`` grid, each owning an ``MXSUB`` by
``MYSUB`` submesh. Residual evaluation exchanges internal-boundary values with
neighboring processes using blocking sends, nonblocking receives, and receive
waits. The received ghost values and local solution are assembled in ``uext``;
``reslocal`` then evaluates diffusion terms and represents the zero Dirichlet
boundary values as algebraic equations.

``PsetupHeat`` and ``PsolveHeat`` implement a diagonal Jacobian
preconditioner, requiring only local calculations. The main program initializes
MPI, creates distributed vectors for states, constraints, variable IDs, and the
preconditioner, imposes nonnegativity, and excludes algebraic components from
the error test with :c:func:`IDASetSuppressAlg`.

.. literalinclude:: ../../../../examples/ida/parallel/idaHeat2D_kry_p.out
   :language: none

.. _IDA.Examples.Parallel.FoodWeb:

An IDABBDPRE example: idaFoodWeb_kry_bbd_p
-------------------------------------------

This program solves the food-web problem from
:ref:`IDA.Examples.Serial.FoodWeb` in parallel with GMRES and
:ref:`IDABBDPRE <IDA.Usage.CC.precond.idabbdpre>`. Each process owns a submesh.
``rescomm`` exchanges ghost data, and ``reslocal`` evaluates the local residual
using an extended array. Homogeneous Neumann conditions are imposed by copying
the first interior mesh line into ghost cells.

The Jacobian block has true half-bandwidth
``NUM_SPECIES * MXSUB``, supplied as ``mudq`` and ``mldq``. To reduce storage
and factorization costs, only half-bandwidths ``mukeep = mlkeep = 2`` are
retained. ``reslocal`` is passed as ``Gres`` to :c:func:`IDABBDPrecInit`; its
``Gcomm`` argument is null because communication is already performed by the
full residual routine. The program calls :c:func:`IDACalcIC` before integrating.

.. literalinclude:: ../../../../examples/ida/parallel/idaFoodWeb_kry_bbd_p.out
   :language: none
