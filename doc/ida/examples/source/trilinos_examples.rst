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

.. _IDA.Examples.Trilinos:

Trilinos example problems
=========================

A shared-memory example: idaHeat2D_kry_tpetra
---------------------------------------------

This C++ example solves the problem from :ref:`IDA.Examples.Serial.Heat2D` with
the Tpetra vector :cite:p:`hoemmen2015tpetra` from Trilinos
:cite:p:`Trilinos-Overview`. Tpetra uses Kokkos :cite:p:`edwards2014kokkos` for
on-node parallelism.

``Tpetra::ScopeGuard`` initializes and finalizes MPI and Tpetra. A Tpetra map
defines zero-based global-to-local indexing, and a Tpetra vector is wrapped
with ``N_VMake_Trilinos``. All other ``N_Vector`` objects are cloned from this
template. The example runs on one MPI rank whether Trilinos was built with or
without MPI. Its residual and preconditioner setup use Kokkos kernels on the
configured default execution space.

.. literalinclude:: ../../../../examples/ida/trilinos/idaHeat2D_kry_tpetra.out
   :language: none

.. _IDA.Examples.Trilinos.MPIPlusX:

An MPI+X example: idaHeat2D_kry_p_tpetra
----------------------------------------

This example parallels :ref:`IDA.Examples.Parallel.Heat2D` using a Tpetra
vector. It requires four MPI ranks and uses the Kokkos default execution space
for on-node parallelism.

Kokkos one-dimensional views serve as MPI buffers for the four subgrid
boundaries. Each buffer has a host mirror. Before MPI calls,
``Kokkos::deep_copy`` updates that mirror; when the buffer already resides in
host memory, the operation introduces no data transfer.

.. literalinclude:: ../../../../examples/ida/trilinos/idaHeat2D_kry_p_tpetra.out
   :language: none
