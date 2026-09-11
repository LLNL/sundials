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

.. _IDA.Examples.PETSc:

PETSc example problems
======================

A PETSc vector example: idaHeat2D_petsc_spgmr
---------------------------------------------

This example solves the same problem as :ref:`IDA.Examples.Parallel.Heat2D`,
but uses the PETSc vector instead of the native SUNDIALS MPI vector. The
numerical output is identical.

PETSc is initialized before use and finalized at shutdown. A two-dimensional
distributed array (DMDA) defines a five-point star stencil, Dirichlet
boundaries, and the mesh partition. ``DMCreateGlobalVector`` creates the PETSc
solution vector, and ``N_VMake_Petsc`` wraps it as an ``N_Vector``. The wrapper
does not own the PETSc vector, so the application destroys the underlying
object after destroying the wrapper. Other vectors are cloned from this
template to preserve its partition and data mapping.

User functions recover the PETSc vector with ``N_VGetVector_Petsc``. The main
PETSc operations used by the residual are:

- ``DMGetLocalVector`` creates local storage including ghost values.
- ``DMGlobalToLocalBegin`` and ``DMGlobalToLocalEnd`` update local and ghost
  data from the global vector.
- ``DMDAVecGetArray`` and ``DMDAVecGetArrayRead`` expose arrays indexed by
  global mesh coordinates.
- ``DMDAGetCorners`` returns the owned region's lower corner and dimensions.
- ``DMDAVecRestoreArray``, ``DMDAVecRestoreArrayRead``, and
  ``DMRestoreLocalVector`` restore borrowed PETSc objects.

The DMDA operations replace the explicit communication and local-indexing
helpers in ``idaHeat2D_kry_p``, separating the model implementation from the
parallel decomposition.

.. literalinclude:: ../../../../examples/ida/petsc/idaHeat2D_petsc_spgmr.out
   :language: none
