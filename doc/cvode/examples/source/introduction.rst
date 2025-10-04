..
   -----------------------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2025, Lawrence Livermore National Security,
   University of Maryland Baltimore County, and the SUNDIALS contributors.
   Copyright (c) 2013-2025, Lawrence Livermore National Security
   and Southern Methodist University.
   Copyright (c) 2002-2013, Lawrence Livermore National Security.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   -----------------------------------------------------------------------------

.. _CVODE.Examples.Intro:

Introduction
============

CVODE includes examples of many types, to illustrate the use of various
integrator options, linear solver options, and vector implementations.

With the exception of "demo"-type example files, the names of the examples are
generally of the form ``[slv][PbName]_[ls]_[prec]_[p]``, where:

- ``[slv]`` identifies the solver (in this case ``cv``).
- ``[PbName]`` identifies the problem.
- ``[ls]`` identifies the linear solver module used (for examples using
  fixed-point iteration for the nonlinear system solver, ``non`` specifies that
  no linear solver was used).
- ``[prec]`` indicates the CVODE preconditioner module used, ``bp`` for
  ``CVBANDPRE`` or ``bbd`` for ``CVBBDPRE`` (only if applicable, for examples
  using a Krylov linear solver);
- ``[p]`` indicates an example using the MPI parallel vector.

The following lists summarize all examples distributed with CVODE. In the
subsequent sections, we give detailed descriptions of some (but not all) of
these examples. We also give our output files for each of these examples, but
users should be cautioned that their results may differ slightly from
these. Differences in solution values may differ within the tolerances, and
differences in cumulative counters, such as numbers of steps or Newton
iterations, may differ from one machine environment to another by as much as 10%
to 20%.

.. note::

   The examples in the CVODE distribution are written in such a way as to
   compile and run for any combination of configuration options during the
   installation of SUNDIALS. As a consequence, they contain portions of code
   that will not be typically present in a user program. For example, programs
   make use of the variables ``SUNDIALS_EXTENDED_PRECISION`` and
   ``SUNDIALS_DOUBLE_PRECISION`` to test if the solver libraries were built in
   extended or double precision, and use the appropriate conversion specifiers
   in ``printf`` functions.

The final section of this report describes a set of tests done with the parallel
version of CVODE, using a problem based on the ``cvDiurnal_kry`` and
``cvDiurnal_kry_p`` examples.

Serial Examples
---------------

Supplied in the ``examples/cvode/serial`` directory are the following serial
examples:

- ``cvRoberts_dns`` solves a chemical kinetics problem consisting of three rate
  equations.

  This program solves the problem with BDF methods and a Newton iteration with
  the dense linear solver and a user-supplied Jacobian. It also uses the
  rootfinding feature of CVODE.

- ``cvRoberts_dns_constraints`` is the same as ``cvRoberts_dns`` but imposes the
  constraint :math:`u \geq 0.0` for all components.

- ``cvRoberts_dnsL`` is the same as ``cvRoberts_dns`` but uses the LAPACK
  dense linear solver.

- ``cvRoberts_dns_uw`` is the same as ``cvRoberts_dns`` but demonstrates the
  user-supplied error weight function feature of CVODE.

- ``cvRoberts_dns_negsol`` is the same as ``cvRoberts_dns`` but demonstrates the
  treatment of negative (unphysical) solution values through the RHS function
  return flag.

- ``cvRocket_dns`` is a simplified model of a rocket, demonstrating the
  preferred way to stop and restart at a root-defined discontinuity.

- ``cvRoberts_klu`` is the same as ``cvRoberts_dns`` but uses the KLU sparse
  direct linear solver.

- ``cvRoberts_block_klu`` solves multiple copies of the ``cvRoberts_dns``
  problem, using the KLU sparse direct linear solver.

- ``cvRoberts_sps`` is the same as ``cvRoberts_dns`` but uses the SuperLU_MT
  sparse direct linear solver.

- ``cvAdvDiff_bnd`` solves the semi-discrete form of an advection-diffusion
  equation in 2D. This program solves the problem with BDF methods and a Newton
  iteration with the banded linear solver and a user-supplied Jacobian.

- ``cvAdvDiff_bndL`` is the same as ``cvAdvDiff_bnd`` but uses the LAPACK
  banded linear solver.

- ``cvDiurnal_kry`` solves the semi-discrete form of a two-species diurnal
  kinetics advection-diffusion PDE system in 2D.

  The problem is solved with BDF methods and a Newton iteration with the GMRES
  linear solver. The block-diagonal part of the Newton matrix is used as a left
  preconditioner. A copy of the block-diagonal part of the Jacobian is saved and
  conditionally reused within the preconditioner setup routine.

- ``cvDiurnal_kry_bp`` solves the same problem as ``cvDiurnal_kry``, but with a
  banded preconditioner generated by difference quotients using the
  ``CVBANDPRE`` module. The problem is solved twice: with preconditioning on the
  left, then on the right.

- ``cvDirectDemo_ls`` is a demonstration program for CVODE with direct linear
  solvers. Two separate problems are solved using both Adams and BDF methods in
  combination with fixed-point and Newton iterations.

  The first problem is the Van der Pol oscillator for which the Newton iteration
  cases use the following types of Jacobian approximations: (1) dense,
  user-supplied, (2) dense, difference-quotient approximation, (3) diagonal,
  with difference-quotient approximation.

  The second problem is a linear ODE with a banded lower triangular matrix
  derived from a 2D advection PDE. In this case, the Newton iteration cases use
  the following types of Jacobian approximation: (1) banded, user-supplied, (2)
  banded, difference-quotient approximation, (3) diagonal, difference-quotient
  approximation.

- ``cvKrylovDemo_ls`` solves the same problem as ``cvDiurnal_kry``, with BDF
  methods, but with three Krylov linear solvers: GMRES, BiCGSTAB, and TFQMR.

- ``cvKrylovDemo_prec`` is a demonstration program using a Newton iteration with
  the GMRES linear solver.

  This program solves a stiff ODE system that arises from a system of partial
  differential equations. The PDE system is a six-species food web population
  model, with predator-prey interaction and diffusion on the unit square in two
  dimensions.

  The preconditioner matrix used is the product of two matrices: (1) a matrix,
  only defined implicitly, based on a fixed number of Gauss-Seidel iterations
  using the diffusion terms only; and (2) a block-diagonal matrix based on the
  partial derivatives of the interaction terms only, using block-grouping.

  Four different runs are made for this problem. The product preconditioner is
  applied on the left and on the right. In each case, both the modified and
  classical Gram-Schmidt options are tested.

- ``cvDisc_dns`` solves two simple problems, one with a discontinuity in the
  solution, and one with a discontinuity in the RHS function.

- ``cvAnalytic_mels`` solves a problem having a simple analytic solution, using
  a custom matrix-imbedded linear solver.

- ``cvParticle_dns`` solves a simple particle motion problem with a solution
  constraint and a user-supplied projection function onto the constraint.

- ``cvPendulum_dns`` solves a pendulum motion problem with and without
  constraints on the solution, and with various tolerances. In the runs with
  constraints, a user-supplied projection function is applied.

OpenMP Examples
---------------

Supplied in the ``examples/cvode/C_openmp`` directory is an example
``cvAdvDiff_bnd_omp``, which solves the same problem as ``cvAdvDiff_bnd`` but
with the OpenMP vector.

MPI Examples
------------

Supplied in the ``examples/cvode/parallel`` directory are the following four
MPI parallel examples:

- ``cvAdvDiff_non_p`` solves the semi-discrete form of a 1D advection-diffusion
  equation.

  This program solves the problem with the option for nonstiff systems,
  i.e. Adams methods with fixed-point iteration.

- ``cvAdvDiff_diag_p`` solves the same problem as ``cvAdvDiff_non_p``, with
  Adams methods, but with a Newton iteration and a diagonal linear solver.

- ``cvDiurnal_kry_p`` is a parallel implementation of ``cvDiurnal_kry``.

- ``cvDiurnal_kry_bbd_p`` solves the same problem as ``cvDiurnal_kry_p``, with
  BDF methods and a Newton iteration with the GMRES linear solver using a
  block-diagonal matrix with banded blocks as a preconditioner, generated by
  difference quotients, using the module ``CVBBDPRE``.

Supplied in the ``examples/cvode/C_mpimanyvector`` directory is an example
``cvDiurnal_kry_mpimanyvec``, which solves the same problem as
``cvDiurnal_kry_p``, but with the MPI ManyVector.

Supplied in the ``examples/cvode/parhyp`` directory is an example
``cvAdvDiff_non_ph``, which solves the same problem as ``cvAdvDiff_non_p`` but
with hypre MPI parallel vectors instead of SUNDIALS MPI parallel vectors.

Supplied in the ``examples/cvode/petsc`` directory are the following examples,
using the PETSc SNES nonlinear solver:

- ``cvAdvDiff_petsc`` solves the same problem as ``cvAdvDiff_non_p``.

- ``cv_petsc_ex7`` solves a problem based on PETSc TS ``ex7.c``, a nonlinear
  system derived from a time-dependent PDE in 2D.

C++ Examples
------------

The following is a list of directories and example names that are written in
C++. Each is based on an example listed earlier, except where noted.

- ``examples/cvode/superludist`` contains an example ``cvAdvDiff_sludist``.

- ``examples/cvode/CXX_serial`` contains two examples -- ``cv_heat2D`` and
  ``cv_kpr``. The latter solves a size-2 system, the Kvaerno-Prothero-Robinson
  test.

- ``examples/cvode/CXX_parallel`` contains an example ``cv_heat2D_p``.

- ``examples/cvode/CXX_parhyp`` contains two examples -- ``cv_heat2D_hypre_ls``
  and ``cv_heat2D_hypre_pfmg``.

GPU Examples
------------

CUDA Examples
^^^^^^^^^^^^^

Supplied in the ``examples/cvode/cuda`` directory are the following examples
using CUDA:

- ``cvAdvDiff_diag_cuda`` solves the same problem as ``cvAdvDiff_non_p``, but
  with the Diagonal linear solver.

- ``cvAdvDiff_kry_cuda`` solves the same problem as ``cvAdvDiff_non_p``, but
  with the GMRES linear solver.

- ``cvAdvDiff_kry_cuda_managed`` is the same as ``cvAdvDiff_kry_cuda`` but uses
  managed memory for the vector data.

HIP Examples
^^^^^^^^^^^^

The ``examples/cvode/hip`` directory contains two examples using HIP,
``cvAdvDiff_diag_hip`` and ``cvAdvDiff_kry_hip``.

SYCL Examples
^^^^^^^^^^^^^

The following examples utilize the SYCL abstraction layer.

- ``examples/cvode/CXX_sycl`` contains an example ``cvAdvDiff_kry_sycl``.

- ``examples/cvode/CXX_onemkl`` contains an example
  ``cvRoberts_blockdiag_onemkl``.

Kokkos Examples
^^^^^^^^^^^^^^^

In the ``examples/cvode/kokkos`` directory are two examples,
``cv_bruss_batched_kokkos`` and ``cv_bruss_batched_kokkos_2D``, using the Kokkos
performance portability layer.

RAJA Examples
^^^^^^^^^^^^^

In the ``examples/cvode/raja`` directory is the ``cvAdvDiff_kry_raja`` example
using the RAJA performance portability layer.

MAGMA Examples
^^^^^^^^^^^^^^

In the ``examples/cvode/magma`` directory is the example
``cv_bruss_batched_magma`` using the MAGMA batched direct linear solver with
CUDA or HIP.

Ginkgo Examples
^^^^^^^^^^^^^^^

In the ``examples/cvode/ginkgo`` directory are two examples,
``cv_heat2D_ginkgo`` and ``cv_kpr_ginkgo``, using linear solvers from the Ginkgo
linear solver. These examples may be run with serial, OpenMP, CUDA, HIP, or SYCL
backends.

Fortran serial Examples
-----------------------

Supplied in the ``examples/cvode/F2003_serial`` directory are the following
examples, all in using the Fortran interface modules:

- ``cv_analytic_fp_f2003`` solves the same problem as ``cvAnalytic_mels``, using
  the fixed-point nonlinear solver.

- ``cv_analytic_sys_dns_f2003`` solves a 3x3 system, also having an analytic
  solution, using the dense linear solver.

- ``cv_analytic_sys_dns_jac_f2003`` solves the same problem as
  ``cv_analytic_sys_dns_f2003``, but with a user-supplied Jacobian.

- ``cv_analytic_sys_klu_f2003`` solves the same problem but with the KLU linear
  solver.

- ``cv_brusselator_dns_f2003`` solves the Brusselator problem, a 3x3 nonlinear
  system, using the dense linear solver.

- ``cv_diurnal_kry_f2003`` solves the same problem as ``cv_Diurnal_kry``

- ``cv_diurnal_kry_bp_f2003`` solves the same problem as ``cv_Diurnal_kry_bp``

- ``cv_advdiff_bnd_f2003`` solves the same problem as ``cv_AdvDiff_bnd``.

- ``cv_roberts_dns_f2003`` solves the same problem as ``cv_roberts_dns``.

- ``cv_roberts_dnsL_f2003`` solves the same problem as ``cv_roberts_dns``, using
  the LAPACK dense solver.

- ``cv_roberts_dns_constraints_f2003`` solves the same problem as
  ``cv_roberts_dns_constraints``.

- ``cv_roberts_klu_f2003`` solves the same problem as ``cv_roberts_klu``.

Fortran MPI Examples
--------------------

Supplied in the ``examples/cvode/F2003_parallel`` directory are the following
MPI parallel examples all using the Fortran interface modules:

- ``cv_diag_non_p_f2003`` solves a simple diagonal nonstiff ODE system.

- ``cv_diag_kry_f2003`` solves a simple diagonal stiff ODE system using the
  SPGMR linear solver.

- ``cv_diag_kry_bbd_f2003`` solves the same problem using the BBD
  preconditioner.
