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

.. _IDAS.Examples.Intro:

Introduction
============

IDAS includes serial and parallel examples of initial value problem (IVP)
integration, forward sensitivity analysis (FSA), and adjoint sensitivity
analysis (ASA), as well as OpenMP and Fortran examples.

With the exception of "demo"-type examples, names generally have the form
``[slv][PbName]_[SA]_[ls]_[prec]_[p]``, where:

- ``[slv]`` identifies the solver (``idas``).
- ``[PbName]`` identifies the problem.
- ``[SA]`` is ``FSA`` for forward sensitivity, ``ASAi`` for adjoint
  sensitivity with an integral output, or ``ASAp`` for adjoint sensitivity
  with a pointwise output.
- ``[ls]`` identifies the linear solver.
- ``[prec]`` identifies a preconditioner (``bbd`` denotes
  :ref:`IDABBDPRE <IDAS.Usage.precond.idabbdpre>`).
- ``[p]`` denotes an example using the MPI parallel vector.

The following lists summarize the examples distributed with IDAS. Detailed
descriptions of selected sensitivity examples follow. The IVP examples are
also described in the :ref:`IDA examples documentation <IDA.Examples>`.

.. note::

   Example output can differ slightly between machines. Solution values may
   differ within the requested tolerances, while cumulative counters such as
   step and nonlinear iteration counts may differ by 10--20%.

Serial examples
---------------

The ``examples/idas/serial`` directory contains examples using the
:ref:`serial vector <NVectors.NVSerial>`.

IVP integration
^^^^^^^^^^^^^^^

- ``idasRoberts_dns`` solves the Robertson chemical kinetics problem
  :cite:p:`Rob:66`, including rootfinding, with the dense linear solver and a
  user-supplied Jacobian.
- ``idasRoberts_klu`` and ``idasRoberts_sps`` solve the same problem with KLU
  and SuperLU_MT, respectively.
- ``idasAkzoNob_dns`` solves the index-one Akzo-Nobel chemical kinetics DAEs
  with the dense linear solver.
- ``idasHeat2D_bnd`` solves a 2D heat equation discretized as a DAE, uses the
  band linear solver, calls ``IDACalcIC`` to correct boundary values, and
  imposes positivity constraints.
- ``idasHeat2D_kry`` solves the same problem with GMRES and a diagonal
  preconditioner.
- ``idasFoodWeb_bnd`` solves a 2D predator-prey reaction-diffusion problem
  with the band linear solver.
- ``idasSlCrank_dns`` solves index-two DAEs for a planar slider-crank mechanism
  and computes time-averaged kinetic energy as a quadrature.
- ``idasKrylovDemo_ls`` solves the heat problem with GMRES, BiCGSTAB, and
  TFQMR, using a diagonal preconditioner.
- ``idasAnalytic_mels`` solves a problem with a known solution using a custom
  matrix-embedded linear solver.

Forward sensitivity analysis
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

- ``idasRoberts_FSA_dns`` computes sensitivities of the Robertson problem with
  respect to its three reaction-rate constants using a dense linear solver and
  a user-supplied Jacobian.
- ``idasRoberts_FSA_klu`` and ``idasRoberts_FSA_sps`` solve the same problem
  with KLU and SuperLU_MT, respectively.
- ``idasSlCrank_FSA_dns`` computes sensitivities of the slider-crank solution
  and cumulative kinetic energy with respect to its spring and damping
  constants.

Adjoint sensitivity analysis
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

- ``idasRoberts_ASAi_dns`` computes gradients of a Robertson-problem
  functional with respect to its reaction-rate constants using dense linear
  solvers and user-supplied Jacobians for both forward and backward problems.
- ``idasRoberts_ASAi_klu`` and ``idasRoberts_ASAi_sps`` solve the same problem
  with KLU and SuperLU_MT, respectively.
- ``idasAkzoNob_ASAi_dns`` computes gradients of an integrated species
  concentration with respect to the initial conditions.
- ``idasHessian_ASA_FSA`` demonstrates a forward-over-adjoint computation of
  Hessian-vector products.

MPI examples
------------

The ``examples/idas/parallel`` directory contains examples using the
:ref:`MPI parallel vector <NVectors.NVParallel>`.

- ``idasHeat2D_kry_p`` solves the 2D heat problem with GMRES and a
  user-supplied diagonal preconditioner.
- ``idasHeat2D_kry_bbd_p`` uses GMRES and IDABBDPRE for the same problem.
- ``idasFoodWeb_kry_p`` solves the food-web problem with GMRES and a
  user-supplied block-diagonal preconditioner.
- ``idasFoodWeb_kry_bbd_p`` solves the same problem using IDABBDPRE.
- ``idasBruss_kry_bbd_p`` solves the two-species Brusselator PDE using GMRES
  and IDABBDPRE.
- ``idasBruss_FSA_kry_bbd_p`` computes sensitivities of the Brusselator with
  respect to two parameters and the gradient of a final-time spatial integral.
- ``idasHeat2D_FSA_kry_bbd_p`` computes sensitivities of the heat problem with
  respect to two PDE coefficients.
- ``idasBruss_ASAp_kry_bbd_p`` computes the Brusselator output gradient using
  adjoint sensitivity analysis.

OpenMP examples
---------------

The ``examples/idas/C_openmp`` directory contains ``idasFoodWeb_bnd_omp`` and
``idasFoodWeb_kry_omp``. They solve the food-web problem using the
:ref:`OpenMP vector <NVectors.OpenMP>`.

Fortran examples
----------------

The ``examples/idas/F2003_serial`` directory contains
``idasHeat2D_kry_f2003`` and ``idasAkzoNob_ASAi_dns_f2003``, Fortran versions
of the corresponding examples above.

.. note::

   The examples are designed to compile for different SUNDIALS configuration
   options. Consequently, they contain conditional code that a typical user
   program may not need. Forward sensitivity examples also accept command-line
   arguments selecting whether sensitivities are computed, the sensitivity
   method, and the error-control strategy.
