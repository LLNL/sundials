..
   Copyright (c) 2002-2026, Lawrence Livermore National Security and
   the SUNDIALS contributors. SPDX-License-Identifier: BSD-3-Clause

.. _CVODES.Examples.Intro:

Introduction
============

CVODES includes serial and parallel examples of IVP integration, forward
sensitivity analysis (FSA), and adjoint sensitivity analysis (ASA), plus
OpenMP and Fortran examples. Names generally have the form
``[slv][PbName]_[SA]_[ls]_[prec]_[p]``. ``FSA`` denotes forward sensitivity,
``ASAi`` adjoint sensitivity for an integral output, ``ASAp`` adjoint
sensitivity for a pointwise output, ``bp`` CVBANDPRE, ``bbd`` CVBBDPRE, and
``p`` an MPI parallel vector.

CVODES contains counterparts of the CVODE IVP examples described in
:ref:`CVODE.Examples`. Its sensitivity examples include:

- Robertson kinetics FSA and ASA variants with dense, KLU, and SuperLU_MT
  solvers, constraints, custom error weights, and sensitivity switching.
- Serial and MPI advection-diffusion FSA examples for diffusion and advection
  parameters.
- Serial and MPI diurnal kinetics FSA examples with GMRES and user
  preconditioning.
- Adjoint examples for advection-diffusion, food-web, and atmospheric
  dispersion models, including CVBBDPRE and forward-over-adjoint Hessian-vector
  products.

The ``C_openmp`` directory provides an OpenMP IVP example.
``F2003_serial`` contains an analytic fixed-point example and a Fortran version
of the advection-diffusion FSA example.

.. note::

   FSA programs accept command-line choices for enabling sensitivities,
   selecting simultaneous or staggered methods, and including sensitivities in
   error control. Output and cumulative counters may vary slightly by platform.
