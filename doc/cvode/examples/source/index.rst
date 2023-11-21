..
   Programmer(s): Daniel R. Reynolds @ SMU
                  modified by Daniel M. Margolis @ SMU
   ----------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2002-2023, Lawrence Livermore National Security
   and Southern Methodist University.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   ----------------------------------------------------------------

.. CVODE_example documentation master file, created by
   sphinx-quickstart on Wed Jun 14 02:10:03 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

==============================================
Example Programs for CVODE
==============================================

Introduction
=============================


This report is intended to serve as a companion document to the User
Documentation of CVODE :cite:p:`cvode_ug`.  It provides details, with
listings, on many of the example programs supplied with the CVODE
distribution package.

This document provides detailed descriptions of CVODE examples of the
following types: serial C examples, parallel C examples, an OpenMP
example, and a *hypre* example.  Additionally, it provides high level
descriptions of all remaining CVODE examples.

With the exception of "demo"-type example files, the names of all the
examples are of the form ``[slv][PbName]_[ls]_[prec]_[p]``, where

* ``[slv]`` identifies the solver (for CVODE examples this is ``cv``;

* ``[PbName]`` identifies the problem;

* ``[ls]`` identifies the linear solver module used (for examples using
  fixed-point iteration for the nonlinear system solver, ``non`` specifies
  that no linear solver was used);

* ``[prec]`` indicates the CVODE preconditioner module used, ``bp`` for
  CVBandPre ``bbd`` for CVBBDPre (only if applicable, for examples using
  a Krylov linear solver);

* ``[p]`` indicates an example using the parallel vector module NVECTOR_PARALLEL.

The following tables summarize the salient features of each of the
example problems in this document.  Each example is designed to be
relatively self-contained, so that you need only study and/or emulate
the problem that is most closely related to your own.  We group these
examples according to programming language (C, C++, Fortran 2003).



Detailed examples
=================

Deep dives into CVODE example problems written in C are listed below.

.. cssclass:: table-bordered

======================================   ========================
Problem                                  Example Type
======================================   ========================
:ref:`deep_dive.cvRoberts_dns`            serial
:ref:`deep_dive.cvAdvDiff_bnd`            serial
:ref:`deep_dive.cvDiurnal_kry`            serial
:ref:`deep_dive.cvAdvDiff_non_p`          parallel
:ref:`deep_dive.cvDiurnal_kry_p`          parallel
:ref:`deep_dive.cvDiurnal_kry_bbd_p`      parallel
:ref:`deep_dive.cvAdvDiff_non_ph`         HYPRE parallel
======================================   ========================



Deep dives into CVODE example problems written in Fortran 2003 are listed
below.

.. cssclass:: table-bordered

======================================   =================
Problem                                  Example Type
======================================   =================
:ref:`deep_dive.cv_diurnal_kry`           serial
:ref:`deep_dive.cv_diag_kry_bbd_p`        parallel
======================================   =================




Tables of Examples
==================

C Examples
-----------

CVODE example problems written in C are summarized in the table
below.

.. cssclass:: table-bordered

==================================== ============  ============  ==========  =============  =========================================================
Problem                              Integrator    Nonlinear     Linear      Size           Extras
==================================== ============  ============  ==========  =============  =========================================================
:ref:`cvAdvDiff_bnd`                 BDF           Newton        Band        1              user Jacobian
:ref:`cvAdvDiff_bndL`                BDF           Newton        Band        1              LAPACK band solver, user Jacobian
:ref:`cvAnalytic_mels`               BDF           Newton        Custom      1              Matrix-embedded Custom Linear Solver
:ref:`cvDirectDemo_ls`               BDF \\        Newton \\     Dense \\    3 \\           Van der Pol (1) user (2) dif-quot (3) diag + 2 \\
                                     ADAMS         Fixed-point   Band        1              2D advection (1) user (2) dif-quot (3) diag + 2
:ref:`cvDisc_dns`                    BDF           Newton        Dense       2              Solves two separate equations, the second twice
:ref:`cvDiurnal_kry`                 BDF           Newton        SPGMR       200            block-diagonal preconditioner
:ref:`cvDiurnal_kry_bp`              BDF           Newton        SPGMR       200            Solved twice with banded prec. on left and right
:ref:`cvHeat2D_klu`                  BDF           Newton        KLU         1100           sparse matrices
:ref:`cvKrylovDemo_ls`               BDF           Newton        4 Krylov    200            SPGMR, SPBCGS, SPTFQMR used
:ref:`cvKrylovDemo_prec`             BDF           Newton        SPGMR       216            preconditioner matrix (1) implicit GS (2) block-diag
:ref:`cvParticle_dns`                BDF           Newton        Dense       2              user Jacobian and user projection
:ref:`cvPendulum_dns`                BDF           Newton        Dense       4, 2           user projection, solves problem via inexact/exact meth.
:ref:`cvRoberts_dns`                 BDF           Newton        Dense       3              user Jacobian, root-finding
:ref:`cvRoberts_dns_constraints`     BDF           Newton        Dense       3              user Jacobian, root-finding, constraints for all comp.
:ref:`cvRoberts_dns_negsol`          BDF           Newton        Dense       3              user Jacobian, root-finding, negative solution comp.
:ref:`cvRoberts_dns_uw`              BDF           Newton        Dense       3              user Jacobian, root-finding, user error weight
:ref:`cvRoberts_dnsL`                BDF           Newton        Dense       3              user Jacobian, root-finding, LAPACK dense solver
:ref:`cvRoberts_klu`                 BDF           Newton        KLU         :math:`3N`     sparse matrices
:ref:`cvRoberts_block_klu`           BDF           Newton        KLU         :math:`3N`     sparse matrices
:ref:`cvRoberts_sps`                 BDF           Newton        SuperLUMT   :math:`3N`     finite-element, sparse matrices
:ref:`cvRocket_dns`                  BDF           Newton        Dense       2              user Jacobian, root-finding
:ref:`cvAdvDiff_bnd_omp`             BDF           Newton        Band        50             user Jacobian, OpenMP
:ref:`cvAdvDiff_kry_ompdev`          BDF           Newton        Band        50             user Jacobian, OpenMP-dev
:ref:`cvAdvDiff_diag_p`              ADAMS         Newton        Diagonal    10             MPI for user routines
:ref:`cvAdvDiff_non_p`               ADAMS         Fixed-point               10             MPI for user routines
:ref:`cvDiurnal_kry_bbd_p`           BDF           Newton        SPGMR       200            MPI for user routines, solved twice with BBD prec.
:ref:`cvDiurnal_kry_p`               BDF           Newton        SPGMR       200            MPI for user routines, block-diagonal left prec.
:ref:`cvAdvDiff_non_ph`              ADAMS         Fixed-point               10             HYPRE parallel vector with IJ interface
:ref:`cvDiurnal_kry_mpimanyvec`      BDF           Newton        SPGMR       200            MPI for user routines, MPIManyVector module
:ref:`cvAdvDiff_kry_cuda`            BDF           Newton        SPGMR       50             CUDA N_Vector module, RHS, and Jacobian product routines
:ref:`cvAdvDiff_kry_raja`            BDF           Newton        SPGMR       50             RAJA N_Vector module, RHS, and Jacobian product routines
==================================== ============  ============  ==========  =============  =========================================================


C++ Examples
-------------

CVODE example problems written in C++ are summarized in the table
below.

.. cssclass:: table-bordered

=================================  ==========  ===========  ==============  ===============  ===================================
Problem                            Integrator  Nonlinear    Linear          Size             Extras
=================================  ==========  ===========  ==============  ===============  ===================================
:ref:`cv_heat2D`                    BDF         Newton       PCG/SPGMR       :math:`nx*ny`
:ref:`cv_kpr`                       BDF         Newton       Dense           2
:ref:`cv_heat2D_p`                  BDF         Newton       PCG/SPGMR       :math:`nx*ny`    MPI parallel vector and routines
:ref:`cv_heat2D_hypre_ls`           BDF         Newton       PCG/SPGMR       :math:`nx*ny`    HYPRE parallel vector
:ref:`cv_heat2D_hypre_pfmg`         BDF         Newton       PCG/SPGMR       :math:`nx*ny`    HYPRE parallel vector, PFMG prec.
:ref:`cvRoberts_blockdiag_onemkl`   BDF         Newton       oneMKL/SPGMR    300              oneMKL linear solver, user Jac.
:ref:`cvAdvDiff_kry_sycl`           BDF         Newton       SPGMR           50               SYCL parallel vector
=================================  ==========  ===========  ==============  ===============  ===================================

Fortran Examples
-----------------

CVODE example problems written in Fortran 2003 are summarized in the table
below.

.. cssclass:: table-bordered

=================================   ==========  ===========  ======  =============  =================================================
Problem                             Integrator  Nonlinear    Linear  Size           Extras
=================================   ==========  ===========  ======  =============  =================================================
:ref:`cv_advdiff_bnd`               BDF         Newton       Band    50             banded preconditioner, user Jacobian
:ref:`cv_analytic_fp`               ADAMS       Fixed-point          1
:ref:`cv_analytic_sys_dns`          BDF         Newton       Dense   3
:ref:`cv_analytic_sys_dns_jac`      BDF         Newton       Dense   3              user Jacobian matrix
:ref:`cv_analytic_sys_klu`          BDF         Newton       KLU     3              sparse matrices, user Jacobian matrix
:ref:`cv_brusselator_dns`           BDF         Newton       Dense   3              user Jacobian matrix
:ref:`cv_diurnal_kry_bp`            BDF         Newton       SPGMR   200            banded preconditioner
:ref:`cv_diurnal_kry`               BDF         Newton       SPGMR   200            user preconditioner and Jacobian
:ref:`cv_roberts_dns`               BDF         Newton       Dense   3              user Jacobian, root-finding
:ref:`cv_roberts_dns_constraints`   BDF         Newton       Dense   3              user Jacobian, root-finding, constraints set
:ref:`cv_roberts_dnsL`              BDF         Newton       Dense   3              user Jacobian, root-finding, LAPACK dense solver
:ref:`cv_roberts_klu`               BDF         Newton       KLU     3              sparse matrices, root-finding, user Jacobian
:ref:`cv_diag_kry_bbd_p`            BDF         Newton       SPGMR   128            parallel MPI, BBD prec. -- left and right solve
:ref:`cv_diag_kry_p`                BDF         Newton       SPGMR   128            parallel MPI, user prec. -- left and right solve
:ref:`cv_diag_non_p`                ADAMS       Fixed-point          128            parallel MPI
=================================   ==========  ===========  ======  =============  =================================================



Full Section List of Examples and Deep Dives
====================================================

The following sections provide the detailed "deep dive" examples, followed by the
high-level descriptions of all CVODE examples.  For each example we show
our output files, but users should be cautioned that their results may differ
slightly from these.  Differences in solution values may differ within the tolerances,
and differences in cumulative counters, such as numbers of steps or Newton iterations,
may differ from one machine environment to another by as much as 10% to 20%.

In the descriptions below, we make frequent references to the CVODE
User Document :cite:p:`cvode_ug`.  All citations to specific sections
are references to parts of that User Document, unless explicitly stated
otherwise.

.. note::

   The examples in the CVODE distribution are written in such a way as
   to compile and run for any combination of configuration options during
   the installation of SUNDIALS (see the "Install" section in the User Guide).
   As a consequence, they contain portions of code that will not be
   typically present in a user program. For example, all C example programs
   make use of the variables SUNDIALS_EXTENDED_PRECISION and SUNDIALS_DOUBLE_PRECISION
   to test if the solver libraries were built in extended or double precision,
   and use the appropriate conversion specifiers in ``printf`` functions.


.. toctree::
   :maxdepth: 1

   c_serial_deep
   c_parallel_deep
   c_parhyp_deep
   c_tests_deep
   f2003_deep
   c_serial
   c_openmp
   c_openmpdev
   c_parallel
   c_parhyp
   c_mpimanyvec
   c_cuda
   c_raja
   cpp_serial
   cpp_parallel
   cpp_hypre
   cpp_onemkl
   cpp_sycl
   f2003_serial
   f2003_parallel
   References

.. only:: html

   * :ref:`search`
