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

.. CVode_example documentation master file, created by
   sphinx-quickstart on Wed Jun 14 02:10:03 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

==============================================
CVode Example documentation
==============================================

CVode Overview
=============================

This is the documentation for the CVode examples.  CVode is a
Krylov method integration package for stiff, nonstiff and
multi-rate systems of ordinary differential equations (ODEs).
The CVode solver is a component of the `SUNDIALS
<https://computing.llnl.gov/projects/sundials>`_ suite of
nonlinear and differential/algebraic equation solvers. It is designed
to have a similar user experience to the `CVODE
<https://computing.llnl.gov/projects/sundials/cvode>`_
solver, with user modes to allow adaptive integration to specified
output times, return after each internal step and root-finding
capabilities, for calculations both in serial and parallel (via
MPI). The default integration and solver options should apply to most
users, though complete control over all internal parameters and time
adaptivity algorithms is enabled through optional interface routines.

CVode is developed by `Lawrence Livermore National Laboratory
<http://www.llnl.gov>`_, with support by the `US Department of Energy
<http://www.doe.gov>`_ through the `FASTMath
<https://scidac5-fastmath.lbl.gov/>`_ SciDAC-5 Institute, under
subcontract 07NA27344.

Along with the CVode solver, we have created a suite of example
problems demonstrating its usage on applications written in C, C++ and
Fortran 2003.  These examples demonstrate a large variety of CVode
solver options, including explicit, implicit and ImEx solvers, root-
finding, Newton and fixed-point nonlinear solvers, direct and iterative
linear solvers, adaptive resize capabilities, and the Fortran solver
interface.  While these examples are not an exhaustive set of all
possible usage scenarios, they are designed to show a variety of
exemplars, and can be used as templates for new problems using CVode's
solvers.  Further information on the CVode package itself may be found
in the accompanying CVode user guide [R2018]_.

The following tables summarize the salient features of each of the
example problems in this document.  Each example is designed to be
relatively self-contained, so that you need only study and/or emulate
the problem that is most closely related to your own.  We group these
examples according to programming language (C, C++, Fortran 2003).


Tables of Examples
==========================

C Examples
-----------

CVode example problems written in C are summarized in the table
below, and are further described in the chapters :ref:`serial_c`,
:ref:`openmp_c`, :ref:`openmpdev_c`, :ref:`parallel_c`, :ref:`parhyp_c`,
:ref:`mpimanyvec_c`, :ref:`cuda_c`, and :ref:`raja_c`.

.. cssclass:: table-bordered

================================   ============  ============  ==========  =============  =========================================================
Problem                            Integrator    Nonlinear     Linear      Size           Extras
================================   ============  ============  ==========  =============  =========================================================
:ref:`cvAdvDiff_bnd`                BDF           Newton        Band        1              user Jacobian
:ref:`cvAdvDiff_bndL`               BDF           Newton        Band        1              LAPACK band solver, user Jacobian
:ref:`cvAnalytic_mels`              BDF           Newton        Custom      1              Matrix-embedded Custom Linear Solver
:ref:`cvDirectDemo_ls`              BDF \\        Newton \\     Dense \\    3 \\           Van der Pol (1) user (2) dif-quot (3) diag + 2 \\
                                    ADAMS         Fixed-point   Band        1              2D advection (1) user (2) dif-quot (3) diag + 2
:ref:`cvDisc_dns`                   BDF           Newton        Dense       2              Solves two separate equations, the second twice
:ref:`cvDiurnal_kry`                BDF           Newton        SPGMR       200            block-diagonal preconditioner
:ref:`cvDiurnal_kry_bp`             BDF           Newton        SPGMR       200            Solved twice with banded prec. on left and right
:ref:`cvHeat2D_klu`                 BDF           Newton        KLU         1100           sparse matrices
:ref:`cvKrylovDemo_ls`              BDF           Newton        4 Krylov    200            SPGMR, SPBCGS, SPTFQMR used
:ref:`cvKrylovDemo_prec`            BDF           Newton        SPGMR       216            preconditioner matrix (1) implicit GS (2) block-diag
:ref:`cvParticle_dns`               BDF           Newton        Dense       2              user Jacobian and user projection
:ref:`cvPendulum_dns`               BDF           Newton        Dense       4, 2           user projection, solves problem via inexact/exact meth.
:ref:`cvRoberts_dns`                BDF           Newton        Dense       3              user Jacobian, root-finding
:ref:`cvRoberts_dns_constraints`    BDF           Newton        Dense       3              user Jacobian, root-finding, constraints for all comp.
:ref:`cvRoberts_dns_negsol`         BDF           Newton        Dense       3              user Jacobian, root-finding, negative solution comp.
:ref:`cvRoberts_dns_uw`             BDF           Newton        Dense       3              user Jacobian, root-finding, user error weight
:ref:`cvRoberts_dnsL`               BDF           Newton        Dense       3              user Jacobian, root-finding, LAPACK dense solver
:ref:`cvRoberts_klu`                BDF           Newton        KLU         :math:`3N`     sparse matrices
:ref:`cvRoberts_block_klu`          BDF           Newton        KLU         :math:`3N`     sparse matrices
:ref:`cvRoberts_sps`                BDF           Newton        SuperLUMT   :math:`3N`     finite-element, sparse matrices
:ref:`cvRocket_dns`                 BDF           Newton        Dense       2              user Jacobian, root-finding
:ref:`cvAdvDiff_bnd_omp`            BDF           Newton        Band        50             user Jacobian, OpenMP
:ref:`cvAdvDiff_kry_ompdev`         BDF           Newton        Band        50             user Jacobian, OpenMP-dev
:ref:`cvAdvDiff_diag_p`             ADAMS         Newton        Diagonal    10             MPI for user routines
:ref:`cvAdvDiff_non_p`              ADAMS         Fixed-point               10             MPI for user routines
:ref:`cvDiurnal_kry_bbd_p`          BDF           Newton        SPGMR       200            MPI for user routines, solved twice with BBD prec.
:ref:`cvDiurnal_kry_p`              BDF           Newton        SPGMR       200            MPI for user routines, block-diagonal left prec.
:ref:`cvAdvDiff_non_ph`             ADAMS         Fixed-point               10             HYPRE parallel vector with IJ interface
:ref:`cvDiurnal_kry_mpimanyvec`     BDF           Newton        SPGMR       200            MPI for user routines, MPIManyVector module
:ref:`cvAdvDiff_kry_cuda`
:ref:`cvAdvDiff_kry_raja`
================================   ============  ============  ==========  =============  =========================================================

C Examples Deep Dives
----------------------

Deep dives into CVode example problems written in C are listed below,
and are further described in great length in the chapters
:ref:`serial_deep_c`, :ref:`parallel_deep_c`, and :ref:`parhyp_deep_c`.
Further testing is discussed in :ref:`tests_deep_c`.

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

C++ Examples
-------------

CVode example problems written in C++ are summarized in the table
below, and are further described in the chapters :ref:`serial_cpp`,
:ref:`parallel_cpp`, :ref:`hypre_cpp`, :ref:`onemkl_cpp`, and :ref:`sycl_cpp`.

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

CVode example problems written in Fortran 2003 are summarized in the table
below, and are further described in the chapters :ref:`serial_f2003` and
:ref:`parallel_f2003`.

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

Fortran Deep Dives
-------------------

Deep dives into CVode example problems written in Fortran 2003 are listed
below, and are further described in great length in the chapter :ref:`deep_f2003`.

.. cssclass:: table-bordered

======================================   =================
Problem                                  Example Type
======================================   =================
:ref:`deep_dive.cv_diurnal_kry`           serial
:ref:`deep_dive.cv_diag_kry_bbd_p`        parallel
======================================   =================


Chapter List of Examples and Deep Dives
====================================================

Further details on many of the above-listed examples are provided
in the following chapters:

.. toctree::
   :maxdepth: 1

   c_serial
   c_openmp
   c_openmpdev
   c_parallel
   c_parhyp
   c_mpimanyvec
   c_cuda
   c_raja
   c_serial_deep
   c_parallel_deep
   c_parhyp_deep
   c_tests_deep
   cpp_serial
   cpp_parallel
   cpp_hypre
   cpp_onemkl
   cpp_sycl
   f2003_serial
   f2003_parallel
   f2003_deep
   References

.. only:: html

   * :ref:`search`
