..
   Programmer(s): Daniel R. Reynolds @ SMU
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


CVode example problems written in C are summarized in the table
below, and are further described in the chapters :ref:`serial_c`,
:ref:`openmp_c`, :ref:`parallel_c` and :ref:`parhyp_c`.

.. cssclass:: table-bordered

================================   ============  ============  ==========  =============  =====================================================
Problem                            Integrator    Nonlinear     Linear      Size           Extras
================================   ============  ============  ==========  =============  =====================================================
:ref:`cvAdvDiff_bnd`               BDF           Newton        Band        1              user Jacobian
:ref:`cvAdvDiff_bndL`              BDF           Newton        Band        1              LAPACK band solver, user Jacobian
:ref:`cvAnalytic_mels`             ----          Newton        Dense       3
:ref:`cvDirectDemo_ls`             BDF           Newton        Dense       3              Van der Pol (1) user (2) dif-quot (3) diag + 2
                                   Adams         fixed-point   Band        1              2D advection (1) user (2) dif-quot (3) diag + 2
:ref:`cvDisc_dns`                  ----          Newton        Dense       3
:ref:`cvDirunal_kry`               BDF           Newton        SPGMR       200            block-diagonal preconditioner
:ref:`cvDirunal_kry_bp`            BDF/BandPre   Newton        SPGMR       200            Solved twice with prec. on left and right
:ref:`cvHeat2D_klu`                BDF           Newton        KLU
:ref:`cvKrylovDemo_ls`             BDF           Newton        3 Krylov    200            SPGMR, SPBCGS, SPTFQMR used
:ref:`cvKrylovDemo_prec`           GMRES         Newton        SPGMR       12             preconditioner matrix (1) implicit GS (2) block-diag
:ref:`cvParticle_dns`              ----          Newton        PCG         N
:ref:`cvPendulum_dns`              ----          Newton        PCG         (dynamic)      adaptive vector resizing
:ref:`cvRoberts_block_klu`         ----          Newton        SPGMR       216            multiple preconditioners
:ref:`cvRoberts_dns`               BDF           Newton        Dense       3              user rootfinding
:ref:`cvRoberts_dns_constraints`   BDF           Newton        Dense       3              user rootfinding
:ref:`cvRoberts_dns_negsol`        ----          Newton        SPGMR       200            HYPRE parallel vector
:ref:`cvRoberts_dns_uw`            BDF           Newton        Dense       3              user error weight
:ref:`cvRoberts_dnsL`              BDF           Newton        Dense       3              LAPACK dense solver, user rootfinding
:ref:`cvRoberts_klu`               BDF           Newton        KLU         :math:`3N`     sparse matrices
:ref:`cvRoberts_sps`               BDF           Newton        SuperLU_MT  :math:`3N`     finite-element, sparse matrices
:ref:`cvRocket_dns`                ----          Newton        SPGMR       200            HYPRE parallel vector
:ref:`cvAdvDiff_bnd_omp`           DIRK          Newton        SPGMR       200            HYPRE parallel vector
:ref:`cvAdvDiff_kry_ompdev`        DIRK          Newton        SPGMR       200            HYPRE parallel vector
:ref:`cvAdvDiff_diag_p`            DIRK          Newton        SPGMR       200            HYPRE parallel vector
:ref:`cvAdvDiff_non_p`             DIRK          Newton        SPGMR       200            HYPRE parallel vector
:ref:`cvAdvDiff_kry_bbd_p`         DIRK          Newton        SPGMR       200            HYPRE parallel vector
:ref:`cvAdvDiff_kry_p`             DIRK          Newton        SPGMR       200            HYPRE parallel vector
:ref:`cvAdvDiff_non_ph`            DIRK          Newton        SPGMR       200            HYPRE parallel vector
:ref:`cvDirunal_kry_mpimanyvec`    DIRK          Newton        SPGMR       200            HYPRE parallel vector
================================   ============  ============  ==========  =============  =====================================================


CVode example problems written in C++ are summarized in the table
below, and are further described in the chapters :ref:`serial_cpp`,
:ref:`parallel_cpp`, :ref:`parhyp_cpp`, :ref:`onemkl_cpp`, and
:ref:`sycl_cpp`.

.. cssclass:: table-bordered

=================================  ==========  ===========  ======  =============  =================================
Problem                            Integrator  Nonlinear    Linear  Size           Extras
=================================  ==========  ===========  ======  =============  =================================
:ref:`cv_heat2D`                    DIRK        Newton       PCG     :math:`nx*ny`  parallel
:ref:`cv_kpr`                       DIRK        Newton       Dense   3
:ref:`cv_heat2D_p`                  DIRK        Newton       Dense   3
:ref:`cv_heat2D_hypre_ls`           DIRK        Newton       Dense   3
:ref:`cv_heat2D_hypre_pfmg`         DIRK        Newton       Dense   3
:ref:`cvRoberts_blockdiag_onemkl`   DIRK        Newton       Dense   3
:ref:`cvAdvDiff_kry_sycl`           DIRK        Newton       Dense   3
=================================  ==========  ===========  ======  =============  =================================


CVode example problems written in Fortran 2003 are summarized in the table
below, and are further described in the chapters :ref:`serial_f2003` and
:ref:`parallel_f2003`.

.. cssclass:: table-bordered

=================================   ==========  ===========  ======  =============  =================================================
Problem                             Integrator  Nonlinear    Linear  Size           Extras
=================================   ==========  ===========  ======  =============  =================================================
:ref:`cv_analytic_fp`               DIRK        Newton       SPGMR   10             banded preconditioner
:ref:`cv_analytic_sys_dns`          DIRK        Newton       Dense   3              LAPACK dense solver, rootfinding
:ref:`cv_analytic_sys_dns_jac`      CV          Newton       Dense   3
:ref:`cv_analytic_sys_klu`          DIRK        Newton       KLU     3N             finite-element, :math:`M\ne I`, sparse matrices
:ref:`cv_brusselator_dns`           DIRK        Newton       SPGMR   10*NProcs      parallel BBD preconditioner
:ref:`cv_roberts_dns_constraints`   ERK         N.A.         N.A.    10*NProcs      parallel
:ref:`cv_roberts_dns`               DIRK        Newton       PCG     :math:`nx*ny`  parallel
:ref:`cv_roberts_dnsL`              DIRK        Newton       PCG     :math:`nx*ny`  parallel
:ref:`cv_roberts_klu`               DIRK        Newton       PCG     :math:`nx*ny`  parallel
:ref:`cv_roberts_sps`               DIRK        Newton       PCG     :math:`nx*ny`  parallel
:ref:`cv_diag_kry_bbd_p`            DIRK        Newton       PCG     :math:`nx*ny`  parallel
:ref:`cv_diag_kry_p`                DIRK        Newton       PCG     :math:`nx*ny`  parallel
:ref:`cv_diag_non_p`                DIRK        Newton       PCG     :math:`nx*ny`  parallel
=================================   ==========  ===========  ======  =============  =================================================





.. only:: html

   Further details on many of the above-listed examples are provided
   in the following chapters:

.. toctree::
   :maxdepth: 1

   c_serial
   c_openmp
   c_parallel
   c_parhyp
   cpp_serial
   cpp_parallel
   f2003_serial
   f2003_parallel
   References

.. only:: html

   * :ref:`search`
