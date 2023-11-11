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

.. ARKode_example documentation master file, created by
   sphinx-quickstart on Sat Dec 22 20:38:03 2012.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

==============================================
ARKode Example documentation
==============================================

This is the documentation for the ARKode examples.  ARKode is an
adaptive step time integration package for stiff, nonstiff and
multi-rate systems of ordinary differential equations (ODEs).  
The ARKode solver is a component of the `SUNDIALS
<https://computing.llnl.gov/projects/sundials>`_ suite of
nonlinear and differential/algebraic equation solvers. It is designed
to have a similar user experience to the `CVODE
<https://computing.llnl.gov/casc/sundials/description/description.html#descr_cvode>`_
solver, with user modes to allow adaptive integration to specified
output times, return after each internal step and root-finding
capabilities, for calculations both in serial and parallel (via
MPI). The default integration and solver options should apply to most
users, though complete control over all internal parameters and time
adaptivity algorithms is enabled through optional interface routines.  

ARKode is developed by `Southern Methodist University
<http://www.smu.edu>`_, with support by the `US Department of Energy
<http://www.doe.gov>`_ through the `FASTMath
<http://www.fastmath-scidac.org/>`_ SciDAC Institute, under subcontract
B598130 from `Lawrence Livermore National Laboratory
<http://www.llnl.gov>`_. 

Along with the ARKode solver, we have created a suite of example
problems demonstrating its usage on applications written in C, C++ and
Fortran.  These examples demonstrate a large variety
of ARKode solver options, including explicit, implicit and ImEx
solvers, root-finding, Newton and fixed-point nonlinear solvers,
direct and iterative linear solvers, adaptive resize capabilities, and
the Fortran solver interface.  While these examples are not an
exhaustive set of all possible usage scenarios, they are designed to
show a variety of exemplars, and can be used as templates for new
problems using ARKode's solvers.  Further information on the ARKode
package itself may be found in the accompanying ARKode user guide
[R2018]_.

The following tables summarize the salient features of each of the
example problems in this document.  Each example is designed to be
relatively self-contained, so that you need only study and/or emulate
the problem that is most closely related to your own.  We group these
examples according to programming language (C, C++, Fortran).


ARKode example problems written in C are summarized in the table
below, and are further described in the chapters :ref:`serial_c`,
:ref:`openmp_c`, :ref:`parallel_c` and :ref:`parhyp_c`. 

.. cssclass:: table-bordered

================================  ==========  ===========  ==========  =============  =====================================================
Problem                           Integrator  Nonlinear    Linear      Size           Extras
================================  ==========  ===========  ==========  =============  =====================================================
:ref:`ark_analytic`               DIRK        Newton       Dense       1              
:ref:`ark_analytic_nonlin`        ERK         N.A.         N.A.        1              ERKStep timestepping module
:ref:`ark_brusselator`            DIRK        Newton       Dense       3              
:ref:`ark_brusselator_fp`         ARK         Fixed-point  N.A.        3              
:ref:`ark_robertson`              DIRK        Newton       Dense       3              
:ref:`ark_robertson_root`         DIRK        Newton       Dense       3              rootfinding
:ref:`ark_brusselator1D`          DIRK        Newton       Band        3N             
:ref:`ark_brusselator1D_omp`      DIRK        Newton       Band        3N             OpenMP-enabled
:ref:`ark_brusselator1D_klu`      DIRK        Newton       KLU         3N             sparse matrices
:ref:`ark_brusselator1D_FEM_slu`  DIRK        Newton       SuperLU_MT  3N             finite-element, :math:`M\ne I`, sparse matrices
:ref:`ark_heat1D`                 DIRK        Newton       PCG         N              
:ref:`ark_heat1D_adapt`           DIRK        Newton       PCG         (dynamic)      adaptive vector resizing
:ref:`ark_KrylovDemo_prec`        DIRK        Newton       SPGMR       216            multiple preconditioners
:ref:`ark_diurnal_kry_bbd_p`      DIRK        Newton       SPGMR       200            parallel, BBD preconditioner
:ref:`ark_diurnal_kry_p`          DIRK        Newton       SPGMR       200            parallel, block-diagonal precond.
:ref:`ark_diurnal_kry_ph`         DIRK        Newton       SPGMR       200            HYPRE parallel vector
================================  ==========  ===========  ==========  =============  =====================================================


ARKode example problems written in C++ are summarized in the table
below, and are further described in the chapters :ref:`serial_cpp` and
:ref:`parallel_cpp`. 

.. cssclass:: table-bordered

=======================  ==========  ===========  ======  =============  =================================
Problem                  Integrator  Nonlinear    Linear  Size           Extras
=======================  ==========  ===========  ======  =============  =================================
:ref:`ark_analytic_sys`  DIRK        Newton       Dense   3              
:ref:`ark_heat2D`        DIRK        Newton       PCG     :math:`nx*ny`  parallel
=======================  ==========  ===========  ======  =============  =================================


..
   ARKode example problems written in Fortran 77 are summarized in the table
   below, and are further described in the chapters :ref:`serial_f77` and
   :ref:`parallel_f77`.

   .. cssclass:: table-bordered

   ==========================   ==========  ===========  ======  =============  =================================
   Problem                      Integrator  Nonlinear    Linear  Size           Extras
   ==========================   ==========  ===========  ======  =============  =================================
   :ref:`fark_diurnal_kry_bp`   DIRK        Newton       SPGMR   10             banded preconditioner
   :ref:`fark_roberts_dnsL`     DIRK        Newton       Dense   3              LAPACK dense solver, rootfinding
   :ref:`fark_diag_kry_bbd_p`   DIRK        Newton       SPGMR   10*NProcs      parallel BBD preconditioner
   :ref:`fark_diag_non_p`       ERK         N.A.         N.A.    10*NProcs      parallel
   ==========================   ==========  ===========  ======  =============  =================================


   ARKode example problems written in Fortran 90 are summarized in the table
   below, and are further described in the chapters :ref:`serial_f90` and
   :ref:`parallel_f90`.

   .. cssclass:: table-bordered

   ==========================  ==========  =========  ======  =============  ===============================================
   Problem                     Integrator  Nonlinear  Linear  Size           Extras
   ==========================  ==========  =========  ======  =============  ===============================================
   :ref:`ark_bruss`            ARK         Newton     Dense   3
   :ref:`ark_bruss1D_FEM_klu`  DIRK        Newton     KLU     3N             finite-element, :math:`M\ne I`, sparse matrices
   :ref:`fark_heat2D`          DIRK        Newton     PCG     :math:`nx*ny`  parallel
   ==========================  ==========  =========  ======  =============  ===============================================





.. only:: html

   Further details on each of the above-listed examples are provided
   in the following chapters:

.. toctree::
   :maxdepth: 1

   c_serial
   c_openmp
   c_parallel
   c_parhyp
   cpp_serial
   cpp_parallel
   References
..
  Remove F77 interface examples
  f77_serial
  f77_parallel
  f90_serial
  f90_parallel

.. only:: html

   * :ref:`search`
