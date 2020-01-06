..
   Programmer(s): Daniel R. Reynolds @ SMU
   ----------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2002-2020, Lawrence Livermore National Security
   and Southern Methodist University.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   ----------------------------------------------------------------
   ARKode documentation master file, created by
   sphinx-quickstart on Sat Dec 22 20:38:03 2012.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

====================================
ARKode Documentation
====================================

This is the documentation for ARKode, an adaptive step time
integration package for stiff, nonstiff and mixed stiff/nonstiff
systems of ordinary differential equations (ODEs) using Runge-Kutta
(i.e. one-step, multi-stage) methods.  The ARKode solver is a
component of the `SUNDIALS
<https://computing.llnl.gov/casc/sundials/main.html>`_ suite of
nonlinear and differential/algebraic equation solvers. It is designed
to have a similar user experience to the `CVODE
<https://computing.llnl.gov/casc/sundials/description/description.html#descr_cvode>`_
solver, including user modes to allow adaptive integration to specified
output times, return after each internal step and root-finding
capabilities, and for calculations in serial, using shared-memory
parallelism (via OpenMP, Pthreads, CUDA, Raja) or distributed-memory
parallelism (via MPI).  The default integration and solver options
should apply to most users, though control over nearly all internal
parameters and time adaptivity algorithms is enabled through optional
interface routines.

ARKode is written in C, with C++ and Fortran interfaces.

ARKode is developed by `Southern Methodist University
<http://www.smu.edu>`_, with support by the `US Department of Energy
<http://www.doe.gov>`_ through the `FASTMath
<http://www.fastmath-scidac.org/>`_ SciDAC Institute, under subcontract
B598130 from `Lawrence Livermore National Laboratory
<http://www.llnl.gov>`_.



.. only:: html

   Documentation sections:

.. toctree::
   :maxdepth: 1

   Introduction
   Mathematics
   Organization
   ARKStep_c_interface/index.rst
   ERKStep_c_interface/index.rst
   MRIStep_c_interface/index.rst
   ARKode_f_interface/index.rst
   ARKodeButcherTable
   nvectors/index.rst
   sunmatrix/index.rst
   sunlinsol/index.rst
   sunnonlinsol/index.rst
   Install
   Constants
   Butcher
   History
   References

.. only:: html

   * :ref:`genindex`
   * :ref:`search`
