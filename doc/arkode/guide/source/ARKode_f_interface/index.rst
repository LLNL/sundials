..
   Programmer(s): Daniel R. Reynolds @ SMU and
                  Cody J. Balos @ LLNL
   ----------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2002-2021, Lawrence Livermore National Security
   and Southern Methodist University.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   ----------------------------------------------------------------

:tocdepth: 3

.. _FortranInterfaces:

=====================================
Using ARKode for Fortran Applications
=====================================

Fortran 2003 interfaces to each of the time-stepping modules as well as a
Fortran 77 style interface to the ARKStep time-stepping module are
provided to support the use of ARKode, for the solution of ODE systems,
in a mixed Fortran/C setting. While ARKode is written in C, it is assumed
here that the user's calling program and user-supplied problem-definining
rotuines are written in Fortran.

.. toctree::
   :maxdepth: 1

   F2003Module
   FARKODE

