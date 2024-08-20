.. ----------------------------------------------------------------
   Programmer(s): Steven B. Roberts @ LLNL
   ----------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2002-2024, Lawrence Livermore National Security
   and Southern Methodist University.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   ----------------------------------------------------------------

.. _ARKODE.Usage.SplittingStep.ARKodeSplittingCoefficients:

=========================================
Operator Splitting Coefficients Structure
=========================================

To store the coefficients representing an operator splitting method, ARKODE
provides the :c:type:`ARKodeSplittingCoefficients` type and several related
utility routines. We use the following notation

