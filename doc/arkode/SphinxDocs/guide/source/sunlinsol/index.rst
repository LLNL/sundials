..
   Programmer(s): Daniel R. Reynolds @ SMU
   ----------------------------------------------------------------
   Copyright (c) 2017, Southern Methodist University.
   All rights reserved.
   For details, see the LICENSE file.
   ----------------------------------------------------------------

:tocdepth: 3

.. _SUNLinSol:

=================================
Linear Solver Data Structures
=================================

The SUNDIALS library comes packaged with a variety of
``SUNLinearSolver`` implementations, designed for simulations
requiring either direct or matrix-free iterative linear solvers for
problems in serial or shared-memory parallel environments.  SUNDIALS
additionally provides a simple interface for generic linear solvers
(akin to a C++ *abstract base class*).  All of the major SUNDIALS
packages (CVODE(s), IDA(s), KINSOL, ARKODE), are constructed to only
depend on these generic linear solver operations, making them
immediately extensible to new user-defined linear solver objects.



.. toctree::
   :maxdepth: 1

   SUNLinSol_API
   SUNLinSol_Dense
   SUNLinSol_Band
   SUNLinSol_LapackDense
   SUNLinSol_LapackBand
   SUNLinSol_KLU
   SUNLinSol_SuperLUMT
   SUNLinSol_SPGMR
   SUNLinSol_SPFGMR
   SUNLinSol_SPBCGS
   SUNLinSol_SPTFQMR
   SUNLinSol_PCG
   SUNLinSol_Examples
   ARKode_requirements
