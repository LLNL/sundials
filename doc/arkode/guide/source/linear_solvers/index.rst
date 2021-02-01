..
   Programmer(s): Daniel R. Reynolds @ SMU
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

.. _LinearSolvers:

========================
Linear Solvers in ARKode
========================

In this chapter, we describe linear solver code components that are
included in SUNDIALS, but which are of potential use as generic
packages in themselves, either in conjunction with the use of SUNDIALS
or separately.

These generic linear solver modules are organized in three families of
solvers, the DLS family, which includes direct linear solvers
appropriate for sequential computations; the SLS family, which
includes sparse matrix solvers; and the SPILS family, which includes
scaled preconditioned iterative (Krylov) linear solvers. The solvers
in each family share common data structures and functions. 


The :ref:`DLS <LinearSolvers.DLS>` family contains the following two
generic linear solvers: 

* The DENSE package, a linear solver for dense matrices either
  specified through a matrix type (defined below) or as simple
  arrays. 
* The BAND package, a linear solver for banded matrices either
  specified through a matrix type (defined below) or as simple
  arrays. 

We further note that this family also includes the BLAS/LAPACK linear
solvers (dense and band) available to the SUNDIALS solvers, but these
are not discussed here. 


The :ref:`SLS <LinearSolvers.SLS>` family contains a sparse matrix
package and interfaces between it and two sparse direct solver packages: 

* The KLU package, a linear solver for compressed-sparse-column (CSC)
  or compressed-sparse-row (CSR) matrices, [KLU]_, [DP2010]_.
* The SUPERLUMT package, a threaded linear solver for
  CSC matrices, [SuperLUMT]_, [L2005]_, [DGL1999]_.


The :ref:`SPILS <LinearSolvers.SPILS>` family contains the following
generic linear solvers: 

* The SPGMR package, a solver for the scaled preconditioned GMRES
  method. 
* The SPFGMR package, a solver for the scaled preconditioned Flexible
  GMRES method. 
* The SPBCG package, a solver for the scaled preconditioned Bi-CGStab
  method. 
* The SPTFQMR package, a solver for the scaled preconditioned TFQMR
  method. 
* The PCG package, a solver for the preconditioned conjugate gradient
  method. 

For reasons related to installation, the names of the files involved
in these packages begin with teh prefix ``sundials_``.  But despite
this, each of the DLS and SPILS solvers is in fact generic, in that
it is usable completely independently of SUNDIALS. 

For the sake of space, the functions for the DENSE, BAND modules
that work with a matrix type, and the functions in the SPGMR, SPFGMR, 
SPBCG, SPTFQMR and PCG modules are only summarized briefly, since
they are less likely to be of direct use in connection with a SUNDIALS
solver.  However, the functions for dense matrices treated as simple
arrays and sparse matrices are fully described, because we anticipate
that they will be useful in the implementation of preconditioners used
with the combination of one of the SUNDIALS solvers and one of the
SPILS linear solvers. 

Lastly, it is possible to supply customized linear solvers to ARKode,
in that the ARKode solvers only require the existence of a minimal set
of generic routines.  Through attaching user-supplied routines for
these function pointers, it is possible to use arbitrary approaches
for solution to the implicit linear systems arising during an ARKode
solve. 

Specifics of these built-in linear solver packages, as well as the
generic linear solver interface, are provided in the following
sub-sections:

.. toctree::
   :maxdepth: 1

   DLS
   SLS
   SPILS
   custom
