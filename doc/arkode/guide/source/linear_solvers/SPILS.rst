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


.. _LinearSolvers.SPILS:

The SPILS modules: SPGMR, SPFGMR, SPBCG, SPTFQMR and PCG
==========================================================

The SPILS modules contain implementations of some of the most commonly
use scaled preconditioned Krylov solvers.  A linear solver module from
the SPILS family can be used in conjunction with any NVECTOR
implementation library.

In the subsections below, we discuss the iterative linear solvers
accessible to ARKode: SPGMR, SPFGMR, SPBCG, SPTFQMR and PCG.  Due to
the similarities between these modules, we provide a more complete
description of only the SPGMR interface, and for the remaining solvers
only discuss the salient differences.



The SPGMR module
-----------------------------------------

The SPGMR package, in the files ``sundials_spgmr.h`` and
``sundials_spgmr.c``, includes an implementation of the scaled
preconditioned GMRES method.  A separate code module, implemented in
``sundials_iterative.h`` and ``sundials_iterative.c``, contains
auxiliary functions that support SPGMR, as well as the other Krylov
solvers in SUNDIALS (SPFGMR, SPBCG, SPTFQMR and PCG).  For full
details, including usage instructions, see the header files
``sundials_spgmr.h`` and ``sundials_iterative.h``. 

The files comprising the SPGMR generic linear solver, and their
locations in the SUNDIALS ``srcdir``, are as follows:

* header files (located in ``srcdir/include/sundials``):

  ``sundials_spgmr.h``, ``sundials_iterative.h``,
  ``sundials_nvector.h``, ``sundials_types.h``, ``sundials_math.h``,
  ``sundials_config.h``

* source files (located in ``srcdir/src/sundials``):

  ``sundials_spgmr.c``, ``sundials_iterative.c``, ``sundials_nvector.c``


Only two of the preprocessing directives in the header file
``sundials_config.h`` are required to use the SPGMR package by itself: 

* (required) definition of the precision of the SUNDIALS type
  ``realtype``. One of the following three lines must be present:

  .. code-block:: c

     #define SUNDIALS_DOUBLE_PRECISION 1
     #define SUNDIALS_SINGLE_PRECISION 1
     #define SUNDIALS_EXTENDED_PRECISION 1

* (optional) use of generic math functions:

  .. code-block:: c

     #define SUNDIALS_USE_GENERIC_MATH 1


The ``sundials_types.h`` header file defines the SUNDIALS ``realtype``
and ``booleantype`` types and the macro ``RCONST``, while the
``sundials_math.h`` header file is needed for the macros ``SUNMIN``,
``SUNMAX``, and ``SUNSQR``, and the functions ``SUNRabs`` and ``SUNRsqrt``.

The generic NVECTOR files, ``sundials_nvector.h`` and
``sundials_nvector.c`` are needed for the definition of the generic
``N_Vector`` type and functions.  The NVECTOR functions used by the
SPGMR module are: :c:func:`N_VDotProd()`, :c:func:`N_VLinearSum()`,
:c:func:`N_VScale()`, :c:func:`N_VProd()`, :c:func:`N_VDiv()`,
:c:func:`N_VConst()`, :c:func:`N_VClone()`,
:c:func:`N_VCloneVectorArray()`, :c:func:`N_VDestroy()`, and
:c:func:`N_VDestroyVectorArray()`. 

The nine files listed above can be extracted from the SUNDIALS
``srcdir`` and compiled by themselves into an SPGMR library or into a
separate user code. 

The following functions are available in the SPGMR package:

* ``SpgmrMalloc``: allocates memory for ``SpgmrSolve``;
* ``SpgmrSolve``: solves :math:`Ax = b` using the SPGMR method;
* ``SpgmrFree``: frees memory allocated by ``SpgmrMalloc``.


The following functions are available in the support package
``sundials_iterative.h`` and ``sundials_iterative.c``:

* ``ModifiedGS``: performs the modified Gram-Schmidt orthogonalization
  procedure;
* ``ClassicalGS``: performs the classical Gram-Schmidt
  orthogonalization procedure;
* ``QRfact``: performs the QR factorization of a Hessenberg matrix;
* ``QRsol``: solves a least squares problem with a Hessenberg matrix
  factored by ``QRfact``. 




The SPFGMR module
-----------------------------------------

The SPFGMR package, in the files ``sundials_spfgmr.h`` and
``sundials_spfgmr.c``, includes an implementation of the scaled
preconditioned Flexible Generalized Minimum Residual method
[S1993]_. For full details, including usage instructions, see the file
``sundials_spfgmr.h``.  

The files needed to use the SPFGMR module by itself are the same as for
the SPGMR module, but with ``sundials_spfgmr.(h,c)`` in place of
``sundials_spgmr.(h,c)``.

The following functions are available in the SPFGMR package:

* ``SpfgmrMalloc``: allocates memory for ``SpfgmrSolve``;
* ``SpfgmrSolve``: solves :math:`Ax = b` using the SPFGMR method;
* ``SpfgmrFree``: frees memory allocated by ``SpfgmrMalloc``.



The SPBCG module
-----------------------------------------

The SPBCG package, in the files ``sundials_spbcgs.h`` and
``sundials_spbcgs.c``, includes an implementation of the scaled
preconditioned Bi-CGStab method. For full details, including usage
instructions, see the file ``sundials_spbcgs.h``.

The files needed to use the SPBCG module by itself are the same as for
the SPGMR module, but with ``sundials_spbcgs.(h,c)`` in place of
``sundials_spgmr.(h,c)``.

The following functions are available in the SPBCG package:

* ``SpbcgMalloc``: allocates memory for ``SpbcgSolve``;
* ``SpbcgSolve``: solves :math:`Ax = b` using the SPBCG method;
* ``SpbcgFree``: frees memory allocated by ``SpbcgMalloc``.



The SPTFQMR module
-----------------------------------------


The SPTFQMR package, in the files ``sundials_sptfqmr.h`` and
``sundials_sptfqmr.c``, includes an implementation of the scaled
preconditioned TFQMR method. For full details, including usage
instructions, see the file ``sundials_sptfqmr.h``.

The files needed to use the SPTFQMR module by itself are the same as
for the SPGMR module, but with ``sundials_sptfqmr.(h,c)`` in place of
``sundials_spgmr.(h,c)``.

The following functions are available in the SPTFQMR package:

* ``SptfqmrMalloc``: allocates memory for ``SptfqmrSolve``;
* ``SptfqmrSolve``: solves :math:`Ax = b` using the SPTFQMR method;
* ``SptfqmrFree``: frees memory allocated by ``SptfqmrMalloc``.



The PCG module
-----------------------------------------

The PCG package, in the files ``sundials_pcg.h`` and
``sundials_pcg.c``, includes an implementation of the 
preconditioned conjugate gradient method.  We note that due to the
requirement of symmetric linear systems for the conjugate gradient
method, this solver should only be used for problems with symmetric
linear operators.  Furthermore, aside from allowing a weight vector
for computing weighted convergence norms, no variable or equation
scaling is allowed for systems using this solver.  For full details,
including usage instructions, see the file ``sundials_pcg.h``.

The files needed to use the PCG module by itself are the same as for
the SPGMR module, but with ``sundials_pcg.(h,c)`` in place of
``sundials_spgmr.(h,c)``.  

The following functions are available in the PCG package:

* ``PcgMalloc``: allocates memory for ``PcgSolve``;
* ``PcgSolve``: solves :math:`Ax = b` using the PCG method;
* ``PcgFree``: frees memory allocated by ``PcgMalloc``.
