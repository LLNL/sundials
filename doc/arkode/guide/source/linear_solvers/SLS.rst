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


.. _LinearSolvers.SLS:

The SLS modules
========================================

SUNDIALS provides a sparse matrix type based on the
compressed-sparse-column and compressed-sparse-row formats, along with
sparse matrix support functions.  In addition, SUNDIALS provides
interfaces to the publicly available KLU and SuperLU_MT sparse direct
solver packages.  The files comprising the SLS matrix module, used in
the KLU and SUPERLUMT linear solver packages, and their locations in
the SUNDIALS ``srcdir``, are as follows: 

* header files (located in ``srcdir/include/sundials``):

  ``sundials_sparse.h``, ``sundials_klu_impl.h``,
  ``sundials_superlumt_impl.h``, ``sundials_types.h``,
  ``sundials_math.h``, ``sundials_config.h``

* source files (located in ``srcdir/src/sundials``):

  ``sundials_sparse.c``, ``sundials_math.c``

Only two of the preprocessing directives in the header file
``sundials_config.h`` are relevant to the SLS package by
itself: 

* (required) definition of the precision of the SUNDIALS type
  ``realtype``. One of the following lines must be present:

  .. code-block:: c
 
     #define SUNDIALS_DOUBLE_PRECISION 1
     #define SUNDIALS_SINGLE_PRECISION 1
     #define SUNDIALS_EXTENDED_PRECISION 1

* (optional) use of generic math functions: 

  .. code-block:: c

     #define SUNDIALS_USE_GENERIC_MATH 1

The ``sundials_types.h`` header file defines the SUNDIALS ``realtype``
and ``booleantype`` types and the macro ``RCONST``, while the
``sundials_math.h`` header file is needed for the macros ``SUNMIN`` and
``SUNMAX``, and the function ``SUNRabs``.

The files listed above for either module can be extracted from the
SUNDIALS ``srcdir`` and compiled by themselves into a separate library
or into a larger user code.



SlsMat
--------------------

SUNDIALS supports operations with compressed-sparse-column (CSC) and
compressed-sparse-row (CSR) matrix formats.  For convenience, integer
sparse matrix identifiers are defined as:\\

.. code-block:: c

   #define CSC_MAT 0
   #define CSR_MAT 1

The type :c:type:`SlsMat`, defined in ``sundials_sparse.h`` is a
pointer to a structure defining generic CSC and CSR matrix formats,
and is used with all linear solvers in the SLS family:  

.. c:type:: SlsMat

   .. code-block:: c

      typedef struct _SlsMat {
        int M;
        int N;
        int NNZ;
	int NP;
        realtype *data;
	int sparsetype;
	int *indexvals;
	int *indexptrs;
        int **rowvals;
        int **colptrs;
	int **colvals;
	int **rowptrs;
      } *SlsMat;

The fields of this structure are as follows (see Figure
:ref:`SlsMat Diagram <SLS_figure>`
for a diagram of the underlying compressed-sparse-column
representation in a sparse matrix of type :c:type:`SlsMat`).  Note that a
sparse matrix of type :c:type:`SlsMat` need not be square.

:M: -- number of rows
:N: --  number of columns
:NNZ: -- maximum number of nonzero entries in the matrix (allocated
   length of **data** and **rowvals** arrays)
:NP: -- number of index pointers (e.g. number of column pointers for 
    CSC matrix). For CSC matrices $NP=N$, and for CSR matrices $NP=M$. This 
    value is set automatically based on the input for **sparsetype**.
:data: -- pointer to a contiguous block of ``realtype`` variables (of
   length **NNZ**), containing the values of the nonzero entries in the
   matrix.
:sparsetype: -- type of the sparse matrix (``CSC_MAT`` or ``CSR_MAT``)
:indexvals: -- pointer to a contiguous block of ``int`` variables (of
   length **NNZ**), containing the row indices (if CSC) or column
   indices (if CSR) of each nonzero matrix entry held in **data**.
:indexptrs: -- pointer to a contiguous block of ``int`` variables (of
  length **NP+1**).  For CSC matrices each entry provides the index of
  the first column entry into the **data** and **indexvals** arrays,
  e.g. if **indexptr[3]=7**, then the first nonzero entry in the
  fourth column of the matrix is located in **data[7]**, and is
  located in row **indexvals[7]** of the matrix.  The last entry
  contains the total number of nonzero values in the matrix and hence
  points just past the end of the active data in the **data** and
  **indexvals** arrays.  For CSR matrices, each entry provides the
  index of the first row intry into the **data** and **indexvals**
  arrays.

.. _SLS_figure:

.. figure:: figs/cscmat.png

   Diagram of the storage for a compressed-sparse-column matrix of
   type :c:type:`SlsMat`: Here ``A`` is an :math:`M \times N` sparse
   matrix of type :c:type:`SlsMat` with storage for up to ``NNZ``
   nonzero entries (the allocated length of both ``data`` and
   ``indexvals``).  The entries in ``indexvals`` may assume values from
   ``0`` to ``M-1``, corresponding to the row index (zero-based) of
   each nonzero value.  The entries in ``data`` contain the values of
   the nonzero entries, with the row ``i``, column ``j`` entry of
   ``A`` (again, zero-based) denoted as ``A(i,j)``.  The ``indexptrs``
   array contains ``N+1`` entries; the first ``N`` denote the starting
   index of each column within the ``indexvals`` and ``data`` arrays,
   while the final entry points one past the final nonzero entry.
   Here, although ``NNZ`` values are allocated, only ``nz`` are
   actually filled in; the greyed-out portions of ``data`` and
   ``indexvals`` indicate extra allocated space.


The following pointers are added to the ``SlsMat`` type for user
convenience, to provide a more intuitive interface to the CSC and CSR
sparse matrix data structures. They are set automatically by the
:c:func:`SparseNewMat()` function, based on the sparse matrix storage
type.  

:rowvals: -- pointer to ``indexvals`` when ``sparsetype`` is ``CSC_MAT``,
    otherwise set to ``NULL``.
:colptrs: -- pointer to ``indexptrs`` when ``sparsetype`` is ``CSC_MAT``,
    otherwise set to ``NULL``.
:colvals: -- pointer to ``indexvals`` when ``sparsetype`` is ``CSR_MAT``,
    otherwise set to ``NULL``.
:rowptrs: -- pointer to ``indexptrs`` when ``sparsetype`` is ``CSR_MAT``,
    otherwise set to ``NULL``.

For example, the :math:`5\times 4` CSC matrix

.. math::

   \left[\begin{array}{cccc} 
      0 & 3 & 1 & 0\\
      3 & 0 & 0 & 2\\
      0 & 7 & 0 & 0\\
      1 & 0 & 0 & 9\\
      0 & 0 & 0 & 5
   \end{array}\right]

could be stored in a :c:type:`SlsMat` structure as either

.. code-block:: c

   M = 5;
   N = 4;
   NNZ = 8;
   NP = N;
   data = {3.0, 1.0, 3.0, 7.0, 1.0, 2.0, 9.0, 5.0};
   sparsetype = CSC_MAT;
   indexvals = {1, 3, 0, 2, 0, 1, 3, 4};
   indexptrs = {0, 2, 4, 5, 8};
   rowvals = &indexvals;
   colptrs = &indexptrs;
   colvals = NULL;
   rowptrs = NULL;

or 

.. code-block:: c

   M = 5;
   N = 4;
   NNZ = 10;
   NP = N;
   data = {3.0, 1.0, 3.0, 7.0, 1.0, 2.0, 9.0, 5.0, *, *};
   sparsetype = CSC_MAT;
   indexvals = {1, 3, 0, 2, 0, 1, 3, 4, *, *};
   indexptrs = {0, 2, 4, 5, 8};
   rowvals = &indexvals;
   colptrs = &indexptrs;
   colvals = NULL;
   rowptrs = NULL;

where the first has no unused space, and the second has additional
storage (the entries marked with ``*`` may contain any values).  Note
in both cases that the final value in ``indexptrs`` is ``8``.  The work
associated with operations on the sparse matrix is proportional to
this value and so one should use the best understanding of the number
of nonzeros here.




Functions in the SPARSE module
-------------------------------------------

The SPARSE module defines functions that act on sparse matrices of
type :c:type:`SlsMat`.  For full details, see the header file
``sundials_sparse.h``.


.. c:function:: SlsMat SparseNewMat(int M, int N, int NNZ, int sparsetype)
   
   Allocates a :c:type:`SlsMat` sparse matrix having *M* rows, *N*
   columns, and storage for *NNZ* nonzero entries and *sparsetype*
   storage type (``CSC_MAT`` or ``CSR_MAT``).

.. c:function:: SlsMat SparseFromDenseMat(DlsMat A)

   Converts a dense matrix of type :c:type:`DlsMat` into a CSC
   matrix of type :c:type:`SlsMat` by retaining only the nonzero
   values of the dense matrix.

.. c:function:: void SparseDestroyMat(SlsMat A)

   Frees memory for a :c:type:`SlsMat` matrix.

.. c:function:: void SparseSetMatToZero(SlsMat A)

   Zeros out a :c:type:`SlsMat` matrix (but retains its storage).

.. c:function:: void SparseCopyMat(SlsMat A, SlsMat B)

   Copies one sparse matrix to another.  It is assumed that the
   matrices have the same row/column dimensions and storage type.  If
   *B* has insufficient storage to hold all the nonzero entries of
   *A*, the data and index arrays in *B* are reallocated to match
   those in *A*.

.. c:function:: void SparseScaleMat(realtype c, SlsMat A)

   Scales every element in the sparse matrix A by the by the
   scalar c.

.. c:function:: void SparseAddIdentityMat(SlsMat A)

   Increments a sparse matrix by the identity matrix.  If *A* is not
   square, only the existing diagonal values are incremented.  Resizes
   the ``data`` and ``indexvals`` arrays of *A* to allow for new nonzero
   entries on the diagonal.

.. c:function:: int SparseAddMat(SlsMat A, SlsMat B)

   Adds two sparse matrices: :math:`A = A+B`.  Resizes the data arrays
   of *A* upon completion to exactly match the nonzero storage for
   the result.  Upon successful completion, the return value is zero;
   otherwise -1 is returned.  It is assumed that both matrices have
   the same size and storage type.

.. c:function:: void SparseReallocMat(SlsMat A)

   This function eliminates unused storage in *A* by reallocating
   the internal ``data`` and ``indexvals`` arrays to contain
   ``indexptrs[N]`` nonzeros.

.. c:function:: int SparseMatvec(SlsMat A, realtype *x, realtype *y)

   Computes the sparse matrix-vector product, :math:`y=Ax`.  If the
   ``SlsMat`` *A* is a sparse matrix of dimension :math:`M\times N`,
   then it is assumed that *x* is a ``realtype`` array of  length
   :math:`N`, and *y* is a ``realtype`` array of length
   :math:`M`. Upon successful completion, the return value is zero;
   otherwise -1 is returned. 

.. c:function:: void SparsePrintMat(DlsMat A)

   Prints a :c:type:`SlsMat` matrix to standard output.





The KLU solver
-------------------------------------------

KLU is a sparse matrix factorization and solver library written by Tim
Davis ([KLU]_, [DP2010]_).   KLU has a symbolic factorization routine
that computes the permutation of the linear system matrix to block
triangular form and the permutations that will pre-order the diagonal
blocks (the only ones that need to be factored) to reduce fill-in
(using AMD, COLAMD, CHOLAMD, natural, or an ordering given by the
user).  Note that SUNDIALS uses the COLAMD ordering by default with
KLU. 

KLU breaks the factorization into two separate parts.  The first is a
symbolic factorization and the second is a numeric factorization that
returns the factored matrix along with final pivot information.  KLU
also has a refactor routine that can be called instead of the numeric
factorization.  This routine will reuse the pivot information.  This
routine also returns diagnostic information that a user can examine to
determine if numerical stability is being lost and a full numerical
factorization should be done instead of the refactor.

The KLU interface in SUNDIALS will perform the symbolic factorization
once.  It then calls the numerical factorization once and will call
the refactor routine until estimates of the numerical conditioning
suggest a new factorization should be completed.  The KLU interface
also has a ``ReInit`` routine that can be used to force a full
refactorization at the next solver setup call.

In order to use the SUNDIALS interface to KLU, it is
assumed that KLU has been installed on the system prior to
installation of SUNDIALS, and that SUNDIALS has been configured
appropriately to link with KLU (see :ref:`Installation` for details).

Designed for serial calculations only, KLU is supported for
calculations employing SUNDIALS' serial or shared-memory parallel
``N_Vector`` modules (see :ref:`NVectors.NVSerial`,
:ref:`NVectors.OpenMP` and :ref:`NVectors.Pthreads`).



The SuperLU_MT solver
-------------------------------------------

SuperLU_MT is a threaded sparse matrix factorization and solver
library written by X. Sherry Li ([SuperLUMT]_, [L2005]_, [DGL1999]_).
The package performs matrix factorization using threads to enhance
efficiency in shared memory parallel environments.  It should be noted
that threads are only used in the factorization step.

In order to use the SUNDIALS interface to SuperLU_MT, it is assumed
that SuperLU_MT has been installed on the system prior to installation
of SUNDIALS, and that SUNDIALS has been configured appropriately to
link with SuperLU_MT (see :ref:`Installation` for details).

Designed for serial and threaded calculations only, SuperLU_MT is
supported for calculations employing SUNDIALS' serial or shared-memory
parallel ``N_Vector`` modules (see :ref:`NVectors.NVSerial`,
:ref:`NVectors.OpenMP` and :ref:`NVectors.Pthreads`).
