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

:tocdepth: 3


.. _LinearSolvers.DLS:

The DLS modules: DENSE and BAND
========================================

The files comprising the DENSE generic linear solver, and their
locations in the SUNDIALS ``srcdir``, are as follows:

* header files (located in ``srcdir/include/sundials``):

  ``sundials_direct.h``, ``sundials_dense.h``, ``sundials_types.h``,
  ``sundials_math.h``, ``sundials_config.h`` 

* source files (located in ``srcdir/src/sundials``):

  ``sundials_direct.c``, ``sundials_dense.c``, ``sundials_math.c``

The files comprising the BAND generic linear solver are as follows: 

* header files (located in ``srcdir/include/sundials``):

  ``sundials_direct.h``, ``sundials_band.h``, ``sundials_types.h``,
  ``sundials_math.h``, ``sundials_config.h`` 

* source files (located in ``srcdir/src/sundials``):

  ``sundials_direct.c``, ``sundials_band.c``, ``sundials_math.c``

Only two of the preprocessing directives in the header file
``sundials_config.h`` are relevant to the DENSE and BAND packages by
themselves.

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
``sundials_math.h`` header file is needed for the macros ``SUNMIN``
and ``SUNMAX``, and the function ``SUNRabs``.

The files listed above for either module can be extracted from the
SUNDIALS ``srcdir`` and compiled by themselves into a separate library
or into a larger user code.



DlsMat
--------------------

The type :c:type:`DlsMat`, defined in ``sundials_direct.h`` is a
pointer to a structure defining a generic matrix, and is used with all
linear solvers in the DLS family: 

.. c:type:: DlsMat

   .. code-block:: c

      typedef struct _DlsMat {
        int type;
        sunindextype M;
        sunindextype N;
        sunindextype ldim;
        sunindextype mu;
        sunindextype ml;
        sunindextype s_mu;
        realtype *data;
        sunindextype ldata;
        realtype **cols;
      } *DlsMat;

For the DENSE module, the relevant fields of this structure are as
follows. Note that a dense matrix of type :c:type:`DlsMat` need not be
square. 

:type: -- ``SUNDIALS_DENSE`` (=1)
:M: -- number of rows
:N: --  number of columns
:ldim: -- leading dimension (:math:`\ge M`)
:data: -- pointer to a contiguous block of ``realtype`` variables 
:ldata: -- length of the data array (:math:`= ldim*N`). The
  ``(i,j)`` element of a dense matrix ``A`` of type ``DlsMat`` (with
  :math:`0 \le i < M` and :math:`0 \le j < N`) is given by the
  expression ``(A->data)[0][j*M+i]`` 
:cols: -- array of pointers. ``cols[j]`` points to the first element
  of the ``j``-th column of the matrix in the array data. The
  ``(i,j)`` element of a dense matrix ``A`` of type ``DlsMat`` (with
  :math:`0 \le i < M` and :math:`0 \le j < N`) is given by the
  expression ``(A->cols)[j][i]`` 

For the BAND module, the relevant fields of this structure are as
follows (see Figure :ref:`DLS Diagram <DLS_figure>` for a diagram of
the underlying data representation in a banded matrix of type
:c:type:`DlsMat`). Note that only square band matrices are allowed.

:type: -- ``SUNDIALS_BAND`` (=2)
:M: -- number of rows
:N: -- number of columns (:math:`N = M`)
:mu: -- upper half-bandwidth, :math:`0 \le mu < min(M,N)`
:ml: -- lower half-bandwidth, :math:`0 \le ml < min(M,N)`
:s_mu: -- storage upper bandwidth, :math:`mu \le s\_mu < N`. The LU
   decomposition routine writes the LU factors into the storage for
   :math:`A`. The upper triangular factor :math:`U`, however, may
   have an upper bandwidth as big as :math:`min(N-1,mu+ml)` because
   of partial pivoting. The ``s_mu`` field holds the upper
   half-bandwidth allocated for :math:`A`. 
:ldim: -- leading dimension (:math:`ldim \ge s\_mu`)
:data: -- pointer to a contiguous block of ``realtype``
   variables. The elements of a banded matrix of type
   :c:type:`DlsMat` are stored columnwise (i.e. columns are stored
   one on top of the other in memory). Only elements within the
   specified half-bandwidths are stored. ``data`` is a pointer to
   ``ldata`` contiguous locations which hold the elements within the
   band of :math:`A`. 
:ldata: -- length of the ``data`` array (:math:`= ldim*(s\_mu+ml+1)`)
:cols: -- array of pointers. ``cols[j]`` is a pointer to the
   uppermost element within the band in the ``j``-th column. This
   pointer may be treated as an array indexed from ``s_mu-mu`` (to
   access the uppermost element within the band in the ``j``-th
   column) to ``s_mu+ml`` (to access the lowest element within the
   band in the ``j``-th column). Indices from 0 to ``s_mu-mu-1`` give
   access to extra storage elements required by the LU decomposition
   function. Finally, ``cols[j][i-j+s_mu]`` is the ``(i,j)``-th
   element, :math:`j-mu \le i \le j+ml`.


.. _DLS_figure:

.. figure:: figs/bandmat.png

   DLS Diagram: Storage for a banded matrix of type :c:type:`DlsMat`. Here
   ``A`` is an :math:`N \times N` band matrix of type :c:type:`DlsMat`
   with upper and lower half-bandwidths ``mu`` and ``ml``,
   respectively. The rows and columns of ``A`` are numbered from
   :math:`0` to :math:`N-1` and the ``(i,j)``-th element of ``A`` is
   denoted ``A(i,j)``. The greyed out areas of the underlying
   component storage are used by the BandGBTRF and BandGBTRS routines.





Accessor macros for the DLS modules
-------------------------------------------

The macros below allow a user to efficiently access individual matrix
elements without writing out explicit data structure references and
without knowing too much about the underlying element storage.  The
only storage assumption needed is that elements are stored columnwise
and that a pointer to the j-th column of elements can be obtained via
the :c:macro:`DENSE_COL` or :c:macro:`BAND_COL` macros. Users should use these
macros whenever possible. 

The following two macros are defined by the DENSE module to provide
access to data in the :c:type:`DlsMat` type:

.. c:macro:: DENSE_ELEM

   **Usage:** ``DENSE_ELEM(A,i,j) = a_ij;``  or  ``a_ij = DENSE_ELEM(A,i,j);``

   This macro references the :math:`(i,j)`-th element of the :math:`M \times N`
   :c:type:`DlsMat` :math:`A`, :math:`0 \le i < M` , :math:`0 \le j < N`.


.. c:macro:: DENSE_COL

   **Usage:** ``col_j = DENSE_COL(A,j);``

   This macro references the :math:`j`-th column of the :math:`M \times N`
   :c:type:`DlsMat` :math:`A`, :math:`0 \le j < N`. The type of the
   expression ``DENSE_COL(A,j)`` is ``realtype *`` . After the 
   assignment in the usage above, ``col_j`` may be treated as an
   array indexed from 0 to :math:`M-1`. The :math:`(i,j)`-th
   element of :math:`A` is referenced by ``col_j[i]``.



The following three macros are defined by the BAND module to provide
access to data in the :c:type:`DlsMat` type:

.. c:macro:: BAND_ELEM

   **Usage:** ``BAND_ELEM(A,i,j) = a_ij;``  or  ``a_ij =
   BAND_ELEM(A,i,j);``

   This macro references the :math:`(i,j)`-th element of the :math:`N \times N`
   band matrix :math:`A`, where :math:`0 \le i`, :math:`j \le N-1`.
   The location :math:`(i,j)` should further satisfy :math:`j-`
   ``(A->mu)`` :math:`\le i \le j+` ``(A->ml)``.

.. c:macro:: BAND_COL

   **Usage:** ``col_j = BAND_COL(A,j);``

   This macro references the diagonal element of the :math:`j`-th column of the
   :math:`N \times N` band matrix :math:`A`, :math:`0 \le j \le
   N-1`. The type of the expression ``BAND_COL(A,j)`` is
   ``realtype *``. The pointer returned by the call ``BAND_COL(A,j)``
   can be treated as an array which is indexed from ``-(A->mu)`` to
   ``(A->ml)``. 

.. c:macro:: BAND_COL_ELEM

   **Usage:** ``BAND_COL_ELEM(col_j,i,j) = a_ij;``  or  ``a_ij =
   BAND_COL_ELEM(col_j,i,j);`` 

   This macro references the :math:`(i,j)`-th entry of the band matrix
   :math:`A` when used in conjunction with :c:macro:`BAND_COL` to reference
   the :math:`j`-th column through ``col_j``. The index :math:`(i,j)`
   should satisfy :math:`j-` ``(A->mu)`` :math:`\le i \le j+` ``(A->ml)``.




Functions in the DENSE module
-------------------------------------------

The DENSE module defines two sets of functions with corresponding
names. The first set contains functions (with names starting with a
capital letter) that act on dense matrices of type :c:type:`DlsMat`. The
second set contains functions (with names starting with a lower case
letter) that act on matrices represented as simple arrays.

The following functions for DlsMat dense matrices are available in the
DENSE package. For full details, see the header files
``sundials_direct.h`` and ``sundials_dense.h``.


.. c:function:: DlsMat NewDenseMat(sunindextype M, sunindextype N)
   
   Allocates a :c:type:`DlsMat` dense matrix.

.. c:function:: void DestroyMat(DlsMat A)

   Frees memory for a :c:type:`DlsMat` matrix

.. c:function:: void PrintMat(DlsMat A)

   Prints a :c:type:`DlsMat` matrix to standard output.

.. c:function:: sunindextype* NewIndexArray(sunindextype N) 
   
   Allocates an array of ``sunindextype`` integers for use as pivots with
   :c:func:`DenseGETRF()` and :c:func:`DenseGETRS()`. 

.. c:function:: int* NewIntArray(int N)

   Allocates an array of ``int`` integers for use as pivots with the
   LAPACK dense solvers.

.. c:function:: realtype* NewRealArray(sunindextype N)
   
   Allocates an array of type ``realtype`` for use as right-hand side
   with :c:func:`DenseGETRS()`.

.. c:function:: void DestroyArray(void* p)

   Frees memory for an array.

.. c:function:: void SetToZero(DlsMat A)

   Loads a matrix with zeros.

.. c:function:: void AddIdentity(DlsMat A)

   Increments a square matrix by the identity matrix.

.. c:function:: void DenseCopy(DlsMat A, DlsMat B)

   Copies one dense matrix to another.

.. c:function:: void DenseScale(realtype c, DlsMat A)

   Scales a dense matrix by a scalar.

.. c:function:: sunindextype DenseGETRF(DlsMat A, sunindextype* p)

   LU factorization with partial pivoting of a dense matrix.

.. c:function:: sunindextype denseGETRF(realtype** a, sunindextype m, sunindextype n, sunindextype* p)

   Solves :math:`Ax = b` using LU factorization (for square matrices
   :math:`A`), using the factorization resulting from :c:func:`DenseGETRF()`.

.. c:function:: sunindextype DensePOTRF(DlsMat A)

   Cholesky factorization of a real symmetric positive definite dense matrix.

.. c:function:: void DensePOTRS(DlsMat A, realtype* b)

   Solves :math:`Ax = b` using the Cholesky factorization of :math:`A`
   resulting from a call to :c:func:`DensePOTRF()`.

.. c:function:: int DenseGEQRF(DlsMat A, realtype* beta, realtype* wrk)

   QR factorization of an :math:`m \times n` dense matrix, with :math:`m \ge n`.

.. c:function:: int DenseORMQR(DlsMat A, realtype* beta, realtype* vn, realtype* vm, realtype* wrk)

   Computes the product :math:`w = Qv`, with :math:`Q` calculated
   using :c:func:`DenseGEQRF()`.  

.. c:function:: int DenseMatvec(DlsMat A, realtype* x, realtype* y)

   Computes the product :math:`y = Ax`, where it is assumed that
   :math:`x` has length equal to the number of columns in the matrix
   :math:`A`, and :math:`y` has length equal to the number of rows in
   the matrix :math:`A`.



The following functions for small dense matrices are available in the
DENSE package.  These functions primarily replicate those defined above
for :c:type:`DlsMat` dense matrices, but act on the individual data
arrays outside of the :c:type:`DlsMat` structure:

.. c:function:: realtype** newDenseMat(sunindextype m, sunindextype n)

   Allocates storage for an :math:`m \times n` dense matrix. It
   returns a pointer to the newly allocated storage if successful. If
   the memory request cannot be satisfied, then the function returns
   ``NULL``.  The underlying type of the dense matrix returned is
   ``realtype**``. If we allocate a dense matrix ``realtype** a`` by
   ``a = newDenseMat(m,n)``, then ``a[j][i]`` references the row ``i``,
   column ``j`` element of the matrix ``a``, :math:`0 \le i < m`,
   :math:`0 \le j < n`, and ``a[j]`` is a pointer to the first element
   in the :math:`j`-th column of ``a``. The location ``a[0]`` contains
   a pointer to :math:`m \times n` contiguous locations which contain
   the elements of ``a``.

.. c:function:: void destroyMat(realtype** a)

   Frees the dense matrix *a* allocated by :c:func:`newDenseMat()`.

.. c:function:: sunindextype* newIndexArray(sunindextype n)

   Allocates an array of *n* integers of ``sunindextype`` type.  It
   returns a pointer to the first element in the array if
   successful. It returns ``NULL`` if the memory request could not be
   satisfied.  

.. c:function:: int* newIntArray(int n)

   Allocates an array of *n* integers of type ``int``.  It returns a
   pointer to the first element in the array if successful. It returns
   ``NULL`` if the memory request could not be satisfied. 

.. c:function:: realtype* newRealArray(sunindextype m)

   Allocates an array of *n* ``realtype`` values. It returns a pointer
   to the first element in the array if successful. It returns
   ``NULL`` if the memory request could not be satisfied. 

.. c:function:: void destroyArray(void* v)

   Frees the array *v* allocated by :c:func:`newIndexArray()`,
   :c:func:`newIntArray()`, or :c:func:`newRealArray()`. 

.. c:function:: void denseCopy(realtype** a, realtype** b, sunindextype m, sunindextype n)

   Copies the :math:`m \times n` dense matrix *a* into the :math:`m
   \times n` dense matrix *b*. 

.. c:function:: void denseScale(realtype c, realtype** a, sunindextype m, sunindextype n)

   Scales every element in the :math:`m \times n` dense matrix *a* by
   the scalar *c*. 

.. c:function:: void denseAddIdentity(realtype** a, sunindextype n)

   Increments the square :math:`n \times n` dense matrix *a* by the
   identity matrix :math:`I_n`.

.. c:function:: sunindextype denseGETRF(realtype** a, sunindextype m, sunindextype n, sunindextype* p)

   Factors the :math:`m \times n` dense matrix *a*, using Gaussian
   elimination with row pivoting. It overwrites the elements of *a*
   with its LU factors and keeps track of the pivot rows chosen in the
   pivot array *p*.

   A successful LU factorization leaves the matrix *a* and the pivot
   array *p* with the following information:

   1. ``p[k]`` contains the row number of the pivot element chosen at
      the beginning of elimination step :math:`k, k = 0, 1, \ldots,
      n-1`.

   2. If the unique LU factorization of *a* is given by :math:`P a =
      LU`, where :math:`P` is a permutation matrix, :math:`L` is a
      :math:`m \times n` lower trapezoidal matrix with all diagonal
      elements equal to 1, and :math:`U` is a :math:`n \times n` upper
      triangular matrix, then the upper triangular part of *a*
      (including its diagonal) contains :math:`U` and the strictly
      lower trapezoidal part of *a* contains the multipliers,
      :math:`I-L`. If *a* is square, :math:`L` is a unit lower
      triangular matrix. 

      :c:func:`denseGETRF()` returns 0 if successful. Otherwise it
      encountered a zero diagonal element during the factorization,
      indicating that the matrix a does not have full column rank. In
      this case it returns the column index (numbered from one) at
      which it encountered the zero.

.. c:function:: void denseGETRS(realtype** a, sunindextype n, sunindextype* p, realtype* b)

   Solves the :math:`n \times n` linear system :math:`ax = b`. It
   assumes that *a* (of size :math:`n \times n`) has been LU-factored
   and the pivot array *p* has been set by a successful call to
   :c:func:`denseGETRF()`. The solution *x* is written into the *b*
   array. 

.. c:function:: sunindextype densePOTRF(realtype** a, sunindextype m)

   Calculates the Cholesky decomposition of the :math:`m \times m`
   dense matrix *a*, assumed to be symmetric positive definite.  Only
   the lower triangle of *a* is accessed and overwritten with the
   Cholesky factor.

.. c:function:: void densePOTRS(realtype** a, sunindextype m, realtype* b)

   Solves the :math:`m \times m` linear system :math:`ax = b`.  It
   assumes that the Cholesky factorization of *a* has been calculated
   in the lower triangular part of *a* by a successful call to
   :c:func:`densePOTRF(m)`. 

.. c:function:: int denseGEQRF(realtype** a, sunindextype m, sunindextype n, realtype* beta, realtype* v)

   Calculates the QR decomposition of the :math:`m \times n` matrix
   *a* (:math:`m \ge n`) using Householder reflections.  On exit, the
   elements on and above the diagonal of *a* contain the :math:`n
   \times n` upper triangular matrix :math:`R`; the elements below the
   diagonal, with the array *beta*, represent the orthogonal matrix
   :math:`Q` as a product of elementary reflectors. The real array
   *wrk*, of length *m*, must be provided as temporary workspace. 

.. c:function:: int denseORMQR(realtype** a, sunindextype m, sunindextype n, realtype* beta, realtype* v, realtype* w, realtype* wrk)

   Calculates the product :math:`w = Qv` for a given vector *v* of
   length *n*, where the orthogonal matrix :math:`Q` is encoded in the
   :math:`m \times n` matrix *a* and the vector *beta* of length *n*,
   after a successful call to :c:func:`denseGEQRF()`. The real array
   *wrk*, of length *m*, must be provided as temporary workspace.

.. c:function:: int denseMatvec(realtype **a, realtype* x, realtype* y, sunindextype m, sunindextype n)

   Computes the product :math:`y = ax`, for an :math:`m\times n`
   matrix :math:`a`, where it is assumed that the vector :math:`x` has
   length :math:`n` and the vector :math:`y` has length :math:`m`.



Functions in the BAND module
-------------------------------------------

The BAND module defines two sets of functions with corresponding
names. The first set contains functions (with names starting with a
capital letter) that act on band matrices of type :c:type:`DlsMat`. The
second set contains functions (with names starting with a lower case
letter) that act on matrices represented as simple arrays.

The following functions for :c:type:`DlsMat` banded matrices are
available in the BAND package. For full details, see the header files
``sundials_direct.h`` and ``sundials_band.h``.  A number of these are
shared with routines from the DENSE package, but are listed again here
for completeness.


.. c:function:: DlsMat NewBandMat(sunindextype N, sunindextype mu, sunindextype ml, sunindextype smu)

   Allocates a :c:type:`DlsMat` band matrix

.. c:function:: void DestroyMat(DlsMat A)

   Frees memory for a :c:type:`DlsMat` matrix

.. c:function:: void PrintMat(DlsMat A)

   Prints a :c:type:`DlsMat` matrix to standard output.

.. c:function:: sunindextype* NewIndexArray(sunindextype N) 
   
   Allocates an array of ``sunindextype`` integers for use as pivots with
   :c:func:`BandGBRF()` and :c:func:`BandGBRS()`. 

.. c:function:: int* NewIntArray(int N)

   Allocates an array of ``int`` integers for use as pivots with the
   LAPACK band solvers.

.. c:function:: realtype* NewRealArray(sunindextype N)
   
   Allocates an array of type ``realtype`` for use as right-hand side
   with :c:func:`BandGBRS()`.

.. c:function:: void DestroyArray(void* p)

   Frees memory for an array.

.. c:function:: void SetToZero(DlsMat A)

   Loads a matrix with zeros.

.. c:function:: void AddIdentity(DlsMat A)

   Increments a square matrix by the identity matrix.

.. c:function:: void BandCopy(DlsMat A, DlsMat B, sunindextype copymu, sunindextype copyml)

   Copies one band matrix to another.

.. c:function:: void BandScale(realtype c, DlsMat A)

   Scales a band matrix by a scalar.

.. c:function:: sunindextype BandGBTRF(DlsMat A, sunindextype* p)

   LU factorization with partial pivoting.

.. c:function:: void BandGBTRS(DlsMat A, sunindextype* p, realtype* b)

   Solves :math:`Ax = b` using LU factorization resulting from
   :c:func:`BandGBTRF()`. 

.. c:function:: int BandMatvec(DlsMat A, realtype* x, realtype* y)

   Computes the product :math:`y = Ax`, where it is assumed that
   :math:`x` and :math:`y` have length equal to the number of rows in
   the square band matrix :math:`A`.



The following functions for small band matrices are available in the
BAND package.  These functions primarily replicate those defined above
for :c:type:`DlsMat` banded matrices, but act on the individual data arrays
outside of the :c:type:`DlsMat` structure:

.. c:function:: realtype** newBandMat(sunindextype n, sunindextype smu, sunindextype ml)

   Allocates storage for a :math:`n \times n` band matrix with lower
   half-bandwidth *ml*. 

.. c:function:: void destroyMat(realtype** a)
 
   Frees the band matrix *a* allocated by :c:func:`newBandMat()`.

.. c:function:: sunindextype* newIndexArray(sunindextype n)

   Allocates an array of *n* integers of type ``sunindextype``. It returns
   a pointer to the first element in the array if successful.  It
   returns ``NULL`` if the memory request could not be satisfied. 

.. c:function:: int* newIntArray(int n)

   Allocates an array of *n* integers of type ``int``. It returns a
   pointer to the first element in the array if successful. It returns
   ``NULL`` if the memory request could not be satisfied. 

.. c:function:: realtype* newRealArray(sunindextype m)

   Allocates an array of *n* ``realtype`` values. It returns a pointer
   to the first element in the array if successful. It returns
   ``NULL`` if the memory request could not be satisfied. 

.. c:function:: void destroyArray(void* v)

   Frees the array *v* allocated by :c:func:`newIndexArray()`,
   :c:func:`newIntArray()`, or :c:func:`newRealArray()`. 

.. c:function:: void bandCopy(realtype** a, realtype** b, sunindextype n, sunindextype a_smu, sunindextype b_smu, sunindextype copymu, sunindextype copyml)

   Copies the :math:`n \times n` band matrix *a* into the :math:`n
   \times n` band matrix *b*. 

.. c:function:: void bandScale(realtype c, realtype** a, sunindextype n, sunindextype mu, sunindextype ml, sunindextype smu)

   Scales every element in the :math:`n \times n` band matrix *a* by
   *c*. 

.. c:function:: void bandAddIdentity(realtype** a, sunindextype n, sunindextype smu)

   Increments the :math:`n \times n` band matrix *a* by the identity
   matrix. 

.. c:function:: sunindextype bandGBTRF(realtype** a, sunindextype n, sunindextype mu, sunindextype ml, sunindextype smu, sunindextype* p)
 
   Factors the :math:`n \times n` band matrix *a*, using Gaussian
   elimination with row pivoting. It overwrites the elements of *a*
   with its LU factors and keeps track of the pivot rows chosen in the
   pivot array *p*. 

.. c:function:: void bandGBTRS(realtype** a, sunindextype n, sunindextype smu, sunindextype ml, sunindextype* p, realtype* b)

   Solves the :math:`n \times n` linear system :math:`ax = b`. It
   assumes that *a* (of size :math:`n \times n`) has been LU-factored
   and the pivot array *p* has been set by a successful call to
   :c:func:`bandGETRF()`. The solution *x* is written into the *b*
   array. 

.. c:function:: int bandMatvec(realtype **a, realtype* x, realtype* y, sunindextype n, sunindextype mu, sunindextype ml, sunindextype smu)

   Computes the product :math:`y = ax`, for an :math:`n\times n`
   square band matrix :math:`a`, having band structure as allocated by
   the parameters *mu*, *ml* and *smu*, and where it is assumed that
   :math:`x` and :math:`y` have length :math:`n`.

