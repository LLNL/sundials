..
   Programmer(s): Mustafa Aggul @ SMU
   ----------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2002-2025, Lawrence Livermore National Security
   and Southern Methodist University.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   ----------------------------------------------------------------

.. _SUNDomEigEst.ARNI:

The SUNDomEigEst_ARNI Module
======================================

The SUNDomEigEst_ARNI implementation of the ``SUNDomEigEstimator`` class performs
the Arnoldi Iteration (ArnI) method :cite:p:`arnoldi51`; this is an iterative dominant
eigenvalue estimator that is designed to be compatible with any ``N_Vector``
implementation that supports a minimal subset of operations (:c:func:`N_VClone()`,
:c:func:`N_VCloneVectorArray()`, :c:func:`N_VDotProd()`, :c:func:`N_VLinearSum()`,
:c:func:`N_VScale()`, :c:func:`N_VDestroy()`, and :c:func:`N_VDestroyVectorArray()`).

ArnI is particularly effective for large, sparse matrices where only the dominant
eigenvalue is needed.  It constructs an orthonormal basis of the Krylov subspace

.. math::

   \mathcal{K}_m(A, \mathbf{v}) = \text{span}\{\mathbf{v}, A \mathbf{v}, A^2 \mathbf{v}, \dots, A^{m-1} \mathbf{v}\}

using the Gram-Schmidt process.  The matrix :math:`A` is projected onto this subspace
to form a small upper Hessenberg matrix :math:`H_m`.  The eigenvalues of :math:`H_m`
approximate some of the eigenvalues of :math:`A`; the dominant eigenvalue of :math:`A` is
well-approximated by the dominant eigenvalue of :math:`H_m`.

ArnI works for the matrices that have a **complex** dominant eigenvalue.  It supports
estimations with a user specified fixed dimension of Krylov subspaces, at least 3.  While
these choice requires a prefixed amount of memory (depending on the dimension), it strictly
determines how good an estimation is.  To improve the estimation accuracy, we found preprocessing
with :c:func:`SUNDomEigEstPreProcess` particularly useful.  This operation is free from any
additional memory requirement and explained below.

The matrix :math:`A` is not required explicitly; only routines that provide :math:`A`
as operator is required.


.. _SUNDomEigEst.ARNI.Usage:

SUNDomEigEst_ARNI Usage
-----------------------

The header file to be included when using this module is ``sundomeigest/sundomeigest_pi.h``.
The SUNDomEigEst_PI module is accessible from all SUNDIALS solvers *without* linking to the
``libsundials_sundomeigestpi`` module library.

The header file to be included when using this module is ``sundomeigest/sundomeigest_arni.h``.
The SUNDomEigEst_ARNI module is accessible from all SUNDIALS solvers *without* linking to the
``libsundials_sundomeigestarni`` module library after enabling LAPACK package.
This LAPACK dependence is limited to the eigenvalue estimation of the Hessenberg matrix with
``dgeev_/sgeev_`` functions.

The module SUNDomEigEst_ARNI provides the following user-callable routines:


.. c:function:: SUNDomEigEstimator SUNDomEigEst_ArnI(N_Vector q, sunindextype krydim, SUNContext sunctx)

   This constructor function creates and allocates memory for an ARNI
   ``SUNDomEigEstimator``.

   **Arguments:**
      * *q* -- a template vector.
      * *krydim* -- the dimension of the Krylov subspaces.
      * *sunctx* -- the :c:type:`SUNContext` object (see :numref:`SUNDIALS.SUNContext`)

   **Return value:**
      If successful, a ``SUNDomEigEstimator`` object.  If *q* is
      incompatible, this routine will return ``NULL``.

   **Notes:**
      This routine will perform consistency checks to ensure that it is
      called with a consistent ``N_Vector`` implementation (i.e.  that it
      supplies the requisite vector operations).

      A ``krydim`` argument that is :math:`\leq 2` will result in the default
      value (3).  This default value is particularly chosen to minimize the memory
      footprint.

.. _SUNDomEigEst.ARNI.Description:

SUNDomEigEst_ARNI Description
-----------------------------


The SUNDomEigEst_ARNI module defines the *content* field of a
``SUNDomEigEstimator`` to be the following structure:

.. code-block:: c

   struct _SUNDomEigEstimatorContent_ArnI {
     SUNATimesFn ATimes;
     void* ATdata;
     N_Vector* V;
     N_Vector q;
     sunindextype krydim;
     sunindextype numwarmups;
     sunrealtype* LAPACK_A;
     sunrealtype* LAPACK_wr;
     sunrealtype* LAPACK_wi;
     sunrealtype* LAPACK_work;
     sunrealtype** LAPACK_arr;
     sunrealtype** Hes;
   };


These entries of the *content* field contain the following
information:

* ``ATimes`` - function pointer to perform :math:`Av` product,

* ``ATData`` - pointer to structure for ``ATimes``,

* ``V, q``   - ``N_Vector`` used for workspace by the ARNI algorithm.

* ``krydim`` - dimension of Krylov subspaces (default is 3),

* ``numwarmups`` - number of preprocessing warmups (default is 0),

* ``LAPACK_A, LAPACK_wr, LAPACK_wi, LAPACK_work`` - ``sunrealtype`` used for workspace by LAPACK,

* ``LAPACK_arr`` - storage for the estimated dominant eigenvalues,

* ``Hes`` - Hessenberg matrix,


This estimator is constructed to perform the following operations:

* During construction all ``N_Vector`` estimator data is allocated, with
  vectors cloned from a template ``N_Vector`` that is input, and
  default generic estimator parameters are set.

* User-facing "set" routines may be called to modify default
  estimator parameters.

* An additional "set" routine must be called by the SUNDIALS estimator
  that interfaces with SUNDomEigEst_ARNI to supply the ``ATimes``
  function pointers and the related data ``ATData``.

* In the "initialize" call, the estimator parameters are checked
  for validity and ARNI estimator memory is allocated.

* In the "preprocess" call, the initial nonzero vector :math:`q_0` is warmed up
  :math:`k=` ``numwarmups`` times as

.. math::

    q_1 = \frac{Aq_0}{||Aq_0||} \quad \cdots \quad q_k = \frac{Aq_{k-1}}{||Aq_{k-1}||}.

* In the "estimate" call the ARNI estimator is performed.

The SUNDomEigEst_ARNI module defines implementations of all
dominant eigenvalue estimator operations listed in
:numref:`SUNDomEigEst.API`:

* ``SUNDomEigEst_ArnIGetID``

* ``SUNDomEigEstSetATimes_ArnI``

* ``SUNDomEigEstInitialize_ArnI``

* ``SUNDomEigEstSetNumPreProcess_ArnI``

* ``SUNDomEigEstPreProcess_ArnI``

* ``SUNDomEigEstComputeHess_ArnI``

* ``SUNDomEigEstimate_ArnI``

* ``SUNDomEigEstFree_ArnI``