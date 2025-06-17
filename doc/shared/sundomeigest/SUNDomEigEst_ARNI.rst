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
the Arnoldi Iteration :cite:p:`arnoldi51` method; this is an iterative dominant
eigenvalue estimator that is designed to be compatible with any ``N_Vector``
implementation that supports a minimal subset of operations

TODO:Check the list.
(:c:func:`N_VClone()`, :c:func:`N_VDotProd()`, :c:func:`N_VScale()`,
:c:func:`N_VLinearSum()`, :c:func:`N_VProd()`, and
:c:func:`N_VDestroy()`).  ARNI requires a prefixed amount of
memory that depends on the user provided dimension of Krylov subspaces.

The matrix :math:`A` is not required explicitly; only routines
that provide :math:`A` as operator is required.


.. _SUNDomEigEst.ARNI.Usage:

SUNDomEigEst_ARNI Usage
-----------------------

The header file to be included when using this module is
``sundomeigest/sundomeigest_arni.h``.  The SUNDomEigEst_ARNI module is accessible from all SUNDIALS solvers
*without* linking to the ``libsundials_sundomeigestarni`` module library after enabling LAPACK package.
This LAPACK dependence is limited to the ``dgeev_`` function after computing Hessenberg matrix internally.

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
      called with a consistent ``N_Vector`` implementation (i.e. that it
      supplies the requisite vector operations).

      A ``krydim`` argument that is :math:`\leq 2` will result in the default
      value (3).


.. c:function:: SUNErrCode SUNDomEigEstComputeHess_ArnI(SUNDomEigEstimator DEE)

   This function computes the Hessenberg matrix.

   **Return value:**

      A :c:type:`SUNErrCode`.

   **Notes:**
      This routine can be called with (:c:func:`SUNDomEigEstPreProcess`) as well.
      It must be called after initialization (:c:func:`SUNDomEigEstInitialize`),
      and also preprocess (if requested) should be performed (:c:func:`SUNDomEigEstPreProcess`)
      before calling this function.


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
     suncomplextype* LAPACK_arr;
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

* ``LAPACK_arr`` - ``suncomplextype`` used for workspace by LAPACK,

* ``Hes`` - Hessenberg matrix,


This estimator is constructed to perform the following operations:

* During construction all ``N_Vector`` estimator data is allocated, with
  vectors cloned from a template ``N_Vector`` that is input, and
  default generic estimator parameters are set.

* User-facing "set" routines may be called to modify default
  estimator parameters.

* Additional "set" routines are called by the SUNDIALS estimator
  that interfaces with SUNDomEigEst_ARNI to supply the
  ``ATimes`` function pointers and the related data ``ATData``.

* In the "initialize" call, the estimator parameters are checked
  for validity and ARNI estimator memory is allocated.

* In the "preprocess" call, the initial random vector :math:`q_0` is warmed up
  :math:`k=` ``numwarmups`` times as :math:`q_1 = \frac{Aq_0}{||Aq_0||} \cdots q_k = \frac{Aq_{k-1}}{||Aq_{k-1}||}`.

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
