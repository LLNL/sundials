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

.. _SUNDomEigEst.PI:

The SUNDomEigEst_PI Module
======================================

The SUNDomEigEst_PI implementation of the ``SUNDomEigEstimator`` class performs
the Power Iteration (PI) method :cite:p:`vonmises29`; this is an iterative dominant
eigenvalue estimator that is designed to be compatible with any ``N_Vector``
implementation that supports a minimal subset of operations (:c:func:`N_VClone()`,
:c:func:`N_VDotProd()`,  :c:func:`N_VScale()`, and :c:func:`N_VDestroy()`).

PI is useful for large, sparse matrices whose dominant eigenvalue  is real-valued 
and has algebraic multiplicity one. The algorithm starts with a non-zero vector 
:math:`\mathbf{v}_{0}`.  It then  iteratively updates this via

.. math::

    \mathbf{v}_{k+1} = \frac{A \mathbf{v}_k}{\|A \mathbf{v}_k\|},

where :math:`\| \cdot \|` denotes the Euclidean norm.  Over successive iterations,
:math:`\mathbf{v}_k` converges to the eigenvector corresponding to
the dominant eigenvalue of :math:`A`.  At each step, the corresponding eigenvalue
can be approximated using the Rayleigh quotient

.. math::

    \lambda_k = \frac{\mathbf{v}_k^T A \mathbf{v}_k}{\|\mathbf{v}_k\|^2}.
The iteration continues until the two successive eigenvalue approximations are
relatively close enough to one another.  That is, for some relative tolerance
:math:`\tau`,
enough to one another.  That is, for some relative tolerance :math:`\tau`,

.. math::

    \frac{\left|\lambda_k - \lambda_{k-1}\right|}{\left|\lambda_k \right|} < \tau.

PI works for the matrices that have a **real** dominant eigenvalue.  If it is strictly
greater than all others (in magnitude), convergence is guaranteed.  The speed of convergence
depends on the ratios of the magnitude of the first two dominant eigenvalues.

The matrix :math:`A` is not required explicitly; only a routine that provides  
the matrix-vector product :math:`Av` is required.  Also, PI requires a fixed 
amount of memory regardless of the number of iterations.  


.. _SUNDomEigEst.PI.Usage:

SUNDomEigEst_PI Usage
---------------------

The header file to be included when using this module is ``sundomeigest/sundomeigest_pi.h``.
The SUNDomEigEst_PI module is accessible from all SUNDIALS solvers *without* linking to the
``libsundials_sundomeigestpi`` module library.

The module SUNDomEigEst_PI provides the following user-callable routines:


.. c:function:: SUNDomEigEstimator SUNDomEigEst_PI(N_Vector q, int max_iters, SUNContext sunctx)

   This constructor function creates and allocates memory for a PI
   ``SUNDomEigEstimator``.

   **Arguments:**
      * *q* -- a template vector.
      * *max_iters* -- maximum number of iterations.
      * *sunctx* -- the :c:type:`SUNContext` object (see :numref:`SUNDIALS.SUNContext`)

   **Return value:**
      If successful, a ``SUNDomEigEstimator`` object.  If *q* is
      incompatible, this routine will return ``NULL``.

   **Notes:**
      This routine will perform consistency checks to ensure that it is
      called with a consistent ``N_Vector`` implementation (i.e.  that it
      supplies the requisite vector operations).

      A ``max_iters`` argument that is :math:`\le0` will result in the default
      value (100).

      Although this default number is not high for large martices,
      it is reasonable since

      1.  most solvers do not need too tight tolerances and consider a safety factor,

      2.  an early (less costly) termination will be a good indicator if PI is compatible.


.. _SUNDomEigEst.PI.Description:

SUNDomEigEst_PI Description
---------------------------


The SUNDomEigEst_PI module defines the *content* field of a
``SUNDomEigEstimator`` to be the following structure:

.. code-block:: c

   struct _SUNDomEigEstimatorContent_PI {
     SUNATimesFn ATimes;
     void* ATdata;
     N_Vector* V;
     N_Vector q;
     int numwarmups;
     int max_iters;
     int curnumiters;
     int maxnumiters;
     int minnumiters;
     long int nATimes;
     sunrealtype powiter_tol;
     sunrealtype curres;
   };


These entries of the *content* field contain the following
information:

* ``ATimes`` - function pointer to perform the product :math:`Av`,

* ``ATData`` - pointer to structure for ``ATimes``,

* ``V, q``   - ``N_Vector`` used for workspace by the PI algorithm.

* ``numwarmups`` - number of preprocessing warmups (default is 0),

* ``max_iters`` - maximum number of iterations (default is 100),

* ``curnumiters`` - current number of power iterations,

* ``maxnumiters`` - maximum number of power iterations so far,

* ``minnumiters`` - minimum number of power iterations so far,

* ``nATimes`` - number of calls to the ``ATimes`` function,

* ``powiter_tol`` - convergence criteria for the power iteration (default is 0.01),

* ``curres`` - current residual of power iterations.


This estimator is constructed to perform the following operations:

* During construction all ``N_Vector`` estimator data is allocated, with
  vectors cloned from a template ``N_Vector`` that is input, and
  default generic estimator parameters are set.

* User-facing "set" routines may be called to modify default
  estimator parameters.

* An additional "set" routine must be called by the SUNDIALS estimator
  that interfaces with SUNDomEigEst_PI to supply the ``ATimes``
  function pointer and the related data ``ATData``.

* In the "initialize" call, the estimator parameters are checked
  for validity and PI estimator memory is allocated.

* In the "preprocess" call, the initial vector :math:`q_0` is warmed up
  :math:`k=` ``numwarmups`` times as

.. math::

    q_1 = \frac{Aq_0}{||Aq_0||} \quad \cdots \quad q_k = \frac{Aq_{k-1}}{||Aq_{k-1}||}.

* In the "estimate" call the PI estimator is performed.

The SUNDomEigEst_PI module defines implementations of all dominant
eigenvalue estimator operations listed in
:numref:`SUNDomEigEst.API`:

* ``SUNDomEigEst_SetATimes_PI``

* ``SUNDomEigEst_Initialize_PI``

* ``SUNDomEigEst_SetNumPreProcess_PI``

* ``SUNDomEigEst_SetTol_PI``

* ``SUNDomEigEst_SetMaxIters_PI``

* ``SUNDomEigEst_PreProcess_PI``

* ``SUNDomEig_Estimate_PI``

* ``SUNDomEigEst_GetCurRes_PI``

* ``SUNDomEigEst_GetCurNumIters_PI``

* ``SUNDomEigEst_GetMaxNumIters_PI``

* ``SUNDomEigEst_GetMinNumIters_PI``

* ``SUNDomEigEst_GetNumATimesCalls_PI``

* ``SUNDomEigEst_PrintStats_PI``

* ``SUNDomEigEstFree_PI``
