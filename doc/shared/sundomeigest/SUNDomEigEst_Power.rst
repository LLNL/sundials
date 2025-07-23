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

.. _SUNDomEigEst.POWER:

The SUNDomEigEst_Power Module
======================================

The SUNDomEigEst_Power implementation of the ``SUNDomEigEstimator`` class performs
the Power Iteration (PI) method :cite:p:`vonmises29`; this is an iterative dominant
eigenvalue estimator that is designed to be compatible with any ``N_Vector``
implementation that supports a minimal subset of operations (:c:func:`N_VClone()`,
:c:func:`N_VDotProd()`,  :c:func:`N_VScale()`, and :c:func:`N_VDestroy()`).

Power iteration is useful for large, sparse matrices whose dominant eigenvalue  
is real-valued and has algebraic multiplicity one. The algorithm starts with a non-zero 
vector :math:`\mathbf{v}_{0}`.  It then  iteratively updates this via

.. math::

    \mathbf{v}_{k+1} = \frac{A \mathbf{v}_k}{\|A \mathbf{v}_k\|},

where :math:`\| \cdot \|` denotes the Euclidean norm.  Over successive iterations,
:math:`\mathbf{v}_k` converges to the eigenvector corresponding to
the dominant eigenvalue of :math:`A`.  At each step, the corresponding eigenvalue
can be approximated using the Rayleigh quotient

.. math::

    \lambda_k = \frac{\mathbf{v}_k^T A \mathbf{v}_k}{\|\mathbf{v}_k\|^2}.

The iteration continues until the two successive eigenvalue approximations are
relatively close enough to one another.  That is, for some :ref:`relative tolerance <pi_rel_tol>`.

Power iteration works for the matrices that have a **real** dominant eigenvalue.  
If it is strictly greater than all others (in magnitude), convergence is guaranteed.  
The speed of convergence depends on the ratios of the magnitude of the first two dominant eigenvalues.

The matrix :math:`A` is not required explicitly; only a routine that provides  
the matrix-vector product :math:`Av` is required.  Also, PI requires a fixed 
amount of memory regardless of the number of iterations.  


.. _SUNDomEigEst.POWER.Usage:

SUNDomEigEst_Power Usage
-----------------------------

The header file to be included when using this module is ``sundomeigest/sundomeigest_power.h``.
The SUNDomEigEst_Power module is accessible from all SUNDIALS solvers *without* linking to the
``libsundials_sundomeigestpi`` module library.

The module SUNDomEigEst_Power provides the following user-callable routines:


.. c:function:: SUNDomEigEstimator SUNDomEigEst_Power(N_Vector q, long int max_iters, int num_warmups, sunrealtype rel_tol, SUNContext sunctx)

   This constructor function creates and allocates memory for a PI
   ``SUNDomEigEstimator``.

   **Arguments:**
      * *q* -- the initial guess for the dominant eigenvector; this should not be a non-dominant eigenvector of the Jacobian.
      * *max_iters* -- maximum number of iterations.
      * *num_warmups* -- number of preprocessing warmups.
      * *rel_tol* -- relative tolerance for convergence check.
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

      Although this default number is not high for large matrices,
      it is reasonable since

      1.  most solvers do not need too tight tolerances and consider a safety factor,

      2.  an early (less costly) termination will be a good indicator whether the power iteration is compatible.

      "Warmup" iterations correspond to power iterations that do not check for 
      convergence.  They can help reduce some computational overhead, and may be useful 
      if the initial guess ``q`` is not a good approximation of the dominant 
      eigenvector.  A ``num_warmups`` argument that is :math:` < 0` will result in the 
      default value (0).  This default is chosen to minimize complexity for the general 
      user.

      When the DEE is used in a time-dependent context, however, it is likely that the
      most-recent ``q`` will provide a suitable initial guess for the subsequent call to 
      :c:func:`SUNDomEig_Estimate`.  Thus, when the DEE is used by LSRKStep (see
      :c:func:`LSRKStepSetDomEigEstimator`), the initial value of ``num_warmups`` will 
      be overwritten after the first :c:func:`SUNDomEig_Estimate` call (see
      :c:func:`LSRKStepSetNumSucceedingWarmups`).

      A ``rel_tol`` argument that is :math:` < 0` will result in the default
      value (0.01).  This default is found particularly small enough for many internal applications.


.. _SUNDomEigEst.POWER.Description:

SUNDomEigEst_Power Description
--------------------------------


The SUNDomEigEst_Power module defines the *content* field of a
``SUNDomEigEstimator`` to be the following structure:

.. code-block:: c

   struct _SUNDomEigEstimatorContent_Power {
     SUNATimesFn ATimes;
     void* ATdata;
     N_Vector* V;
     N_Vector q;
     int num_warmups;
     long int max_iters;
     long int cur_num_iters;
     long int max_num_iters;
     long int min_num_iters;
     long int num_ATimes;
     sunrealtype powiter_tol;
     sunrealtype cur_res;
   };


These entries of the *content* field contain the following
information:

* ``ATimes`` - function pointer to perform the product :math:`Av`,

* ``ATData`` - pointer to structure for ``ATimes``,

* ``V, q``   - ``N_Vector`` used for workspace by the PI algorithm.

* ``num_warmups`` - number of preprocessing warmups (default is 0),

* ``max_iters`` - maximum number of iterations (default is 100),

* ``cur_num_iters`` - current number of power iterations,

* ``max_num_iters`` - maximum number of power iterations in any single estimate so far,

* ``min_num_iters`` - minimum number of power iterations in any single estimate so far,

* ``num_ATimes`` - number of calls to the ``ATimes`` function,

* ``powiter_tol`` - convergence criteria for the power iteration (default is 0.01),

* ``cur_res`` - current residual of power iterations.


This estimator is constructed to perform the following operations:

* During construction all ``N_Vector`` estimator data is allocated, with
  vectors cloned from a template ``N_Vector`` that is input, and
  default generic estimator parameters are set.

* User-facing "set" routines may be called to modify default
  estimator parameters.

* An additional "set" routine must be called by the SUNDIALS estimator
  that interfaces with SUNDomEigEst_Power to supply the ``ATimes``
  function pointer and the related data ``ATData``.

* In the "initialize" call, the estimator parameters are checked
  for validity and the initial eigenvector is normalized.

* In the "estimate" call, the initial nonzero vector :math:`q_0` is warmed up
  :math:`k=` ``num_warmups`` times as follows unless otherwise is set by an
  integrator such as by calling :c:func:`LSRKStepSetNumSucceedingWarmups`. 
  Then, the PI estimator is performed.

.. math::

    q_1 = \frac{Aq_0}{||Aq_0||} \quad \cdots \quad q_k = \frac{Aq_{k-1}}{||Aq_{k-1}||}.

The SUNDomEigEst_Power module defines implementations of all dominant
eigenvalue estimator operations listed in
:numref:`SUNDomEigEst.API`:

* ``SUNDomEigEst_SetATimes_Power``

* ``SUNDomEigEst_SetMaxIters_Power``

* ``SUNDomEigEst_SetNumPreProcess_Power``

* ``SUNDomEigEst_SetTol_Power``

* ``SUNDomEigEst_Initialize_Power``

* ``SUNDomEig_Estimate_Power``

* ``SUNDomEigEst_GetCurRes_Power``

* ``SUNDomEigEst_GetCurNumIters_Power``

* ``SUNDomEigEst_GetMaxNumIters_Power``

* ``SUNDomEigEst_GetMinNumIters_Power``

* ``SUNDomEigEst_GetNumATimesCalls_Power``

* ``SUNDomEigEst_PrintStats_Power``

* ``SUNDomEigEst_Destroy_Power``
