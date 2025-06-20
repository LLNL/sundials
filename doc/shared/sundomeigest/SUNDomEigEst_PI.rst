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
the Power Iteration :cite:p:`vonmises29` method; this is an iterative dominant
eigenvalue estimator that is designed to be compatible with any ``N_Vector``
implementation that supports a minimal subset of operations

TODO:Check the list.
(:c:func:`N_VClone()`, :c:func:`N_VDotProd()`, :c:func:`N_VScale()`,
:c:func:`N_VLinearSum()`, :c:func:`N_VProd()`, and
:c:func:`N_VDestroy()`).  PI requires a fixed amount of memory regardless of
the number of Iterations.

The matrix :math:`A` is not required explicitly; only routines
that provide :math:`A` as operator is required.


.. _SUNDomEigEst.PI.Usage:

SUNDomEigEst_PI Usage
---------------------

The header file to be included when using this module is
``sundomeigest/sundomeigest_pi.h``.  The SUNDomEigEst_PI module is accessible from all SUNDIALS solvers
*without* linking to the ``libsundials_sundomeigestpi`` module library.

The module SUNDomEigEst_PI provides the following user-callable routines:


.. c:function:: SUNDomEigEstimator SUNDomEigEst_PI(N_Vector q, sunindextype max_powiter, SUNContext sunctx)

   This constructor function creates and allocates memory for a PI
   ``SUNDomEigEstimator``.

   **Arguments:**
      * *q* -- a template vector.
      * *max_powiter* -- maximum number of power iterations.
      * *sunctx* -- the :c:type:`SUNContext` object (see :numref:`SUNDIALS.SUNContext`)

   **Return value:**
      If successful, a ``SUNDomEigEstimator`` object.  If *q* is
      incompatible, this routine will return ``NULL``.

   **Notes:**
      This routine will perform consistency checks to ensure that it is
      called with a consistent ``N_Vector`` implementation (i.e. that it
      supplies the requisite vector operations).

      A ``max_powiter`` argument that is :math:`\le0` will result in the default
      value (100).


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
     sunindextype max_powiter;
     sunindextype numwarmups;
     sunrealtype powiter_tol;
     sunrealtype res;
     sunindextype numiters;
   };


These entries of the *content* field contain the following
information:

* ``ATimes`` - function pointer to perform :math:`Av` product,

* ``ATData`` - pointer to structure for ``ATimes``,

* ``V, q``   - ``N_Vector`` used for workspace by the PI algorithm.

* ``max_powiter`` - maximum number of power iterations (default is 100),

* ``numwarmups`` - number of preprocessing warmups (default is 0),

* ``powiter_tol`` - convergence criteria for the power iteration (default is 0.01),

* ``res`` - current residual of power iterations,

* ``numiters`` - current number of power iterations.


This estimator is constructed to perform the following operations:

* During construction all ``N_Vector`` estimator data is allocated, with
  vectors cloned from a template ``N_Vector`` that is input, and
  default generic estimator parameters are set.

* User-facing "set" routines may be called to modify default
  estimator parameters.

* Additional "set" routines are called by the SUNDIALS estimator
  that interfaces with SUNDomEigEst_PI to supply the
  ``ATimes`` function pointers and the related data ``ATData``.

* In the "initialize" call, the estimator parameters are checked
  for validity and PI estimator memory is allocated.

* In the "preprocess" call, the initial random vector :math:`q_0` is warmed up
  :math:`k=` ``numwarmups`` times as :math:`q_1 = \frac{Aq_0}{||Aq_0||} \cdots q_k = \frac{Aq_{k-1}}{||Aq_{k-1}||}`.

* In the "estimate" call the PI estimator is performed.

The SUNDomEigEst_PI module defines implementations of all
dominant eigenvalue estimator operations listed in
:numref:`SUNDomEigEst.API`:

* ``SUNDomEigEst_PIGetID``

* ``SUNDomEigEstSetATimes_PI``

* ``SUNDomEigEstInitialize_PI``

* ``SUNDomEigEstSetNumPreProcess_PI``

* ``SUNDomEigEstSetMaxPowerIter_PI``

* ``SUNDomEigEstPreProcess_PI``

* ``SUNDomEigEstimate_PI``

* ``SUNDomEigEstNumIters_PI``

* ``SUNDomEigEstRes_PI``

* ``SUNDomEigEstFree_PI``
