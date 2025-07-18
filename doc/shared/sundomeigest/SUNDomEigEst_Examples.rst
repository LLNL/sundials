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

.. _SUNDomEigEst.Examples:

SUNDomEigEstimator Examples
======================================

There are ``SUNDomEigEstimator`` examples that may be installed for each
implementation; these make use of the functions in ``test_sundomeigest.c``.
These example functions show simple usage of the ``SUNDomEigEstimator`` family
of modules.  The command-line inputs to the examples depend on the estimator type,
and are output to ``stdout`` if the example is run without the
appropriate number of command-line arguments.

The following is a list of the example functions in ``test_sundomeigest.c``:

* ``Test_SUNDomEigEst_SetATimes`` Verifies that ``SUNDomEigEst_SetATimes`` can
  be called and returns successfully.

* ``Test_SUNDomEigEst_SetMaxIters`` Verifies that
  ``SUNDomEigEst_SetMaxIters`` can be called and returns successfully.


  ``SUNDomEigEst_SetMaxIters`` is not an option for some estimators, e.g.,
  Arnoldi iteration.  In this case, it should return with SUN_SUCCESS.

* ``Test_SUNDomEigEst_SetNumPreProcess`` Verifies that
  ``SUNDomEigEst_SetNumPreProcess`` can be called and returns successfully.

* ``Test_SUNDomEigEst_SetTol`` Verifies that
  ``SUNDomEigEst_SetTol`` can be called and returns successfully.


  ``SUNDomEigEst_SetTol`` is not an option for some estimators, e.g.,
  Arnoldi iteration.  In this case, it should return with SUN_SUCCESS.

* ``Test_SUNDomEigEst_Initialize``: Verifies that ``SUNDomEigEst_Initialize``
  can be called and returns successfully.

* ``Test_SUNDomEig_Estimate``: Verifies that ``SUNDomEig_Estimate``
  can be called and returns successfully.  The estimated dominant eigenvalue is
  :math:`\lambda_{\max} = \lambda_R + \lambda_I i` such that
  :math:`|\lambda| = \max\{|\lambda_i| : A \vec{v_i} = \lambda_i \vec{v_i}, \ \vec{v_i} \neq \vec{0} \}`.
  This test compares the estimated dominant eigenvalue to the known value
  and returns a passing flag if the estimation is within a specified relative
  tolerance; otherwise, it returns a failure flag.

* ``Test_SUNDomEigEst_GetCurRes`` Verifies that
  ``SUNDomEigEst_GetCurRes`` can be called and returns successfully.


  ``SUNDomEigEst_GetCurRes`` is not an option for some estimators, e.g.,
  Arnoldi iteration.  In this case, it should return with SUN_SUCCESS
  and `cur_res = 0`.  A failure flag returns otherwise.

* ``Test_SUNDomEigEst_GetCurNumIters`` Verifies that
  ``SUNDomEigEst_GetCurNumIters`` can be called and returns successfully.


  ``SUNDomEigEst_GetCurNumIters`` is not an option for some estimators, e.g.,
  Arnoldi iteration.  In this case, it should return with SUN_SUCCESS
  and `curniter = 0`.  A failure flag returns otherwise.

* ``Test_SUNDomEigEst_GetMaxNumIters`` Verifies that
  ``SUNDomEigEst_GetMaxNumIters`` can be called and returns successfully.

  ``SUNDomEigEst_GetMaxNumIters`` is not an option for some estimators, e.g.,
  Arnoldi iteration.  In this case, it should return with SUN_SUCCESS
  and `max_niter = 0`.  A failure flag returns otherwise.

* ``Test_SUNDomEigEst_GetMinNumIters`` Verifies that
  ``SUNDomEigEst_GetMinNumIters`` can be called and returns successfully.

  ``SUNDomEigEst_GetMinNumIters`` is not an option for some estimators, e.g.,
  Arnoldi iteration.  In this case, it should return with SUN_SUCCESS
  and `min_niter = 0`.  A failure flag returns otherwise.

* ``Test_SUNDomEigEst_GetNumATimesCalls`` Verifies that
  ``SUNDomEigEst_GetNumATimesCalls`` can be called and returns successfully.

* ``Test_SUNDomEigEst_PrintStats`` Verifies that
  ``SUNDomEigEst_PrintStats`` can be called and returns successfully.

We'll note that these tests should be performed in a particular
order.  For all estimators,
``SUNDomEigEst_SetATimes`` must be called before other set routines e.g., 
``SUNDomEigEst_SetMaxIters``, ``SUNDomEigEst_SetNumPreProcess``, 
``SUNDomEigEst_SetTol`` (if applicable).
Then, ``Test_SUNDomEigEst_Initialize`` must be called before
``Test_SUNDomEig_Estimate``.
For the estimator stats ``Test_SUNDomEigEst_GetCurRes``, ``Test_SUNDomEigEst_GetCurNumIters``,
``Test_SUNDomEigEst_GetMaxNumIters``, ``Test_SUNDomEigEst_GetMinNumIters``, 
``SUNDomEigEst_PrintStats`` and ``SUNDomEigEst_GetNumATimesCalls`` should be called after 
``Test_SUNDomEig_Estimate``.
These are called in the appropriate order in all of the example problems.
