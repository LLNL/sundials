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
of modules.  The inputs to the examples depend on the estimator type,
and are output to ``stdout`` if the example is run without the
appropriate number of command-line arguments.

The following is a list of the example functions in ``test_sundomeigest.c``:

* ``Test_SUNDomEigEstGetID``: Verifies the returned estimator identifier against
  the value that should be returned.

* ``Test_SUNDomEigEstSetATimes`` Verifies that ``SUNDomEigEstSetATimes`` can
  be called and returns successfully.

* ``Test_SUNDomEigEstSetMaxPowerIter`` Verifies that
  ``SUNDomEigEstSetMaxPowerIter`` can be called and returns successfully.


  ``SUNDomEigEstSetMaxPowerIter`` is not an option for some estimators, e.g.,
  Arnoldi iteration.  In this case, it should return with SUN_SUCCESS.

* ``Test_SUNDomEigEstSetNumPreProcess`` Verifies that
  ``SUNDomEigEstSetNumPreProcess`` can be called and returns successfully.

* ``Test_SUNDomEigEstSetTol`` Verifies that
  ``SUNDomEigEstSetTol`` can be called and returns successfully.


  ``SUNDomEigEstSetTol`` is not an option for some estimators, e.g.,
  Arnoldi iteration.  In this case, it should return with SUN_SUCCESS.

* ``Test_SUNDomEigEstInitialize``: Verifies that ``SUNDomEigEstInitialize``
  can be called and returns successfully.

* ``Test_SUNDomEigEstPreProcess``: Verifies that ``SUNDomEigEstPreProcess``
  can be called and returns successfully.

* ``Test_SUNDomEigEstComputeHess``: Verifies that ``SUNDomEigEstComputeHess``
  can be called and returns successfully.


  ``SUNDomEigEstComputeHess`` is not an option for some estimators, e.g.,
  Power iteration.  In this case, it should return with SUN_SUCCESS.
  A failure flag returns otherwise.

* ``Test_SUNDomEigEstimate``: Verifies that ``SUNDomEigEstimate``
  can be called and returns successfully.  Estimated dominant eigenvalue is
  :math:`\lambda_{\max} = \lambda` such that
  :math:`|\lambda| = \max\{|\lambda_i| : A \vec{v_i} = \lambda_i \vec{v_i}, \ \vec{v_i} \neq \vec{0} \}`.
  This test compares the estimated dominant eigenvalue to the known value
  and returns a passing flag if the estimation is within a specified relative
  tolerance; otherwise, it returns a failure flag.

* ``Test_SUNDomEigEstNumIters`` Verifies that
  ``SUNDomEigEstNumIters`` can be called and returns successfully.


  ``SUNDomEigEstNumIters`` is not an option for some estimators, e.g.,
  Arnoldi iteration.  In this case, it should return with SUN_SUCCESS
  and `niter = 0`.  A failure flag returns otherwise.

* ``Test_SUNDomEigEstRes`` Verifies that
  ``SUNDomEigEstRes`` can be called and returns successfully.


  ``SUNDomEigEstRes`` is not an option for some estimators, e.g.,
  Arnoldi iteration.  In this case, it should return with SUN_SUCCESS
  and `res = 0`.  A failure flag returns otherwise.

We'll note that these tests should be performed in a particular
order.  For all estimators,
``SUNDomEigEstSetATimes`` must be called
before other set routines e.g., ``SUNDomEigEstSetMaxPowerIter``,
``SUNDomEigEstSetNumPreProcess``, ``SUNDomEigEstSetTol`` (if applicable).
Then, ``Test_SUNDomEigEstInitialize`` must be called before
``Test_SUNDomEigEstPreProcess`` (if applicable).
``Test_SUNDomEigEstComputeHess`` (if the estimator requires)
must be called next and before ``Test_SUNDomEigEstimate``.
For the estimator stats ``Test_SUNDomEigEstNumIters`` and ``Test_SUNDomEigEstRes``
should be called after ``Test_SUNDomEigEstimate``.
These are called in the appropriate order in all of the example problems.
