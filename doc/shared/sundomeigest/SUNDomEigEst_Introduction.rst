.. ----------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2002-2025, Lawrence Livermore National Security
   and Southern Methodist University.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   ----------------------------------------------------------------

.. _SUNDomEigEst:

##############################
Dominant Eigenvalue Estimators
##############################

For problems that require the dominant eigenvalue of a matrix (Jacobian),
the SUNDIALS packages operate using generic dominant eigenvalue estimator modules
defined through the :c:type:`SUNDomEigEstimator`, or "SUNDomEigEst", API.
This allows SUNDIALS packages to utilize any valid SUNDomEigEst
implementation that provides a set of required functions.  These
functions can be divided into three categories.  The first are the core
estimator functions.  The second group consists of "set" routines
to supply the dominant eigenvalue estimator object with functions provided by the
SUNDIALS package, or for modification of estimator parameters.  The last
group consists of "get" routines for retrieving artifacts (statistics,
residual, etc.) from the estimator.  All of these functions
are defined in the header file ``sundials/sundials_domeigestimator.h``.

The implementations provided with SUNDIALS work in coordination
with the SUNDIALS :c:type:`N_Vector`, and optionally :c:type:`SUNMatrix`,
modules to provide a set of compatible data structures for the estimator.
Moreover, advanced users can provide a customized ``SUNDomEigEstimator``
implementation to any SUNDIALS package, particularly in cases where they
provide their own ``N_Vector`` and/or ``SUNMatrix`` modules.

While Krylov-based estimators preset the number of Krylov subspace
dimensions, resulting in a tolerance-free estimation, SUNDIALS requires
that iterative estimators stop when the residual meets a prescribed
tolerance, i.e.,

.. math::

   ||\lambda_{k+1} - \lambda_k|| < \text{tol}.

For users interested in providing their own SUNDomEigEst module, the
following section presents the SUNDomEigEst API and its implementation
beginning with the definition of SUNDomEigEst functions in
:numref:`SUNDomEigEst.CoreFn` -- :numref:`SUNDomEigEst.GetFn`. This is followed by
the definition of functions supplied to an estimator implementation in
:numref:`SUNDomEigEst.SUNSuppliedFn`. The estimator return codes are described
in :numref:`SUNDomEigEst.ErrorCodes`. The ``SUNDomEigEstimator`` type and the
generic SUNDomEigEst module are defined in :numref:`SUNDomEigEst.Generic`.
:numref:`SUNDomEigEst.API.Custom` lists the requirements for supplying a custom
SUNDomEigEst module and discusses some intended use cases. Users wishing to
supply their own SUNDomEigEst module are encouraged to use the SUNDomEigEst
implementations provided with SUNDIALS as a template for supplying custom
dominant eigenvalue estimator modules. The section that then follows describes
the SUNDomEigEst functions required by this SUNDIALS package, and provides
additional package specific details. Then the remaining sections of this
chapter present the SUNDomEigEst modules provided with SUNDIALS.
