.. ----------------------------------------------------------------
   Programmer(s): Daniel R. Reynolds @ SMU
                  David J. Gardner @ LLNL
   ----------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2002-2024, Lawrence Livermore National Security
   and Southern Methodist University.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   ----------------------------------------------------------------

.. _LSRKSTEP.Usage.UserSupplied:

User-supplied functions
=============================

Currently, there is only one user-supplied function for LSRKStep (in addition to the :c:type:`ARKRhsFn` arguments that define the IVP).




.. _LSRKStep.Usage.DomEig:

The dominant eigenvalue estimation
----------------------------------

When running LSRKStep with either the RKC or RKL methods, the user must supply a dominant eigenvalue estimation function of type :c:type:`ARKDomEigFn`:

.. c:type:: int (*ARKDomEigFn)(sunrealtype* t, N_Vector y, sunrealtype* lambdaR, sunrealtype* lambdaI, void* user_data)

   These functions compute the dominant eigenvalue of the Jacobian of the ODE right-hand side for a given
   value of the independent variable :math:`t` and state vector :math:`y`.

   :param t: the current value of the independent variable.
   :param y: the current value of the dependent variable vector.
   :param lambdaR: The real part of the dominant eigenvalue.
   :param lambdaI: The imaginary part of the dominant eigenvalue.   
   :param user_data: the `user_data` pointer that was passed to
                     :c:func:`ARKodeSetUserData`.

   :return: An *ARKDomEigFn* should return 0 if successful and any nonzero for a failure.