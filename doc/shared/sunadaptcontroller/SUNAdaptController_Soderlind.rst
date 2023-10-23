..
   Programmer(s): Daniel R. Reynolds @ SMU
   ----------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2002-2023, Lawrence Livermore National Security
   and Southern Methodist University.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   ----------------------------------------------------------------

.. _SUNAdaptController.Soderlind:

The SUNAdaptController_Soderlind Module
======================================

The Soderlind implementation of the SUNAdaptController class,
SUNAdaptController_Soderlind, implements a general structure for temporal
control proposed by G. Söderlind in :cite:p:`Sod:98`, :cite:p:`Sod:03`
and :cite:p:`Sod:06`.  This controller has the form

.. math::
   h' = h_n \varepsilon_n^{-k_1/(p+1)} \varepsilon_{n-1}^{-k_2/(p+1)} \varepsilon_{n-2}^{-k_3/(p+1)} \left(\dfrac{h_n}{h_{n-1}}\right)^{k_4} \left(\dfrac{h_{n-1}}{h_{n-2}}\right)^{k_5}

with default parameter values :math:`k_1 = 1.25`, :math:`k_2 = 0.5`,
:math:`k_3 = -0.75`, :math:`k_4 = 0.25` and :math:`k_5 = 0.75`, where
:math:`p` is the global order of the time integration method.  In this estimate,
a floor of :math:`\varepsilon > 10^{-10}` is enforced to avoid division-by-zero
errors.  During the first two steps (when :math:`\varepsilon_{n-2}`,
:math:`\varepsilon_{n-1}`, :math:`h_{n-2}` and :math:`h_{n-2}` may be unavailable),
the corresponding terms are merely omitted during estimation of :math:`h'`.

The SUNAdaptController_Soderlind controller is implemented as a derived
SUNAdaptController class, and defines its *content* field as:

.. code-block:: c

   struct _SUNAdaptControllerContent_Soderlind {
     sunrealtype k1;
     sunrealtype k2;
     sunrealtype k3;
     sunrealtype k4;
     sunrealtype k5;
     sunrealtype bias;
     sunrealtype ep;
     sunrealtype epp;
     sunrealtype hp;
     sunrealtype hpp;
     int firststeps;
   };

These entries of the *content* field contain the following information:

* ``k1, k2, k3, k4, k5`` - controller parameters above.

* ``bias`` - error bias factor, that converts from an input temporal error
  estimate via :math:`\varepsilon = \text{bias}*\text{dsm}`.

* ``ep, epp`` - storage for the two previous error estimates,
  :math:`\varepsilon_{n-1}` and :math:`varepsilon_{n-2}`.

* ``hp, hpp`` - storage for the previous two step sizes, :math:`h_{n-1}`
  and :math:`h_{n-2}`.

* ``firststeps`` - counter to handle first two steps (where previous
  step sizes and errors are unavailable).

The header file to be included when using this module is
``sunadaptcontroller/sunadaptcontroller_soderlind.h``.

We note that through appropriate selection of the parameters :math:`k_1 - k_5`,
this controller may replicate the behavior of each of the PID, PI and I
SUNController implementations.

The SUNAdaptController_Soderlind class provides implementations of all operations
relevant to a `SUN_ADAPTCONTROLLER_H` controller listed in
:numref:`SUNAdaptController.Description.operations`. This class
also provides the following additional user-callable routines:


.. c:function:: SUNAdaptController SUNAdaptController_Soderlind(SUNContext sunctx)

   This constructor function creates and allocates memory for a SUNAdaptController_Soderlind
   object, and inserts its default parameters.

   :param sunctx: the current :c:type:`SUNContext` object.
   :return: if successful, a usable :c:type:`SUNAdaptController` object;
            otherwise it will return ``NULL``.

   Usage:

   .. code-block:: c

      SUNAdaptController C = SUNAdaptController_Soderlind(sunctx);

.. c:function:: int SUNAdaptController_SetParams_Soderlind(SUNAdaptController C, sunrealtype k1, sunrealtype k2, sunrealtype k3, sunrealtype k4, sunrealtype k5)

   This user-callable function provides control over the relevant parameters
   above.  This should be called *before* the time integrator is called to evolve
   the problem.

   :param C: the SUNAdaptController_Soderlind object.
   :param k1: parameter used within the controller time step estimate.
   :param k2: parameter used within the controller time step estimate.
   :param k3: parameter used within the controller time step estimate.
   :param k4: parameter used within the controller time step estimate.
   :param k5: parameter used within the controller time step estimate.
   :return: error code indication success or failure (see :numref:`SUNAdaptController.Description.errorCodes`).

   Usage:

   .. code-block:: c

      /* Specify parameters for Soderlind's H_{0}312 controller */
      retval = SUNAdaptController_SetParams_Soderlind(C, 0.25, 0.5, 0.25, -0.75, -0.25);
