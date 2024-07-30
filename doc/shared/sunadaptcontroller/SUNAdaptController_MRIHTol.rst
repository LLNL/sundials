..
   Programmer(s): Daniel R. Reynolds @ SMU
   ----------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2002-2024, Lawrence Livermore National Security
   and Southern Methodist University.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   ----------------------------------------------------------------

.. _SUNAdaptController.MRIHTol:

The SUNAdaptController_MRIHTol Module
=======================================

.. versionadded:: x.y.z

The MRIHTol implementation of the SUNAdaptController class,
SUNAdaptController_MRIHTol, implements a general structure for telescopic
multirate temporal control.  This controller has the form

.. math::
   h' = h_n \varepsilon_n^{-k_1/(p+1)} \varepsilon_{n-1}^{-k_2/(p+1)} \varepsilon_{n-2}^{-k_3/(p+1)} \left(\dfrac{h_n}{h_{n-1}}\right)^{k_4} \left(\dfrac{h_{n-1}}{h_{n-2}}\right)^{k_5}

with default parameter values :math:`k_1 = 1.25`, :math:`k_2 = 0.5`,
:math:`k_3 = -0.75`, :math:`k_4 = 0.25`, and :math:`k_5 = 0.75`, where
:math:`p` is the global order of the time integration method.  In this estimate,
a floor of :math:`\varepsilon_* > 10^{-10}` is enforced to avoid division-by-zero
errors.  During the first two steps (when :math:`\varepsilon_{n-2}`,
:math:`\varepsilon_{n-1}`, :math:`h_{n-2}`, and :math:`h_{n-2}` may be unavailable),
the corresponding terms are merely omitted during estimation of :math:`h'`.

The SUNAdaptController_Soderlind controller is implemented as a derived
SUNAdaptController class, and defines its *content* field as:

.. code-block:: c

   struct _SUNAdaptControllerContent_Soderlind {
     SUNAdaptController HControl;
     SUNAdaptController TolControl;
     sunrealtype inner_max_relch;
     sunrealtype inner_min_tolfac;
     sunrealtype inner_max_tolfac;
   };

These entries of the *content* field contain the following information:

* ``HControl`` - single time-scale SUNAdaptController object to adapt
  the current step size, :math:`H`.

* ``TolControl`` - single time-scale SUNAdaptController object to adapt
  the fast time scale relative tolerance factor, :math:`reltol`.

* ``inner_max_relch`` - bound on the relative change in :math:`tolfac` from one
  step to the next.

* ``inner_min_tolfac`` - absolute lower bound on :math:`tolfac`, to avoid calling
  the inner integrator with an arbitrarily small tolerance.

* ``inner_max_tolfac`` - absolute upper bound on :math:`tolfac`, to avoid calling
  the inner integrator with too loose of a tolerance.

The header file to be included when using this module is
``sunadaptcontroller/sunadaptcontroller_mrihtol.h``.

The SUNAdaptController_Soderlind class provides implementations of all operations
relevant to a ``SUN_ADAPTCONTROLLER_MRI_TOL`` controller listed in
:numref:`SUNAdaptController.Description.operations`. This class
also provides the following additional user-callable routines:


.. c:function:: SUNAdaptController SUNAdaptController_MRIHTol(SUNContext sunctx, SUNAdaptController HControl, SUNAdaptController TolControl)

   This constructor creates and allocates memory for a SUNAdaptController_MRIHTol
   object, and inserts its default parameters.

   :param sunctx: the current :c:type:`SUNContext` object.
   :param HControl: the slow time step adaptivity controller object.
   :param TolControl: the fast time scale tolerance factor adaptivity controller object.
   :return: if successful, a usable :c:type:`SUNAdaptController` object;
            otherwise it will return ``NULL``.

   Usage:

   .. code-block:: c

      SUNAdaptController C = SUNAdaptController_MRIHTol(sunctx, HControl, TolControl);

.. c:function:: SUNErrCode SUNAdaptController_SetParams_MRIHTol(SUNAdaptController C, sunrealtype inner_max_relch, sunrealtype inner_min_tolfac, sunrealtype inner_max_tolfac)

   This user-callable function provides control over the relevant parameters
   above.  This should be called *before* the time integrator is called to evolve
   the problem.

   :param C: the SUNAdaptController_MRIHTol object.
   :param inner_max_relch: parameter used within the controller time step estimate.
   :param inner_min_tolfac: parameter used within the controller time step estimate.
   :param inner_max_tolfac: parameter used within the controller time step estimate.
   :return: :c:type:`SUNErrCode` indicating success or failure.

   Usage:

   .. code-block:: c

      retval = SUNAdaptController_SetParams_MRIHTol(C, 20.0, 1e-5, 1.0);
