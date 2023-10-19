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

.. _SUNAdaptController.PI:

The SUNAdaptController_PI Module
======================================

The PI implementation of the SUNAdaptController class, SUNAdaptController_PI, implements a
standard PI temporal adaptivity controller.  Like with SUNAdaptController_PID, the PI
controller derives from those found in :cite:p:`KenCarp:03`, :cite:p:`Sod:98`,
:cite:p:`Sod:03` and :cite:p:`Sod:06`, but it differs in that it only uses the
two most recent step sizes in its adaptivity algorithm,

.. math::
   h' \;=\; h_n\; \varepsilon_n^{-k_1/ord}\; \varepsilon_{n-1}^{k_2/ord}.

where the constants :math:`k_1` and :math:`k_2` default to 0.8 and 0.31,
respectively, and :math:`ord = p+1+adj`, where both :math:`p` and :math:`adj` are
described below. In this estimate, a floor of :math:`\varepsilon > 10^{-10}` is enforced to
avoid division-by-zero errors.

This is implemented as a derived SUNAdaptController class, and defines its *content*
field as:

.. code-block:: c

   struct _SUNAdaptControllerContent_PI {
     sunrealtype k1;
     sunrealtype k2;
     sunrealtype bias;
     sunrealtype ep;
     int p;
     int adj;
   };

These entries of the *content* field contain the following information:

* ``k1, k2`` - controller parameters above.

* ``bias`` - error bias factor, that converts from an input temporal error
  estimate via :math:`\varepsilon = \text{bias}*\text{dsm}`.

* ``ep`` - storage for the previous error estimate, :math:`\varepsilon_{n-1}`.

* ``p`` - asymptotic order to use in error control.  This is provided by
  the time integrator, corresponding to the order of accuracy for the time
  integration method, the embedding, or the minimum of the two.

* ``adj`` - order of accuracy adjustment to use within the controller [default ``-1``].


The header file to be included when using this module is
``sunadaptcontroller/sunadaptcontroller_pi.h``.

The SUNAdaptController_PI class provides implementations of all operations
relevant to a `SUN_ADAPTCONTROLLER_H` controller listed in
:numref:`SUNAdaptController.Description.operations`. The SUNAdaptController_PI class
also provides the following additional user-callable routines:


.. c:function:: SUNAdaptController SUNAdaptController_PI(SUNContext sunctx)

   This constructor function creates and allocates memory for a SUNAdaptController_PI
   object, and inserts its default parameters.

   :param sunctx: the current :c:type:`SUNContext` object.
   :return: if successful, a usable :c:type:`SUNAdaptController` object; otherwise it will return ``NULL``.

   Usage:

   .. code-block:: c

      SUNAdaptController C = SUNAdaptController_PI(sunctx);

.. c:function:: int SUNAdaptController_SetParams_PI(SUNAdaptController C, sunrealtype k1, sunrealtype k2)

   This user-callable function provides control over the relevant parameters
   above.  This should be called *before* the time integrator is called to evolve
   the problem.

   :param C: the SUNAdaptController_PI object.
   :param k1: parameter used within the controller time step estimate (only stored if non-negative).
   :param k2: parameter used within the controller time step estimate (only stored if non-negative).
   :return: error code indication success or failure (see :numref:`SUNAdaptController.Description.errorCodes`).

   Usage:

   .. code-block:: c

      retval = SUNAdaptController_SetParams_PI(C, 0.9, 0.3);
