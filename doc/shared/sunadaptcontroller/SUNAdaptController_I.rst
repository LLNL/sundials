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

.. _SUNAdaptController.I:

The SUNAdaptController_I Module
======================================

The I implementation of the SUNAdaptController class, SUNAdaptController_I, implements a
standard I temporal adaptivity controller.  This is the standard time adaptivity
control algorithm in use by most publicly-available ODE solver codes.  It bases
the prospective time step estimate entirely off of the current local error
estimate,

.. math::
   h' \;=\; h_n\; \varepsilon_n^{-k_1/ord},

where :math:`ord = p+1`, where :math:`p` is described below. In this estimate, a floor of
:math:`\varepsilon > 10^{-10}` is enforced to avoid division-by-zero errors. By default
the constant :math:`k_1=1`.

This is implemented as a derived SUNAdaptController class, and defines its *content*
field as:

.. code-block:: c

   struct _SUNAdaptControllerContent_I {
     sunrealtype k1;
     sunrealtype bias;
     int p;
   };

These entries of the *content* field contain the following information:

* ``k1`` - controller parameter above.

* ``bias`` - error bias factor, that converts from an input temporal error
  estimate via :math:`\varepsilon = \text{bias}*\text{dsm}`.

* ``p`` - asymptotic order to use in error control (provided by the time integrator).


The header file to be included when using this module is
``sunadaptcontroller/sunadaptcontroller_i.h``.

The SUNAdaptController_I class provides implementations of all operations
relevant to a `SUN_ADAPTCONTROLLER_H` controller listed in
:numref:`SUNAdaptController.Description.operations`. The SUNAdaptController_I class
also provides the following additional user-callable routines:


.. c:function:: SUNAdaptController SUNAdaptController_I(SUNContext sunctx)

   This constructor function creates and allocates memory for a SUNAdaptController_I
   object, and inserts its default parameters.

   :param sunctx: the current :c:type:`SUNContext` object.
   :return: if successful, a usable :c:type:`SUNAdaptController` object; otherwise it will return ``NULL``.

   Usage:

   .. code-block:: c

      SUNAdaptController C = SUNAdaptController_I(sunctx);

.. c:function:: int SUNAdaptController_SetParams_I(SUNAdaptController C, sunrealtype k1)

   This user-callable function provides control over the relevant parameters
   above.  This should be called *before* the time integrator is called to evolve
   the problem.

   :param C: the SUNAdaptController_I object.
   :param k1: parameter used within the controller time step estimate (only stored if non-negative).
   :return: error code indication success or failure (see :numref:`SUNAdaptController.Description.errorCodes`).

   Usage:

   .. code-block:: c

      retval = SUNAdaptController_SetParams_I(C, 0.95);
