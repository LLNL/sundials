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

.. _SUNAdaptController.PID:

The SUNAdaptController_PID Module
======================================

The PID implementation of the SUNAdaptController class, SUNAdaptController_PID, implements a
standard PID temporal adaptivity controller.  It derives from those found in
:cite:p:`KenCarp:03`, :cite:p:`Sod:98`, :cite:p:`Sod:03` and :cite:p:`Sod:06`,
and uses three temporal error estimates, :math:`\varepsilon_n`,
:math:`\varepsilon_{n-1}` and :math:`\varepsilon_{n-2}` in determination of a
prospective step size,

.. math::
   h' \;=\; h_n\; \varepsilon_n^{-k_1/(p+1)}\; \varepsilon_{n-1}^{k_2/(p+1)}\;
        \varepsilon_{n-2}^{-k_3/(p+1)},

where the constants :math:`k_1`, :math:`k_2` and :math:`k_3` default to 0.58,
0.21 and 0.1, respectively, and :math:`p` is described below.
In this estimate, a floor of :math:`\varepsilon > 10^{-10}` is enforced to avoid
division-by-zero errors.

This is implemented as a derived SUNAdaptController class, and defines its *content*
field as:

.. code-block:: c

   struct _SUNAdaptControllerContent_PID {
     sunrealtype k1;
     sunrealtype k2;
     sunrealtype k3;
     sunrealtype bias;
     sunrealtype ep;
     sunrealtype epp;
     int p;
   };

These entries of the *content* field contain the following information:

* ``k1, k2, k3`` - controller parameters above.

* ``bias`` - error bias factor, that converts from an input temporal error
  estimate via :math:`\varepsilon = \text{bias}*\text{dsm}`.

* ``ep, epp`` - storage for the two previous error estimates,
  :math:`\varepsilon_{n-1}` and :math:`varepsilon_{n-2}`.

* ``p`` - asymptotic order to use in error control (provided by the time integrator).


The header file to be included when using this module is
``sunadaptcontroller/sunadaptcontroller_pid.h``.


The SUNAdaptController_PID class provides implementations of all operations
relevant to a `SUN_ADAPTCONTROLLER_H` controller listed in
:numref:`SUNAdaptController.Description.operations`. The SUNAdaptController_PID class
also provides the following additional user-callable routines:


.. c:function:: SUNAdaptController SUNAdaptController_PID(SUNContext sunctx)

   This constructor function creates and allocates memory for a SUNAdaptController_PID
   object, and inserts its default parameters.

   :param sunctx: the current :c:type:`SUNContext` object.
   :return: if successful, a usable :c:type:`SUNAdaptController` object; otherwise it will return ``NULL``.

   Usage:

   .. code-block:: c

      SUNAdaptController C = SUNAdaptController_PID(sunctx);

.. c:function:: int SUNAdaptController_SetParams_PID(SUNAdaptController C, sunrealtype k1, sunrealtype k2, sunrealtype k3)

   This user-callable function provides control over the relevant parameters
   above.  This should be called *before* the time integrator is called to evolve
   the problem.

   :param C: the SUNAdaptController_PID object.
   :param k1: parameter used within the controller time step estimate (only stored if non-negative).
   :param k2: parameter used within the controller time step estimate (only stored if non-negative).
   :param k3: parameter used within the controller time step estimate (only stored if non-negative).
   :return: error code indication success or failure (see :numref:`SUNAdaptController.Description.errorCodes`).

   Usage:

   .. code-block:: c

      retval = SUNAdaptController_SetParams_PID(C, 0.6, 0.2, -1.0);
