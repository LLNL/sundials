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

.. _SUNAdaptController.ExpGus:

The SUNAdaptController_ExpGus Module
======================================

The explicit Gustafsson implementation of the SUNAdaptController class,
SUNAdaptController_ExpGus, implements a controller proposed by K. Gustafsson for
explicit Runge--Kutta methods in :cite:p:`Gust:91`.  This controller has the
form

.. math::
   h' \;=\; \begin{cases}
      h_1\; \varepsilon_1^{-1/(p+1)}, &\quad\text{on the first step}, \\
      h_n\; \varepsilon_n^{-k_1/(p+1)}\;
        \left(\dfrac{\varepsilon_n}{\varepsilon_{n-1}}\right)^{k_2/(p+1)}, &
      \quad\text{on subsequent steps},
   \end{cases}

with default values :math:`k_1=0.367` and :math:`k_2=0.268`, and where :math:`p`
is the global order of the time integration method. In this estimate, a floor of
:math:`\varepsilon > 10^{-10}` is enforced to avoid division-by-zero errors.

The SUNAdaptController_ExpGus controller is implemented as a derived SUNAdaptController class,
and defines its *content* field as:

.. code-block:: c

   struct _SUNAdaptControllerContent_ExpGus {
     sunrealtype k1;
     sunrealtype k2;
     sunrealtype bias;
     sunrealtype ep;
     sunbooleantype firststep;
   };

These entries of the *content* field contain the following information:

* ``k1, k2`` - controller parameters above.

* ``bias`` - error bias factor, that converts from an input temporal error
  estimate via :math:`\varepsilon = \text{bias}*\text{dsm}`.

* ``ep`` - storage for the previous error estimate, :math:`\varepsilon_{n-1}`.

* ``firststep`` - flag indicating whether a step has successfully completed, in which
  case the formula above transitions from :math:`h_1` to :math:`h_n`.

The header file to be included when using this module is
``sunadaptcontroller/sunadaptcontroller_expgus.h``.


The SUNAdaptController_ExpGus class provides implementations of all operations
relevant to a `SUN_ADAPTCONTROLLER_H` controller listed in
:numref:`SUNAdaptController.Description.operations`. The
SUNAdaptController_ExpGus class also provides the following additional user-callable
routines:


.. c:function:: SUNAdaptController SUNAdaptController_ExpGus(SUNContext sunctx)

   This constructor function creates and allocates memory for a SUNAdaptController_ExpGus
   object, and inserts its default parameters.

   :param sunctx: the current :c:type:`SUNContext` object.
   :return: if successful, a usable :c:type:`SUNAdaptController` object; otherwise it will return ``NULL``.

   Usage:

   .. code-block:: c

      SUNAdaptController C = SUNAdaptController_ExpGus(sunctx);

.. c:function:: int SUNAdaptController_SetParams_ExpGus(SUNAdaptController C, sunrealtype k1, sunrealtype k2)

   :param C: the SUNAdaptController_ExpGus object.
   :param k1: parameter used within the controller time step estimate (only stored if non-negative).
   :param k2: parameter used within the controller time step estimate (only stored if non-negative).
   :return: error code indication success or failure (see :numref:`SUNAdaptController.Description.errorCodes`).

   Usage:

   .. code-block:: c

      retval = SUNAdaptController_SetParams_ExpGus(C, 0.4, 0.25);
