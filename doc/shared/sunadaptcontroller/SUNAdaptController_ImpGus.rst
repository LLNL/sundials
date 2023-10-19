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

.. _SUNAdaptController.ImpGus:

The SUNAdaptController_ImpGus Module
======================================

The implicit Gustafsson implementation of the SUNAdaptController class,
SUNAdaptController_ImpGus, implements a controller proposed by K. Gustafsson for
implicit Runge--Kutta methods in :cite:p:`Gust:94`.  This controller has the
form

.. math::
   h' = \begin{cases}
      h_1 \varepsilon_1^{-1/ord}, &\quad\text{on the first step}, \\
      h_n \left(\dfrac{h_n}{h_{n-1}}\right) \varepsilon_n^{-k_1/ord}
        \left(\dfrac{\varepsilon_n}{\varepsilon_{n-1}}\right)^{-k_2/ord}, &
      \quad\text{on subsequent steps},
   \end{cases}

with default parameter values :math:`k_1 = 0.98` and :math:`k_2 = 0.95`, and where
:math:`ord = p+1+adj`, where both :math:`p` and :math:`adj` are described below. In
this estimate, a floor of :math:`\varepsilon > 10^{-10}` is enforced to avoid
division-by-zero errors.

The SUNAdaptController_ImpGus controller is implemented as a derived SUNAdaptController class,
and defines its *content* field as:

.. code-block:: c

   struct _SUNAdaptControllerContent_ImpGus {
     sunrealtype k1;
     sunrealtype k2;
     sunrealtype bias;
     sunrealtype ep;
     sunrealtype hp;
     int p;
     int adj;
     sunbooleantype firststep;
   };

These entries of the *content* field contain the following information:

* ``k1, k2`` - controller parameters above.

* ``bias`` - error bias factor, that converts from an input temporal error
  estimate via :math:`\varepsilon = \text{bias}*\text{dsm}`.

* ``ep`` - storage for the previous error estimate, :math:`\varepsilon_{n-1}`.

* ``hp`` - storage for the previous step size, :math:`h_{n-1}`.

* ``p`` - asymptotic order to use in error control.  This is provided by
  the time integrator, corresponding to the order of accuracy for the time
  integration method, the embedding, or the minimum of the two.

* ``adj`` - order of accuracy adjustment to use within the controller [default ``-1``].

* ``firststep`` - flag indicating whether any time steps have completed
  successfully (and thus to transition from :math:`h_1` to :math:`h_n` in
  the formula above).

The header file to be included when using this module is
``sunadaptcontroller/sunadaptcontroller_impgus.h``.


The SUNAdaptController_ImpGus class provides implementations of all operations
relevant to a `SUN_ADAPTCONTROLLER_H` controller listed in
:numref:`SUNAdaptController.Description.operations`. This class
also provides the following additional user-callable routines:


.. c:function:: SUNAdaptController SUNAdaptController_ImpGus(SUNContext sunctx)

   This constructor function creates and allocates memory for a SUNAdaptController_ImpGus
   object, and inserts its default parameters.

   :param sunctx: the current :c:type:`SUNContext` object.
   :return: if successful, a usable :c:type:`SUNAdaptController` object; otherwise it will return ``NULL``.

   Usage:

   .. code-block:: c

      SUNAdaptController C = SUNAdaptController_ImpGus(sunctx);

.. c:function:: int SUNAdaptController_SetParams_ImpGus(SUNAdaptController C, sunrealtype k1, sunrealtype k2)

   This user-callable function provides control over the relevant parameters
   above.  This should be called *before* the time integrator is called to evolve
   the problem.

   :param C: the SUNAdaptController_ImpGus object.
   :param k1: parameter used within the controller time step estimate (only stored if non-negative).
   :param k2: parameter used within the controller time step estimate (only stored if non-negative).
   :return: error code indication success or failure (see :numref:`SUNAdaptController.Description.errorCodes`).

   Usage:

   .. code-block:: c

      retval = SUNAdaptController_SetParams_ImpGus(C, 1.0, 0.9);
