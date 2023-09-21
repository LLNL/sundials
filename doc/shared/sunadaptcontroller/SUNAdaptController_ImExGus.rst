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

.. _SUNAdaptController.ImExGus:

The SUNAdaptController_ImExGus Module
======================================

The ImEx Gustafsson implementation of the SUNAdaptController class, SUNAdaptController_ImExGus,
implements a combination of two adaptivity controllers proposed
by K. Gustafsson.  His "explicit" controller was proposed in :cite:p:`Gust:91`,
is primarily useful with explicit Runge--Kutta methods, and has the form

.. math::
   h' \;=\; \begin{cases}
      h_1\; \varepsilon_1^{-1/p}, &\quad\text{on the first step}, \\
      h_n\; \varepsilon_n^{-k_1/p}\;
        \left(\dfrac{\varepsilon_n}{\varepsilon_{n-1}}\right)^{k_2/p}, &
      \quad\text{on subsequent steps}.
   \end{cases}
   :label: expGusController

The default values of :math:`k_1` and :math:`k_2` are 0.367 and 0.268,
respectively.  Gustafsson proposed a controller for implicit Runge--Kutta
methods in :cite:p:`Gust:94`, with the form

.. math::
   h' = \begin{cases}
      h_1 \varepsilon_1^{-1/p}, &\quad\text{on the first step}, \\
      h_n \left(\dfrac{h_n}{h_{n-1}}\right) \varepsilon_n^{-k_1/p}
        \left(\dfrac{\varepsilon_n}{\varepsilon_{n-1}}\right)^{-k_2/p}, &
      \quad\text{on subsequent steps}.
   \end{cases}
   :label: impGusController

The default parameter values are :math:`k_1 = 0.98` and :math:`k_2 = 0.95`.

The SUNAdaptController_ImExGus controller implements both formulas
:eq:`expGusController` and :eq:`impGusController`, and sets its recommended step
size as the minimum of these two.  It is implemented as a derived SUNAdaptController
class, and defines its *content* field as:

.. code-block:: c

   struct _SUNAdaptControllerContent_ImExGus {
     realtype k1e;
     realtype k2e;
     realtype k1i;
     realtype k2i;
     realtype bias;
     realtype ep;
     realtype hp;
     int p;
     sunbooleantype pq;
     sunbooleantype firststep;
   };

These entries of the *content* field contain the following information:

* ``k1e, k2e`` - explicit controller parameters used in :eq:`expGusController`.

* ``k1i, k2i`` - implicit controller parameters used in :eq:`impGusController`.

* ``bias`` - error bias factor, that converts from an input temporal error
  estimate via :math:`\varepsilon = \text{bias}*\text{dsm}`.

* ``ep`` - storage for the previous error estimate, :math:`\varepsilon_{n-1}`.

* ``hp`` - storage for the previous step size, :math:`h_{n-1}`.

* ``p`` - asymptotic order to use in error control.

* ``pq`` - flag indicating whether ``p`` corresponds to the order of accuracy
  for the time integration method (``SUNTRUE``) or the embedding (``SUNFALSE``).


The header file to be included when using this module is
``sunadaptcontroller/sunadaptcontroller_imexgus.h``.


The SUNAdaptController_ImExGus class provides implementations of all operations
relevant to a `SUN_ADAPTCONTROLLER_H` controller listed in
:numref:`SUNAdaptController.Description.operations`. The
SUNAdaptController_ImExGus class also provides the following additional user-callable
routines:


.. c:function:: SUNAdaptController SUNAdaptController_ImExGus(SUNContext sunctx)

   This constructor function creates and allocates memory for a
   SUNAdaptController_ImExGus object, and inserts its default parameters.  The only
   argument is the SUNDIALS context object.  Upon successful completion it will
   return a :c:type:`SUNAdaptController` object; otherwise it will return ``NULL``.


.. c:function:: int SUNAdaptController_SetParams_ImExGus(SUNAdaptController C, sunbooleantype pq, realtype k1e, realtype k2e, realtype k1i, realtype k2i)

   This user-callable function provides control over the relevant parameters
   above.  The *pq* input is stored directly.  The *k1e*, *k2e*, *k1i* and *k2i*
   parameters are only stored if the corresponding input is non-negative.  Upon
   completion, this returns ``SUNADAPTCONTROLLER_SUCCESS``.
