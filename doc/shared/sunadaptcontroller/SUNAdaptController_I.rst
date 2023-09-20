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
   h' \;=\; h_n\; \varepsilon_n^{-k_1/p}.

By default the constant :math:`k_1=1`.

This is implemented as a derived SUNAdaptController class, and defines its *content*
field as:

.. code-block:: c

   struct SUNAdaptControllerContent_PI_ {
     realtype k1;
     realtype bias;
     int p;
     sunbooleantype pq;
   };

These entries of the *content* field contain the following information:

* ``k1`` - controller parameter above.

* ``bias`` - error bias factor, that converts from an input temporal error
  estimate via :math:`\varepsilon = \text{bias}*\text{dsm}`.

* ``p`` - asymptotic order to use in error control.

* ``pq`` - flag indicating whether ``p`` corresponds to the order of accuracy
  for the time integration method (``SUNTRUE``) or the embedding (``SUNFALSE``).


The header file to be included when using this module is
``sunadaptcontroller/sunadaptcontroller_i.h``.

The SUNAdaptController_I class provides implementations of all operations
relevant to a `SUN_ADAPTCONTROLLER_H` controller listed in
:numref:`SUNAdaptController.Description.operations`. The SUNAdaptController_I class
also provides the following additional user-callable routines:


.. c:function:: SUNAdaptController SUNAdaptControllerI(SUNContext sunctx)

   This constructor function creates and allocates memory for a SUNAdaptController_I
   object, and inserts its default parameters.  The only argument is the
   SUNDIALS context object.  Upon successful completion it will return a
   :c:type:`SUNAdaptController` object; otherwise it will return ``NULL``.


.. c:function:: int SUNAdaptControllerI_SetParams(SUNAdaptController C, sunbooleantype pq, realtype k1)

   This user-callable function provides control over the relevant parameters
   above.  The *pq* input is stored directly.  The *k1* parameter is only stored
   if the input is non-negative.  Upon completion, this returns
   ``SUNADAPTCONTROLLER_SUCCESS``.
