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

.. _SUNControl.I:

The SUNControl_I Module
======================================

The I implementation of the SUNControl class, SUNControl_I, implements a
standard I temporal adaptivity controller.  This is the standard time adaptivity
control algorithm in use by most publicly-available ODE solver codes.  It bases
the prospective time step estimate entirely off of the current local error
estimate,

.. math::
   h' \;=\; h_n\; \varepsilon_n^{-k_1/p}.

By default the constant :math:`k_1=1`.

This is implemented as a derived SUNControl class, and defines its *content*
field as:

.. code-block:: c

   struct _SUNControlContent_PI {
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
``suncontrol/suncontrol_i.h``.

The SUNControl_I class provides implementations of all controller operations
listed in :numref:`SUNControl.Description.operations`. The SUNControl_I class
also provides the following additional user-callable routines:


.. c:function:: SUNControl SUNControlI(SUNContext sunctx)

   This constructor function creates and allocates memory for a SUNControl_I
   object, and inserts its default parameters.  The only argument is the
   SUNDIALS context object.  Upon successful completion it will return a
   :c:type:`SUNControl` object; otherwise it will return ``NULL``.


.. c:function:: int SUNControlI_SetParams(SUNControl C, sunbooleantype pq, realtype k1)

   This user-callable function provides control over the relevant parameters
   above.  The *pq* input is stored directly.  The *k1* parameter is only stored
   if the input is non-negative.  Upon completion, this returns
   ``SUNCONTROL_SUCCESS``.
