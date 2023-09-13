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

.. _SUNControl.PID:

The SUNControl_PID Module
======================================

The PID implementation of the SUNControl class, SUNControl_PID, implements a
standard PID temporal adaptivity controller.  It derives from those found in
:cite:p:`KenCarp:03`, :cite:p:`Sod:98`, :cite:p:`Sod:03` and :cite:p:`Sod:06`,
and uses three temporal error estimates, :math:`\varepsilon_n`,
:math:`\varepsilon_{n-1}` and :math:`\varepsilon_{n-2}` in determination of a
prospective step size,

.. math::
   h' \;=\; h_n\; \varepsilon_n^{-k_1/p}\; \varepsilon_{n-1}^{k_2/p}\;
        \varepsilon_{n-2}^{-k_3/p},

where the constants :math:`k_1`, :math:`k_2` and :math:`k_3` default to 0.58,
0.21 and 0.1, respectively. In this estimate, a floor of :math:`\varepsilon >
10^{-10}` is enforced to avoid division-by-zero errors.

This is implemented as a derived SUNControl class, and defines its *content*
field as:

.. code-block:: c

   struct _SUNControlContent_PID {
     realtype k1;
     realtype k2;
     realtype k3;
     realtype bias;
     realtype ep;
     realtype epp;
     int p;
     sunbooleantype pq;
   };

These entries of the *content* field contain the following information:

* ``k1, k2, k3`` - controller parameters above.

* ``bias`` - error bias factor, that converts from an input temporal error
  estimate via :math:`\varepsilon = \text{bias}*\text{dsm}`.

* ``ep, epp`` - storage for the two previous error estimates,
  :math:`\varepsilon_{n-1}` and :math:`varepsilon_{n-2}`.

* ``p`` - asymptotic order to use in error control.

* ``pq`` - flag indicating whether ``p`` corresponds to the order of accuracy
  for the time integration method (``SUNTRUE``) or the embedding (``SUNFALSE``).


The header file to be included when using this module is
``suncontrol/suncontrol_pid.h``.


The SUNControl_PID class provides implementations of all controller operations
listed in :numref:`SUNControl.Description.operations`. The SUNControl_PID class
also provides the following additional user-callable routines:


.. c:function:: SUNControl SUNControlPID(SUNContext sunctx)

   This constructor function creates and allocates memory for a SUNControl_PID
   object, and inserts its default parameters.  The only argument is the
   SUNDIALS context object.  Upon successful completion it will return a
   :c:type:`SUNControl` object; otherwise it will return ``NULL``.


.. c:function:: int SUNControlPID_SetParams(SUNControl C, sunbooleantype pq, realtype k1, realtype k2, realtype k3)

   This user-callable function provides control over the relevant parameters
   above.  The *pq* input is stored directly.  The *k1*, *k2* and *k3* are only
   stored if the corresponding input is non-negative.  Upon completion, this
   returns ``SUNCONTROL_SUCCESS``.
