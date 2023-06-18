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

.. _SUNControl.PI:

The SUNControl_PI Module
======================================

The PI implementation of the SUNControl class, SUNControl_PI, implements a standard PI temporal adaptivity controller.  Like with SUNControl_PID, the PI controller derives from those found in :cite:p:`KenCarp:03`, :cite:p:`Sod:98`, :cite:p:`Sod:03` and :cite:p:`Sod:06`, but it differs in that it only uses the two most recent step sizes in its adaptivity algorithm,

.. math::
   h' \;=\; h_n\; \varepsilon_n^{-k_1/p}\; \varepsilon_{n-1}^{k_2/p}.

where the constants :math:`k_1` and :math:`k_2` default to 0.8 and 0.31, respectively.

This is implemented as a derived SUNControl class, and defines its *content* field as:

.. code-block:: c

   struct _SUNControlContent_PI {
     realtype k1;
     realtype k2;
     realtype bias;
     realtype ep;
     int p;
     sunbooleantype pq;
   };

These entries of the *content* field contain the following information:

* ``k1, k2`` - controller parameters above.

* ``bias`` - error bias factor, that converts from an input temporal error estimate via :math:`\varepsilon = \text{bias}*\text{dsm}`.

* ``ep`` - storage for the previous error estimate, :math:`\varepsilon_{n-1}`.

* ``p`` - asymptotic order to use in error control.

* ``pq`` - flag indicating whether ``p`` corresponds to the order of accuracy for the time integration method (``SUNTRUE``) or the embedding (``SUNFALSE``).


The header file to be included when using this module is ``suncontrol/suncontrol_pi.h``.

The SUNControl_PI class provides implementations of all controller operations listed in :numref:`SUNControl.Description.operations`. The SUNControl_PI class also provides the following additional user-callable routines:


.. c:function:: SUNControl SUNControlPI(SUNContext sunctx)

   This constructor function creates and allocates memory for a SUNControl_PI object, and inserts its default parameters.  The only argument is the SUNDIALS context object.  Upon successful completion it will return a :c:type:`SUNControl` object; otherwise it will return ``NULL``.


.. c:function:: int SUNControlPI_SetParams(SUNControl C, sunbooleantype pq, realtype k1, realtype k2)

   This user-callable function provides control over the relevant parameters above.  The *pq* input is stored directly.  The *k1* and *k2* are only stored if the corresponding input is non-negative.  Upon completion, this returns ``SUNCONTROL_SUCCESS``.
