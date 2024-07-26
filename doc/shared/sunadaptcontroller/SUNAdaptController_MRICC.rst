..
   Programmer(s): Daniel R. Reynolds @ SMU
   ----------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2002-2024, Lawrence Livermore National Security
   and Southern Methodist University.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   ----------------------------------------------------------------

.. _SUNAdaptController.MRICC:

The SUNAdaptController_MRICC Module
=======================================

.. versionadded:: x.y.z

The MRICC implementation of the SUNAdaptController class,
SUNAdaptController_MRICC, implements the "constant-constant" multirate temporal
controller proposed in :cite:p:`Fish:23`.  This controller has the form

.. math::
   H' &= H_n \left(\varepsilon^s_n\right)^{\alpha},\\
   M' &= \left\lceil\frac{H_n}{h_n}\right\rceil \left(\varepsilon^s_n\right)^{\beta_1} \left(\varepsilon^f_n\right)^{\beta_2},\\
   h' &= H'/M'

where :math:`\alpha = \frac{k_1}{P}`, :math:`\beta_1 = \frac{(p+1)k_1}{Pp}`, and
:math:`\beta_2 = -\frac{k_2}{p}`, and where :math:`P` and :math:`p` are the global
orders of accuracy for the slow and fast time integration methods, respectively.
The default parameter values are :math:`k_1 = 0.42` and :math:`k_2 = 0.44`.  In
these estimate, a floor of :math:`\varepsilon^*_* > 10\epsilon_{mach}` is enforced
to avoid division-by-zero errors, where :math:`\epsilon_{mach}` is floating point
roundoff for the current working precision.

The SUNAdaptController_MRICC controller is implemented as a derived
SUNAdaptController class having type ``SUN_ADAPTCONTROLLER_MRI_H``, and its
*content* field is:

.. code-block:: c

   struct _SUNAdaptControllerContent_MRICC {
     sunrealtype k1;
     sunrealtype k2;
     sunrealtype bias;
     int p;
   };

These entries of the *content* field contain the following information:

* ``k1, k2`` - controller parameters above.

* ``bias`` - error bias factor, that converts from input temporal error
  estimates via :math:`\varepsilon^s = \text{bias}*\text{DSM}` and
  :math:`\varepsilon^f = \text{bias}*\text{dsm}`.  The default bias value is 1.5.

* ``p`` - global order of accuracy for the fast time scale solver.

The header file to be included when using this module is
``sunadaptcontroller/sunadaptcontroller_mricc.h``.

The SUNAdaptController_MRICC class provides implementations of all operations
relevant to a ``SUN_ADAPTCONTROLLER_MRI_H`` controller listed in
:numref:`SUNAdaptController.Description.operations`. This class
also provides the following additional user-callable routines:


.. c:function:: SUNAdaptController SUNAdaptController_MRICC(SUNContext sunctx, int p)

   This constructor creates and allocates memory for a SUNAdaptController_MRICC
   object, and inserts its default parameters.

   :param sunctx: the current :c:type:`SUNContext` object.
   :param p: the global order of accuracy for the fast time scale solver.
   :return: if successful, a usable :c:type:`SUNAdaptController` object;
            otherwise it will return ``NULL``.

   Usage:

   .. code-block:: c

      SUNAdaptController C = SUNAdaptController_MRICC(sunctx, 3);

.. c:function:: SUNErrCode SUNAdaptController_SetParams_MRICC(SUNAdaptController C, sunrealtype k1, sunrealtype k2)

   This user-callable function provides control over the relevant parameters
   above.  This should be called *before* the time integrator is called to evolve
   the problem.

   :param C: the SUNAdaptController_MRICC object.
   :param k1: parameter used within the controller time step estimate.
   :param k2: parameter used within the controller time step estimate.
   :return: :c:type:`SUNErrCode` indicating success or failure.

   Usage:

   .. code-block:: c

      retval = SUNAdaptController_SetParams_MRICC(C, 0.42, 0.44);
