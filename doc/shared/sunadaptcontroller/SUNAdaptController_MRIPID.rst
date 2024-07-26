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

.. _SUNAdaptController.MRIPID:

The SUNAdaptController_MRIPID Module
=======================================

.. versionadded:: x.y.z

The MRIPID implementation of the SUNAdaptController class,
SUNAdaptController_MRIPID, implements the "PIDMR" multirate temporal
controller proposed in :cite:p:`Fish:23`.  This controller has the form

.. math::
   H' &= H_n \left(\varepsilon^s_n\right)^{\alpha_1} \left(\varepsilon^s_{n-1}\right)^{\alpha_2} \left(\varepsilon^s_{n-2}\right)^{\alpha_3},\\
   M' &= M_n \left(\varepsilon^s_n\right)^{\beta_{11}} \left(\varepsilon^s_{n-1}\right)^{\beta_{12}} \left(\varepsilon^s_{n-2}\right)^{\beta_{13}} \left(\varepsilon^f_n\right)^{\beta_{21}} \left(\varepsilon^f_{n-1}\right)^{\beta_{22}} \left(\varepsilon^f_{n-2}\right)^{\beta_{23}},\\
   h' &= H'/M'

where :math:`M_n = \left\lceil\frac{H_n}{h_n}\right\rceil`,
:math:`\alpha_1 = \frac{k_{11}+k_{12}+k_{13}}{3P}`,
:math:`\alpha_2 = -\frac{k_{11}+k_{12}}{3P}`,
:math:`\alpha_3 = \frac{k_{11}}{3P}`,
:math:`\beta_{11} = \frac{(p+1)(k_{11}+k_{12}+k_{13})}{3Pp}`,
:math:`\beta_{12} = -\frac{(p+1)(k_{11}+k_{12})}{3Pp}`,
:math:`\beta_{13} = \frac{(p+1)k_{11}}{3Pp}`,
:math:`\beta_{21} = -\frac{k_{21}+k_{22}+k_{23}}{3p}`,
:math:`\beta_{22} = \frac{k_{21}+k_{22}}{3p}`,
and :math:`\beta_{23} = -\frac{k_{21}}{3p}`, and where :math:`P` and :math:`p`
are the global orders of accuracy for the slow and fast time integration methods,
respectively. The default parameter values are :math:`k_{11} = 0.34`, :math:`k_{12} = 0.1`,
:math:`k_{13} = 0.78`, :math:`k_{21} = 0.46`, :math:`k_{22} = 0.42`, and :math:`k_{23} = 0.74`.
In these estimate, a floor of :math:`\varepsilon^*_* > 10\epsilon_{mach}` is enforced to avoid
division-by-zero errors, where :math:`\epsilon_{mach}` is floating point roundoff for the
current working precision. During the first steps (when :math:`\varepsilon^s_{n-1}`,
:math:`\varepsilon^s_{n-2}`, :math:`\varepsilon^f_{n-1}`, and :math:`\varepsilon^f_{n-2}`
may be unavailable), those corresponding values are set to 1, effectively removing them from
the estimation formula.

The SUNAdaptController_MRIPID controller is implemented as a derived
SUNAdaptController class having type ``SUN_ADAPTCONTROLLER_MRI_H``, and its
*content* field is:

.. code-block:: c

   struct _SUNAdaptControllerContent_MRIPID {
     sunrealtype k11;
     sunrealtype k12;
     sunrealtype k13;
     sunrealtype k21;
     sunrealtype k22;
     sunrealtype k23;
     sunrealtype bias;
     sunrealtype esp;
     sunrealtype efp;
     sunrealtype espp;
     sunrealtype efpp;
     int p;
   };

These entries of the *content* field contain the following information:

* ``k11, k12, k13, k21, k22, k23`` - controller parameters above.

* ``bias`` - error bias factor, that converts from input temporal error
  estimates via :math:`\varepsilon^s = \text{bias}*\text{DSM}` and
  :math:`\varepsilon^f = \text{bias}*\text{dsm}`.  The default bias value is 1.5.

* ``esp, efp, espp, efpp`` - storage for the two previous slow and fast error
  estimates, :math:`\varepsilon^s_{n-1}`, :math:`\varepsilon^f_{n-1}`,
  :math:`\varepsilon^s_{n-2}`, and  :math:`\varepsilon^f_{n-2}`.

* ``p`` - global order of accuracy for the fast time scale solver.

The header file to be included when using this module is
``sunadaptcontroller/sunadaptcontroller_mripid.h``.

The SUNAdaptController_MRIPID class provides implementations of all operations
relevant to a ``SUN_ADAPTCONTROLLER_MRI_H`` controller listed in
:numref:`SUNAdaptController.Description.operations`. This class
also provides the following additional user-callable routines:


.. c:function:: SUNAdaptController SUNAdaptController_MRIPID(SUNContext sunctx, int p)

   This constructor creates and allocates memory for a SUNAdaptController_MRIPID
   object, and inserts its default parameters.

   :param sunctx: the current :c:type:`SUNContext` object.
   :param p: the global order of accuracy for the fast time scale solver.
   :return: if successful, a usable :c:type:`SUNAdaptController` object;
            otherwise it will return ``NULL``.

   Usage:

   .. code-block:: c

      SUNAdaptController C = SUNAdaptController_MRIPID(sunctx, 3);

.. c:function:: SUNErrCode SUNAdaptController_SetParams_MRIPID(SUNAdaptController C, sunrealtype k11, sunrealtype k12, sunrealtype k13, sunrealtype k21, sunrealtype k22, sunrealtype k23)

   This user-callable function provides control over the relevant parameters
   above.  This should be called *before* the time integrator is called to evolve
   the problem.

   :param C: the SUNAdaptController_MRIPID object.
   :param k11: parameter used within the controller time step estimate.
   :param k12: parameter used within the controller time step estimate.
   :param k13: parameter used within the controller time step estimate.
   :param k21: parameter used within the controller time step estimate.
   :param k22: parameter used within the controller time step estimate.
   :param k23: parameter used within the controller time step estimate.
   :return: :c:type:`SUNErrCode` indicating success or failure.

   Usage:

   .. code-block:: c

      retval = SUNAdaptController_SetParams_MRIPID(C, 0.34 0.1, 0.78, 0.46, 0.42, 0.74);
