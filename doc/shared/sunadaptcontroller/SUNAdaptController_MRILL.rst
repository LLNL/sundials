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

.. _SUNAdaptController.MRILL:

The SUNAdaptController_MRILL Module
=======================================

.. versionadded:: x.y.z

The MRILL implementation of the SUNAdaptController class,
SUNAdaptController_MRILL, implements the "linear-linear" multirate temporal
controller proposed in :cite:p:`Fish:23`.  This controller has the form

.. math::
   H' &= H_n \left(\frac{H_n}{H_{n-1}}\right) \left(\varepsilon^s_n\right)^{\alpha_1} \left(\varepsilon^s_{n-1}\right)^{\alpha_2},\\
   M' &= M_n \left(\frac{M_n}{M_{n-1}}\right) \left(\varepsilon^s_n\right)^{\beta_{11}} \left(\varepsilon^s_{n-1}\right)^{\beta_{12}} \left(\varepsilon^f_n\right)^{\beta_{21}} \left(\varepsilon^f_{n-1}\right)^{\beta_{22}},\\
   h' &= H'/M'

where :math:`M_n = \left\lceil\frac{H_n}{h_n}\right\rceil`,
:math:`M_{n-1} = \left\lceil\frac{H_{n-1}}{h_{n-1}}\right\rceil`,
:math:`\alpha_1 = \frac{k_{11}+k_{12}}{2P}`, :math:`\alpha_2 = -\frac{k_{11}}{2P}`,
:math:`\beta_{11} = \frac{(p+1)(k_{11}+k_{12})}{2Pp}`,
:math:`\beta_{12} = -\frac{(p+1)k_{11}}{2Pp}`,
:math:`\beta_{21} = -\frac{k_{21}+k_{22}}{2p}`,  and
:math:`\beta_{22} = \frac{k_{21}}{2p}`, and where :math:`P` and :math:`p` are the global
orders of accuracy for the slow and fast time integration methods, respectively.
The default parameter values are :math:`k_{11} = 0.82`, :math:`k_{12} = 0.54`,
:math:`k_{21} = 0.94` and :math:`k_{22} = 0.9`.  In these estimate, a floor of
:math:`\varepsilon^*_* > 10\epsilon_{mach}` is enforced to avoid division-by-zero errors,
where :math:`\epsilon_{mach}` is floating point roundoff for the current working precision.
During the first step (when :math:`\varepsilon^s_{n-1}`,
:math:`\varepsilon^f_{n-1}`, :math:`H_{n-1}`, and :math:`h_{n-1}` are unavailable),
the corresponding temporal errors and step size and multirate ratios,
:math:`\varepsilon^s_{n-1}`, :math:`\varepsilon^f_{n-1}`,
:math:`\left(\frac{H_n}{H_{n-1}}\right)` and :math:`\left(\frac{M_n}{M_{n-1}}\right)`, are
set to 1, effectively removing them from the estimation formula.

The SUNAdaptController_MRILL controller is implemented as a derived
SUNAdaptController class having type ``SUN_ADAPTCONTROLLER_MRI_H``, and its
*content* field is:

.. code-block:: c

   struct _SUNAdaptControllerContent_MRILL {
     sunrealtype k11;
     sunrealtype k12;
     sunrealtype k21;
     sunrealtype k22;
     sunrealtype bias;
     sunrealtype esp;
     sunrealtype efp;
     sunrealtype hsp;
     sunrealtype hfp;
     int p;
     sunbooleantype firststep;
   };

These entries of the *content* field contain the following information:

* ``k11, k12, k21, k22`` - controller parameters above.

* ``bias`` - error bias factor, that converts from input temporal error
  estimates via :math:`\varepsilon^s = \text{bias}*\text{DSM}` and
  :math:`\varepsilon^f = \text{bias}*\text{dsm}`.  The default bias value is 1.5.

* ``esp, efp`` - storage for the previous slow and fast error estimates,
  :math:`\varepsilon^s_{n-1}` and :math:`\varepsilon^f_{n-1}`.

* ``hsp, hfp`` - storage for the previous slow and fast step sizes, :math:`H_{n-1}`
  and :math:`h_{n-1}`.

* ``p`` - global order of accuracy for the fast time scale solver.

* ``firststep`` - counter to handle first (where the previous step
  sizes and errors are unavailable).

The header file to be included when using this module is
``sunadaptcontroller/sunadaptcontroller_mrill.h``.

The SUNAdaptController_MRILL class provides implementations of all operations
relevant to a ``SUN_ADAPTCONTROLLER_MRI_H`` controller listed in
:numref:`SUNAdaptController.Description.operations`. This class
also provides the following additional user-callable routines:


.. c:function:: SUNAdaptController SUNAdaptController_MRILL(SUNContext sunctx, int p)

   This constructor creates and allocates memory for a SUNAdaptController_MRILL
   object, and inserts its default parameters.

   :param sunctx: the current :c:type:`SUNContext` object.
   :param p: the global order of accuracy for the fast time scale solver.
   :return: if successful, a usable :c:type:`SUNAdaptController` object;
            otherwise it will return ``NULL``.

   Usage:

   .. code-block:: c

      SUNAdaptController C = SUNAdaptController_MRILL(sunctx, 3);

.. c:function:: SUNErrCode SUNAdaptController_SetParams_MRILL(SUNAdaptController C, sunrealtype k11, sunrealtype k12, sunrealtype k21, sunrealtype k22)

   This user-callable function provides control over the relevant parameters
   above.  This should be called *before* the time integrator is called to evolve
   the problem.

   :param C: the SUNAdaptController_MRILL object.
   :param k11: parameter used within the controller time step estimate.
   :param k12: parameter used within the controller time step estimate.
   :param k21: parameter used within the controller time step estimate.
   :param k22: parameter used within the controller time step estimate.
   :return: :c:type:`SUNErrCode` indicating success or failure.

   Usage:

   .. code-block:: c

      retval = SUNAdaptController_SetParams_MRILL(C, 0.82, 0.54, 0.94, 0.9);
