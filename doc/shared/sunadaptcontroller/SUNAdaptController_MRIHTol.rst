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

.. _SUNAdaptController.MRIHTol:

The SUNAdaptController_MRIHTol Module
=======================================

.. versionadded:: x.y.z

Mathematical motivation
-----------------------

The MRIHTol implementation of the SUNAdaptController class,
SUNAdaptController_MRIHTol, implements a general structure for telescopic
multirate temporal control.  A SUNAdaptController_MRIHTol object is constructed
using two single-rate controller objects, *HControl* and *TolControl*.  The
MRIHTol controller assumes that overall solution error at a given time scale
results from two types of error:

#. "Slow" temporal errors introduced at the current time scale,

   .. math::
      \varepsilon^s_{n} = C(t_n) \left(h_n^s\right)^{P+1},
      :label: slow_error_assumption

   where :math:`C(t)` is independent of the current time scale step size :math:`h^s`
   but may vary in time.

#. "Fast" errors introduced through calls to the next-fastest ("inner") solver,
   :math:`\varepsilon^f_{n}`.  If this inner solver is called to evolve IVPs over
   time intervals :math:`[t_{0,i}, t_{F,i}]` with a relative tolerance
   :math:`\text{reltol}_n^f`, then it will result in accumulated errors over these
   intervals of the form

   .. math::
      \varepsilon^f_{n} = c(t_n) \left(\text{reltol}_n^f\right) \left(t_{F,i}-t_{0,i}\right),

   where :math:`c(t)` is independent of the tolerance or subinterval width but may vary in
   time, or equivalently,

   .. math::
      \varepsilon^f_{n} = \kappa(t_n) \left(\text{tolfac}_n^f\right),
      :label: inner_solver_assumption

   where :math:`\text{reltol}_n^f = \text{reltol}^s \text{tolfac}_n^f`,
   :math:`\text{reltol}^s` is the relative tolerance that was supplied to the
   current time scale solver, and where
   :math:`\kappa(t_n) = c(t_n) \text{reltol}^s \left(t_{F,i}-t_{0,i}\right)` is
   independent of the relative tolerance factor, :math:`\text{tolfac}_n^f`.

Single-rate controllers are constructed to adapt a single parameter, e.g.,
:math:`\delta`, under an assumption that solution error :math:`\varepsilon` depends
asymptotically on this parameter via the form

.. math::
   \varepsilon = \mathcal{O}(\delta^{q+1}).

Both :eq:`slow_error_assumption` and :eq:`inner_solver_assumption` fit this form,
with control parameters :math:`h^s` and :math:`\text{tolfac}^f_n`, and "orders"
:math:`P` and :math:`0`, respectively.  Thus an MRIHTol controller employs
*HControl* to adapt :math:`h_n^s` to control the current time scale error
:math:`\varepsilon^s_n`, and it employs *TolControl* to adapt
:math:`\text{tolfac}_n^f` to control the accumulated inner solver error
:math:`\varepsilon^f_n`.

To avoid overly large changes in calls to the inner solver, we apply bounds on the
results from *TolControl*.  If *TolControl* predicts a control parameter
:math:`\text{tolfac}'`, we obtain the eventual tolerance factor via
enforcing the following bounds:

.. math::
   \frac{\text{tolfac}_{n}^f}{\text{tolfac}'} &\le \text{max}_{relch},\\
   \frac{\text{tolfac}'}{\text{tolfac}_{n}^f} &\le \text{max}_{relch},\\
   \text{tolfac}_{min} &\le \text{tolfac}' \le \text{tolfac}_{max}.

The default values for these bounds are :math:`\text{max}_{relch} = 20`,
:math:`\text{tolfac}_{min} = 10^{-5}`, and :math:`\text{tolfac}_{max} = 1`.


Implementation
--------------

The SUNAdaptController_MRIHTol controller is implemented as a derived
SUNAdaptController class, and defines its *content* field as:

.. code-block:: c

   struct SUNAdaptControllerContent_MRIHTol
   {
     SUNAdaptController HControl;
     SUNAdaptController TolControl;
     sunrealtype inner_max_relch;
     sunrealtype inner_min_tolfac;
     sunrealtype inner_max_tolfac;
   };

These entries of the *content* field contain the following information:

* ``HControl`` - single time-scale SUNAdaptController object to adapt
  the current step size, :math:`h^s_n`.

* ``TolControl`` - single time-scale SUNAdaptController object to adapt
  the inner solver relative tolerance factor, :math:`\text{reltol}^f_n`.

* ``inner_max_relch`` - the parameter :math:`\text{max}_{relch}` above.

* ``inner_min_tolfac`` - the parameter :math:`\text{tolfac}_{min}` above.

* ``inner_max_tolfac`` - the parameter :math:`\text{tolfac}_{max}` above.

The header file to be included when using this module is
``sunadaptcontroller/sunadaptcontroller_mrihtol.h``.

The SUNAdaptController_MRIHTol class provides implementations of all operations
relevant to a ``SUN_ADAPTCONTROLLER_MRI_TOL`` controller listed in
:numref:`SUNAdaptController.Description.operations`. This class
also provides the following additional user-callable routines:


.. c:function:: SUNAdaptController SUNAdaptController_MRIHTol(SUNContext sunctx, SUNAdaptController HControl, SUNAdaptController TolControl)

   This constructor creates and allocates memory for a SUNAdaptController_MRIHTol
   object, and inserts its default parameters.

   :param sunctx: the current :c:type:`SUNContext` object.
   :param HControl: the slow time step adaptivity controller object.
   :param TolControl: the inner solver tolerance factor adaptivity controller object.
   :return: if successful, a usable :c:type:`SUNAdaptController` object;
            otherwise it will return ``NULL``.


.. c:function:: SUNErrCode SUNAdaptController_SetParams_MRIHTol(SUNAdaptController C, sunrealtype inner_max_relch, sunrealtype inner_min_tolfac, sunrealtype inner_max_tolfac)

   This user-callable function provides control over the relevant parameters
   above.  This should be called *before* the time integrator is called to evolve
   the problem.

   :param C: the SUNAdaptController_MRIHTol object.
   :param inner_max_relch: the parameter :math:`\text{max}_{relch}`.
   :param inner_min_tolfac: the parameter :math:`\text{tolfac}_{min}`.
   :param inner_max_tolfac: the parameter :math:`\text{tolfac}_{max}`.
   :return: :c:type:`SUNErrCode` indicating success or failure.
