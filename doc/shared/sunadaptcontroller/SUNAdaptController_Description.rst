..
   ----------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2002-2024, Lawrence Livermore National Security
   and Southern Methodist University.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   ----------------------------------------------------------------

.. _SUNAdaptController.Description:

The SUNAdaptController API
==========================

.. versionadded:: 6.7.0

.. versionchanged:: x.y.z

The SUNAdaptController base class provides a common API for accuracy-based adaptivity
controllers to be used by SUNDIALS integrators. These controllers estimate step
sizes (among other things) such that the next step solution satisfies a desired
temporal accuracy, while striving to maximize computational efficiency. We note
that in the descriptions below, we frequently use *dsm* to represent
temporal error. This is **not** the raw temporal error estimate; instead, it is
a norm of the temporal error estimate after scaling by the user-supplied
accuracy tolerances (see :eq:`ARKODE_WRMS_NORM`),

.. math::
   \text{dsm} = \left( \frac{1}{N} \sum_{i=1}^N
   \left(\frac{\text{error}_i}{\text{rtol}\cdot |y_{n-1,i}| + \text{atol}_i}\right)^2\right)^{1/2}.

Thus *dsm* values below one represent errors estimated to be more accurate than
needed, whereas errors above one are considered to be larger than allowable.

The :c:type:`SUNAdaptController` class is modeled after SUNDIALS' other object-oriented
classes, in that this class contains a pointer to an implementation-specific
*content*, an *ops* structure with generic controller operations, and a
:c:type:`SUNContext` object.

A :c:type:`SUNAdaptController` is a pointer to the
:c:struct:`_generic_SUNAdaptController` structure:

.. c:type:: struct _generic_SUNAdaptController *SUNAdaptController

.. c:struct:: _generic_SUNAdaptController

   .. c:member:: void* content

      Pointer to the controller-specific member data

   .. c:member:: SUNAdaptController_Ops ops;

      A virtual table of controller operations provided by a specific
      implementation

   .. c:member:: SUNContext sunctx

      The SUNDIALS simulation context

The virtual table structure is defined as

.. c:type:: struct _generic_SUNAdaptController_Ops *SUNAdaptController_Ops

.. c:struct:: _generic_SUNAdaptController_Ops

   The structure defining :c:type:`SUNAdaptController` operations.

   .. c:member:: SUNAdaptController_Type (*gettype)(SUNAdaptController C)

      The function implementing :c:func:`SUNAdaptController_GetType`

   .. c:member:: SUNErrCode (*destroy)(SUNAdaptController C)

      The function implementing :c:func:`SUNAdaptController_Destroy`

   .. c:member:: SUNErrCode (*estimatestep)(SUNAdaptController C, sunrealtype h, int p, sunrealtype dsm, sunrealtype* hnew)

      The function implementing :c:func:`SUNAdaptController_EstimateStep`

   .. c:member:: SUNErrCode (*estimatemristeps)(SUNAdaptController C, sunrealtype H, sunrealtype h, int P, sunrealtype DSM, sunrealtype dsm, sunrealtype* Hnew, sunrealtype* hnew)

      The function implementing :c:func:`SUNAdaptController_EstimateMRISteps`

   .. c:member:: SUNErrCode (*estimatesteptol)(SUNAdaptController C, sunrealtype H, sunrealtype tolfac, int P, sunrealtype DSM, sunrealtype dsm, sunrealtype* Hnew, sunrealtype* tolfacnew)

      The function implementing :c:func:`SUNAdaptController_EstimateStepTol`

   .. c:member:: SUNErrCode (*reset)(SUNAdaptController C)

      The function implementing :c:func:`SUNAdaptController_Reset`

   .. c:member:: SUNErrCode (*setdefaults)(SUNAdaptController C)

      The function implementing :c:func:`SUNAdaptController_SetDefaults`

   .. c:member:: SUNErrCode (*write)(SUNAdaptController C, FILE* fptr)

      The function implementing :c:func:`SUNAdaptController_Write`

   .. c:member:: SUNErrCode (*seterrorbias)(SUNAdaptController C, sunrealtype bias)

      The function implementing :c:func:`SUNAdaptController_SetErrorBias`

   .. c:member:: SUNErrCode (*updateh)(SUNAdaptController C, sunrealtype h, sunrealtype dsm)

      The function implementing :c:func:`SUNAdaptController_UpdateH`

   .. c:member:: SUNErrCode (*updatemrih)(SUNAdaptController C, sunrealtype H, sunrealtype h, sunrealtype DSM, sunrealtype dsm)

      The function implementing :c:func:`SUNAdaptController_UpdateMRIH`

   .. c:member:: SUNErrCode (*updatemritol)(SUNAdaptController C, sunrealtype H, sunrealtype tolfac, sunrealtype DSM, sunrealtype dsm)

      The function implementing :c:func:`SUNAdaptController_UpdateMRITol`

   .. c:member:: SUNErrCode (*space)(SUNAdaptController C, long int *lenrw, long int *leniw)

      The function implementing :c:func:`SUNAdaptController_Space`


.. _SUNAdaptController.Description.controllerTypes:

SUNAdaptController Types
------------------------

The time integrators in SUNDIALS adapt a variety of parameters to achieve
accurate and efficient computations. To this end, each SUNAdaptController implementation
should note its type, so that integrators will understand the types of
adaptivity that the controller is designed to perform. These are encoded in the
following set of SUNAdaptController types:

.. c:enum:: SUNAdaptController_Type

   The enumerated type :c:type:`SUNAdaptController_Type` defines the enumeration
   constants for SUNDIALS error controller types

.. c:enumerator:: SUN_ADAPTCONTROLLER_NONE

   Empty object that performs no control.

.. c:enumerator:: SUN_ADAPTCONTROLLER_H

   Controls a single-rate step size.

.. c:enumerator:: SUN_ADAPTCONTROLLER_MRI_H

   Controls both slow and fast time steps within a multirate simulation that has
   two time scales.

.. c:enumerator:: SUN_ADAPTCONTROLLER_MRI_TOL

   Controls both a slow time step and a tolerance factor to apply on the next-faster
   time scale within a multirate simulation that has an arbitrary number of time scales.



.. _SUNAdaptController.Description.operations:

SUNAdaptController Operations
-----------------------------

The base SUNAdaptController class defines and implements all SUNAdaptController functions. Most
of these routines are merely wrappers for the operations defined by a particular
SUNAdaptController implementation, which are accessed through the *ops* field of the
``SUNAdaptController`` structure. The base SUNAdaptController class provides the
constructor

.. c:function:: SUNAdaptController SUNAdaptController_NewEmpty(SUNContext sunctx)

  This function allocates a new generic ``SUNAdaptController`` object and initializes
  its content pointer and the function pointers in the operations structure to
  ``NULL``.

  :param sunctx: the :c:type:`SUNContext` object (see :numref:`SUNDIALS.SUNContext`)

  :returns: If successful, a generic :c:type:`SUNAdaptController` object. If
            unsuccessful, a ``NULL`` pointer will be returned.

Each of the following methods are *optional* for any specific SUNAdaptController
implementation, however some may be required based on the implementation's
:c:type:`SUNAdaptController_Type` (see Section :numref:`SUNAdaptController.Description.controllerTypes`). We
note these requirements below. Additionally, we note the behavior of the base SUNAdaptController methods when they perform an action other than only a successful return.

.. c:function:: void SUNAdaptController_DestroyEmpty(SUNAdaptController C)

  This routine frees the generic ``SUNAdaptController`` object, under the
  assumption that any implementation-specific data that was allocated within the
  underlying content structure has already been freed. It will additionally test
  whether the ops pointer is ``NULL``, and, if it is not, it will free it as
  well.

  :param C: the :c:type:`SUNAdaptController` object.
  :return: :c:type:`SUNErrCode` indicating success or failure.

  Usage:

  .. code-block:: c

    retval = SUNAdaptController_DestroyEmpty(C);

.. c:function:: SUNAdaptController_Type SUNAdaptController_GetType(SUNAdaptController C)

   Returns the type identifier for the controller *C*. Returned values
   are given in Section :numref:`SUNAdaptController.Description.controllerTypes`

   :param C: the :c:type:`SUNAdaptController` object.
   :return: :c:type:`SUNAdaptController_Type` type identifier.

   Usage:

   .. code-block:: c

      SUNAdaptController_Type id = SUNAdaptController_GetType(C);

.. c:function:: SUNErrCode SUNAdaptController_Destroy(SUNAdaptController C)

   Deallocates the controller *C*. If this method is not provided by the
   implementation, the base class method will free both the *content* and
   *ops* objects -- this should be sufficient unless a controller implementation
   performs dynamic memory allocation of its own (note that the
   SUNDIALS-provided SUNAdaptController implementations do not need to supply this
   routine).

   :param C: the :c:type:`SUNAdaptController` object.
   :return: :c:type:`SUNErrCode` indicating success or failure.

   Usage:

   .. code-block:: c

      retval = SUNAdaptController_Destroy(C);

.. c:function:: SUNErrCode SUNAdaptController_EstimateStep(SUNAdaptController C, sunrealtype h, int p, sunrealtype dsm, sunrealtype* hnew)

   Estimates a single-rate step size. This routine is required for controllers
   of type ``SUN_ADAPTCONTROLLER_H``.  If this is not provided by the
   implementation, the base class method will set ``*hnew = h`` and return.

   :param C: the :c:type:`SUNAdaptController` object.
   :param h: the step size from the previous step attempt.
   :param p: the current order of accuracy for the time integration method.
   :param dsm: the local temporal estimate from the previous step attempt.
   :param hnew: (output) the estimated step size.
   :return: :c:type:`SUNErrCode` indicating success or failure.

   Usage:

   .. code-block:: c

      retval = SUNAdaptController_EstimateStep(C, hcur, p, dsm, &hnew);

.. c:function:: SUNErrCode SUNAdaptController_EstimateMRISteps(SUNAdaptController C, sunrealtype H, sunrealtype h, int P, sunrealtype DSM, sunrealtype dsm, sunrealtype* Hnew, sunrealtype* hnew)

   Estimates slow and fast step sizes within a two-time-scale multirate application.
   This routine is required for controllers of type ``SUN_ADAPTCONTROLLER_MRI_H``.
   If this is not provided by the implementation, the base class method will set
   ``*Hnew = H`` and ``*hnew = h`` and return.

   :param C: the :c:type:`SUNAdaptController` object.
   :param H: the slow step size from the previous step attempt.
   :param h: the fast size from the previous step attempt.
   :param P: the current order of accuracy for the slow time scale integration method.
   :param DSM: the slow time scale local temporal estimate from the previous step attempt.
   :param dsm: the fast time scale local temporal estimate from the previous step attempt.
   :param Hnew: (output) the estimated slow step size.
   :param hnew: (output) the estimated fast step size.
   :return: :c:type:`SUNErrCode` indicating success or failure.

   Usage:

   .. code-block:: c

      retval = SUNAdaptController_EstimateMRISteps(C, Hcur, hcur, P, DSM, dsm, &Hnew, &hnew);

.. c:function:: SUNErrCode SUNAdaptController_EstimateStepTol(SUNAdaptController C, sunrealtype H, sunrealtype tolfac, int P, sunrealtype DSM, sunrealtype dsm, sunrealtype* Hnew, sunrealtype* tolfacnew)

   Estimates a slow step size and a fast tolerance multiplication factor
   for two adjacent time scales within a multirate application.  This
   routine is required for controllers of type ``SUN_ADAPTCONTROLLER_MRI_TOL``.
   If this is not provided by the implementation, the base class method will set
   ``*Hnew = H`` and ``*tolfacnew = tolfac`` and return.

   :param C: the :c:type:`SUNAdaptController` object.
   :param H: the slow step size from the previous step attempt.
   :param tolfac: the current relative tolerance factor for the next-faster time scale.
   :param P: the current order of accuracy for the slow time scale integration method.
   :param DSM: the slow time scale local temporal estimate from the previous step attempt.
   :param dsm: the fast time scale local temporal estimate from the previous step attempt.
   :param Hnew: (output) the estimated slow step size.
   :param tolfacnew: (output) the estimated relative tolerance factor.
   :return: :c:type:`SUNErrCode` indicating success or failure.

   .. note::

      If the current [slow] time scale has relative tolerance :math:`rtol`, then the
      next-faster time scale will be called with relative tolerance :math:`tolfac * rtol`.

   Usage:

   .. code-block:: c

      retval = SUNAdaptController_EstimateStepTol(C, Hcur, tolfac, P, DSM, dsm, &Hnew, &tolfacnew);

.. c:function:: SUNErrCode SUNAdaptController_Reset(SUNAdaptController C)

   Resets the controller to its initial state, e.g., if it stores a small number
   of previous *dsm* or *h* values.

   :param C:  the :c:type:`SUNAdaptController` object.
   :return: :c:type:`SUNErrCode` indicating success or failure.

   Usage:

   .. code-block:: c

      retval = SUNAdaptController_Reset(C);

.. c:function:: SUNErrCode SUNAdaptController_SetDefaults(SUNAdaptController C)

   Sets the controller parameters to their default values.

   :param C:  the :c:type:`SUNAdaptController` object.
   :return: :c:type:`SUNErrCode` indicating success or failure.

   Usage:

   .. code-block:: c

      retval = SUNAdaptController_SetDefaults(C);

.. c:function:: SUNErrCode SUNAdaptController_Write(SUNAdaptController C, FILE* fptr)

   Writes all controller parameters to the indicated file pointer.

   :param C:  the :c:type:`SUNAdaptController` object.
   :param fptr:  the output stream to write the parameters to.
   :return: :c:type:`SUNErrCode` indicating success or failure.

   Usage:

   .. code-block:: c

      retval = SUNAdaptController_Write(C, stdout);

.. c:function:: SUNErrCode SUNAdaptController_SetErrorBias(SUNAdaptController C, sunrealtype bias)

   Sets an error bias factor for scaling the local error factors. This is
   typically used to slightly exaggerate the temporal error during the
   estimation process, leading to a more conservative estimated step size.

   :param C:  the :c:type:`SUNAdaptController` object.
   :param bias:  the error bias factor -- an input :math:`\leq 0` indicates to use
                 the default value for the controller.
   :return: :c:type:`SUNErrCode` indicating success or failure.

   Usage:

   .. code-block:: c

      retval = SUNAdaptController_SetErrorBias(C, 1.2);

.. c:function:: SUNErrCode SUNAdaptController_UpdateH(SUNAdaptController C, sunrealtype h, sunrealtype dsm)

   Notifies a controller of type ``SUN_ADAPTCONTROLLER_H`` that a successful time step
   was taken with stepsize *h* and local error factor *dsm*, indicating that these
   can be saved for subsequent controller functions. This is typically relevant for
   controllers that store a history of either step sizes or error estimates for
   performing the estimation process.

   :param C:  the :c:type:`SUNAdaptController` object.
   :param h:  the successful step size.
   :param dsm:  the successful temporal error estimate.
   :return: :c:type:`SUNErrCode` indicating success or failure.

   Usage:

   .. code-block:: c

      retval = SUNAdaptController_UpdateH(C, h, dsm);

.. c:function:: SUNErrCode SUNAdaptController_UpdateMRIH(SUNAdaptController C, sunrealtype H, sunrealtype h, sunrealtype DSM, sunrealtype dsm)

   Notifies a controller of type ``SUN_ADAPTCONTROLLER_MRI_H`` that a successful time step
   was taken with slow stepsize *H* and fast stepsize *h*, and that the step had slow and
   fast local error factors *DSM* and *dsm*, indicating that these can be saved for
   subsequent controller functions. This is typically relevant for controllers that store a
   history of either step sizes or error estimates for performing the estimation process.

   :param C:  the :c:type:`SUNAdaptController` object.
   :param H:  the successful slow step size.
   :param h:  the successful fast step size.
   :param DSM:  the successful slow temporal error estimate.
   :param dsm:  the successful fast temporal error estimate.
   :return: :c:type:`SUNErrCode` indicating success or failure.

   Usage:

   .. code-block:: c

      retval = SUNAdaptController_UpdateMRIH(C, H, h, DSM, dsm);

.. c:function:: SUNErrCode SUNAdaptController_UpdateMRITol(SUNAdaptController C, sunrealtype H, sunrealtype tolfac, sunrealtype DSM, sunrealtype dsm)

   Notifies a controller of type ``SUN_ADAPTCONTROLLER_MRI_TOL`` that a successful time step
   was taken with slow stepsize *H* and fast relative tolerance factor *tolfac*, and that the
   step had slow and fast local error factors *DSM* and *dsm*, indicating that these can be
   saved for subsequent controller functions. This is typically relevant for controllers that
   store a history of either step sizes or error estimates for performing the estimation process.

   :param C:  the :c:type:`SUNAdaptController` object.
   :param H:  the successful slow step size.
   :param tolfac:  the successful fast time scale relative tolerance factor.
   :param DSM:  the successful slow temporal error estimate.
   :param dsm:  the successful fast temporal error estimate.
   :return: :c:type:`SUNErrCode` indicating success or failure.

   Usage:

   .. code-block:: c

      retval = SUNAdaptController_UpdateMRITol(C, H, tolfac, DSM, dsm);

.. c:function:: SUNErrCode SUNAdaptController_Space(SUNAdaptController C, long int *lenrw, long int *leniw)

   Informative routine that returns the memory requirements of the
   :c:type:`SUNAdaptController` object.

   :param C:  the :c:type:`SUNAdaptController` object..
   :param lenrw: (output)  number of ``sunsunrealtype`` words stored in the
                 controller.
   :param leniw: (output)  number of ``sunindextype`` words stored in the
                 controller. This may also include pointers, `int` and
                 `long int` words.
   :return: :c:type:`SUNErrCode` indicating success or failure.

   Usage:

   .. code-block:: c

      retval = SUNAdaptController_Space(C, &lenrw, &leniw);



C/C++ API Usage
---------------

Specific SUNDIALS adaptivity controller modules can be used in C and C++ programs by including
the corresponding header file for that module, e.g. ``sunadaptcontroller/sunadaptcontroller_XYZ.h``.

Example usage (here ``SUNAdaptController_XYZ`` is a placeholder for an actual SUNAdaptController
constructor):

.. code-block:: c

    #include <stdio.h>
    #include <stdlib.h>
    #include <sundials/sundials_context.h>
    #include <sundials/sundials_types.h>
    #include <sunadaptcontroller/sunadaptcontroller_XYZ.h>

    int main()
    {
        /* Create a SUNContext object */
        SUNContext sunctx = ...;

        /* Create a SUNAdaptController object */
        SUNAdaptController C = SUNAdaptController_XYZ(sunctx);

        /* Use the control object */

        /* Destroy the control object */
        retval = SUNAdaptController_Destroy(C);

        return 0;
    }
