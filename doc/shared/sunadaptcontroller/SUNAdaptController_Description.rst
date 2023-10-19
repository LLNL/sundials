..
   ----------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2002-2023, Lawrence Livermore National Security
   and Southern Methodist University.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   ----------------------------------------------------------------

.. _SUNAdaptController.Description:

The SUNAdaptController API
==========================

.. versionadded:: x.x.x

The SUNAdaptController base class provides a common API for accuracy-based adaptivity
controllers to be used by SUNDIALS integrators. These controllers estimate step
sizes (among other things) such that the next step solution satisfies a desired
temporal accuracy, while striving to maximize computational efficiency. We note
that in the descriptions below, we frequently use the object *dsm* to represent
temporal error. This is **not** the raw temporal error estimate; instead, it is
a norm of the temporal error estimate after scaling by the user-supplied
accuracy tolerances (the SUNDIALS WRMS-norm),

.. math::
   \text{dsm} = \left( \frac{1}{N} \sum_{i=1}^N
   \left(\frac{\text{error}_i}{\text{rtol}\cdot |y_{n-1,i}| + \text{atol}_i}\right)^2\right)^{1/2}.

Thus *dsm* values below one represent errors estimated to be more accurate than
needed, whereas errors above one are considered to be larger than allowable.

The base ``SUNAdaptController`` class is modeled after SUNDIALS' other object-oriented
classes, in that this class contains a pointer to an implementation-specific
*content*, an *ops* structure with generic controller operations, and a
:c:type:`SUNContext` object. Specifically, the type ``SUNAdaptController`` is defined
as:

.. c:type:: struct _generic_SUNAdaptController *SUNAdaptController

and the base class structure is defined as

.. code-block:: C

   struct _generic_SUNAdaptController {
        void* content;
        generic_SUNAdaptController_Ops_* ops;
        SUNContext sunctx;
    };

Here, ``_generic_SUNAdaptController_Ops`` is the pointer to a structure containing
function pointers to the various controller operations, and is defined as

.. code-block:: c

    struct _generic_SUNAdaptController_Ops {
        SUNAdaptController_Type (*getid)(SUNAdaptController C);
        int (*destroy)(SUNAdaptController C);
        int (*estimatestep)(SUNAdaptController C, sunrealtype h, sunrealtype dsm, sunrealtype* hnew);
        int (*estimatestepandorder)(SUNAdaptController C, sunrealtype h, int q, sunrealtype dsm, sunrealtype* hnew, int *qnew);
        int (*estimatemristeps)(SUNAdaptController C, sunrealtype H, sunrealtype h, sunrealtype DSM, sunrealtype dsm, sunrealtype* Hnew, sunrealtype *hnew);
        int (*estimatesteptol)(SUNAdaptController C, sunrealtype H, sunrealtype tolfac, sunrealtype DSM, sunrealtype dsm, sunrealtype *Hnew, sunrealtype* tolfacnew);
        int (*reset)(SUNAdaptController C);
        int (*setdefaults)(SUNAdaptController C);
        int (*write)(SUNAdaptController C, FILE* fptr);
        int (*setmethodorder)(SUNAdaptController C, int p);
        int (*adjustcontrollerorder)(SUNAdaptController C, int adj);
        int (*seterrorbias)(SUNAdaptController C, sunrealtype bias);
        int (*update)(SUNAdaptController C, sunrealtype h, sunrealtype dsm);
        int (*updatemrih)(SUNAdaptController C, sunrealtype H, sunrealtype h, sunrealtype DSM, sunrealtype dsm);
        int (*updatemritol)(SUNAdaptController C, sunrealtype H, sunrealtype tolfac, sunrealtype DSM, sunrealtype dsm);
        int (*space)(SUNAdaptController C, long int *lenrw, long int *leniw);
    };


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

.. c:enumerator:: SUN_ADAPTCONTROLLER_HQ

   Controls a single-rate step size and method order.

.. c:enumerator:: SUN_ADAPTCONTROLLER_MRI_H

   Controls two multirate step sizes.

.. c:enumerator:: SUN_ADAPTCONTROLLER_MRI_TOL

   Controls slow multirate step size and fast tolerance.



.. _SUNAdaptController.Description.operations:

SUNAdaptController Operations
-----------------------------

The base SUNAdaptController class defines and implements all SUNAdaptController functions. Most
of these routines are merely wrappers for the operations defined by a particular
SUNAdaptController implementation, which are accessed through the *ops* field of the
``SUNAdaptController`` structure. However, the base SUNAdaptController class also provides the
convenience routine

.. c:function:: SUNAdaptController SUNAdaptController_NewEmpty(SUNContext sunctx)

  This function allocates a new generic ``SUNAdaptController`` object and initializes
  its content pointer and the function pointers in the operations structure to
  ``NULL``.

  :param sunctx: the :c:type:`SUNContext` object (see :numref:`SUNDIALS.SUNContext`)

  :returns: If successful, a generic :c:type:`SUNAdaptController` object. If
            unsuccessful, a ``NULL`` pointer will be returned.

Each of the following routines are *optional* for any specific SUNAdaptController
implementation, however some may be required based on the implementation's
:c:type:`SUNAdaptController_Type` (see Section :numref:`SUNAdaptController.Description.controllerTypes`). We
note these requirements, as well as the behavior of the base SUNAdaptController wrapper
routine, below.

.. c:function:: SUNAdaptController_Type SUNAdaptController_GetType(SUNAdaptController C)

   Returns the type identifier for the controller *C*. Returned values
   are given in Section :numref:`SUNAdaptController.Description.controllerTypes`

   :param C: the :c:type:`SUNAdaptController` object.
   :return: :c:type:`SUNAdaptController_Type` type identifier.

   Usage:

   .. code-block:: c

      SUNAdaptController_Type id = SUNAdaptController_GetType(C);

.. c:function:: int SUNAdaptController_Destroy(SUNAdaptController C)

   Deallocates the controller *C*. If this is not provided by the
   implementation, the base wrapper routine will free both the *content* and
   *ops* objects -- this should be sufficient unless a controller implementation
   performs dynamic memory allocation of its own (note that the
   SUNDIALS-provided SUNAdaptController implementations do not need to supply this
   routine).

   :param C: the :c:type:`SUNAdaptController` object.
   :return: error code indicating success failure
            (see :numref:`SUNAdaptController.Description.errorCodes`).

   Usage:

   .. code-block:: c

      retval = SUNAdaptController_Destroy(C);

.. c:function:: int SUNAdaptController_EstimateStep(SUNAdaptController C, sunrealtype h, sunrealtype dsm, sunrealtype* hnew)

   Estimates a single-rate step size. This routine is required for controllers
   of type ``SUN_ADAPTCONTROLLER_H``.

   :param C: the :c:type:`SUNAdaptController` object.
   :param h: the step size from the previous step attempt.
   :param dsm: the local temporal estimate from the previous step attempt.
   :param hnew: (output) pointer to the estimated step size.
   :return: error code indicating success failure
            (see :numref:`SUNAdaptController.Description.errorCodes`).

   Usage:

   .. code-block:: c

      retval = SUNAdaptController_EstimateStep(C, hcur, dsm, &hnew);

.. c:function:: int SUNAdaptController_EstimateStepAndOrder(SUNAdaptController C, sunrealtype h, int q, sunrealtype dsm, sunrealtype* hnew, int* qnew)

   Estimates a single-rate step size and corresponding method order. This
   routine is required for controllers of type ``SUN_ADAPTCONTROLLER_HQ``.

   :param C: the :c:type:`SUNAdaptController` object.
   :param h: the step size from the previous step attempt.
   :param q: the method order from the previous step attempt.
   :param dsm: the local temporal estimate from the previous step attempt.
   :param hnew: (output)  pointer to the estimated step size.
   :param qnew: (output)  pointer to the estimated method order.
   :return: error code indicating success failure
            (see :numref:`SUNAdaptController.Description.errorCodes`).

   Usage:

   .. code-block:: c

      retval = SUNAdaptController_EstimateStepAndOrder(C, hcur, qcur, dsm, &hnew, &qnew);

.. c:function:: int SUNAdaptController_EstimateMRISteps(SUNAdaptController C, sunrealtype H, sunrealtype h, sunrealtype DSM, sunrealtype dsm, sunrealtype* Hnew, sunrealtype *hnew)

   Estimates the slow and fast multirate step sizes. This routine is required
   for controllers of type ``SUN_ADAPTCONTROLLER_MRI_H``.

   :param C: the :c:type:`SUNAdaptController` object.
   :param H: the slow step size from the previous multirate step attempt.
   :param h: the fast step size from the previous multirate step attempt.
   :param DSM: the local slow temporal error estimate from the previous step
               attempt.
   :param dsm: the local fast temporal error estimate from the previous step
               attempt.
   :param Hnew: (output) pointer to the estimated slow step size.
   :param hnew: (output) pointer to the estimated fast step size.
   :return: error code indicating success failure
            (see :numref:`SUNAdaptController.Description.errorCodes`).

   Usage:

   .. code-block:: c

      retval = SUNAdaptController_EstimateMRISteps(C, Hcur, hcur, DSM, &Hnew, &hnew);

.. c:function:: int SUNAdaptController_EstimateStepTol(SUNAdaptController C, sunrealtype H, sunrealtype tolfac, sunrealtype DSM, sunrealtype *Hnew, sunrealtype* tolfacnew)

   Estimates the slow step size and recommended fast relative tolerance factor
   for a multirate step. This routine is required for controllers of type
   ``SUN_ADAPTCONTROLLER_MRI_TOL``.

   :param C: the :c:type:`SUNAdaptController` object.
   :param H: the slow step size from the previous multirate step attempt.
   :param tolfac: the ratio of fast/slow relative tolerances,
                  :math:`\text{reltol}/\text{RELTOL}`, from the previous
                  multirate step attempt.
   :param DSM: the local slow temporal error estimate from the previous step
               attempt.
   :param dsm: the local fast temporal error estimate from the previous step
               attempt.
   :param Hnew: (output) pointer to the estimated slow step size.
   :param tolfacnew: (output) pointer to the estimated relative tolerance
                     ratio.
   :return: error code indicating success failure
            (see :numref:`SUNAdaptController.Description.errorCodes`).

   Usage:

   .. code-block:: c

      retval = SUNAdaptController_EstimateStepTol(C, Hcur, tolfaccur, DSM, &Hnew, &tolfacnew);

.. c:function:: int SUNAdaptController_Reset(SUNAdaptController C)

   Resets the controller to its initial state, e.g., if it stores a small number
   of previous *dsm* or *h* values. The return value is an integer flag denoting
   success/failure of the routine (see
   :numref:`SUNAdaptController.Description.errorCodes`).

   :param C:  the :c:type:`SUNAdaptController` object.
   :return: error code indicating success failure
            (see :numref:`SUNAdaptController.Description.errorCodes`).

   Usage:

   .. code-block:: c

      retval = SUNAdaptController_Reset(C);

.. c:function:: int SUNAdaptController_SetDefaults(SUNAdaptController C)

   Sets the controller parameters to their default values.

   :param C:  the :c:type:`SUNAdaptController` object..
   :return: error code indicating success failure
            (see :numref:`SUNAdaptController.Description.errorCodes`).

   Usage:

   .. code-block:: c

      retval = SUNAdaptController_SetDefaults(C);

.. c:function:: int SUNAdaptController_Write(SUNAdaptController C, FILE* fptr)

   Writes all controller parameters to the indicated file pointer.

   :param C:  the :c:type:`SUNAdaptController` object.
   :param fptr:  the output stream to write the parameters.
   :return: error code indicating success failure
            (see :numref:`SUNAdaptController.Description.errorCodes`).

   Usage:

   .. code-block:: c

      retval = SUNAdaptController_Write(C, stdout);

.. c:function:: int SUNAdaptController_SetMethodOrder(SUNAdaptController C, int p)

   Called by the time integrator to inform the controller of the asymptotic
   order of accuracy for the method.

   :param C:  the :c:type:`SUNAdaptController` object.
   :param p:  the asymptotic method order of accuracy to use.
   :return: error code indicating success failure
            (see :numref:`SUNAdaptController.Description.errorCodes`).

   Usage:

   .. code-block:: c

      retval = SUNAdaptController_SetMethodOrder(C, 3, 2);

.. c:function:: int SUNAdaptController_AdjustControllerOrder(SUNAdaptController C, int adj)

   Called by a user to request that the controller adjust the order specified by the time
   integration method when performing temporal adaptivity, e.g., if the user expects order
   reduction due to problem stiffness, they may request that the controller assume a
   reduced order of accuracy for the method by specifying a value :math:`adj < 0`.  This
   adjustment will apply to all subsequent time step adaptivity estimates used by the
   controller, and may be undone by a call with :math:`adj = 0`.

   :param C:  the :c:type:`SUNAdaptController` object.
   :param p:  the adjustment that will be applied to the values *p* and *q* from
              :c:func:`SUNAdaptController_SetMethodOrder`.
   :return:  error code indicating success failure
             (see :numref:`SUNAdaptController.Description.errorCodes`).

   Usage:

   .. code-block:: c

      retval = SUNAdaptController_AdjustControllerOrder(C, -1);

.. c:function:: int SUNAdaptController_SetErrorBias(SUNAdaptController C, sunrealtype bias)

   Sets an error bias factor for scaling the local error factors. This is
   typically used to slightly exaggerate the temporal error during the
   estimation process, leading to a more conservative estimated step size.

   :param C:  the :c:type:`SUNAdaptController` object.
   :param bias:  the error bias factor -- an input :math:`\leq 0` indicates to use
                 the default value for the controller.
   :return: error code indicating success failure
            (see :numref:`SUNAdaptController.Description.errorCodes`).

   Usage:

   .. code-block:: c

      retval = SUNAdaptController_SetErrorBias(C, 1.2);

.. c:function:: int SUNAdaptController_Update(SUNAdaptController C, sunrealtype h, sunrealtype dsm)

   Notifies the controller of a successful time step of size *h* and with
   temporal error estimate *dsm*. This is typically used for controllers that
   store a history of either step sizes or error estimates for performing the
   estimation process.

   :param C:  the :c:type:`SUNAdaptController` object.
   :param h:  the successful step size.
   :param dsm:  the successful temporal error estimate.
   :return: error code indicating success failure
            (see :numref:`SUNAdaptController.Description.errorCodes`).

   Usage:

   .. code-block:: c

      retval = SUNAdaptController_Update(C, h, dsm);

.. c:function:: int SUNAdaptController_UpdateMRIH(SUNAdaptController C, sunrealtype H, sunrealtype h, sunrealtype DSM, sunrealtype dsm)

   Notifies the controller of a successful multirate time step of sizes *H* and
   *h*, and with temporal error estimates *DSM* and *dsm*. This is used for
   controllers of type *SUN_ADAPTCONTROLLER_MRI_H* that store a history of either
   step size inputs or resulting error estimates for performing the estimation
   process.

   :param C:  the :c:type:`SUNAdaptController` object.
   :param H:  the successful slow step size.
   :param h:  the successful fast step size.
   :param DSM:  the successful slow temporal error estimate.
   :param dsm:  the successful fast temporal error estimate.
   :return: error code indicating success failure
            (see :numref:`SUNAdaptController.Description.errorCodes`).

   Usage:

   .. code-block:: c

      retval = SUNAdaptController_UpdateMRIH(C, H, h, DSM, dsm);

.. c:function:: int SUNAdaptController_UpdateMRITol(SUNAdaptController C, sunrealtype H, sunrealtype tolfac, sunrealtype DSM, sunrealtype dsm)

   Notifies the controller of a successful multirate time step of size *H* and
   fast tolerance factor *tolfac*, that resulted in temporal error estimates
   *DSM* and *dsm*. This is typically used for controllers of type
   *SUN_ADAPTCONTROLLER_MRI_TOL* that store a history of either control inputs or
   resulting error estimates for performing the estimation process.

   :param C:  the :c:type:`SUNAdaptController` object.
   :param H:  the successful slow step size.
   :param tolfac:  the successful fast relative tolerance factor.
   :param DSM:  the successful slow temporal error estimate.
   :param dsm:  the successful fast temporal error estimate.
   :return: error code indicating success failure
            (see :numref:`SUNAdaptController.Description.errorCodes`).

   Usage:

   .. code-block:: c

      retval = SUNAdaptController_Update(C, h, dsm);

.. c:function:: int SUNAdaptController_Space(SUNAdaptController C, long int *lenrw, long int *leniw)

   Informative routine that returns the memory requirements of the
   :c:type:`SUNAdaptController` object.

   :param C:  the :c:type:`SUNAdaptController` object..
   :param lenrw: (output)  number of ``sunsunrealtype`` words stored in the
                 controller.
   :param leniw: (output)  number of ``sunindextype`` words stored in the
                 controller. This may also include pointers, `int` and
                 `long int` words.
   :return: error code indicating success failure
            (see :numref:`SUNAdaptController.Description.errorCodes`).

   Usage:

   .. code-block:: c

      retval = SUNAdaptController_Space(C, &lenrw, &leniw);



.. _SUNAdaptController.Description.errorCodes:

SUNAdaptController Error Codes
------------------------------

SUNAdaptController functions return one of the following set of error codes:

* ``SUNADAPTCONTROLLER_SUCCESS`` (0) -- successful call.

* ``SUNADAPTCONTROLLER_ILL_INPUT`` (-1001) -- an illegal input has been provided to the function.

* ``SUNADAPTCONTROLLER_MEM_FAIL`` (-1002) -- a memory access or allocation failed.

* ``SUNADAPTCONTROLLER_USER_FCN_FAIL`` (-1003) -- a user-supplied function returned a nonzero [error] value.

* ``SUNADAPTCONTROLLER_OPERATION_FAIL`` (-1004) -- catch-all for errors not in the above list.

.. note::
   The SUNDIALS time integrators do not rely on these specific return values (only
   on whether the returned values are 0 (successful) or negative (failure).  Thus,
   user-defined implementations are not required to use these specific error codes,
   so long as the zero/negative structure is followed.


C/C++ API Usage
---------------

Specific SUNDIALS adaptivity controller modules can be used in C and C++ programs by including
the corresponding header file for that module, e.g. ``sunadaptcontroller/sunadaptcontrollerXYZ.h``.

Example usage (here ``SUNAdaptControllerXYZ`` is a placeholder for an actual SUNAdaptController
constructor):

.. code-block:: c

    #include <stdio.h>
    #include <stdlib.h>
    #include <sundials/sundials_context.h>
    #include <sundials/sundials_types.h>
    #include <sunadaptcontroller/sunadaptcontrollerXYZ.h>

    int main()
    {
        /* Create a SUNContext object */
        SUNContext sunctx = ...;

        /* Create a SUNAdaptController object */
        SUNAdaptController C = SUNAdaptControllerXYZ(sunctx);

        /* Use the control object */

        /* Destroy the control object */
        retval = SUNAdaptController_Destroy(C);

        return 0;
    }
