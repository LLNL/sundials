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

.. _SUNControl.Description:

The SUNControl API
==================

The SUNControl base class provides a common API for accuracy-based adaptivity
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

The base ``SUNControl`` class is modeled after SUNDIALS' other object-oriented
classes, in that this class contains a pointer to an implementation-specific
*content*, an *ops* structure with generic controller operations, and a
:c:type:`SUNContext` object. Specifically, the type ``SUNControl`` is defined
as:

.. c:type:: struct generic_SUNControl_ *SUNControl

and the base class structure is defined as

.. code-block:: C

   struct generic_SUNControl_ {
        void* content;
        generic_SUNControl_Ops_* ops;
        SUNContext sunctx;
    };

Here, ``generic_SUNControl_Ops_`` is the pointer to a structure containing
function pointers to the various controller operations, and is defined as

.. code-block:: c

    struct generic_SUNControl_Ops_ {
        SUNControl_Type (*getid)(SUNControl C);
        void            (*destroy)(SUNControl C);
        int             (*estimatestep)(SUNControl C, realtype h, realtype dsm, realtype* hnew);
        int             (*estimatestepandorder)(SUNControl C, realtype h, int q, realtype dsm, realtype* hnew, int *qnew);
        int             (*estimatemristeps)(SUNControl C, realtype H, realtype h, realtype DSM, realtype dsm, realtype* Hnew, realtype *hnew);
        int             (*estimatesteptol)(SUNControl C, realtype H, realtype tolfac, realtype DSM, realtype dsm, realtype *Hnew, realtype* tolfacnew);
        int             (*reset)(SUNControl C);
        int             (*setdefaults)(SUNControl C);
        int             (*write)(SUNControl C, FILE* fptr);
        int             (*setmethodorder)(SUNControl C, int q);
        int             (*setembeddingorder)(SUNControl C, int p);
        int             (*seterrorbias)(SUNControl C, realtype bias);
        int             (*update)(SUNControl C, realtype h, realtype dsm);
        int             (*updatemrih)(SUNControl C, realtype H, realtype h, realtype DSM, realtype dsm);
        int             (*updatemritol)(SUNControl C, realtype H, realtype tolfac, realtype DSM, realtype dsm);
        int             (*space)(SUNControl C, long int *lenrw, long int *leniw);
    };


.. _SUNControl.Description.controllerTypes:

SUNControl Types
----------------

The time integrators in SUNDIALS adapt a variety of parameters to achieve
accurate and efficient computations. To this end, each SUNControl implementation
should note its type, so that integrators will understand the types of
adaptivity that the controller is designed to perform. These are encoded in the
following set of SUNControl types:

.. c:enum:: SUNControl_Type

   The enumerated type :c:type:`SUNControl_Type` defines the enumeration
   constants for SUNDIALS error controller types

.. c:enumerator:: SUNDIALS_CONTROL_NONE

   Empty object that performs no control.

.. c:enumerator:: SUNDIALS_CONTROL_H

   Controls a single-rate step size.

.. c:enumerator:: SUNDIALS_CONTROL_HQ

   Controls a single-rate step size and method order.

.. c:enumerator:: SUNDIALS_CONTROL_MRI_H

   Controls two multirate step sizes.

.. c:enumerator:: SUNDIALS_CONTROL_MRI_TOL

   Controls slow multirate step size and fast tolerance.



.. _SUNControl.Description.operations:

SUNControl Operations
---------------------

The base SUNControl class defines and implements all SUNControl functions. Most
of these routines are merely wrappers for the operations defined by a particular
SUNControl implementation, which are accessed through the *ops* field of the
``SUNControl`` structure. However, the base SUNControl class also provides the
convenience routine

.. c:function:: SUNControl SUNControl_NewEmpty(SUNContext sunctx)

  This function allocates a new generic ``SUNControl`` object and initializes
  its content pointer and the function pointers in the operations structure to
  ``NULL``.

  :param sunctx: the :c:type:`SUNContext` object (see :numref:`SUNDIALS.SUNContext`)

  :returns: If successful, a generic :c:type:`SUNControl` object. If
            unsuccessful, a ``NULL`` pointer will be returned.

Each of the following routines are *optional* for any specific SUNControl
implementation, however some may be required based on the implementation's
SUNControl_Type (see Section :numref:`SUNControl.Description.controllerTypes`). We
note these requirements, as well as the behavior of the base SUNControl wrapper
routine, below.

.. c:function:: SUNControl_Type SUNControl_GetType(SUNControl C)

   Returns the type identifier for the controller *C*. Returned values
   are given in Section :numref:`SUNControl.Description.controllerTypes`

   :param C: the :c:type:`SUNControl` object.
   :return: :c:type:`SUNControl_Type` type identifier.

   Usage:

   .. code-block:: c

      SUNControl_Type id = SUNControl_GetType(C);

.. c:function:: void SUNControl_Destroy(SUNControl C)

   Deallocates the controller *C*. If this is not provided by the
   implementation, the base wrapper routine will free both the *content* and
   *ops* objects -- this should be sufficient unless a controller implementation
   performs dynamic memory allocation of its own (note that the
   SUNDIALS-provided SUNControll implementations do not need to supply this
   routine).

   :param C: the :c:type:`SUNControl` object.

   Usage:

   .. code-block:: c

      SUNControl_Destroy(C);

.. c:function:: int SUNControl_EstimateStep(SUNControl C, realtype h, realtype dsm, realtype* hnew)

   Estimates a single-rate step size. This routine is required for controllers
   of type ``SUNDIALS_CONTROL_H``.

   :param C: the :c:type:`SUNControl` object.
   :param h: the step size from the previous step attempt.
   :param dsm: the local temporal estimate from the previous step attempt.
   :param hnew: (output) pointer to the estimated step size.
   :return: error code indicating success failure
            (see :numref:`SUNControl.Description.errorCodes`).

   Usage:

   .. code-block:: c

      retval = SUNControl_EstimateStep(C, hcur, dsm, &hnew);

.. c:function:: int SUNControl_EstimateStepAndOrder(SUNControl C, realtype h, int q, realtype dsm, realtype* hnew, int* qnew)

   Estimates a single-rate step size and corresponding method order. This
   routine is required for controllers of type ``SUNDIALS_CONTROL_HQ``.

   :param C: the :c:type:`SUNControl` object.
   :param h: the step size from the previous step attempt.
   :param q: the method order from the previous step attempt.
   :param dsm: the local temporal estimate from the previous step attempt.
   :param hnew: (output)  pointer to the estimated step size.
   :param qnew: (output)  pointer to the estimated method order.
   :return: error code indicating success failure
            (see :numref:`SUNControl.Description.errorCodes`).

   Usage:

   .. code-block:: c

      retval = SUNControl_EstimateStepAndOrder(C, hcur, qcur, dsm, &hnew, &qnew);

.. c:function:: int SUNControl_EstimateMRISteps(SUNControl C, realtype H, realtype h, realtype DSM, realtype dsm, realtype* Hnew, realtype *hnew)

   Estimates the slow and fast multirate step sizes. This routine is required
   for controllers of type ``SUNDIALS_CONTROL_MRI_H``.

   :param C: the :c:type:`SUNControl` object.
   :param H: the slow step size from the previous multirate step attempt.
   :param h: the fast step size from the previous multirate step attempt.
   :param DSM: the local slow temporal error estimate from the previous step
               attempt.
   :param dsm: the local fast temporal error estimate from the previous step
               attempt.
   :param Hnew: (output) pointer to the estimated slow step size.
   :param hnew: (output) pointer to the estimated fast step size.
   :return: error code indicating success failure
            (see :numref:`SUNControl.Description.errorCodes`).

   Usage:

   .. code-block:: c

      retval = SUNControl_EstimateMRISteps(C, Hcur, hcur, DSM, &Hnew, &hnew);

.. c:function:: int SUNControl_EstimateStepTol(SUNControl C, realtype H, realtype tolfac, realtype DSM, realtype *Hnew, realtype* tolfacnew)

   Estimates the slow step size and recommended fast relative tolerance factor
   for a multirate step. This routine is required for controllers of type
   ``SUNDIALS_CONTROL_MRI_TOL``.

   :param C: the :c:type:`SUNControl` object.
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
            (see :numref:`SUNControl.Description.errorCodes`).

   Usage:

   .. code-block:: c

      retval = SUNControl_EstimateStepTol(C, Hcur, tolfaccur, DSM, &Hnew, &tolfacnew);

.. c:function:: int SUNControl_Reset(SUNControl C)

   Resets the controller to its initial state, e.g., if it stores a small number
   of previous *dsm* or *h* values. The return value is an integer flag denoting
   success/failure of the routine (see
   :numref:`SUNControl.Description.errorCodes`).

   :param C:  the :c:type:`SUNControl` object.
   :return: error code indicating success failure
            (see :numref:`SUNControl.Description.errorCodes`).

   Usage:

   .. code-block:: c

      retval = SUNControl_Reset(C);

.. c:function:: int SUNControl_SetDefaults(SUNControl C)

   Sets the controller parameters to their default values.

   :param C:  the :c:type:`SUNControl` object..
   :return: error code indicating success failure
            (see :numref:`SUNControl.Description.errorCodes`).

   Usage:

   .. code-block:: c

      retval = SUNControl_SetDefaults(C);

.. c:function:: int SUNControl_Write(SUNControl C, FILE* fptr)

   Writes all controller parameters to the indicated file pointer.

   :param C:  the :c:type:`SUNControl` object.
   :param fptr:  the output stream to write the parameters.
   :return: error code indicating success failure
            (see :numref:`SUNControl.Description.errorCodes`).

   Usage:

   .. code-block:: c

      retval = SUNControl_Write(C, stdout);

.. c:function:: int SUNControl_SetMethodOrder(SUNControl C, int q)

   Called by the time integrator to inform the controller of the asymptotic
   order of accuracy for the method.

   :param C:  the :c:type:`SUNControl` object.
   :param q:  the asymptotic order of accuracy for the time integration method.
   :return: error code indicating success failure
            (see :numref:`SUNControl.Description.errorCodes`).

   Usage:

   .. code-block:: c

      retval = SUNControl_SetMethodOrder(C, 3);

.. c:function:: int SUNControl_SetEmbeddingOrder(SUNControl C, int p)

   Called by the time integrator to inform the controller of the asymptotic
   order of accuracy for the method embedding.

   :param C:  the :c:type:`SUNControl` object.
   :param p:  the asymptotic order of accuracy for the time integration method
              embedding.
   :return: error code indicating success failure
            (see :numref:`SUNControl.Description.errorCodes`).

   Usage:

   .. code-block:: c

      retval = SUNControl_SetEmbeddingOrder(C, 2);

.. c:function:: int SUNControl_SetErrorBias(SUNControl C, realtype bias)

   Sets an error bias factor for scaling the local error factors. This is
   typically used to slightly exaggerate the temporal error during the
   estimation process, leading to a more conservative estimated step size.

   :param C:  the :c:type:`SUNControl` object.
   :param bias:  the error bias factor.
   :return: error code indicating success failure
            (see :numref:`SUNControl.Description.errorCodes`).

   Usage:

   .. code-block:: c

      retval = SUNControl_SetErrorBias(C, 1.2);

.. c:function:: int SUNControl_Update(SUNControl C, realtype h, realtype dsm)

   Notifies the controller of a successful time step of size *h* and with
   temporal error estimate *dsm*. This is typically used for controllers that
   store a history of either step sizes or error estimates for performing the
   estimation process.

   :param C:  the :c:type:`SUNControl` object.
   :param h:  the successful step size.
   :param dsm:  the successful temporal error estimate.
   :return: error code indicating success failure
            (see :numref:`SUNControl.Description.errorCodes`).

   Usage:

   .. code-block:: c

      retval = SUNControl_Update(C, h, dsm);

.. c:function:: int SUNControl_UpdateMRIH(SUNControl C, realtype H, realtype h, realtype DSM, realtype dsm)

   Notifies the controller of a successful multirate time step of sizes *H* and
   *h*, and with temporal error estimates *DSM* and *dsm*. This is used for
   controllers of type *SUNDIALS_CONTROL_MRI_H* that store a history of either
   step size inputs or resulting error estimates for performing the estimation
   process.

   :param C:  the :c:type:`SUNControl` object.
   :param H:  the successful slow step size.
   :param h:  the successful fast step size.
   :param DSM:  the successful slow temporal error estimate.
   :param dsm:  the successful fast temporal error estimate.
   :return: error code indicating success failure
            (see :numref:`SUNControl.Description.errorCodes`).

   Usage:

   .. code-block:: c

      retval = SUNControl_UpdateMRIH(C, H, h, DSM, dsm);

.. c:function:: int SUNControl_UpdateMRITol(SUNControl C, realtype H, realtype tolfac, realtype DSM, realtype dsm)

   Notifies the controller of a successful multirate time step of size *H* and
   fast tolerance factor *tolfac*, that resulted in temporal error estimates
   *DSM* and *dsm*. This is typically used for controllers of type
   *SUNDIALS_CONTROL_MRI_TOL* that store a history of either control inputs or
   resulting error estimates for performing the estimation process.

   :param C:  the :c:type:`SUNControl` object.
   :param H:  the successful slow step size.
   :param tolfac:  the successful fast relative tolerance factor.
   :param DSM:  the successful slow temporal error estimate.
   :param dsm:  the successful fast temporal error estimate.
   :return: error code indicating success failure
            (see :numref:`SUNControl.Description.errorCodes`).

   Usage:

   .. code-block:: c

      retval = SUNControl_Update(C, h, dsm);

.. c:function:: int SUNControl_Space(SUNControl C, long int *lenrw, long int *leniw)

   Informative routine that returns the memory requirements of the
   :c:type:`SUNControl` object.

   :param C:  the :c:type:`SUNControl` object..
   :param lenrw: (output)  number of ``sunrealtype`` words stored in the
                 controller.
   :param leniw: (output)  number of ``sunindextype`` words stored in the
                 controller. This may also include pointers, `int` and
                 `long int` words.
   :return: error code indicating success failure
            (see :numref:`SUNControl.Description.errorCodes`).

   Usage:

   .. code-block:: c

      retval = SUNControl_Space(C, &lenrw, &leniw);



.. _SUNControl.Description.errorCodes:

SUNControl Error Codes
----------------------

SUNControl functions return one of the following set of error codes:

* ``SUNCONTROL_SUCCESS`` (0) -- successful call.

* ``SUNCONTROL_ILL_INPUT`` (-1001) -- an illegal input has been provided to the function.

* ``SUNCONTROL_MEM_FAIL`` (-1002) -- a memory access or allocation failed.

* ``SUNCONTROL_USER_FCN_FAIL`` (-1003) -- a user-supplied function returned a nonzero [error] value.

* ``SUNCONTROL_OPERATION_FAIL`` (-1004) -- catch-all for errors not in the above list.


C/C++ API Usage
---------------

The SUNDIALS Controller module can be used in C and C++ programs by including
the header file ``sundials/sundials_controller.h``.

Example usage (here ``SUNControlXYZ`` is a placeholder for an actual SUNControl
implementation constructor):

.. code-block:: c

    #include <stdio.h>
    #include <stdlib.h>
    #include <sundials/sundials_context.h>
    #include <sundials/sundials_types.h>
    #include <sundials/sundials_controller.h>

    int main()
    {
        /* Create a SUNContext object */
        SUNContext sunctx = ...;

        /* Create a SUNControl object */
        SUNControl C = SUNControlXYZ(sunctx);

        /* Use the control object */

        /* Destroy the control object */
        SUNControl_Destroy(C);

        return 0;
    }
