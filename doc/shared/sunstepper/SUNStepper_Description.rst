.. ----------------------------------------------------------------
   Programmer(s): Steven B. Roberts @LLNL
                  David J. Gardner @ LLNL
   ----------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2002-2024, Lawrence Livermore National Security
   and Southern Methodist University.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   ----------------------------------------------------------------

.. _SUNStepper.Description:

The SUNStepper API
==================

.. versionadded:: x.y.z

As with other SUNDIALS classes, the :c:type:`SUNStepper` abstract base class is
implemented using a C structure containing a ``content`` pointer to the derived
class member data and a structure of function pointers to the derived class
implementations of the virtual methods.

.. c:type:: SUNStepper

   An object for solving the IVP :eq:`SUNStepper_IVP`.

   The actual definition of the ``SUNStepper`` structure is kept private to
   allow for the object internals to change without impacting user code. The
   following sections describe the base class methods and the virtual methods
   that a must be provided by a derived class.

.. _SUNStepper.Description.BaseMethods:

Base Class Methods
------------------

This section describes methods provided by the :c:type:`SUNStepper` abstract
base class that aid the user in implementing derived classes. This includes
functions for creating and destroying a generic base class object, attaching and
retrieving the derived class ``content`` pointer, and setting function pointers
to derived class method implementations.

.. _SUNStepper.Description.BaseMethods.CreateDestroy:

Creating and Destroying an Object
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In addition to creating an empty :c:type:`SUNStepper` using
:c:func:`SUNStepper_Create` described below, there is the
:c:func:`ARKodeCreateSUNStepper` function to construct a :c:type:`SUNStepper`
from an ARKODE integrator.

.. c:function:: SUNErrCode SUNStepper_Create(SUNContext sunctx, SUNStepper *stepper)

   This function creates a :c:type:`SUNStepper` object to which a user should
   attach the member data (content) pointer and method function pointers.

   :param sunctx: the SUNDIALS simulation context.
   :param stepper: a pointer to a stepper object.
   :return: A :c:type:`SUNErrCode` indicating success or failure.

   **Example usage:**

   .. code-block:: C

      /* create an instance of the base class */
      SUNStepper stepper = NULL;
      SUNErrCode err = SUNStepper_Create(sunctx, &stepper);

   .. note::

      See :numref:`SUNStepper.Description.BaseMethods.Content` and
      :numref:`SUNStepper.Description.BaseMethods.AttachFunctions`
      for details on how to attach member data and method function pointers.


.. c:function:: SUNErrCode SUNStepper_Destroy(SUNStepper *stepper)

   This function destroys a :c:type:`SUNStepper` object.

   :param stepper: a pointer to a stepper object.
   :return: A :c:type:`SUNErrCode` indicating success or failure.

   .. note::

      This function only frees memory allocated within the base class and the
      base class structure itself. The user is responsible for freeing any
      memory allocated for the member data (content).


.. _SUNStepper.Description.BaseMethods.SteppingFunctions:

Stepping Functions
^^^^^^^^^^^^^^^^^^

.. c:function:: SUNErrCode SUNStepper_Evolve(SUNStepper stepper, sunrealtype tout, N_Vector vret, sunrealtype* tret)

   This function evolves the ODE :eq:`SUNStepper_IVP` towards the time ``tout``
   and stores the solution at time ``tret`` in ``vret``.

   :param stepper: the stepper object.
   :param tout: the time to evolve towards.
   :param vret: on output, the state at time ``tret``.
   :param tret: the time corresponding to the output value ``vret``.
   :return: A :c:type:`SUNErrCode` indicating success or failure.


.. c:function:: SUNErrCode SUNStepper_FullRhs(SUNStepper stepper, sunrealtype t, N_Vector v, N_Vector f)

   This function computes the full right-hand side function of the ODE,
   :math:`f(t, v) + r(t)` in :eq:`SUNStepper_IVP` for a given value of the
   independent variable ``t`` and state vector ``v``.

   :param stepper: the stepper object.
   :param t: the current value of the independent variable.
   :param v: the current value of the dependent variable vector.
   :param f: the output vector for the ODE right-hand side,
      :math:`f(t, v) + r(t)`, in :eq:`SUNStepper_IVP`.
   :return: A :c:type:`SUNErrCode` indicating success or failure.


.. c:function:: SUNErrCode SUNStepper_Reset(SUNStepper stepper, sunrealtype tR, N_Vector vR)

   This function resets the stepper state to the provided independent variable
   value and dependent variable vector.

   :param stepper: the stepper object.
   :param tR: the value of the independent variable :math:`t_R`.
   :param vR: the value of the dependent variable vector :math:`v(t_R)`.
   :return: A :c:type:`SUNErrCode` indicating success or failure.


.. c:function:: SUNErrCode SUNStepper_SetStopTime(SUNStepper stepper, sunrealtype tstop)

   This function specifies the value of the independent variable :math:`t` past
   which the solution is not to proceed.

   :param stepper: the stepper object.
   :param tstop: stopping time for the stepper.
   :return: A :c:type:`SUNErrCode` indicating success or failure.


.. c:function:: SUNErrCode SUNStepper_SetForcing(SUNStepper stepper, sunrealtype tshift, sunrealtype tscale, N_Vector* forcing, int nforcing)
   
   This function sets the data necessary to compute the forcing term
   :eq:`SUNStepper_forcing`. This includes the shift and scaling factors for the
   normalized time :math:`\frac{t - t_{\text{shift}}}{t_{\text{scale}}}` and the
   array of polynomial coefficient vectors :math:`\widehat{f}_k`.

   :param stepper: a stepper object.
   :param tshift: the time shift to apply to the current time when computing
      the forcing, :math:`t_{\text{shift}}`.
   :param tscale: the time scaling to apply to the current time when computing
      the forcing, :math:`t_{\text{scale}}`.
   :param forcing: a pointer to an array of forcing vectors,
      :math:`\widehat{f}_k`.
   :param nforcing: the number of forcing vectors, :math:`n_{\text{forcing}}`. A
      value of 0 effectively eliminates the forcing term.
   :return: A :c:type:`SUNErrCode` indicating success or failure.

   .. note::

      When integrating the ODE :eq:`SUNStepper_IVP` the :c:type:`SUNStepper` is
      responsible for evaluating ODE right-hand side function :math:`f(t, v)` as
      well as computing and applying the forcing term :eq:`SUNStepper_forcing`
      to obtain the full right-hand side of the ODE :eq:`SUNStepper_IVP`.

.. _SUNStepper.Description.BaseMethods.Content:

Attaching and Accessing the Content Pointer
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. c:function:: SUNErrCode SUNStepper_SetContent(SUNStepper stepper, void *content)

   This function attaches a member data (content) pointer to a
   :c:type:`SUNStepper` object.

   :param stepper: a stepper object.
   :param content: a pointer to the stepper member data.
   :return: A :c:type:`SUNErrCode` indicating success or failure.


.. c:function:: SUNErrCode SUNStepper_GetContent(SUNStepper stepper, void **content)

   This function retrieves the member data (content) pointer from a
   :c:type:`SUNStepper` object.

   :param stepper: a stepper object.
   :param content: a pointer to set to the stepper member data pointer.
   :return: A :c:type:`SUNErrCode` indicating success or failure.


Handling Warnings and Errors
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

An implementation of a :c:type:`SUNStepper` may have a system of warning and error
handling that cannot be encoded as a :c:type:`SUNErrCode` which is the return
type of all :c:type:`SUNStepper` functions. Therefore, we provide the following
function to get and set a separate flag associated with a stepper.

.. c:function:: SUNErrCode SUNStepper_SetLastFlag(SUNStepper stepper, int last_flag)

   This function sets a flag that can be used by :c:type:`SUNStepper` implementations to
   indicate warnings or errors that occurred during an operation, e.g.,
   :c:func:`SUNStepper_Evolve`.

   :param stepper: the stepper object.
   :param last_flag: the flag value.
   :return: A :c:type:`SUNErrCode` indicating success or failure.

.. c:function:: SUNErrCode SUNStepper_GetLastFlag(SUNStepper stepper, int *last_flag)

   This function provides the last value of the flag used by the :c:type:`SUNStepper`
   implementation to indicate warnings or errors that occurred during an
   operation, e.g., :c:func:`SUNStepper_Evolve`.

   :param stepper: the stepper object.
   :param last_flag: A pointer to where the flag value will be written.
   :return: A :c:type:`SUNErrCode` indicating success or failure.


.. _SUNStepper.Description.BaseMethods.AttachFunctions:

Setting Member Functions
^^^^^^^^^^^^^^^^^^^^^^^^

The functions in this section are used to specify how each operation on a
:c:type:`SUNStepper` implementation is performed. Technically, all of these
functions are optional to call; the functions that need to be attached are
determined by the "consumer" of the :c:type:`SUNStepper`.

.. c:function:: SUNErrCode SUNStepper_SetEvolveFn(SUNStepper stepper, SUNStepperEvolveFn fn)

   This function attaches a :c:type:`SUNStepperEvolveFn` function to a
   :c:type:`SUNStepper` object.

   :param stepper: a stepper object.
   :param fn: the :c:type:`SUNStepperEvolveFn` function to attach.
   :return: A :c:type:`SUNErrCode` indicating success or failure.


.. c:function:: SUNErrCode SUNStepper_SetFullRhsFn(SUNStepper stepper, SUNStepperFullRhsFn fn)

   This function attaches a :c:type:`SUNStepperFullRhsFn` function to a
   :c:type:`SUNStepper` object.

   :param stepper: a stepper object.
   :param fn: the :c:type:`SUNStepperFullRhsFn` function to attach.
   :return: A :c:type:`SUNErrCode` indicating success or failure.


.. c:function:: SUNErrCode SUNStepper_SetResetFn(SUNStepper stepper, SUNStepperResetFn fn)

   This function attaches a :c:type:`SUNStepperResetFn` function to a
   :c:type:`SUNStepper` object.

   :param stepper: a stepper object.
   :param fn: the :c:type:`SUNStepperResetFn` function to attach.
   :return: A :c:type:`SUNErrCode` indicating success or failure.


.. c:function:: SUNErrCode SUNStepper_SetStopTimeFn(SUNStepper stepper, SUNStepperSetStopTimeFn fn)

   This function attaches a :c:type:`SUNStepperSetStopTimeFn` function to a
   :c:type:`SUNStepper` object.

   :param stepper: a stepper object.
   :param fn: the :c:type:`SUNStepperSetStopTimeFn` function to attach.
   :return: A :c:type:`SUNErrCode` indicating success or failure.


.. c:function:: SUNErrCode SUNStepper_SetForcingFn(SUNStepper stepper, SUNStepperSetForcingFn fn)

   This function attaches a :c:type:`SUNStepperSetForcingFn` function to a
   :c:type:`SUNStepper` object.

   :param stepper: a stepper object.
   :param fn: the :c:type:`SUNStepperSetForcingFn` function to attach.
   :return: A :c:type:`SUNErrCode` indicating success or failure.


.. _SUNStepper.Description.ImplMethods:

Implementation Specific Methods
-------------------------------

This section describes the virtual methods defined by the :c:type:`SUNStepper`
abstract base class.


.. c:type:: SUNErrCode (*SUNStepperEvolveFn)(SUNStepper stepper, sunrealtype tout, N_Vector v, sunrealtype* tret, int* stop_reason)

   This type represents a function with the signature of
   :c:func:`SUNStepper_Evolve`.


.. c:type:: SUNErrCode (*SUNStepperFullRhsFn)(SUNStepper stepper, sunrealtype t, N_Vector v, N_Vector f)

   This type represents a function with the signature of
   :c:func:`SUNStepper_FullRhs`.


.. c:type:: SUNErrCode (*SUNStepperResetFn)(SUNStepper stepper, sunrealtype tR, N_Vector vR)

   This type represents a function with the signature of
   :c:func:`SUNStepper_Reset`.


.. c:type:: SUNErrCode (*SUNStepperSetStopTimeFn)(SUNStepper stepper, sunrealtype tstop)

   This type represents a function with the signature of
   :c:func:`SUNStepper_SetStopTime`.


.. c:type:: SUNErrCode (*SUNStepperSetForcingFn)(SUNStepper stepper, sunrealtype tshift, sunrealtype tscale, N_Vector* forcing, int nforcing)

   This type represents a function with the signature of
   :c:func:`SUNStepper_SetForcing`.

