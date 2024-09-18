.. ----------------------------------------------------------------
   Programmer(s): David J. Gardner @ LLNL
                  Steven B. Roberts @LLNL
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
class member data and a structure of function pointers the derived class
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
retrieving the derived class ``content`` pointer, setting function pointers to
derived class method implementations, and accessing base class data e.g., for
computing the forcing term :eq:`SUNStepper_forcing`.

.. _SUNStepper.Description.BaseMethods.CreateDestroy:

Creating and Destroying an Object
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. c:function:: SUNErrCode SUNStepper_Create(SUNContext sunctx, SUNStepper *stepper)

   This function creates a :c:type:`SUNStepper` object to which a user should
   attach the member data (content) pointer and method function pointers.

   **Arguments:**
      * ``sunctx`` -- the SUNDIALS simulation context.
      * ``stepper`` -- a pointer to a stepper object.

   **Return value:**
      * A :c:type:`SUNErrCode` indicating success or failure.

   **Example usage:**

   .. code-block:: C

      /* create an instance of the base class */
      SUNStepper stepper = NULL;
      err = SUNStepper_Create(&stepper);

   .. note::

      See :numref:`SUNStepper.Description.BaseMethods.Content` and
      :numref:`SUNStepper.Description.BaseMethods.AttachFunctions`
      for details on how to attach member data and method function pointers.


.. c:function:: SUNErrCode SUNStepper_Destroy(SUNStepper *stepper)

   This function destroys a :c:type:`SUNStepper` object.

   **Arguments:**
      * *stepper* -- a pointer to a stepper object.

   **Return value:**
      * A :c:type:`SUNErrCode` indicating success or failure.

   **Example usage:**

   .. code-block:: C

      /* destroy an instance of the base class */
      err = SUNStepper_Destroy(&stepper);

   .. note::

      This function only frees memory allocated within the base class and the
      base class structure itself. The user is responsible for freeing any
      memory allocated for the member data (content).

.. _SUNStepper.Description.BaseMethods.Content:

Attaching and Accessing the Content Pointer
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. c:function:: SUNErrCode SUNStepper_SetContent(SUNStepper stepper, void *content)

   This function attaches a member data (content) pointer to a
   :c:type:`SUNStepper` object.

   **Arguments:**
      * *stepper* -- a stepper object.
      * *content* -- a pointer to the stepper member data.

   **Return value:**
      * A :c:type:`SUNErrCode` indicating success or failure.

   **Example usage:**

   .. code-block:: C

      /* set the stepper content pointer */
      MyStepperContent my_object_data;
      err = SUNStepper_SetContent(stepper, &my_object_data);


.. c:function:: SUNErrCode SUNStepper_GetContent(SUNStepper stepper, void **content)

   This function retrieves the member data (content) pointer from a
   :c:type:`SUNStepper` object.

   **Arguments:**
      * *stepper* -- a stepper object.
      * *content* -- a pointer to set to the stepper member data pointer.

   **Return value:**
      * A :c:type:`SUNErrCode` indicating success or failure.

   **Example usage:**

   .. code-block:: C

      /* get the stepper content pointer */
      void             *content;
      MyStepperContent *my_object_data;

      err = SUNStepper_GetContent(stepper, &content);
      my_object_data = (MyStepperContent*) content;


.. _SUNStepper.Description.BaseMethods.AttachFunctions:

Setting Member Functions
^^^^^^^^^^^^^^^^^^^^^^^^

.. c:function:: SUNErrCode SUNStepper_SetEvolveFn(SUNStepper stepper, SUNStepperEvolveFn fn)

   This function attaches a :c:type:`SUNStepperEvolveFn` function to a
   :c:type:`SUNStepper` object.

   **Arguments:**
      * *stepper* -- a stepper object.
      * *fn* -- the :c:type:`SUNStepperEvolveFn` function to attach.

   **Return value:**
      * A :c:type:`SUNErrCode` indicating success or failure.

   **Example usage:**

   .. code-block:: C

      /* set the stepper evolve function */
      err = SUNStepper_SetEvolveFn(stepper, MyEvolve);


.. c:function:: SUNErrCode SUNStepper_SetFullRhsFn(SUNStepper stepper, SUNStepperFullRhsFn fn)

   This function attaches a :c:type:`SUNStepperFullRhsFn` function to a
   :c:type:`SUNStepper` object.

   **Arguments:**
      * *stepper* -- a stepper object.
      * *fn* -- the :c:type:`SUNStepperFullRhsFn` function to attach.

   **Return value:**
      * A :c:type:`SUNErrCode` indicating success or failure.

   **Example usage:**

   .. code-block:: C

      /* set the stepper full right-hand side function */
      err = SUNStepper_SetFullRhsFn(stepper, MyFullRHS);


.. c:function:: SUNErrCode SUNStepper_SetResetFn(SUNStepper stepper, SUNStepperResetFn fn)

   This function attaches a :c:type:`SUNStepperResetFn` function to a
   :c:type:`SUNStepper` object.

   **Arguments:**
      * *stepper* -- a stepper object.
      * *fn* -- the :c:type:`SUNStepperResetFn` function to attach.

   **Return value:**
      * A :c:type:`SUNErrCode` indicating success or failure.

   **Example usage:**

   .. code-block:: C

      /* set the stepper reset function */
      err = SUNStepper_SetResetFn(stepper, MyReset);


.. _SUNStepper.Description.BaseMethods.Forcing:

Applying and Accessing Forcing Data
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

When integrating the ODE :eq:`SUNStepper_IVP` the :c:type:`SUNStepper` is
responsible for evaluating ODE right-hand side function :math:`f(t, v)` as well
as computing and applying the forcing term :eq:`SUNStepper_forcing` to obtain
the full right-hand side of the ODE :eq:`SUNStepper_IVP`. The functions in this
section can be used to either apply the forcing or access the data necessary to
construct the forcing polynomial.


.. c:function:: SUNErrCode SUNStepper_AddForcing(SUNStepper stepper, sunrealtype t, N_Vector f)

   This function computes the forcing term :eq:`SUNStepper_forcing` at the input
   time *t* and adds it to input vector *f*, i.e., the right-hand side vector.

   **Arguments:**
      * *stepper* -- a stepper object.
      * *t* -- the time at which the forcing should be evaluated.
      * *f* -- the vector to which the forcing should be applied.

   **Return value:**
      * A :c:type:`SUNErrCode` indicating success or failure.

   **Example usage:**

   .. code-block:: C

      /* compute the forcing term and add it the fast RHS vector */
      err = SUNStepper_AddForcing(stepper, t, f);


.. c:function:: SUNErrCode SUNStepper_GetForcingData(SUNStepper stepper, sunrealtype *tshift, sunrealtype *tscale, N_Vector **forcing, int *nforcing)

   This function provides access to data necessary to compute the forcing term
   :eq:`SUNStepper_forcing`. This includes the shift and scaling factors for the
   normalized time :math:`\frac{t - t_{\text{shift}}}{t_{\text{scale}}}` and the
   array of polynomial coefficient vectors :math:`\widehat{f}_k`.

   **Arguments:**
      * *stepper* -- a stepper object.
      * *tshift* -- the time shift to apply to the current time when computing
        the forcing, :math:`t_{\text{shift}}`.
      * *tscale* -- the time scaling to apply to the current time when computing
        the forcing, :math:`t_{\text{scale}}`.
      * *forcing* -- a pointer to an array of forcing vectors,
        :math:`\widehat{f}_k`.
      * *nforcing* -- the number of forcing vectors, :math:`n_{\text{forcing}}`.

   **Return value:**
      * A :c:type:`SUNErrCode` indicating success or failure.

   **Example usage:**

   .. code-block:: C

      SUNErrCode err
      int        k;
      int        nforcing_vecs;   /* number of forcing vectors */
      double     tshift, tscale;  /* time normalization values */
      double     tau;             /* normalized time           */
      double     tau_k;           /* tau raised to the power k */
      N_Vector   *forcing_vecs;   /* array of forcing vectors  */

      /* get the forcing data from the stepper */
      err = SUNStepper_GetForcingData(stepper, &tshift, &tscale,
                                      &forcing_vecs, &nforcing_vecs);

      /* compute the normalized time, initialize tau^k */
      tau   = (t - tshift) / tscale;
      tau_k = 1.0;

      /* compute the polynomial forcing terms and add them to fast RHS vector */
      for (k = 0; k < nforcing_vecs; k++)
      {
        N_VLinearSum(1.0, f_fast, tau_k, forcing_vecs[k], f_fast);
        tau_k *= tau;
      }


.. _SUNStepper.Description.ImplMethods:

Implementation Specific Methods
-------------------------------

This section describes the required and optional virtual methods defined by the
:c:type:`SUNStepper` abstract base class.

Required Member Functions
^^^^^^^^^^^^^^^^^^^^^^^^^

An :c:type:`SUNStepper` *must* provide implementations of the following
member functions:


.. c:type:: SUNErrCode (*SUNStepperEvolveFn)(SUNStepper stepper, sunrealtype t0, sunrealtype tout, N_Vector v)

   This function advances the state vector *v* for the ODE system from time *t0*
   to time *tout*.

   **Arguments:**
      * *stepper* -- the stepper object.
      * *t0* -- the initial time for the integration.
      * *tout* -- the final time for the integration.
      * *v* -- on input the state at time *t0* and, on output, the state at time
        *tout*.

   **Return value:**
      * A :c:type:`SUNErrCode` indicating success or failure.


.. c:type:: SUNErrCode (*SUNStepperResetFn)(SUNStepper stepper, sunrealtype tR, N_Vector vR)

   This function resets the stepper state to the provided independent variable
   value and dependent variable vector.

   **Arguments:**
      * *stepper* -- the stepper object.
      * *tR* -- the value of the independent variable :math:`t_R`.
      * *vR* -- the value of the dependent variable vector :math:`v(t_R)`.

   **Return value:**
      * A :c:type:`SUNErrCode` indicating success or failure.


Optional Member Functions
^^^^^^^^^^^^^^^^^^^^^^^^^

An :c:type:`SUNStepper` *may* provide implementations of any of the following
member functions:

.. c:type:: SUNErrCode (*SUNStepperFullRhsFn)(SUNStepper stepper, sunrealtype t, N_Vector v, N_Vector f, int mode)

   This function computes the full right-hand side function of the ODE,
   :math:`f(t, v)` in :eq:`SUNStepper_IVP` for a given value of the independent
   variable *t* and state vector *y*.

   **Arguments:**
      * *stepper* -- the stepper object.
      * *t* -- the current value of the independent variable.
      * *v* -- the current value of the dependent variable vector.
      * *f* -- the output vector for the ODE right-hand side, :math:`f(t, v)`,
        in :eq:`SUNStepper_IVP`.
      * *mode* -- a flag indicating the purpose for which the right-hand side
        function evaluation is called.

        * ``ARK_FULLRHS_START`` -- called at the beginning of the simulation
        * ``ARK_FULLRHS_END``   -- called at the end of a successful step
        * ``ARK_FULLRHS_OTHER`` -- called elsewhere e.g., for dense output

   **Return value:**
      * A :c:type:`SUNErrCode` indicating success or failure.
