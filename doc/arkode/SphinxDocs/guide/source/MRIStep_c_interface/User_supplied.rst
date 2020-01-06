..
   Programmer(s): David J. Gardner @ LLNL
   ----------------------------------------------------------------
   Based on ERKStep by Daniel R. Reynolds @ SMU
   ----------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2002-2020, Lawrence Livermore National Security
   and Southern Methodist University.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   ----------------------------------------------------------------

:tocdepth: 3



.. _MRIStep_CInterface.UserSupplied:

User-supplied functions
=============================

The user-supplied functions for MRIStep consist of:

* functions that defines the ODE (required),

* a function that handles error and warning messages (optional),

..
   * a function that provides the error weight vector (optional),

..
   * a function that handles adaptive time step error control (optional),

..
   * a function that handles explicit time step stability (optional),

* a function that defines the root-finding problem(s) to solve
  (optional),

* a function that handles vector resizing operations, if the
  underlying vector structure supports resizing (as opposed to
  deletion/recreation), and if the user plans to call
  :c:func:`MRIStepResize()` (optional).




.. _MRIStep_CInterface.ODERHS:

ODE right-hand side
-----------------------------

The user must supply two functions of type :c:type:`ARKRhsFn` to
specify the right-hand side of the ODE system:


.. c:type:: typedef int (*ARKRhsFn)(realtype t, N_Vector y, N_Vector ydot, void* user_data)

   This function computes a portion of the ODE right-hand side for a given
   value of the independent variable :math:`t` and state vector :math:`y`.

   **Arguments:**
      * *t* -- the current value of the independent variable.
      * *y* -- the current value of the dependent variable vector.
      * *ydot* -- the output vector that forms a portion the ODE RHS :math:`f(t,y)`.
      * *user_data* -- the `user_data` pointer that was passed to :c:func:`MRIStepSetUserData()`.

   **Return value:**
   An *ARKRhsFn* should return 0 if successful, a positive value if a
   recoverable error occurred, or a negative value if it failed unrecoverably.
   As the MRIStep module only supports fixed step sizes at this time any
   non-zero return value will halt the integration.

   **Notes:** Allocation of memory for `ydot` is handled within the
   MRIStep module.  A recoverable failure error return from the
   *ARKRhsFn* is typically used to flag a value of the dependent
   variable :math:`y` that is "illegal" in some way (e.g., negative
   where only a non-negative value is physically meaningful).

   ..
      If such a return is made, MRIStep will attempt to recover by reducing the
      step size in order to avoid this recoverable error return.  There are some
      situations in which recovery is not possible even if the right-hand side
      function returns a recoverable error flag. One is when this occurs at the
      very first call to the *ARKRhsFn* (in which case MRIStep returns
      *ARK_FIRST_RHSFUNC_ERR*).




.. _MRIStep_CInterface.ErrorHandler:

Error message handler function
--------------------------------------

As an alternative to the default behavior of directing error and
warning messages to the file pointed to by `errfp` (see
:c:func:`MRIStepSetErrFile()`), the user may provide a function of type
:c:type:`ARKErrHandlerFn` to process any such messages.



.. c:type:: typedef void (*ARKErrHandlerFn)(int error_code, const char* module, const char* function, char* msg, void* user_data)

   This function processes error and warning messages from
   MRIStep and its sub-modules.

   **Arguments:**
      * *error_code* -- the error code.
      * *module* -- the name of the MRIStep module reporting the error.
      * *function* -- the name of the function in which the error occurred.
      * *msg* -- the error message.
      * *user_data* -- a pointer to user data, the same as the
        *eh_data* parameter that was passed to :c:func:`MRIStepSetErrHandlerFn()`.

   **Return value:**
   An *ARKErrHandlerFn* function has no return value.

   **Notes:** *error_code* is negative for errors and positive
   (*ARK_WARNING*) for warnings.  If a function that returns a
   pointer to memory encounters an error, it sets *error_code* to
   0.




..
   .. _MRIStep_CInterface.ErrorWeight:

   Error weight function
   --------------------------------------

   As an alternative to providing the relative and absolute tolerances,
   the user may provide a function of type :c:type:`ARKEwtFn` to compute a
   vector *ewt* containing the weights in the WRMS norm
   :math:`\|v\|_{WRMS} = \left(\frac{1}{n} \sum_{i=1}^n \left(ewt_i\; v_i\right)^2
   \right)^{1/2}`.  These weights will be used in place of those defined
   in the section :ref:`Mathematics.Error.Norm`.



   .. c:type:: typedef int (*ARKEwtFn)(N_Vector y, N_Vector ewt, void* user_data)

      This function computes the WRMS error weights for the vector
      :math:`y`.

      **Arguments:**
         * *y* -- the dependent variable vector at which the
           weight vector is to be computed.
         * *ewt* -- the output vector containing the error weights.
         * *user_data* -- a pointer to user data, the same as the
           *user_data* parameter that was passed to :c:func:`MRIStepSetUserData()`.

      **Return value:**
      An *ARKEwtFn* function must return 0 if it
      successfully set the error weights, and -1 otherwise.

      **Notes:** Allocation of memory for *ewt* is handled within MRIStep.

      The error weight vector must have all components positive.  It is
      the user's responsibility to perform this test and return -1 if it
      is not satisfied.



   .. _MRIStep_CInterface.AdaptivityFn:

   Time step adaptivity function
   --------------------------------------

   As an alternative to using one of the built-in time step adaptivity
   methods for controlling solution error, the user may provide a
   function of type :c:type:`ARKAdaptFn` to compute a target step size
   :math:`h` for the next integration step.  These steps should be chosen
   as the maximum value such that the error estimates remain below 1.



   .. c:type:: typedef int (*ARKAdaptFn)(N_Vector y, realtype t, realtype h1, realtype h2, realtype h3, realtype e1, realtype e2, realtype e3, int q, int p, realtype* hnew, void* user_data)

      This function implements a time step adaptivity algorithm
      that chooses :math:`h` satisfying the error tolerances.

      **Arguments:**
         * *y* -- the current value of the dependent variable vector.
         * *t* -- the current value of the independent variable.
         * *h1* -- the current step size, :math:`t_n - t_{n-1}`.
         * *h2* -- the previous step size, :math:`t_{n-1} - t_{n-2}`.
         * *h3* -- the step size :math:`t_{n-2}-t_{n-3}`.
         * *e1* -- the error estimate from the current step, :math:`n`.
         * *e2* -- the error estimate from the previous step, :math:`n-1`.
         * *e3* -- the error estimate from the step :math:`n-2`.
         * *q* -- the global order of accuracy for the method.
         * *p* -- the global order of accuracy for the embedded method.
         * *hnew* -- the output value of the next step size.
         * *user_data* -- a pointer to user data, the same as the
           *h_data* parameter that was passed to :c:func:`MRIStepSetAdaptivityFn()`.

      **Return value:**
      An *ARKAdaptFn* function should return 0 if it
      successfully set the next step size, and a non-zero value otherwise.




   .. _MRIStep_CInterface.StabilityFn:

   Explicit stability function
   --------------------------------------

   A user may supply a function to predict the maximum stable step size
   for the explicit Runge Kutta method on this problem.  While the
   accuracy-based time step adaptivity algorithms may be sufficient
   for retaining a stable solution to the ODE system, these may be
   inefficient if :math:`f(t,y)` contains moderately stiff terms.  In
   this scenario, a user may provide a function of type :c:type:`ARKExpStabFn`
   to provide this stability information to MRIStep.  This function
   must set the scalar step size satisfying the stability restriction for
   the upcoming time step.  This value will subsequently be bounded by
   the user-supplied values for the minimum and maximum allowed time
   step, and the accuracy-based time step.



   .. c:type:: typedef int (*ARKExpStabFn)(N_Vector y, realtype t, realtype* hstab, void* user_data)

      This function predicts the maximum stable step size for the ODE system.

      **Arguments:**
         * *y* -- the current value of the dependent variable vector.
         * *t* -- the current value of the independent variable.
         * *hstab* -- the output value with the absolute value of the
           maximum stable step size.
         * *user_data* -- a pointer to user data, the same as the
           *estab_data* parameter that was passed to :c:func:`MRIStepSetStabilityFn()`.

      **Return value:**
      An *ARKExpStabFn* function should return 0 if it
      successfully set the upcoming stable step size, and a non-zero
      value otherwise.

      **Notes:**  If this function is not supplied, or if it returns
      *hstab* :math:`\le 0.0`, then MRIStep will assume that there is no explicit
      stability restriction on the time step size.



.. _MRIStep_CInterface.RootfindingFn:

Rootfinding function
--------------------------------------

If a rootfinding problem is to be solved during the integration of the
ODE system, the user must supply a function of type :c:type:`ARKRootFn`.



.. c:type:: typedef int (*ARKRootFn)(realtype t, N_Vector y, realtype* gout, void* user_data)

   This function implements a vector-valued function
   :math:`g(t,y)` such that the roots of the *nrtfn* components
   :math:`g_i(t,y)` are sought.

   **Arguments:**
      * *t* -- the current value of the independent variable.
      * *y* -- the current value of the dependent variable vector.
      * *gout* -- the output array, of length *nrtfn*, with components :math:`g_i(t,y)`.
      * *user_data* -- a pointer to user data, the same as the
        *user_data* parameter that was passed to :c:func:`MRIStepSetUserData()`.

   **Return value:**
   An *ARKRootFn* function should return 0 if successful
   or a non-zero value if an error occurred (in which case the
   integration is halted and MRIStep returns *ARK_RTFUNC_FAIL*).

   **Notes:** Allocation of memory for *gout* is handled within MRIStep.



.. _MRIStep_CInterface.VecResizeFn:

Vector resize function
--------------------------------------

For simulations involving changes to the number of equations and
unknowns in the ODE system (e.g. when using spatial adaptivity in a
PDE simulation), the MRIStep integrator may be "resized" between
integration steps, through calls to the :c:func:`MRIStepResize()`
function. Typically, when performing adaptive simulations the solution
is stored in a customized user-supplied data structure, to enable
adaptivity without repeated allocation/deallocation of memory.  In
these scenarios, it is recommended that the user supply a customized
vector kernel to interface between SUNDIALS and their problem-specific
data structure.  If this vector kernel includes a function of type
:c:type:`ARKVecResizeFn` to resize a given vector implementation, then
this function may be supplied to :c:func:`MRIStepResize()` so that all
internal MRIStep vectors may be resized, instead of deleting and
re-creating them at each call.  This resize function should have the
following form:


.. c:type:: typedef int (*ARKVecResizeFn)(N_Vector y, N_Vector ytemplate, void* user_data)

   This function resizes the vector *y* to match the dimensions of the
   supplied vector, *ytemplate*.

   **Arguments:**
      * *y* -- the vector to resize.
      * *ytemplate* -- a vector of the desired size.
      * *user_data* -- a pointer to user data, the same as the
        *resize_data* parameter that was passed to :c:func:`MRIStepResize()`.

   **Return value:**
   An *ARKVecResizeFn* function should return 0 if it successfully
   resizes the vector *y*, and a non-zero value otherwise.

   **Notes:**  If this function is not supplied, then MRIStep will
   instead destroy the vector *y* and clone a new vector *y* off of
   *ytemplate*.

.. _MRIStep_CInterface.PreInnerFn:

Pre inner integrator communication function
--------------------------------------------

The user may supply a function of type :c:type:`MRIStepPreInnerFn` that will be
called *before* each inner integration to perform any communication or
memory transfers of forcing data supplied by the the outer integrator to inner
integrator for the inner integration.


.. c:type:: typedef int (*MRIStepPreInnerFn)(realtype t, N_Vector* f, int num_vecs, void* user_data)

   **Arguments:**
      * *t* -- the current value of the independent variable.
      * *f* -- an ``N_Vector`` array of outer forcing vectors.
      * *num_vecs* -- the number of vectors in the ``N_Vector`` array.
      * *user_data* -- the `user_data` pointer that was passed to
        :c:func:`MRIStepSetUserData()`.

   **Return value:**
   An *MRIStepPreInnerFn* function should return 0 if successful, a positive value
   if a recoverable error occurred, or a negative value if an unrecoverable
   error occurred. As the MRIStep module only supports fixed step sizes at this
   time any non-zero return value will halt the integration.

   **Notes:**
   In a heterogeneous computing environment if any data copies between the host
   and device vector data are necessary, this is where that should occur.


.. _MRIStep_CInterface.PostInnerFn:

Post inner integrator communication function
---------------------------------------------

The user may supply a function of type :c:type:`MRIStepPostInnerFn` that will be
called *after* each inner integration to perform any communication or
memory transfers of state data supplied by the inner integrator to the
outer integrator for the outer integration.


.. c:type:: typedef int (*MRIStepPostInnerFn)(realtype t, N_Vector y, void* user_data)

   **Arguments:**
      * *t* -- the current value of the independent variable.
      * *y* -- the current value of the dependent variable vector.
      * *user_data* -- the `user_data` pointer that was passed to
        :c:func:`MRIStepSetUserData()`.

   **Return value:**
   An *MRIStepPostInnerFn* function should return 0 if successful, a positive value
   if a recoverable error occurred, or a negative value if an unrecoverable
   error occurred. As the MRIStep module only supports fixed step sizes at this
   time any non-zero return value will halt the integration.

   **Notes:**
   In a heterogeneous computing environment if any data copies between the host
   and device vector data are necessary, this is where that should occur.
