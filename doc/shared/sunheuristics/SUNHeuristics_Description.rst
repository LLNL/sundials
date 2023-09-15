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

.. _SUNHeuristics.Description:

The SUNHeuristics API
=====================

The SUNHeuristics base class provides a common API for time step heuristic
constraints that may be applied to step sizes by SUNDIALS integrators.  These
heuristics can include a variety of things, including bounds on how much the
time step size can grow or shrink in a single step, absolute maximum and minimum
allowable step sizes, stability-limitations (e.g., CFL) on explicit steps, etc.

The base ``SUNHeuristics`` class is modeled after SUNDIALS' other
object-oriented classes, in that this class contains a pointer to an
implementation-specific *content*, an *ops* structure with generic heuristics
operations, and a :c:type:`SUNContext` object.  Specifically, the type
``SUNHeuristics`` is defined as:

.. c:type:: struct generic_SUNHeuristics_ *SUNHeuristics

and the base class structure is defined as

.. code-block:: C

   struct generic_SUNHeuristics_ {
        void* content;
        generic_SUNHeuristics_Ops_* ops;
        SUNContext sunctx;
    };

Here, ``generic_SUNHeuristics_Ops_`` is the pointer to a structure containing
function pointers to the various heuristics operations, and is defined as

.. code-block:: c

    struct generic_SUNHeuristics_Ops_ {

      SUNHeuristics_ID (*getid)(SUNHeuristics H);
      void             (*destroy)(SUNHeuristics H);
      int              (*constrainstep)(SUNHeuristics H, realtype hcur,
                           realtype hnew, realtype *hconstr);
      int              (*etestfail)(SUNHeuristics H, realtype hcur,
                                    realtype hnew, int nef, realtype *hconstr);
      int              (*convfail)(SUNHeuristics H, realtype hcur, realtype *hconstr);
      int              (*boundreduction)(SUNHeuristics H, realtype hcur,
                                         realtype hnew, realtype *hconstr);
      int              (*boundfirststep)(SUNHeuristics H, realtype h0,
                                         realtype *h0constr);
      int              (*reset)(SUNHeuristics H);
      int              (*update)(SUNHeuristics H);
      int              (*setdefaults)(SUNHeuristics H);
      int              (*write)(SUNHeuristics H, FILE* fptr);
      int              (*setmaxstep)(SUNHeuristics H, realtype hmax);
      int              (*setminstep)(SUNHeuristics H, realtype hmin);
      int              (*setexpstabfn)(SUNHeuristics H, SUNExpStabFn EStab,
                                       void* estab_data);
      int              (*setcflfraction)(SUNHeuristics H, realtype cfl_frac);
      int              (*setsafetyfactor)(SUNHeuristics H, realtype safety);
      int              (*setmaxgrowth)(SUNHeuristics H, realtype mx_growth);
      int              (*setminreduction)(SUNHeuristics H, realtype eta_min);
      int              (*setfixedstepbounds)(SUNHeuristics H, realtype lb, realtype ub);
      int              (*setmaxfirstgrowth)(SUNHeuristics H, realtype etamx1);
      int              (*setmaxefailgrowth)(SUNHeuristics H, realtype etamxf);
      int              (*setsmallnumefails)(SUNHeuristics H, int small_nef);
      int              (*setmaxcfailgrowth)(SUNHeuristics H, realtype etacf);
      int              (*getnumexpsteps)(SUNHeuristics H, long int* expsteps);
      int              (*getnumaccsteps)(SUNHeuristics H, long int* accsteps);
      int              (*space)(SUNHeuristics H, long int *lenrw, long int *leniw);
    };

.. _SUNHeuristics.Description.heuristicsIDs:

SUNHeuristics Type Definitions
------------------------------

The time integrators in SUNDIALS leverage a myriad of heuristics to achieve
accurate and efficient computations; however, most follow a rather standard
pattern, or request no heuristic control whatsoever.  Therefore, each
SUNHeuristics implementation should indicate its type, informing the integrators
of whether they perform heuristic control or not:

.. c:enum:: SUNHeuristics_ID

   The enumerated type :c:type:`SUNHeuristics_ID` defines the enumeration
   constants for SUNDIALS heuristics implementations

.. c:enumerator:: SUNDIALS_HEURISTICS_STD

   Performs "standard" heuristic stepsize constraints.

.. c:enumerator:: SUNDIALS_HEURISTICS_NULL

   Performs no heuristic constraints.



Additionally, if the heuristics object can be used to limit stepsizes based on
explicit stability, then we define the following user-supplied function type:

.. c:type:: int (*SUNExpStabFn)(realtype *hstab, void *user_data)

   This function predicts a maximum stable step size.

   :param hstab: (output) the absolute value of the maximum stable step size.
   :param user_data: a pointer to user data for evaluation.
   :return: a *SUNExpStabFn* function should return 0 if it is successful, and a
            non-zero value otherwise.




.. _SUNHeuristics.Description.operations:

SUNHeuristics Operations
------------------------

The base SUNHeuristics class defines and implements all SUNHeuristics functions.
Most of these routines are merely wrappers for the operations defined by a
particular SUNHeuristics implementation, which are accessed through the *ops*
field of the ``SUNHeuristics`` structure.  However, the base SUNHeuristics class
also provides the convenience routine

.. c:function:: SUNHeuristics SUNHeuristics_NewEmpty(SUNContext sunctx)

  This function allocates a new generic ``SUNHeuristics`` object and initializes
  its content pointer and the function pointers in the operations structure to
  ``NULL``.

  :param sunctx: the :c:type:`SUNContext` object (see
                 :numref:`SUNDIALS.SUNContext`)
  :returns: If successful, a generic :c:type:`SUNHeuristics` object.  If
            unsuccessful, a ``NULL`` pointer will be returned.


Each of the following routines are *optional* for any specific SUNHeuristics implementation.


.. c:function:: SUNHeuristics_ID SUNHeuristics_GetID(SUNHeuristics H)

   Returns the type identifier for the heuristics object *H*.  Returned values
   are given in Section :numref:`SUNHeuristics.Description.heuristicsIDs`

   :param H: the :c:type:`SUNHeuristics` object.
   :return: :c:type:`SUNHeuristics_ID` type identifier.

   Usage:

   .. code-block:: c

      SUNHeuristics_ID id = SUNHeuristics_GetID(H);

.. c:function:: void SUNHeuristics_Destroy(SUNHeuristics H)

   Deallocates the heuristics object *H*.  If this is not provided by the
   implementation, the base wrapper routine will free both the *content* and
   *ops* objects -- this should be sufficient unless an heuristic implementation
   performs dynamic memory allocation of its own (note that the
   SUNDIALS-provided SUNHeuristics implementations do not need to supply this
   routine).

   :param H: the :c:type:`SUNHeuristics` object.

   Usage:

   .. code-block:: c

      SUNHeuristics_Destroy(H);

.. c:function:: int SUNHeuristics_ConstrainStep(SUNHeuristics H, realtype hcur, realtype hnew, realtype* hconstr)

   Main constraint-application function.  This will attempt to change the step
   *hcur* to *hnew*, applying any heuristic bounds on the step size adjustments.

   :param H: the :c:type:`SUNHeuristics` object.
   :param H: the heuristics object.
   :param hcur: the step size from the previous step attempt.
   :param hnew: the requested step size for the upcoming step attempt.
   :param hconstr: (output) pointer to the constrained step size.
   :return: error code indicating success failure (see
            :numref:`SUNHeuristics.Description.errorCodes`).

   Usage:

   .. code-block:: c

      retval = SUNHeuristics_ConstrainStep(H, hcur, hnew, &hnew);

.. c:function:: int SUNHeuristics_ETestFail(SUNHeuristics H, realtype hcur, realtype hnew, int nef, realtype* hconstr)

   Function to apply constraints following a step with unacceptable temporal
   error.

   :param H: the heuristics object.
   :param hcur: the step size that led to the error test failure.
   :param hnew: the requested step size for the upcoming step attempt (e.g.,
                from a :c:type:`SUNControl` object).
   :param nef: the integrator-provided counter of how many temporal error test
               failures have occurred on this time step.
   :param hconstr: (output) pointer to the constrained step size.
   :return: error code indicating success failure (see
            :numref:`SUNHeuristics.Description.errorCodes`).

   Usage:

   .. code-block:: c

      retval = SUNHeuristics_ETestFail(H, hcur, hnew, nef, &hnew);

.. c:function:: int SUNHeuristics_BoundReduction(SUNHeuristics H, realtype hcur, realtype hnew, realtype *hconstr)

   This ensures that a step size reduction is within user-prescribed bounds.

   :param H: the heuristics object.
   :param hcur: the step size from the previous step attempt.
   :param hnew: the requested step size for the upcoming step attempt (e.g.,
                from a :c:type:`SUNControl` object).
   :param hconstr: (output) pointer to the constrained step size.
   :return:
      * *SUNHEURISTICS_SUCCESS* if successful
      * *SUNHEURISTICS_CANNOT_DECREASE* if a reduction is requested but no
        reduction is possible

   Usage:

   .. code-block:: c

      retval = SUNHeuristics_BoundReduction(H, hcur, hnew, &hnew);

.. c:function:: int SUNHeuristics_BoundFirstStep(SUNHeuristics H, realtype h0, realtype *h0constr)

   This bounds the initial step by user-provided min/max step values.

   :param H: the heuristics object.
   :param h0: the requested initial step size.
   :param h0constr: (output) pointer to the constrained initial step size.
   :return: error code indicating success failure (see
            :numref:`SUNHeuristics.Description.errorCodes`).

   Usage:

   .. code-block:: c

      retval = SUNHeuristics_BoundFirstStep(H, h0, &h0);

.. c:function:: int SUNHeuristics_ConvFail(SUNHeuristics H, realtype hcur, realtype *hconstr)

   Function to apply constraints following a step with an algebraic solver
   convergence failure.

   :param H: the heuristics object.
   :param hcur: the step size that led to the convergence failure.
   :param hconstr: (output) pointer to the constrained step size.
   :return: error code indicating success failure (see
            :numref:`SUNHeuristics.Description.errorCodes`).

   Usage:

   .. code-block:: c

      retval = SUNHeuristics_ConvFail(H, hcur, &hnew);

.. c:function::int SUNHeuristics_Reset(SUNHeuristics H)

   Function to reset the heuristics to its initial state: zeros any internal
   counters, and resets any stepsize growth factor bounds.

   :param H: the heuristics object.
   :return: error code indicating success failure (see
            :numref:`SUNHeuristics.Description.errorCodes`).

   Usage:

   .. code-block:: c

      retval = SUNHeuristics_Reset(H);

.. c:function::int SUNHeuristics_Update(SUNHeuristics H)

   Function to notify the heuristics object that a time step has succeeded,
   indicating e.g. that the stepsize growh factor should should be set to its
   "default" state.

   :param H: the heuristics object.
   :return: error code indicating success failure (see
            :numref:`SUNHeuristics.Description.errorCodes`).

   Usage:

   .. code-block:: c

      retval = SUNHeuristics_Update(H);

.. c:function::int SUNHeuristics_SetDefaults(SUNHeuristics H)

   Function to set the heuristics parameters to their default values.

   :param H: the heuristics object.
   :return: error code indicating success failure (see
            :numref:`SUNHeuristics.Description.errorCodes`).

   Usage:

   .. code-block:: c

      retval = SUNHeuristics_SetDefaults(H);

.. c:function::int SUNHeuristics_Write(SUNHeuristics H, FILE* fptr)

   Writes all controller parameters to the indicated file pointer.

   :param H: the heuristics object.
   :param fptr: the output stream to write the parameters.
   :return: error code indicating success failure (see
            :numref:`SUNHeuristics.Description.errorCodes`).

   Usage:

   .. code-block:: c

      retval = SUNHeuristics_Write(H, stdout);

.. c:function::int SUNHeuristics_SetMaxStep(SUNHeuristics H, realtype hmax)

   Function to inform the heuristics object about a maximum allowed absolute
   step size.

   :param H: the heuristics object.
   :param hmax: maximum absolute step size allowed (:math:`\text{hmax} \le 0`
                implies :math:`\text{hmax}=\infty`).
   :return: error code indicating success failure (see
            :numref:`SUNHeuristics.Description.errorCodes`).

   Usage:

   .. code-block:: c

      retval = SUNHeuristics_SetMaxStep(H, 1.0e-2);

.. c:function::int SUNHeuristics_SetMinStep(SUNHeuristics H, realtype hmin)

   Function to inform the heuristics object about a minimum allowed absolute
   step size.

   :param H: the heuristics object.
   :param hmin: minimum absolute step size allowed (:math:`\text{hmin} \le 0`
                implies no minimum).
   :return: error code indicating success failure (see
            :numref:`SUNHeuristics.Description.errorCodes`).

   Usage:

   .. code-block:: c

      retval = SUNHeuristics_SetMinStep(H, 1.0e-5);

.. c:function::int SUNHeuristics_SetExpStabFn(SUNHeuristics H, SUNExpStabFn EStab, void* estab_data)

   Function to provide a user-supplied function for the maximum stable step
   size.

   :param H: the heuristics object.
   :param EStab: user-supplied function specifying the maximum stable step size
                 (``EStab == NULL`` disables).
   :param estab_data: user-supplied data pointer that should be provided on all
                      calls to *EStab*.
   :return: error code indicating success failure (see
            :numref:`SUNHeuristics.Description.errorCodes`).

   Usage:

   .. code-block:: c

      retval = SUNHeuristics_SetExpStabFn(H, CFLStabilityBound, (void*) (&mystruct));

.. c:function::int SUNHeuristics_SetCFLFraction(SUNHeuristics H, realtype cfl_frac)

   Function to set a CFL safety factor that should be applied to the outputs of
   *EStab*, above.

   :param H: the heuristics object.
   :param cfl_frac: safety factor.  Allowable values must be within the open
                    interval (0,1), :math:`\text{cfl\_frac} \le 0` implies a
                    reset to the default value.
   :return: error code indicating success failure (see
            :numref:`SUNHeuristics.Description.errorCodes`).

   Usage:

   .. code-block:: c

      retval = SUNHeuristics_SetCFLFraction(H, 0.5);

.. c:function::int SUNHeuristics_SetSafetyFactor(SUNHeuristics C, realtype safety)

   Function to set a step size safety factor that should be used to constrain an
   error-controller-recommended step size.

   :param H: the heuristics object.
   :param safety: safety factor.  Allowable values must be within the open
                  interval (0,1), :math:`\text{safety} \le 0` implies a reset
                  to the default value.
   :return: error code indicating success failure (see
            :numref:`SUNHeuristics.Description.errorCodes`).

   Usage:

   .. code-block:: c

      retval = SUNHeuristics_SetSafetyFactor(H, 0.95);

.. c:function::int SUNHeuristics_SetMaxGrowth(SUNHeuristics H, realtype mx_growth)

   Function to set maximum stepsize growth factor for general steps.

   :param H: the heuristics object.
   :param mx_growth: maximum step size growth factor.  Allowable values must be
                     greater than 1; illegal values imply a reset to the
                     default.
   :return: error code indicating success failure (see
            :numref:`SUNHeuristics.Description.errorCodes`).

   **Notes:**
      The first step uses a separate maximum growth factor.

   .. code-block:: c

      retval = SUNHeuristics_SetMaxGrowth(H, 20.0);

.. c:function::int SUNHeuristics_SetMaxFirstGrowth(SUNHeuristics H, realtype etamx1)

   Function to set maximum stepsize growth factor for the first internal time
   step.

   :param H: the heuristics object.
   :param etamx1: maximum first step size growth factor.  Allowable values must
                  be greater than 1; illegal values imply a reset to the
                  default.
   :return: error code indicating success failure (see
            :numref:`SUNHeuristics.Description.errorCodes`).

   Usage:

   .. code-block:: c

      retval = SUNHeuristics_SetMaxFirstGrowth(H, 10000.0);

.. c:function::int SUNHeuristics_SetFixedStepBounds(SUNHeuristics H, realtype lb, realtype ub)

   Function to specify the step size growth interval within which the step size
   will remain unchanged.  Allowable values must enclose the value 1.0.  Any
   illegal interval implies a reset to the default.

   :param H: the heuristics object.
   :param lb: interval lower bound.
   :param ub: interval upper bound.
   :return: error code indicating success failure (see
            :numref:`SUNHeuristics.Description.errorCodes`).

   Usage:

   .. code-block:: c

      retval = SUNHeuristics_SetFixedStepBounds(H, 1.0, 1.5);

.. c:function::int SUNHeuristics_SetMinReduction(SUNHeuristics H, realtype eta_min)

   Function to set a lower bound for the stepsize adjustment factor following a
   temporal error failure.

   :param H: the heuristics object.
   :param eta_min: minimum step size reduction factor.  Allowable values must be
                   in the open interval (0,1); illegal values imply a reset to
                   the default.
   :return: error code indicating success failure (see
            :numref:`SUNHeuristics.Description.errorCodes`).

   Usage:

   .. code-block:: c

      retval = SUNHeuristics_SetMinReduction(H, 1.0e-2);

.. c:function::int SUNHeuristics_SetMaxEFailGrowth(SUNHeuristics H, realtype etamxf)

   Function to set an upper bound for the stepsize adjustment factor following a
   temporal error failure (once at least *small_nef* failures have occurred, see
   :c:func:`SUNHeuristics_SetSmallNumEFails`).

   :param H: the heuristics object.
   :param etamxf: temporal error failure step size reduction factor.  Allowable
                  values must be in the interval (0,1]; illegal values imply a
                  reset to the default.
   :return: error code indicating success failure (see
            :numref:`SUNHeuristics.Description.errorCodes`).

   Usage:

   .. code-block:: c

      retval = SUNHeuristics_SetMaxEFailGrowth(H, 0.8);

.. c:function::int SUNHeuristics_SetSmallNumEFails(SUNHeuristics H, int small_nef)

   Function to specify the step adaptivity constant *small_nef* (see
   :c:func:`SUNHeuristics_SetMaxEFailGrowth`).

   :param H: the heuristics object.
   :param small_nef: minimum number of error test failures before *etamxf* takes
                     effect. Allowable values are greater than 0. Illegal
                     values imply a reset to the default.
   :return: error code indicating success failure (see
            :numref:`SUNHeuristics.Description.errorCodes`).

   Usage:

   .. code-block:: c

      retval = SUNHeuristics_SetSmallNumEFails(H, 3);

.. c:function::int SUNHeuristics_SetMaxCFailGrowth(SUNHeuristics H, realtype etacf)

   Function to specify an upper bound for the stepsize adjustment factor
   following an algebraic solver convergence failure.

   :param H: the heuristics object.
   :param etacf: convergence failure step adjustment factor. Allowable values
                 must be in the interval (0,1]; illegal values imply a reset to
                 the default.
   :return: error code indicating success failure (see
            :numref:`SUNHeuristics.Description.errorCodes`).

   Usage:

   .. code-block:: c

      retval = SUNHeuristics_SetMaxCFailGrowth(H, 0.3);

.. c:function::int SUNHeuristics_GetNumExpSteps(SUNHeuristics H, long int* expsteps)

   Function to return the current number of stability-limited steps.

   :param H: the heuristics object.
   :param expsteps: (output) pointer to hold the stability-limited step counter.
   :return: error code indicating success failure (see
            :numref:`SUNHeuristics.Description.errorCodes`).

   Usage:

   .. code-block:: c

      retval = SUNHeuristics_GetNumExpSteps(H, &expsteps);

.. c:function::int SUNHeuristics_GetNumAccSteps(SUNHeuristics H, long int* accsteps)

   Function to return the current number of accuracy-limited steps.

   :param H: the heuristics object.
   :param expsteps: (output) pointer to hold the accuracy-limited step counter.
   :return: error code indicating success failure (see
            :numref:`SUNHeuristics.Description.errorCodes`).

   Usage:

   .. code-block:: c

      retval = SUNHeuristics_GetNumAccSteps(H, &accsteps);

.. c:function:: int SUNHeuristics_Space(SUNHeuristics H, long int *lenrw, long int *leniw)

   Informative routine that returns the memory requirements of the heuristics object.

   :param H: the heuristics object.
   :param lenrw: (output) number of ``sunrealtype`` words stored in the
                 heuristics object.
   :param leniw: (output) number of ``sunindextype`` words stored in the
                 heuristics object.  This may also include pointers, `int` and
                 `long int` words.
   :return: error code indicating success failure (see
            :numref:`SUNHeuristics.Description.errorCodes`).

   Usage:

   .. code-block:: c

      retval = SUNHeuristics_Space(H, &lenrw, &leniw);



.. _SUNHeuristics.Description.errorCodes:

SUNHeuristics Error Codes
-------------------------

SUNHeuristics functions return one of the following set of error codes:

* ``SUNHEURISTICS_SUCCESS`` (0) -- successful call.

* ``SUNHEURISTICS_ILL_INPUT`` (-1101) -- an illegal input has been provided to
  the function.

* ``SUNHEURISTICS_CANNOT_DECREASE`` (-1102) -- a step reduction was required but
  was not possible (already at minimum).

* ``SUNHEURISTICS_USER_FCN_FAIL`` (-1103) -- a user-supplied function returned a
  nonzero [error] value.


C/C++ API Usage
---------------

The SUNDIALS Heuristics module can be used in C and C++ programs by including
the header file ``sundials/sundials_heuristisc.h``.

Example usage (here ``SUNHeuristicsXYZ`` is a placeholder for an actual
SUNHeuristics implementation constructor):

.. code-block:: c

    #include <stdio.h>
    #include <stdlib.h>
    #include <sundials/sundials_context.h>
    #include <sundials/sundials_types.h>
    #include <sundials/sundials_heuristics.h>

    int main()
    {
        /* Create a SUNContext object */
        SUNContext sunctx = ...;

        /* Create a SUNHeuristics object */
        SUNHeuristics H = SUNHeuristicsXYZ(sunctx);

        /* Use the heuristics object */

        /* Destroy the heuristics object */
        SUNHeuristics_Destroy(C);

        return 0;
    }
