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

.. _SUNTimestepHeuristics.Description:

The SUNTimestepHeuristics API
=============================

.. versionadded:: x.x.x

The SUNTimestepHeuristics base class provides a common API for time step heuristic
constraints that may be applied to step sizes by SUNDIALS integrators.  These
heuristics can include a variety of things, including bounds on how much the
time step size can grow or shrink in a single step, absolute maximum and minimum
allowable step sizes, stability-limitations (e.g., CFL) on explicit steps, etc.

The base ``SUNTimestepHeuristics`` class is modeled after SUNDIALS' other
object-oriented classes, in that this class contains a pointer to an
implementation-specific *content*, an *ops* structure with generic heuristics
operations, and a :c:type:`SUNContext` object.  Specifically, the type
``SUNTimestepHeuristics`` is defined as:

.. c:type:: struct generic_SUNTimestepHeuristics_ *SUNTimestepHeuristics

and the base class structure is defined as

.. code-block:: C

   struct generic_SUNTimestepHeuristics_ {
        void* content;
        generic_SUNTimestepHeuristics_Ops_* ops;
        SUNContext sunctx;
    };

Here, ``generic_SUNTimestepHeuristics_Ops_`` is the pointer to a structure containing
function pointers to the various heuristics operations, and is defined as

.. code-block:: c

    struct generic_SUNTimestepHeuristics_Ops_ {

      SUNTimestepHeuristics_ID (*getid)(SUNTimestepHeuristics H);
      int (*destroy)(SUNTimestepHeuristics H);
      int (*constrainstep)(SUNTimestepHeuristics H, realtype hcur,
                           realtype hnew, realtype *hconstr);
      int (*etestfail)(SUNTimestepHeuristics H, realtype hcur,
                       realtype hnew, int nef, realtype *hconstr);
      int (*convfail)(SUNTimestepHeuristics H, realtype hcur, realtype *hconstr);
      int (*boundreduction)(SUNTimestepHeuristics H, realtype hcur,
                            realtype hnew, realtype *hconstr);
      int (*boundfirststep)(SUNTimestepHeuristics H, realtype h0,
                            realtype *h0constr);
      int (*reset)(SUNTimestepHeuristics H);
      int (*update)(SUNTimestepHeuristics H);
      int (*setdefaults)(SUNTimestepHeuristics H);
      int (*write)(SUNTimestepHeuristics H, FILE* fptr);
      int (*setmaxstep)(SUNTimestepHeuristics H, realtype hmax);
      int (*setminstep)(SUNTimestepHeuristics H, realtype hmin);
      int (*setexpstabfn)(SUNTimestepHeuristics H, SUNExpStabFn EStab,
                          void* estab_data);
      int (*setcflfraction)(SUNTimestepHeuristics H, realtype cfl_frac);
      int (*setsafetyfactor)(SUNTimestepHeuristics H, realtype safety);
      int (*setmaxgrowth)(SUNTimestepHeuristics H, realtype mx_growth);
      int (*setminreduction)(SUNTimestepHeuristics H, realtype eta_min);
      int (*setfixedstepbounds)(SUNTimestepHeuristics H, realtype lb, realtype ub);
      int (*setmaxfirstgrowth)(SUNTimestepHeuristics H, realtype etamx1);
      int (*setmaxefailgrowth)(SUNTimestepHeuristics H, realtype etamxf);
      int (*setsmallnumefails)(SUNTimestepHeuristics H, int small_nef);
      int (*setmaxcfailgrowth)(SUNTimestepHeuristics H, realtype etacf);
      int (*getnumexpsteps)(SUNTimestepHeuristics H, long int* expsteps);
      int (*getnumaccsteps)(SUNTimestepHeuristics H, long int* accsteps);
      int (*space)(SUNTimestepHeuristics H, long int *lenrw, long int *leniw);
    };

.. _SUNTimestepHeuristics.Description.heuristicsIDs:

SUNTimestepHeuristics Type Definitions
--------------------------------------

The time integrators in SUNDIALS leverage a myriad of heuristics to achieve
accurate and efficient computations; however, most follow a rather standard
pattern, or request no heuristic control whatsoever.  Therefore, each
SUNTimestepHeuristics implementation should indicate its type, informing the integrators
of whether they perform heuristic control or not:

.. c:enum:: SUNTimestepHeuristics_ID

   The enumerated type :c:type:`SUNTimestepHeuristics_ID` defines the enumeration
   constants for SUNDIALS heuristics implementations

.. c:enumerator:: SUN_TIMESTEPHEURISTICS_STD

   Performs "standard" heuristic stepsize constraints.

.. c:enumerator:: SUN_TIMESTEPHEURISTICS_NULL

   Performs no heuristic constraints.



Additionally, if the heuristics object can be used to limit stepsizes based on
explicit stability, then we define the following user-supplied function type:

.. c:type:: int (*SUNExpStabFn)(realtype *hstab, void *user_data)

   This function predicts a maximum stable step size.

   :param hstab: (output) the absolute value of the maximum stable step size.
   :param user_data: a pointer to user data for evaluation.
   :return: a *SUNExpStabFn* function should return 0 if it is successful, and a
            non-zero value otherwise.




.. _SUNTimestepHeuristics.Description.operations:

SUNTimestepHeuristics Operations
--------------------------------

The base SUNTimestepHeuristics class defines and implements all SUNTimestepHeuristics functions.
Most of these routines are merely wrappers for the operations defined by a
particular SUNTimestepHeuristics implementation, which are accessed through the *ops*
field of the ``SUNTimestepHeuristics`` structure.  However, the base SUNTimestepHeuristics class
also provides the convenience routine

.. c:function:: SUNTimestepHeuristics SUNTimestepHeuristics_NewEmpty(SUNContext sunctx)

  This function allocates a new generic ``SUNTimestepHeuristics`` object and initializes
  its content pointer and the function pointers in the operations structure to
  ``NULL``.

  :param sunctx: the :c:type:`SUNContext` object (see
                 :numref:`SUNDIALS.SUNContext`)
  :returns: If successful, a generic :c:type:`SUNTimestepHeuristics` object.  If
            unsuccessful, a ``NULL`` pointer will be returned.


Each of the following routines are *optional* for any specific SUNTimestepHeuristics implementation.


.. c:function:: SUNTimestepHeuristics_ID SUNTimestepHeuristics_GetID(SUNTimestepHeuristics H)

   Returns the type identifier for the heuristics object *H*.  Returned values
   are given in Section :numref:`SUNTimestepHeuristics.Description.heuristicsIDs`

   :param H: the :c:type:`SUNTimestepHeuristics` object.
   :return: :c:type:`SUNTimestepHeuristics_ID` type identifier.

   Usage:

   .. code-block:: c

      SUNTimestepHeuristics_ID id = SUNTimestepHeuristics_GetID(H);

.. c:function:: int SUNTimestepHeuristics_Destroy(SUNTimestepHeuristics H)

   Deallocates the heuristics object *H*.  If this is not provided by the
   implementation, the base wrapper routine will free both the *content* and
   *ops* objects -- this should be sufficient unless an heuristic implementation
   performs dynamic memory allocation of its own (note that the
   SUNDIALS-provided SUNTimestepHeuristics implementations do not need to supply this
   routine).

   :param H: the :c:type:`SUNTimestepHeuristics` object.
   :return: error code indicating success failure (see
            :numref:`SUNTimestepHeuristics.Description.errorCodes`).

   Usage:

   .. code-block:: c

      retval = SUNTimestepHeuristics_Destroy(H);

.. c:function:: int SUNTimestepHeuristics_ConstrainStep(SUNTimestepHeuristics H, realtype hcur, realtype hnew, realtype* hconstr)

   Main constraint-application function.  This will attempt to change the step
   *hcur* to *hnew*, applying any heuristic bounds on the step size adjustments.

   :param H: the :c:type:`SUNTimestepHeuristics` object.
   :param H: the heuristics object.
   :param hcur: the step size from the previous step attempt.
   :param hnew: the requested step size for the upcoming step attempt.
   :param hconstr: (output) pointer to the constrained step size.
   :return: error code indicating success failure (see
            :numref:`SUNTimestepHeuristics.Description.errorCodes`).

   Usage:

   .. code-block:: c

      retval = SUNTimestepHeuristics_ConstrainStep(H, hcur, hnew, &hnew);

.. c:function:: int SUNTimestepHeuristics_ETestFail(SUNTimestepHeuristics H, realtype hcur, realtype hnew, int nef, realtype* hconstr)

   Function to apply constraints following a step with unacceptable temporal
   error.

   :param H: the heuristics object.
   :param hcur: the step size that led to the error test failure.
   :param hnew: the requested step size for the upcoming step attempt (e.g.,
                from a :c:type:`SUNAdaptController` object).
   :param nef: the integrator-provided counter of how many temporal error test
               failures have occurred on this time step.
   :param hconstr: (output) pointer to the constrained step size.
   :return: error code indicating success failure (see
            :numref:`SUNTimestepHeuristics.Description.errorCodes`).

   Usage:

   .. code-block:: c

      retval = SUNTimestepHeuristics_ETestFail(H, hcur, hnew, nef, &hnew);

.. c:function:: int SUNTimestepHeuristics_BoundReduction(SUNTimestepHeuristics H, realtype hcur, realtype hnew, realtype *hconstr)

   This ensures that a step size reduction is within user-prescribed bounds.

   :param H: the heuristics object.
   :param hcur: the step size from the previous step attempt.
   :param hnew: the requested step size for the upcoming step attempt (e.g.,
                from a :c:type:`SUNAdaptController` object).
   :param hconstr: (output) pointer to the constrained step size.
   :return:
      * *SUNTIMESTEPHEURISTICS_SUCCESS* if successful
      * *SUNTIMESTEPHEURISTICS_CANNOT_DECREASE* if a reduction is requested but no
        reduction is possible

   Usage:

   .. code-block:: c

      retval = SUNTimestepHeuristics_BoundReduction(H, hcur, hnew, &hnew);

.. c:function:: int SUNTimestepHeuristics_BoundFirstStep(SUNTimestepHeuristics H, realtype h0, realtype *h0constr)

   This bounds the initial step by user-provided min/max step values.

   :param H: the heuristics object.
   :param h0: the requested initial step size.
   :param h0constr: (output) pointer to the constrained initial step size.
   :return: error code indicating success failure (see
            :numref:`SUNTimestepHeuristics.Description.errorCodes`).

   Usage:

   .. code-block:: c

      retval = SUNTimestepHeuristics_BoundFirstStep(H, h0, &h0);

.. c:function:: int SUNTimestepHeuristics_ConvFail(SUNTimestepHeuristics H, realtype hcur, realtype *hconstr)

   Function to apply constraints following a step with an algebraic solver
   convergence failure.

   :param H: the heuristics object.
   :param hcur: the step size that led to the convergence failure.
   :param hconstr: (output) pointer to the constrained step size.
   :return: error code indicating success failure (see
            :numref:`SUNTimestepHeuristics.Description.errorCodes`).

   Usage:

   .. code-block:: c

      retval = SUNTimestepHeuristics_ConvFail(H, hcur, &hnew);

.. c:function::int SUNTimestepHeuristics_Reset(SUNTimestepHeuristics H)

   Function to reset the heuristics to its initial state: zeros any internal
   counters, and resets any stepsize growth factor bounds.

   :param H: the heuristics object.
   :return: error code indicating success failure (see
            :numref:`SUNTimestepHeuristics.Description.errorCodes`).

   Usage:

   .. code-block:: c

      retval = SUNTimestepHeuristics_Reset(H);

.. c:function::int SUNTimestepHeuristics_Update(SUNTimestepHeuristics H)

   Function to notify the heuristics object that a time step has succeeded,
   indicating e.g. that the stepsize growh factor should should be set to its
   "default" state.

   :param H: the heuristics object.
   :return: error code indicating success failure (see
            :numref:`SUNTimestepHeuristics.Description.errorCodes`).

   Usage:

   .. code-block:: c

      retval = SUNTimestepHeuristics_Update(H);

.. c:function::int SUNTimestepHeuristics_SetDefaults(SUNTimestepHeuristics H)

   Function to set the heuristics parameters to their default values.

   :param H: the heuristics object.
   :return: error code indicating success failure (see
            :numref:`SUNTimestepHeuristics.Description.errorCodes`).

   Usage:

   .. code-block:: c

      retval = SUNTimestepHeuristics_SetDefaults(H);

.. c:function::int SUNTimestepHeuristics_Write(SUNTimestepHeuristics H, FILE* fptr)

   Writes all controller parameters to the indicated file pointer.

   :param H: the heuristics object.
   :param fptr: the output stream to write the parameters.
   :return: error code indicating success failure (see
            :numref:`SUNTimestepHeuristics.Description.errorCodes`).

   Usage:

   .. code-block:: c

      retval = SUNTimestepHeuristics_Write(H, stdout);

.. c:function::int SUNTimestepHeuristics_SetMaxStep(SUNTimestepHeuristics H, realtype hmax)

   Function to inform the heuristics object about a maximum allowed absolute
   step size.

   :param H: the heuristics object.
   :param hmax: maximum absolute step size allowed (:math:`\text{hmax} \le 0`
                implies :math:`\text{hmax}=\infty`).
   :return: error code indicating success failure (see
            :numref:`SUNTimestepHeuristics.Description.errorCodes`).

   Usage:

   .. code-block:: c

      retval = SUNTimestepHeuristics_SetMaxStep(H, 1.0e-2);

.. c:function::int SUNTimestepHeuristics_SetMinStep(SUNTimestepHeuristics H, realtype hmin)

   Function to inform the heuristics object about a minimum allowed absolute
   step size.

   :param H: the heuristics object.
   :param hmin: minimum absolute step size allowed (:math:`\text{hmin} \le 0`
                implies no minimum).
   :return: error code indicating success failure (see
            :numref:`SUNTimestepHeuristics.Description.errorCodes`).

   Usage:

   .. code-block:: c

      retval = SUNTimestepHeuristics_SetMinStep(H, 1.0e-5);

.. c:function::int SUNTimestepHeuristics_SetExpStabFn(SUNTimestepHeuristics H, SUNExpStabFn EStab, void* estab_data)

   Function to provide a user-supplied function for the maximum stable step
   size.

   :param H: the heuristics object.
   :param EStab: user-supplied function specifying the maximum stable step size
                 (``EStab == NULL`` disables).
   :param estab_data: user-supplied data pointer that should be provided on all
                      calls to *EStab*.
   :return: error code indicating success failure (see
            :numref:`SUNTimestepHeuristics.Description.errorCodes`).

   Usage:

   .. code-block:: c

      retval = SUNTimestepHeuristics_SetExpStabFn(H, CFLStabilityBound, (void*) (&mystruct));

.. c:function::int SUNTimestepHeuristics_SetCFLFraction(SUNTimestepHeuristics H, realtype cfl_frac)

   Function to set a CFL safety factor that should be applied to the outputs of
   *EStab*, above.

   :param H: the heuristics object.
   :param cfl_frac: safety factor.  Allowable values must be within the open
                    interval (0,1), :math:`\text{cfl\_frac} \le 0` implies a
                    reset to the default value.
   :return: error code indicating success failure (see
            :numref:`SUNTimestepHeuristics.Description.errorCodes`).

   Usage:

   .. code-block:: c

      retval = SUNTimestepHeuristics_SetCFLFraction(H, 0.5);

.. c:function::int SUNTimestepHeuristics_SetSafetyFactor(SUNTimestepHeuristics C, realtype safety)

   Function to set a step size safety factor that should be used to constrain an
   error-controller-recommended step size.

   :param H: the heuristics object.
   :param safety: safety factor.  Allowable values must be within the open
                  interval (0,1), :math:`\text{safety} \le 0` implies a reset
                  to the default value.
   :return: error code indicating success failure (see
            :numref:`SUNTimestepHeuristics.Description.errorCodes`).

   Usage:

   .. code-block:: c

      retval = SUNTimestepHeuristics_SetSafetyFactor(H, 0.95);

.. c:function::int SUNTimestepHeuristics_SetMaxGrowth(SUNTimestepHeuristics H, realtype mx_growth)

   Function to set maximum stepsize growth factor for general steps.

   :param H: the heuristics object.
   :param mx_growth: maximum step size growth factor.  Allowable values must be
                     greater than 1; illegal values imply a reset to the
                     default.
   :return: error code indicating success failure (see
            :numref:`SUNTimestepHeuristics.Description.errorCodes`).

   **Notes:**
      The first step uses a separate maximum growth factor.

   .. code-block:: c

      retval = SUNTimestepHeuristics_SetMaxGrowth(H, 20.0);

.. c:function::int SUNTimestepHeuristics_SetMaxFirstGrowth(SUNTimestepHeuristics H, realtype etamx1)

   Function to set maximum stepsize growth factor for the first internal time
   step.

   :param H: the heuristics object.
   :param etamx1: maximum first step size growth factor.  Allowable values must
                  be greater than 1; illegal values imply a reset to the
                  default.
   :return: error code indicating success failure (see
            :numref:`SUNTimestepHeuristics.Description.errorCodes`).

   Usage:

   .. code-block:: c

      retval = SUNTimestepHeuristics_SetMaxFirstGrowth(H, 10000.0);

.. c:function::int SUNTimestepHeuristics_SetFixedStepBounds(SUNTimestepHeuristics H, realtype lb, realtype ub)

   Function to specify the step size growth interval within which the step size
   will remain unchanged.  Allowable values must enclose the value 1.0.  Any
   illegal interval implies a reset to the default.

   :param H: the heuristics object.
   :param lb: interval lower bound.
   :param ub: interval upper bound.
   :return: error code indicating success failure (see
            :numref:`SUNTimestepHeuristics.Description.errorCodes`).

   Usage:

   .. code-block:: c

      retval = SUNTimestepHeuristics_SetFixedStepBounds(H, 1.0, 1.5);

.. c:function::int SUNTimestepHeuristics_SetMinReduction(SUNTimestepHeuristics H, realtype eta_min)

   Function to set a lower bound for the stepsize adjustment factor following a
   temporal error failure.

   :param H: the heuristics object.
   :param eta_min: minimum step size reduction factor.  Allowable values must be
                   in the open interval (0,1); illegal values imply a reset to
                   the default.
   :return: error code indicating success failure (see
            :numref:`SUNTimestepHeuristics.Description.errorCodes`).

   Usage:

   .. code-block:: c

      retval = SUNTimestepHeuristics_SetMinReduction(H, 1.0e-2);

.. c:function::int SUNTimestepHeuristics_SetMaxEFailGrowth(SUNTimestepHeuristics H, realtype etamxf)

   Function to set an upper bound for the stepsize adjustment factor following a
   temporal error failure (once at least *small_nef* failures have occurred, see
   :c:func:`SUNTimestepHeuristics_SetSmallNumEFails`).

   :param H: the heuristics object.
   :param etamxf: temporal error failure step size reduction factor.  Allowable
                  values must be in the interval (0,1]; illegal values imply a
                  reset to the default.
   :return: error code indicating success failure (see
            :numref:`SUNTimestepHeuristics.Description.errorCodes`).

   Usage:

   .. code-block:: c

      retval = SUNTimestepHeuristics_SetMaxEFailGrowth(H, 0.8);

.. c:function::int SUNTimestepHeuristics_SetSmallNumEFails(SUNTimestepHeuristics H, int small_nef)

   Function to specify the step adaptivity constant *small_nef* (see
   :c:func:`SUNTimestepHeuristics_SetMaxEFailGrowth`).

   :param H: the heuristics object.
   :param small_nef: minimum number of error test failures before *etamxf* takes
                     effect. Allowable values are greater than 0. Illegal
                     values imply a reset to the default.
   :return: error code indicating success failure (see
            :numref:`SUNTimestepHeuristics.Description.errorCodes`).

   Usage:

   .. code-block:: c

      retval = SUNTimestepHeuristics_SetSmallNumEFails(H, 3);

.. c:function::int SUNTimestepHeuristics_SetMaxCFailGrowth(SUNTimestepHeuristics H, realtype etacf)

   Function to specify an upper bound for the stepsize adjustment factor
   following an algebraic solver convergence failure.

   :param H: the heuristics object.
   :param etacf: convergence failure step adjustment factor. Allowable values
                 must be in the interval (0,1]; illegal values imply a reset to
                 the default.
   :return: error code indicating success failure (see
            :numref:`SUNTimestepHeuristics.Description.errorCodes`).

   Usage:

   .. code-block:: c

      retval = SUNTimestepHeuristics_SetMaxCFailGrowth(H, 0.3);

.. c:function::int SUNTimestepHeuristics_GetNumExpSteps(SUNTimestepHeuristics H, long int* expsteps)

   Function to return the current number of stability-limited steps.

   :param H: the heuristics object.
   :param expsteps: (output) pointer to hold the stability-limited step counter.
   :return: error code indicating success failure (see
            :numref:`SUNTimestepHeuristics.Description.errorCodes`).

   Usage:

   .. code-block:: c

      retval = SUNTimestepHeuristics_GetNumExpSteps(H, &expsteps);

.. c:function::int SUNTimestepHeuristics_GetNumAccSteps(SUNTimestepHeuristics H, long int* accsteps)

   Function to return the current number of accuracy-limited steps.

   :param H: the heuristics object.
   :param expsteps: (output) pointer to hold the accuracy-limited step counter.
   :return: error code indicating success failure (see
            :numref:`SUNTimestepHeuristics.Description.errorCodes`).

   Usage:

   .. code-block:: c

      retval = SUNTimestepHeuristics_GetNumAccSteps(H, &accsteps);

.. c:function:: int SUNTimestepHeuristics_Space(SUNTimestepHeuristics H, long int *lenrw, long int *leniw)

   Informative routine that returns the memory requirements of the heuristics object.

   :param H: the heuristics object.
   :param lenrw: (output) number of ``sunrealtype`` words stored in the
                 heuristics object.
   :param leniw: (output) number of ``sunindextype`` words stored in the
                 heuristics object.  This may also include pointers, `int` and
                 `long int` words.
   :return: error code indicating success failure (see
            :numref:`SUNTimestepHeuristics.Description.errorCodes`).

   Usage:

   .. code-block:: c

      retval = SUNTimestepHeuristics_Space(H, &lenrw, &leniw);



.. _SUNTimestepHeuristics.Description.errorCodes:

SUNTimestepHeuristics Error Codes
-------------------------

SUNTimestepHeuristics functions return one of the following set of error codes:

* ``SUNTIMESTEPHEURISTICS_SUCCESS`` (0) -- successful call.

* ``SUNTIMESTEPHEURISTICS_ILL_INPUT`` (-1101) -- an illegal input has been provided to
  the function.

* ``SUNTIMESTEPHEURISTICS_CANNOT_DECREASE`` (-1102) -- a step reduction was required but
  was not possible (already at minimum).

* ``SUNTIMESTEPHEURISTICS_USER_FCN_FAIL`` (-1103) -- a user-supplied function returned a
  nonzero [error] value.


C/C++ API Usage
---------------

The SUNDIALS Heuristics module can be used in C and C++ programs by including
the header file ``sundials/sundials_timestepheuristics.h``.

Example usage (here ``SUNTimestepHeuristicsXYZ`` is a placeholder for an actual
SUNTimestepHeuristics implementation constructor):

.. code-block:: c

    #include <stdio.h>
    #include <stdlib.h>
    #include <sundials/sundials_context.h>
    #include <sundials/sundials_types.h>
    #include <sundials/sundials_timestepheuristics.h>

    int main()
    {
        /* Create a SUNContext object */
        SUNContext sunctx = ...;

        /* Create a SUNTimestepHeuristics object */
        SUNTimestepHeuristics H = SUNTimestepHeuristicsXYZ(sunctx);

        /* Use the heuristics object */

        /* Destroy the heuristics object */
        retval = SUNTimestepHeuristics_Destroy(C);

        return 0;
    }
