.. ----------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2002-2023, Lawrence Livermore National Security
   and Southern Methodist University.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   ----------------------------------------------------------------

.. _ARKODE.Usage.SPRKStep.UserCallable:

SPRKStep User-callable functions
==================================

This section describes the functions that are called by the
user to setup and then solve an IVP using the SPRKStep time-stepping
module. Some of these are required; however, starting with
:numref:`ARKODE.Usage.SPRKStep.OptionalInputs`, the functions listed involve
optional inputs/outputs or restarting, and those paragraphs may be
skipped for a casual use of ARKODE's SPRKStep module. In any case,
refer to the preceding section, :numref:`ARKODE.Usage.SPRKStep.Skeleton`,
for the correct order of these calls.

On an error, each user-callable function returns a negative value  (or
``NULL`` if the function returns a pointer) and sends an error message
to the error handler routine, which prints the message to ``stderr``
by default. However, the user can set a file as error output or can
provide their own error handler function (see
:numref:`ARKODE.Usage.SPRKStep.OptionalInputs` for details).



.. _ARKODE.Usage.SPRKStep.Initialization:

SPRKStep initialization and deallocation functions
------------------------------------------------------


.. c:function:: void* SPRKStepCreate(ARKRhsFn f1, ARKRhsFn f2, sunrealtype t0,\
                                     N_Vector y0, SUNContext sunctx)

   This function allocates and initializes memory for a problem to
   be solved using the SPRKStep time-stepping module in ARKODE.

   :param f1: the name of the C function (of type :c:func:`ARKRhsFn()`) defining :math:`f_1(t,q) = \frac{\partial V(t,q)}{\partial q}`
   :param f2: the name of the C function (of type :c:func:`ARKRhsFn()`) defining :math:`f_2(t,p) = \frac{\partial T(t,p)}{\partial p}`
   :param t0: the initial value of :math:`t`
   :param y0: the initial condition vector :math:`y(t_0)`
   :param sunctx: the :c:type:`SUNContext` object (see :numref:`SUNDIALS.SUNContext`)

   :returns: If successful, a pointer to initialized problem memory of type
             ``void*``, to be passed to all user-facing SPRKStep routines listed
             below.  If unsuccessful, a ``NULL`` pointer will be returned, and
             an error message will be printed to ``stderr``.

   .. warning::

      SPRKStep requires a partitioned problem where ``f1`` should only modify
      the q variables and ``f2`` should only modify the p variables (or vice
      versa). However, the vector passed to these functions is the full vector
      with both p and q. The ordering of the variables is determined implicitly
      by the user when they set the initial conditions.


.. c:function:: void SPRKStepFree(void** arkode_mem)

   This function frees the problem memory *arkode_mem* created by
   :c:func:`SPRKStepCreate`.

   :param arkode_mem: pointer to the SPRKStep memory block.


.. _ARKODE.Usage.SPRKStep.RootFinding:

Rootfinding initialization function
--------------------------------------

As described in :numref:`ARKODE.Mathematics.Rootfinding`, while
solving the IVP, ARKODE's time-stepping modules have the capability to
find the roots of a set of user-defined functions.  To activate the
root-finding algorithm, call the following function.  This is normally
called only once, prior to the first call to
:c:func:`SPRKStepEvolve()`, but if the rootfinding problem is to be
changed during the solution, :c:func:`SPRKStepRootInit()` can also be
called prior to a continuation call to :c:func:`SPRKStepEvolve()`.

.. note::

   The solution is interpolated to the times at which roots are found.


.. c:function:: int SPRKStepRootInit(void* arkode_mem, int nrtfn, ARKRootFn g)

   Initializes a rootfinding problem to be solved during the
   integration of the ODE system. It *must* be called after
   :c:func:`SPRKStepCreate`, and before :c:func:`SPRKStepEvolve()`.

   To disable the rootfinding feature after it has already been initialized, or
   to free memory associated with SPRKStep's rootfinding module, call
   :c:func:`SPRKStepRootInit` with *nrtfn = 0*.

   Similarly, if a new IVP is to be solved with a call to
   :c:func:`SPRKStepReInit()`, where the new IVP has no rootfinding problem but
   the prior one did, then call :c:func:`SPRKStepRootInit` with *nrtfn = 0*.

   :param arkode_mem: pointer to the SPRKStep memory block.
   :param nrtfn: number of functions :math:`g_i`, an integer :math:`\ge` 0.
   :param g: name of user-supplied function, of type :c:func:`ARKRootFn()`,
      defining the functions :math:`g_i` whose roots are sought.

   :retval ARK_SUCCESS: if successful
   :retval ARK_MEM_NULL: if the SPRKStep memory was ``NULL``
   :retval ARK_MEM_FAIL: if there was a memory allocation failure
   :retval ARK_ILL_INPUT: if *nrtfn* is greater than zero but *g* = ``NULL``.


.. _ARKODE.Usage.SPRKStep.Integration:

SPRKStep solver function
-------------------------

This is the central step in the solution process -- the call to perform
the integration of the IVP.  One of the input arguments (*itask*)
specifies one of two modes as to where SPRKStep is to return a
solution.  These modes are modified if the user has set a stop time
(with a call to the optional input function :c:func:`SPRKStepSetStopTime()`) or
has requested rootfinding.



.. c:function:: int SPRKStepEvolve(void* arkode_mem, sunrealtype tout, N_Vector yout, sunrealtype *tret, int itask)

   Integrates the ODE over an interval in :math:`t`.

   :param arkode_mem: pointer to the SPRKStep memory block.
   :param tout: the next time at which a computed solution is desired.
   :param yout: the computed solution vector.
   :param tret: the time corresponding to *yout* (output).
   :param itask: a flag indicating the job of the solver for the next user step.

                 The *ARK_NORMAL* option causes the solver to take internal
                 steps until it has just overtaken a user-specified output
                 time, *tout*, in the direction of integration,
                 i.e. :math:`t_{n-1} <` *tout* :math:`\le t_{n}` for forward
                 integration, or :math:`t_{n} \le` *tout* :math:`< t_{n-1}` for
                 backward integration.  It will then compute an approximation
                 to the solution :math:`y(tout)` by interpolation (using one
                 of the dense output routines described in
                 :numref:`ARKODE.Mathematics.Interpolation`).

                 The *ARK_ONE_STEP* option tells the solver to only take a
                 single internal step, :math:`y_{n-1} \to y_{n}`, and return the
                 solution at that point, :math:`y_{n}`, in the vector *yout*.

   :retval ARK_SUCCESS: if successful.
   :retval ARK_ROOT_RETURN: if :c:func:`SPRKStepEvolve()` succeeded, and
                            found one or more roots.  If the number of root functions,
                            *nrtfn*, is greater than 1, call
                            :c:func:`SPRKStepGetRootInfo()` to see which :math:`g_i` were
                            found to have a root at (*\*tret*).
   :retval ARK_TSTOP_RETURN: if :c:func:`SPRKStepEvolve()` succeeded and
                             returned at *tstop*.
   :retval ARK_MEM_NULL: if the *arkode_mem* argument was ``NULL``.
   :retval ARK_NO_MALLOC: if *arkode_mem* was not allocated.
   :retval ARK_ILL_INPUT: if one of the inputs to
                          :c:func:`SPRKStepEvolve()` is illegal, or some other
                          input to the solver was either illegal or missing.
                          Details will be provided in the error message. Typical
                          causes of this failure are a root of one of the root
                          functions was found both at a point :math:`t` and also
                          very near :math:`t`.
   :retval ARK_TOO_MUCH_WORK: if the solver took *mxstep* internal steps but
                              could not reach *tout*. The default value for
                              *mxstep* is *MXSTEP_DEFAULT = 500*.
   :retval ARK_ERR_FAILURE: if error test failures occurred either too many
                            times (*ark_maxnef*) during one internal time step
                            or occurred with :math:`|h| = h_{min}`.
   :retval ARK_VECTOROP_ERR: a vector operation error occurred.

   .. note::

      The input vector *yout* can use the same memory as the
      vector *y0* of initial conditions that was passed to
      :c:func:`SPRKStepCreate`.

      In *ARK_ONE_STEP* mode, *tout* is used only on the first call, and
      only to get the direction and a rough scale of the independent
      variable.

      All failure return values are negative and so testing the
      return argument for negative values will trap all
      :c:func:`SPRKStepEvolve()` failures.

      Since interpolation may reduce the accuracy in the reported
      solution, if full method accuracy is desired the user should issue
      a call to :c:func:`SPRKStepSetStopTime()` before the call to
      :c:func:`SPRKStepEvolve()` to specify a fixed stop time to
      end the time step and return to the user.  Upon return from
      :c:func:`SPRKStepEvolve()`, a copy of the internal solution
      :math:`y_{n}` will be returned in the vector *yout*.  Once the
      integrator returns at a *tstop* time, any future testing for
      *tstop* is disabled (and can be re-enabled only though a new call
      to :c:func:`SPRKStepSetStopTime()`). Interpolated outputs may or may not
      conserve the Hamiltonian. Our testing has shown that Lagrange
      interpolation typically performs well in this regard, while Hermite
      interpolation does not. As such, SPRKStep uses the Lagrange interpolation
      module by default.

      On any error return in which one or more internal steps were taken
      by :c:func:`SPRKStepEvolve()`, the returned values of *tret* and
      *yout* correspond to the farthest point reached in the integration.
      On all other error returns, *tret* and *yout* are left unchanged
      from those provided to the routine.




.. _ARKODE.Usage.SPRKStep.OptionalInputs:

Optional input functions
-------------------------

There are numerous optional input parameters that control the behavior
of SPRKStep, each of which may be modified from its default value through
calling an appropriate input function.  The following tables list all
optional input functions, grouped by which aspect of SPRKStep they control.
Detailed information on the calling syntax and arguments for each
function are then provided following each table.

The optional inputs are grouped into the following categories:

* General SPRKStep options (:numref:`ARKODE.Usage.SPRKStep.SPRKStepInputTable`),

* IVP method solver options (:numref:`ARKODE.Usage.SPRKStep.SPRKStepMethodInputTable`),

* Rootfinding options (:numref:`ARKODE.Usage.SPRKStep.SPRKStepRootfindingInputTable`).

For the most casual use of SPRKStep, relying on the default set of
solver parameters, the reader can skip to section on user-supplied
functions, :numref:`ARKODE.Usage.UserSupplied`.

We note that, on an error return, all of the optional input functions send an
error message to the error handler function. All error return values are
negative, so a test on the return arguments for negative values will catch all
errors. Finally, a call to a ``SPRKStepSet***`` function can generally be made
from the user's calling program at any time and, if successful, takes effect
immediately. For ``SPRKStepSet***`` functions that cannot be called at any time,
this is explicitly noted in the function documentation.



.. _ARKODE.Usage.SPRKStep.SPRKStepInput:

Optional inputs for SPRKStep
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. _ARKODE.Usage.SPRKStep.SPRKStepInputTable:
.. table:: Optional inputs for SPRKStep

  +-----------------------------------------------------+------------------------------------------+------------------------+
  | Optional input                                      | Function name                            | Default                |
  +=====================================================+==========================================+========================+
  | Return SPRKStep solver parameters to their defaults | :c:func:`SPRKStepSetDefaults()`          | internal               |
  +-----------------------------------------------------+------------------------------------------+------------------------+
  | Set dense output interpolation type                 | :c:func:`SPRKStepSetInterpolantType()`   | ``ARK_INTERP_LAGRANGE``|
  +-----------------------------------------------------+------------------------------------------+------------------------+
  | Set dense output polynomial degree                  | :c:func:`SPRKStepSetInterpolantDegree()` | 5                      |
  +-----------------------------------------------------+------------------------------------------+------------------------+
  | Supply a pointer to a diagnostics output file       | :c:func:`SPRKStepSetDiagnostics()`       | ``NULL``               |
  +-----------------------------------------------------+------------------------------------------+------------------------+
  | Supply a pointer to an error output file            | :c:func:`SPRKStepSetErrFile()`           | ``stderr``             |
  +-----------------------------------------------------+------------------------------------------+------------------------+
  | Supply a custom error handler function              | :c:func:`SPRKStepSetErrHandlerFn()`      | internal fn            |
  +-----------------------------------------------------+------------------------------------------+------------------------+
  | Set fixed step size (required user input)           | :c:func:`SPRKStepSetFixedStep()`         | user defined           |
  +-----------------------------------------------------+------------------------------------------+------------------------+
  | Maximum no. of internal steps before *tout*         | :c:func:`SPRKStepSetMaxNumSteps()`       | 500                    |
  +-----------------------------------------------------+------------------------------------------+------------------------+
  | Set a value for :math:`t_{stop}`                    | :c:func:`SPRKStepSetStopTime()`          | undefined              |
  +-----------------------------------------------------+------------------------------------------+------------------------+
  | Disable the stop time                               | :c:func:`SPRKStepClearStopTime`          | N/A                    |
  +-----------------------------------------------------+------------------------------------------+------------------------+
  | Supply a pointer for user data                      | :c:func:`SPRKStepSetUserData()`          | ``NULL``               |
  +-----------------------------------------------------+------------------------------------------+------------------------+


.. c:function:: int SPRKStepSetDefaults(void* arkode_mem)

   Resets all optional input parameters to SPRKStep's original
   default values.

   :param arkode_mem: pointer to the SPRKStep memory block.

   :retval ARK_SUCCESS: if successful
   :retval ARK_MEM_NULL: if the SPRKStep memory is ``NULL``
   :retval ARK_ILL_INPUT: if an argument has an illegal value

   .. note::

      Does not change problem-defining function pointer *f*
      or the *user_data* pointer.

      Also leaves alone any data structures or options related to
      root-finding (those can be reset using :c:func:`SPRKStepRootInit()`).


.. c:function:: int SPRKStepSetInterpolantType(void* arkode_mem, int itype)

   Specifies use of the Lagrange or Hermite interpolation modules (used for
   dense output -- interpolation of solution output values and implicit
   method predictors).

   :param arkode_mem: pointer to the SPRKStep memory block.
   :param itype: requested interpolant type (``ARK_INTERP_HERMITE`` or ``ARK_INTERP_LAGRANGE``)

   :retval ARK_SUCCESS: if successful
   :retval ARK_MEM_NULL: if the SPRKStep memory is ``NULL``
   :retval ARK_MEM_FAIL: if the interpolation module cannot be allocated
   :retval ARK_ILL_INPUT: if the *itype* argument is not recognized or the
                          interpolation module has already been initialized

   .. note::

      The Hermite interpolation module is described in
      :numref:`ARKODE.Mathematics.Interpolation.Hermite`, and the Lagrange interpolation module
      is described in :numref:`ARKODE.Mathematics.Interpolation.Lagrange`.

      This routine frees any previously-allocated interpolation module, and re-creates
      one according to the specified argument.  Thus any previous calls to
      :c:func:`SPRKStepSetInterpolantDegree()` will be nullified.

      This routine must be called *after* the call to :c:func:`SPRKStepCreate`.
      After the first call to :c:func:`SPRKStepEvolve()` the interpolation type may
      not be changed without first calling :c:func:`SPRKStepReInit()`.

      If this routine is not called, the Lagrange interpolation module will be used.

      Interpolated outputs may or may not conserve the Hamiltonian. Our testing
      has shown that Lagrange interpolation typically performs well in this
      regard, while Hermite interpolation does not.



.. c:function:: int SPRKStepSetInterpolantDegree(void* arkode_mem, int degree)

   Specifies the degree of the polynomial interpolant
   used for dense output (i.e. interpolation of solution output values).
   Allowed values are between 0 and 5.

   :param arkode_mem: pointer to the SPRKStep memory block.
   :param degree: requested polynomial degree.

   :retval ARK_SUCCESS: if successful
   :retval ARK_MEM_NULL: if the SPRKStep memory or interpolation module are ``NULL``
   :retval ARK_INTERP_FAIL: if this is called after :c:func:`SPRKStepEvolve()`
   :retval ARK_ILL_INPUT: if an argument has an illegal value or the
                          interpolation module has already been initialized

   .. note::

      This routine should be called *after* :c:func:`SPRKStepCreate` and *before*
      :c:func:`SPRKStepEvolve()`. After the first call to :c:func:`SPRKStepEvolve()`
      the interpolation degree may not be changed without first calling
      :c:func:`SPRKStepReInit()`.

      If a user calls both this routine and :c:func:`SPRKStepSetInterpolantType()`, then
      :c:func:`SPRKStepSetInterpolantType()` must be called first.

      Since the accuracy of any polynomial interpolant is limited by the
      accuracy of the time-step solutions on which it is based, the *actual*
      polynomial degree that is used by SPRKStep will be the minimum of
      :math:`q-1` and the input *degree*, for :math:`q > 1` where :math:`q` is
      the order of accuracy for the time integration method.

      When `q = 1`, a linear interpolant is the default to ensure values
      obtained by the integrator are returned at the ends of the time interval.


.. c:function:: int SPRKStepSetErrFile(void* arkode_mem, FILE* errfp)

   Specifies a pointer to the file where all SPRKStep warning and error
   messages will be written if the default internal error handling
   function is used.

   The default value for *errfp* is ``stderr``.

   Passing a ``NULL`` value disables all future error message output
   (except for the case wherein the SPRKStep memory pointer is
   ``NULL``).  This use of the function is strongly discouraged.

   :param arkode_mem: pointer to the SPRKStep memory block.
   :param errfp: pointer to the output file.

   :retval ARK_SUCCESS: if successful
   :retval ARK_MEM_NULL: if the SPRKStep memory is ``NULL``
   :retval ARK_ILL_INPUT: if an argument has an illegal value

   .. note::

      If used, this routine should be called before any other
      optional input functions, in order to take effect for subsequent
      error messages.


.. c:function:: int SPRKStepSetErrHandlerFn(void* arkode_mem, ARKErrHandlerFn ehfun, void* eh_data)

   Specifies the optional user-defined function to be used
   in handling error messages.

   :param arkode_mem: pointer to the SPRKStep memory block.
   :param ehfun: name of user-supplied error handler function.
   :param eh_data: pointer to user data passed to *ehfun* every time it is called.

   :retval ARK_SUCCESS: if successful
   :retval ARK_MEM_NULL: if the SPRKStep memory is ``NULL``
   :retval ARK_ILL_INPUT: if an argument has an illegal value

   .. note::

      Error messages indicating that the SPRKStep solver memory is
      ``NULL`` will always be directed to ``stderr``.


.. c:function:: int SPRKStepSetFixedStep(void* arkode_mem, sunrealtype hfixed)

   Sets the time step size used within SPRKStep.

   :param arkode_mem: pointer to the SPRKStep memory block.
   :param hfixed: value of the fixed step size to use.

   :retval ARK_SUCCESS: if successful
   :retval ARK_MEM_NULL: if the SPRKStep memory is ``NULL``
   :retval ARK_ILL_INPUT: if an argument has an illegal value


.. c:function:: int SPRKStepSetMaxNumSteps(void* arkode_mem, long int mxsteps)

   Specifies the maximum number of steps to be taken by the
   solver in its attempt to reach the next output time, before SPRKStep
   will return with an error.

   Passing *mxsteps* = 0 results in SPRKStep using the
   default value (500).

   Passing *mxsteps* < 0 disables the test (not recommended).

   :param arkode_mem: pointer to the SPRKStep memory block.
   :param mxsteps: maximum allowed number of internal steps.

   :retval ARK_SUCCESS: if successful
   :retval ARK_MEM_NULL: if the SPRKStep memory is ``NULL``
   :retval ARK_ILL_INPUT: if an argument has an illegal value


.. c:function:: int SPRKStepSetStopTime(void* arkode_mem, sunrealtype tstop)

   Specifies the value of the independent variable
   :math:`t` past which the solution is not to proceed.

   The default is that no stop time is imposed.

   Once the integrator returns at a stop time, any future testing for
   ``tstop`` is disabled (and can be reenabled only though a new call to
   :c:func:`SPRKStepSetStopTime`).

   A stop time not reached before a call to :c:func:`SPRKStepReInit` or
   :c:func:`SPRKStepReset` will remain active but can be disabled by calling
   :c:func:`SPRKStepClearStopTime`.

   :param arkode_mem: pointer to the SPRKStep memory block.
   :param tstop: stopping time for the integrator.

   :retval ARK_SUCCESS: if successful
   :retval ARK_MEM_NULL: if the SPRKStep memory is ``NULL``
   :retval ARK_ILL_INPUT: if an argument has an illegal value


.. c:function:: int SPRKStepClearStopTime(void* arkode_mem)

   Disables the stop time set with :c:func:`SPRKStepSetStopTime`.

   The stop time can be reenabled though a new call to
   :c:func:`SPRKStepSetStopTime`.

   :param arkode_mem: pointer to the SPRKStep memory block.

   :retval ARK_SUCCESS: if successful
   :retval ARK_MEM_NULL: if the SPRKStep memory is ``NULL``


.. c:function:: int SPRKStepSetUserData(void* arkode_mem, void* user_data)

   Specifies the user data block *user_data* and
   attaches it to the main SPRKStep memory block.

   If specified, the pointer to *user_data* is passed to all
   user-supplied functions for which it is an argument; otherwise
   ``NULL`` is passed.

   :param arkode_mem: pointer to the SPRKStep memory block.
   :param user_data: pointer to the user data.

   :retval ARK_SUCCESS: if successful
   :retval ARK_MEM_NULL: if the SPRKStep memory is ``NULL``
   :retval ARK_ILL_INPUT: if an argument has an illegal value


.. _ARKODE.Usage.SPRKStep.SPRKStepMethodInput:

Optional inputs for IVP method selection
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. _ARKODE.Usage.SPRKStep.SPRKStepMethodInputTable:
.. table:: Optional inputs for IVP method selection

   +-----------------------------+-------------------------------------------+-------------------------------------+
   |       Optional input        |               Function name               |               Default               |
   +=============================+===========================================+=====================================+
   | Set integrator method order | :c:func:`SPRKStepSetOrder()`              | 4                                   |
   +-----------------------------+-------------------------------------------+-------------------------------------+
   | Set SPRK method             | :c:func:`SPRKStepSetMethod()`             | :c:type:`ARKODE_SPRK_MCLACHLAN_4_4` |
   +-----------------------------+-------------------------------------------+-------------------------------------+
   | Set SPRK method by name     | :c:func:`SPRKStepSetMethodName()`         | "ARKODE_SPRK_MCLACHLAN_4_4"         |
   +-----------------------------+-------------------------------------------+-------------------------------------+
   | Use compensated summation   | :c:func:`SPRKStepSetUseCompensatedSums()` | false                               |
   +-----------------------------+-------------------------------------------+-------------------------------------+


.. c:function:: int SPRKStepSetOrder(void* arkode_mem, int ord)

   Specifies the order of accuracy for the SPRK integration method.

   The allowed values are :math:`1,2,3,4,5,6,8,10`. Any illegal input will
   result in the default value of 4.

   Since *ord* affects the memory requirements for the internal
   SPRKStep memory block, it cannot be changed after the first call to
   :c:func:`SPRKStepEvolve()`, unless :c:func:`SPRKStepReInit()` is called.

   :param arkode_mem: pointer to the SPRKStep memory block.
   :param ord: requested order of accuracy.

   :retval ARK_SUCCESS: if successful
   :retval ARK_MEM_NULL: if the SPRKStep memory is ``NULL``
   :retval ARK_ILL_INPUT: if an argument has an illegal value

   .. warning::

      This overrides any previously set method so it should not be used with
      :c:func:`SPRKStepSetMethod` or :c:func:`SPRKStepMethodByName`.


.. c:function:: int SPRKStepSetMethod(void* arkode_mem, ARKodeSPRKTable sprk_table)

   Specifies the SPRK method.

   :param arkode_mem: pointer to the SPRKStep memory block.
   :param sprk_table: the SPRK method table.

   :retval ARK_SUCCESS: if successful
   :retval ARK_MEM_NULL: if the SPRKStep memory is ``NULL``
   :retval ARK_ILL_INPUT: if an argument has an illegal value

   .. note::

      No error checking is performed on the coefficients contained in the
      table to ensure its declared order of accuracy.


.. c:function:: int SPRKStepSetMethodName(void* arkode_mem, const char* method)

   Specifies the SPRK method by its name.

   :param arkode_mem: pointer to the SPRKStep memory block.
   :param method: the SPRK method name.

   :retval ARK_SUCCESS: if successful
   :retval ARK_MEM_NULL: if the SPRKStep memory is ``NULL``
   :retval ARK_ILL_INPUT: if an argument has an illegal value


.. c:function:: int SPRKStepSetUseCompensatedSums(void* arkode_mem, sunbooleantype onoff)

   Specifies if :ref:`compensated summation (and the incremental form) <ARKODE.Mathematics.SPRKStep.Compensated>`
   should be used where applicable.

   This increases the computational cost by 2 extra vector operations per stage
   and an additional 5 per time step. It also requires one extra vector to be
   stored. However, it is signficantly more robust to roundoff error
   accumulation.

   :param arkode_mem: pointer to the SPRKStep memory block.
   :param onoff: should compensated summation be used (1) or not (0)

   :retval ARK_SUCCESS: if successful
   :retval ARK_MEM_NULL: if the SPRKStep memory is ``NULL``
   :retval ARK_ILL_INPUT: if an argument has an illegal value


.. _ARKODE.Usage.SPRKStep.SPRKStepRootfindingInput:


Rootfinding optional input functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The following functions can be called to set optional inputs to
control the rootfinding algorithm, the mathematics of which are
described in :numref:`ARKODE.Mathematics.Rootfinding`.


.. _ARKODE.Usage.SPRKStep.SPRKStepRootfindingInputTable:
.. table:: Rootfinding optional input functions

   +-----------------------------------------+-------------------------------------------+----------+
   | Optional input                          | Function name                             | Default  |
   +=========================================+===========================================+==========+
   | Direction of zero-crossings to monitor  | :c:func:`SPRKStepSetRootDirection()`      | both     |
   +-----------------------------------------+-------------------------------------------+----------+
   | Disable inactive root warnings          | :c:func:`SPRKStepSetNoInactiveRootWarn()` | enabled  |
   +-----------------------------------------+-------------------------------------------+----------+



.. c:function:: int SPRKStepSetRootDirection(void* arkode_mem, int* rootdir)

   Specifies the direction of zero-crossings to be located and returned.

   The default behavior is to monitor for both zero-crossing directions.

   :param arkode_mem: pointer to the SPRKStep memory block.
   :param rootdir: state array of length *nrtfn*, the number of root
      functions :math:`g_i`  (the value of *nrtfn* was supplied in
      the call to :c:func:`SPRKStepRootInit()`).  If ``rootdir[i] ==
      0`` then crossing in either direction for :math:`g_i` should be
      reported.  A value of +1 or -1 indicates that the solver
      should report only zero-crossings where :math:`g_i` is
      increasing or decreasing, respectively.

   :retval ARK_SUCCESS: if successful
   :retval ARK_MEM_NULL: if the SPRKStep memory is ``NULL``
   :retval ARK_ILL_INPUT: if an argument has an illegal value


.. c:function:: int SPRKStepSetNoInactiveRootWarn(void* arkode_mem)

   Disables issuing a warning if some root function appears
   to be identically zero at the beginning of the integration.

   SPRKStep will not report the initial conditions as a possible zero-crossing
   (assuming that one or more components :math:`g_i` are zero at the initial
   time). However, if it appears that some :math:`g_i` is identically zero at
   the initial time (i.e., :math:`g_i` is zero at the initial time *and* after
   the first step), SPRKStep will issue a warning which can be disabled with
   this optional input function.

   :param arkode_mem: pointer to the SPRKStep memory block.

   :retval ARK_SUCCESS: if successful
   :retval ARK_MEM_NULL: if the SPRKStep memory is ``NULL``


.. _ARKODE.Usage.SPRKStep.InterpolatedOutput:

Interpolated output function
--------------------------------

An optional function :c:func:`SPRKStepGetDky()` is available to obtain
additional values of solution-related quantities.  This function
should only be called after a successful return from
:c:func:`SPRKStepEvolve()`, as it provides interpolated values either of
:math:`y` or of its derivatives.
interpolated to any value of :math:`t` in the last internal step taken
by :c:func:`SPRKStepEvolve()`.



.. c:function:: int SPRKStepGetDky(void* arkode_mem, sunrealtype t, int k, N_Vector dky)

   Computes the *k*-th derivative of the function :math:`y` at the time *t*,
   i.e., :math:`y^{(k)}(t)`, for values of the independent variable satisfying
   :math:`t_n-h_n \le t \le t_n`, with :math:`t_n` as current internal time
   reached, and :math:`h_n` is the last internal step size successfully used by
   the solver. A user may access the values :math:`t_n` and :math:`h_n` via the
   functions :c:func:`SPRKStepGetCurrentTime()` and
   :c:func:`SPRKStepGetLastStep()`, respectively.

   This routine uses an interpolating polynomial of degree *min(degree, 5)*,
   where *degree* is the argument provided to
   :c:func:`SPRKStepSetInterpolantDegree()`.  The user may request *k* in the
   range {0,..., *min(degree, kmax)*} where *kmax* depends on the choice of
   interpolation module. For Hermite interpolants *kmax = 5* and for Lagrange
   interpolants *kmax = 3*.

   :param arkode_mem: pointer to the SPRKStep memory block.
   :param t: the value of the independent variable at which the
        derivative is to be evaluated.
   :param k: the derivative order requested.
   :param dky: output vector (must be allocated by the user).

   :retval ARK_SUCCESS: if successful
   :retval ARK_BAD_K: if *k* is not in the range {0,..., *min(degree, kmax)*}.
   :retval ARK_BAD_T: if *t* is not in the interval :math:`[t_n-h_n, t_n]`
   :retval ARK_BAD_DKY: if the *dky* vector was ``NULL``
   :retval ARK_MEM_NULL: if the SPRKStep memory is ``NULL``

   .. note::

      Dense outputs may or may not conserve the Hamiltonian. Our testing has
      shown that Lagrange interpolation typically performs well in this regard,
      while Hermite interpolation does not.

   .. warning::

      It is only legal to call this function after a successful
      return from :c:func:`SPRKStepEvolve()`.


.. _ARKODE.Usage.SPRKStep.OptionalOutputs:

Optional output functions
------------------------------

SPRKStep provides an extensive set of functions that can be used to
obtain solver performance information.  We organize these into groups:

#. General SPRKStep output routines are in
   :numref:`ARKODE.Usage.SPRKStep.SPRKStepMainOutputs`,

#. Output routines regarding root-finding results are in
   :numref:`ARKODE.Usage.SPRKStep.SPRKStepRootOutputs`,

#. General usability routines (e.g. to print the current SPRKStep
   parameters, or output the current Butcher tables) are in
   :numref:`ARKODE.Usage.SPRKStep.SPRKStepExtraOutputs`.

Following each table, we elaborate on each function.

Some of the optional outputs, especially the various counters, can be
very useful in determining the efficiency of various methods inside
SPRKStep.  For example:

* The counters *nsteps* and *nf_evals* provide a rough measure of the
  overall cost of a given run, and can be compared between runs with
  different solver options to suggest which set of options is the most
  efficient.

.. * The ratio *nsteps/step_attempts* can measure the quality of the
..   time step adaptivity algorithm, since a poor algorithm will result
..   in more failed steps, and hence a lower ratio.

It is therefore recommended that users retrieve and output these
statistics following each run, and take some time to investigate
alternate solver options that will be more optimal for their
particular problem of interest.



.. _ARKODE.Usage.SPRKStep.SPRKStepMainOutputs:

Main solver optional output functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. _ARKODE.Usage.SPRKStep.SPRKStepMainOutputsTable:
.. table:: Main solver optional output functions

  +-----------------------------------------------------+--------------------------------------------+
  | Optional output                                     | Function name                              |
  +=====================================================+============================================+
  | Cumulative number of internal steps                 | :c:func:`SPRKStepGetNumSteps`              |
  +-----------------------------------------------------+--------------------------------------------+
  | Step size used for the last successful step         | :c:func:`SPRKStepGetLastStep`              |
  +-----------------------------------------------------+--------------------------------------------+
  | Step size to be attempted on the next step          | :c:func:`SPRKStepGetCurrentStep`           |
  +-----------------------------------------------------+--------------------------------------------+
  | Current internal time reached by the solver         | :c:func:`SPRKStepGetCurrentTime`           |
  +-----------------------------------------------------+--------------------------------------------+
  | Current internal state reached by the solver        | :c:func:`SPRKStepGetCurrentState`          |
  +-----------------------------------------------------+--------------------------------------------+
  | Single accessor to many statistics at once          | :c:func:`SPRKStepGetStepStats`             |
  +-----------------------------------------------------+--------------------------------------------+
  | Print all statistics                                | :c:func:`SPRKStepPrintAllStats`            |
  +-----------------------------------------------------+--------------------------------------------+
  | Name of constant associated with a return flag      | :c:func:`SPRKStepGetReturnFlagName`        |
  +-----------------------------------------------------+--------------------------------------------+
  | No. of attempted steps                              | :c:func:`SPRKStepGetNumStepAttempts`       |
  +-----------------------------------------------------+--------------------------------------------+
  | No. of calls to right-hand side functions           | :c:func:`SPRKStepGetNumRhsEvals`           |
  +-----------------------------------------------------+--------------------------------------------+
  | Current method table                                | :c:func:`SPRKStepGetCurrentMethod`         |
  +-----------------------------------------------------+--------------------------------------------+
  | Retrieve a pointer for user data                    | :c:func:`SPRKStepGetUserData`              |
  +-----------------------------------------------------+--------------------------------------------+



.. c:function:: int SPRKStepGetNumSteps(void* arkode_mem, long int* nsteps)

   Returns the cumulative number of internal steps taken by
   the solver (so far).

   :param arkode_mem: pointer to the SPRKStep memory block.
   :param nsteps: number of steps taken in the solver.

   :retval ARK_SUCCESS: if successful
   :retval ARK_MEM_NULL: if the SPRKStep memory was ``NULL``


.. c:function:: int SPRKStepGetLastStep(void* arkode_mem, sunrealtype* hlast)

   Returns the integration step size taken on the last successful
   internal step.

   :param arkode_mem: pointer to the SPRKStep memory block.
   :param hlast: step size taken on the last internal step.

   :retval ARK_SUCCESS: if successful
   :retval ARK_MEM_NULL: if the SPRKStep memory was ``NULL``


.. c:function:: int SPRKStepGetCurrentStep(void* arkode_mem, sunrealtype* hcur)

   Returns the integration step size to be attempted on the next internal step.

   :param arkode_mem: pointer to the SPRKStep memory block.
   :param hcur: step size to be attempted on the next internal step.

   :retval ARK_SUCCESS: if successful
   :retval ARK_MEM_NULL: if the SPRKStep memory was ``NULL``


.. c:function:: int SPRKStepGetCurrentTime(void* arkode_mem, sunrealtype* tcur)

   Returns the current internal time reached by the solver.

   :param arkode_mem: pointer to the SPRKStep memory block.
   :param tcur: current internal time reached.

   :retval ARK_SUCCESS: if successful
   :retval ARK_MEM_NULL: if the SPRKStep memory was ``NULL``


.. c:function:: int SPRKStepGetCurrentState(void *arkode_mem, N_Vector *ycur)

   Returns the current internal solution reached by the solver.

   :param arkode_mem: pointer to the SPRKStep memory block.
   :param ycur: current internal solution

   :retval ARK_SUCCESS: if successful
   :retval ARK_MEM_NULL: if the SPRKStep memory was ``NULL``

   .. warning::

      Users should exercise extreme caution when using this function,
      as altering values of *ycur* may lead to undesirable behavior, depending
      on the particular use case and on when this routine is called.


.. c:function:: int SPRKStepGetStepStats(void* arkode_mem, long int* nsteps, sunrealtype* hinused, sunrealtype* hlast, sunrealtype* hcur, sunrealtype* tcur)

   Returns many of the most useful optional outputs in a single call.

   :param arkode_mem: pointer to the SPRKStep memory block.
   :param nsteps: number of steps taken in the solver.
   :param hinused: actual value of initial step size.
   :param hlast: step size taken on the last internal step.
   :param hcur: step size to be attempted on the next internal step.
   :param tcur: current internal time reached.

   :retval ARK_SUCCESS: if successful
   :retval ARK_MEM_NULL: if the SPRKStep memory was ``NULL``


.. c:function:: int SPRKStepPrintAllStats(void* arkode_mem, FILE* outfile, SUNOutputFormat fmt)

   Outputs all of the integrator statistics.

   :param arkode_mem: pointer to the SPRKStep memory block.
   :param outfile: pointer to output file.
   :param fmt: the output format:

       * :c:enumerator:`SUN_OUTPUTFORMAT_TABLE` -- prints a table of values
       * :c:enumerator:`SUN_OUTPUTFORMAT_CSV` -- prints a comma-separated list
         of key and value pairs e.g., ``key1,value1,key2,value2,...``

   :retval ARK_SUCCESS: -- if the output was successfully.
   :retval ARK_MEM_NULL: -- if the SPRKStep memory was ``NULL``.
   :retval ARK_ILL_INPUT: -- if an invalid formatting option was provided.

   .. note::

      The file ``scripts/sundials_csv.py`` provides python utility functions to
      read and output the data from a SUNDIALS CSV output file using the key
      and value pair format.



.. c:function:: char *SPRKStepGetReturnFlagName(long int flag)

   Returns the name of the SPRKStep constant corresponding to *flag*.
   See :ref:`ARKODE.Constants`.

   :param flag: a return flag from an SPRKStep function.

   :returns: The return value is a string containing the name of the
             corresponding constant.


.. c:function:: int SPRKStepGetNumStepAttempts(void* arkode_mem, long int* step_attempts)

   Returns the cumulative number of steps attempted by the solver (so far).

   :param arkode_mem: pointer to the SPRKStep memory block.
   :param step_attempts: number of steps attempted by solver.

   :retval ARK_SUCCESS: if successful
   :retval ARK_MEM_NULL: if the SPRKStep memory was ``NULL``


.. c:function:: int SPRKStepGetNumRhsEvals(void* arkode_mem, long int* nf1, long int* nf2)

   Returns the number of calls to the user's right-hand
   side functions, :math:`f_1` and :math:`f_2` (so far).

   :param arkode_mem: pointer to the SPRKStep memory block.
   :param nf1: number of calls to the user's :math:`f_1(t,p)` function.
   :param nf2: number of calls to the user's :math:`f_2(t,q)` function.

   :retval ARK_SUCCESS: if successful
   :retval ARK_MEM_NULL: if the SPRKStep memory was ``NULL``



.. c:function:: int SPRKStepGetCurrentMethod(void* arkode_mem, ARKodeSPRKTable *sprk_table)

   Returns the SPRK method coefficient table currently in use by the solver.

   :param arkode_mem: pointer to the SPRKStep memory block.
   :param sprk_table: pointer to the SPRK method table.

   :retval ARK_SUCCESS: if successful
   :retval ARK_MEM_NULL: if the SPRKStep memory was ``NULL``


.. c:function:: int SPRKStepGetUserData(void* arkode_mem, void** user_data)

   Returns the user data pointer previously set with
   :c:func:`SPRKStepSetUserData`.

   :param arkode_mem: pointer to the SPRKStep memory block.
   :param user_data: memory reference to a user data pointer

   :retval ARK_SUCCESS: if successful
   :retval ARK_MEM_NULL: if the ARKStep memory was ``NULL``


.. _ARKODE.Usage.SPRKStep.SPRKStepRootOutputs:

Rootfinding optional output functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


.. _ARKODE.Usage.SPRKStep.SPRKStepRootOutputsTable:
.. table:: Rootfinding optional output functions

   +--------------------------------------------------+---------------------------------+
   | Optional output                                  | Function name                   |
   +==================================================+=================================+
   | Array showing roots found                        | :c:func:`SPRKStepGetRootInfo()` |
   +--------------------------------------------------+---------------------------------+
   | No. of calls to user root function               | :c:func:`SPRKStepGetNumGEvals()`|
   +--------------------------------------------------+---------------------------------+



.. c:function:: int SPRKStepGetRootInfo(void* arkode_mem, int* rootsfound)

   Returns an array showing which functions were found to have a root.

   For the components of :math:`g_i` for which a root was found, the sign of
   ``rootsfound[i]`` indicates the direction of zero-crossing. A value of +1
   indicates that :math:`g_i` is increasing, while a value of -1 indicates a
   decreasing :math:`g_i`.

   The user must allocate space for *rootsfound* prior to calling this function.

   :param arkode_mem: pointer to the SPRKStep memory block.
   :param rootsfound: array of length *nrtfn* with the indices of the
        user functions :math:`g_i` found to have a root (the value of
        *nrtfn* was supplied in the call to
        :c:func:`SPRKStepRootInit()`).  For :math:`i = 0 \ldots`
        *nrtfn*-1, ``rootsfound[i]`` is nonzero if :math:`g_i` has a
        root, and 0 if not.

   :retval ARK_SUCCESS: if successful
   :retval ARK_MEM_NULL: if the SPRKStep memory was ``NULL``


.. c:function:: int SPRKStepGetNumGEvals(void* arkode_mem, long int* ngevals)

   Returns the cumulative number of calls made to the
   user's root function :math:`g`.

   :param arkode_mem: pointer to the SPRKStep memory block.
   :param ngevals: number of calls made to :math:`g` so far.

   :retval ARK_SUCCESS: if successful
   :retval ARK_MEM_NULL: if the SPRKStep memory was ``NULL``


.. _ARKODE.Usage.SPRKStep.SPRKStepExtraOutputs:

General usability functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The following optional routine may be called by a user to inquire
about existing solver parameters.  While it would not typically be called
during the course of solving an initial value problem, it may be useful
for users wishing to better understand SPRKStep.


.. _ARKODE.Usage.SPRKStep.SPRKStepExtraOutputsTable:
.. table:: General usability functions

   +----------------------------------------+--------------------------------------+
   | Optional routine                       | Function name                        |
   +----------------------------------------+--------------------------------------+
   | Output all SPRKStep solver parameters  | :c:func:`SPRKStepWriteParameters()`  |
   +----------------------------------------+--------------------------------------+


.. c:function:: int SPRKStepWriteParameters(void* arkode_mem, FILE *fp)

   Outputs all SPRKStep solver parameters to the provided file pointer.

   The *fp* argument can be ``stdout`` or ``stderr``, or it may point to a
   specific file created using ``fopen``.

   When run in parallel, only one process should set a non-NULL value for this
   pointer, since parameters for all processes would be identical.

   :param arkode_mem: pointer to the SPRKStep memory block.
   :param fp: pointer to use for printing the solver parameters.

   :retval ARK_SUCCESS: if successful
   :retval ARK_MEM_NULL: if the SPRKStep memory was ``NULL``


.. _ARKODE.Usage.SPRKStep.Reinitialization:

SPRKStep re-initialization function
-------------------------------------

To reinitialize the SPRKStep module for the solution of a new problem,
where a prior call to :c:func:`SPRKStepCreate` has been made, the
user must call the function :c:func:`SPRKStepReInit()`.  The new
problem must have the same size as the previous one.  This routine
retains the current settings for all SPRKStep module options and
performs the same input checking and initializations that are done in
:c:func:`SPRKStepCreate`, but it performs no memory allocation as it
assumes that the existing internal memory is sufficient for the new
problem.  A call to this re-initialization routine deletes the
solution history that was stored internally during the previous
integration, and deletes any previously-set *tstop* value specified via a
call to :c:func:`SPRKStepSetStopTime()`.  Following a successful call to
:c:func:`SPRKStepReInit()`, call :c:func:`SPRKStepEvolve()` again for the
solution of the new problem.

The use of :c:func:`SPRKStepReInit()` requires that the number of
Runge--Kutta stages, denoted by *s*, be no larger for the new problem than
for the previous problem.  This condition is automatically fulfilled
if the method order *q* is left unchanged.

One potential use of the :c:func:`SPRKStepReInit()` function is in the
treating of jump discontinuities in the RHS function :cite:p:`Tao:22`.
In lieu of including if statements within the RHS function to handle
discontinuities, it may be more computationally efficient to stop at each
point of discontinuity (e.g., through use of tstop or the rootfinding feature)
and restart the integrator with a readjusted ODE model, using a call to
this routine. We note that for the solution to retain temporal accuracy,
the RHS function should not incorporate the discontinuity.


.. c:function:: int SPRKStepReInit(void* arkode_mem, ARKRhsFn f1, ARKRhsFn f2, sunrealtype t0, N_Vector y0)

   Provides required problem specifications and re-initializes the SPRKStep
   time-stepper module.

   All previously set options are retained but may be updated by calling the
   appropriate "Set" functions.

   If an error occurred, :c:func:`SPRKStepReInit()` also sends an error message
   to the error handler function.

   :param arkode_mem: pointer to the SPRKStep memory block.
   :param f1: the name of the C function (of type :c:func:`ARKRhsFn()`) defining :math:`f1(t,q) = \frac{\partial V(t,q)}{\partial q}`
   :param f2: the name of the C function (of type :c:func:`ARKRhsFn()`) defining :math:`f2(t,p) = \frac{\partial T(t,p)}{\partial p}`
   :param t0: the initial value of :math:`t`.
   :param y0: the initial condition vector :math:`y(t_0)`.

   :retval ARK_SUCCESS: if successful
   :retval ARK_MEM_NULL:  if the SPRKStep memory was ``NULL``
   :retval ARK_MEM_FAIL:  if a memory allocation failed
   :retval ARK_ILL_INPUT: if an argument has an illegal value.


.. _ARKODE.Usage.SPRKStep.Reset:

SPRKStep reset function
-----------------------

To reset the SPRKStep module to a particular state :math:`(t_R,y(t_R))` for the
continued solution of a problem, where a prior
call to :c:func:`SPRKStepCreate` has been made, the user must call the function
:c:func:`SPRKStepReset()`.  Like :c:func:`SPRKStepReInit()` this routine retains
the current settings for all SPRKStep module options and performs no memory
allocations but, unlike :c:func:`SPRKStepReInit()`, this routine performs only a
*subset* of the input checking and initializations that are done in
:c:func:`SPRKStepCreate`. In particular this routine retains all internal
counter values. Like :c:func:`SPRKStepReInit()`, a call to
:c:func:`SPRKStepReset()` will delete any previously-set *tstop* value specified
via a call to :c:func:`SPRKStepSetStopTime()`.  Following a successful call to
:c:func:`SPRKStepReset()`, call :c:func:`SPRKStepEvolve()` again to continue
solving the problem. By default the next call to :c:func:`SPRKStepEvolve()` will
use the step size computed by SPRKStep prior to calling :c:func:`SPRKStepReset()`.


.. c:function:: int SPRKStepReset(void* arkode_mem, sunrealtype tR, N_Vector yR)

   Resets the current SPRKStep time-stepper module state to the provided
   independent variable value and dependent variable vector.

   All previously set options are retained but may be updated by calling
   the appropriate "Set" functions.

   If an error occurred, :c:func:`SPRKStepReset()` also sends an error message
   to the error handler function.

   :param arkode_mem: pointer to the SPRKStep memory block.
   :param tR: the value of the independent variable :math:`t`.
   :param yR: the value of the dependent variable vector :math:`y(t_R)`.

   :retval ARK_SUCCESS: if successful
   :retval ARK_MEM_NULL:  if the SPRKStep memory was ``NULL``
   :retval ARK_MEM_FAIL:  if a memory allocation failed
   :retval ARK_ILL_INPUTL: if an argument has an illegal value.

   .. note::

      By default the next call to :c:func:`SPRKStepEvolve()` will use the step
      size computed by SPRKStep prior to calling :c:func:`SPRKStepReset()`.
