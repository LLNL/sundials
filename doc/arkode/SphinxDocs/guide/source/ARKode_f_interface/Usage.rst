..
   Programmer(s): Daniel R. Reynolds @ SMU
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

.. _FInterface.Usage:

Usage of the FARKODE interface module
==========================================

The usage of FARKODE requires calls to a variety of interface
functions, depending on the method options selected, and two or more
user-supplied routines which define the problem to be solved.  These
function calls and user routines are summarized separately below.
Some details are omitted, and the user is referred to the description
of the corresponding C interface ARKStep functions for complete
information on the arguments of any given user-callable interface
routine, or of a given user-supplied function called by an interface
function.  The usage of FARKODE for rootfinding and with
preconditioner modules is described in later subsections.



.. _FInterface.RHS:

Right-hand side specification
--------------------------------------

The user must in all cases supply the following Fortran routines:

.. f:subroutine:: FARKIFUN(T, Y, YDOT, IPAR, RPAR, IER)

   Sets the *YDOT* array to :math:`f^I(t,y)`, the implicit portion of
   the right-hand side of the ODE system, as function of the
   independent variable *T* :math:`=t` and the array of dependent state
   variables *Y* :math:`=y`.

   **Arguments:**
      * *T* (``realtype``, input) -- current value of the independent variable.
      * *Y* (``realtype``, input) -- array containing state variables.
      * *YDOT* (``realtype``, output) -- array containing state derivatives.
      * *IPAR* (``long int``, input) -- array containing integer user
        data that was passed to :f:func:`FARKMALLOC()`.
      * *RPAR* (``realtype``, input) -- array containing real user
        data that was passed to :f:func:`FARKMALLOC()`.
      * *IER* (``int``, output) -- return flag (0 success, >0
        recoverable error, <0 unrecoverable error).


.. f:subroutine:: FARKEFUN(T, Y, YDOT, IPAR, RPAR, IER)

   Sets the *YDOT* array to :math:`f^E(t,y)`, the explicit portion of
   the right-hand side of the ODE system, as function of the
   independent variable *T* :math:`=t` and the array of dependent state
   variables *Y* :math:`=y`.

   **Arguments:**
      * *T* (``realtype``, input) -- current value of the independent variable.
      * *Y* (``realtype``, input) -- array containing state variables.
      * *YDOT* (``realtype``, output) -- array containing state derivatives.
      * *IPAR* (``long int``, input) -- array containing integer user
        data that was passed to :f:func:`FARKMALLOC()`.
      * *RPAR* (``realtype``, input) -- array containing real user
        data that was passed to :f:func:`FARKMALLOC()`.
      * *IER* (``int``, output) -- return flag (0 success, >0
        recoverable error, <0 unrecoverable error).

For purely explicit problems, although the routine
:f:func:`FARKIFUN()` must exist, it will never be called, and may
remain empty.  Similarly, for purely implicit problems,
:f:func:`FARKEFUN()` will never be called and must exist and may
remain empty.



.. _FInterface.NVector:

NVECTOR module initialization
--------------------------------------

If using one of the NVECTOR modules supplied with SUNDIALS, the user
must make a call of the form

.. code::

   CALL FNVINITS(4, NEQ, IER)
   CALL FNVINITP(COMM, 4, NLOCAL, NGLOBAL, IER)
   CALL FNVINITOMP(4, NEQ, NUM_THREADS, IER)
   CALL FNVINITPTS(4, NEQ, NUM_THREADS, IER)
   CALL FNVINITPH(COMM, 4, NLOCAL, NGLOBAL, IER)

in which the specific arguments are as described in the
appropriate section of the Chapter :ref:`NVectors`.



.. _FInterface.SUNMatrix:

SUNMATRIX module initialization
--------------------------------------

In the case of using either an implicit or ImEx method, the solution
of each Runge-Kutta stage may involve the solution of linear systems
related to the Jacobian :math:`J = \frac{\partial f^I}{\partial y}` of
the implicit portion of the ODE system.  If using a Newton iteration
with direct SUNLINSOL linear solver module and one of the SUNMATRIX
modules supplied with SUNDIALS, the user must make a call of the form

.. code::

   CALL FSUNBANDMATINIT(4, N, MU, ML, SMU, IER)
   CALL FSUNDENSEMATINIT(4, M, N, IER)
   CALL FSUNSPARSEMATINIT(4, M, N, NNZ, SPARSETYPE, IER)

in which the specific arguments are as described in the appropriate
section of the Chapter :ref:`SUNMatrix`.  Note that these matrix
options are usable only in a serial or multi-threaded environment.

As described in the section :ref:`Mathematics.MassSolve`, in the case
of using a problem with a non-identity mass matrix (no matter whether
the integrator is implicit, explicit or ImEx), linear systems of the
form :math:`Mx=b` must be solved, where :math:`M` is the system mass
matrix.  If these are to be solved with a direct SUNLINSOL linear
solver module and one of the SUNMATRIX modules supplied with SUNDIALS,
the user must make a call of the form

.. code::

   CALL FSUNBANDMASSMATINIT(N, MU, ML, SMU, IER)
   CALL FSUNDENSEMASSMATINIT(M, N, IER)
   CALL FSUNSPARSEMASSMATINIT(M, N, NNZ, SPARSETYPE, IER)

in which the specific arguments are as described in the appropriate
section of the Chapter :ref:`SUNMatrix`, again noting that these are
only usable in a serial or multi-threaded environment.



.. _FInterface.SUNLinSol:

SUNLINSOL module initialization
--------------------------------------

If using a Newton iteration with one of the SUNLINSOL linear
solver modules supplied with SUNDIALS, the user must make a call
of the form

.. code::

   CALL FSUNBANDLINSOLINIT(4, IER)
   CALL FSUNDENSELINSOLINIT(4, IER)
   CALL FSUNKLUINIT(4, IER)
   CALL FSUNLAPACKBANDINIT(4, IER)
   CALL FSUNLAPACKDENSEINIT(4, IER)
   CALL FSUNPCGINIT(4, PRETYPE, MAXL, IER)
   CALL FSUNSPBCGSINIT(4, PRETYPE, MAXL, IER)
   CALL FSUNSPFGMRINIT(4, PRETYPE, MAXL, IER)
   CALL FSUNSPGMRINIT(4, PRETYPE, MAXL, IER)
   CALL FSUNSPTFQMRINIT(4, PRETYPE, MAXL, IER)
   CALL FSUNSUPERLUMTINIT(4, NUM_THREADS, IER)

in which the specific arguments are as described in the
appropriate section of the Chapter :ref:`SUNLinSol`.  Note that the
dense, band and sparse solvers are usable only in a serial or
multi-threaded environment.

Once one of these has been initialized, its solver parameters may be
modified using a call to the functions

.. code::

   CALL FSUNKLUSETORDERING(4, ORD_CHOICE, IER)
   CALL FSUNSUPERLUMTSETORDERING(4, ORD_CHOICE, IER)
   CALL FSUNPCGSETPRECTYPE(4, PRETYPE, IER)
   CALL FSUNPCGSETMAXL(4, MAXL, IER)
   CALL FSUNSPBCGSSETPRECTYPE(4, PRETYPE, IER)
   CALL FSUNSPBCGSSETMAXL(4, MAXL, IER)
   CALL FSUNSPFGMRSETGSTYPE(4, GSTYPE, IER)
   CALL FSUNSPFGMRSETPRECTYPE(4, PRETYPE, IER)
   CALL FSUNSPGMRSETGSTYPE(4, GSTYPE, IER)
   CALL FSUNSPGMRSETPRECTYPE(4, PRETYPE, IER)
   CALL FSUNSPTFQMRSETPRECTYPE(4, PRETYPE, IER)
   CALL FSUNSPTFQMRSETMAXL(4, MAXL, IER)

where again the call sequences are described in the appropriate
sections of the Chapter :ref:`SUNLinSol`.


Similarly, in the case of using one of the SUNLINSOL linear solver
modules supplied with SUNDIALS to solve a problem with a non-identity
mass matrix, the user must make a call of the form

.. code::

   CALL FSUNMASSBANDLINSOLINIT(IER)
   CALL FSUNMASSDENSELINSOLINIT(IER)
   CALL FSUNMASSKLUINIT(IER)
   CALL FSUNMASSLAPACKBANDINIT(IER)
   CALL FSUNMASSLAPACKDENSEINIT(IER)
   CALL FSUNMASSPCGINIT(PRETYPE, MAXL, IER)
   CALL FSUNMASSSPBCGSINIT(PRETYPE, MAXL, IER)
   CALL FSUNMASSSPFGMRINIT(PRETYPE, MAXL, IER)
   CALL FSUNMASSSPGMRINIT(PRETYPE, MAXL, IER)
   CALL FSUNMASSSPTFQMRINIT(PRETYPE, MAXL, IER)
   CALL FSUNMASSSUPERLUMTINIT(NUM_THREADS, IER)

in which the specific arguments are as described in the
appropriate section of the Chapter :ref:`SUNLinSol`.

Once one of these has been initialized, its solver parameters may be
modified using a call to the functions

.. code::

   CALL FSUNMASSKLUSETORDERING(ORD_CHOICE, IER)
   CALL FSUNMASSSUPERLUMTSETORDERING(ORD_CHOICE, IER)
   CALL FSUNMASSPCGSETPRECTYPE(PRETYPE, IER)
   CALL FSUNMASSPCGSETMAXL(MAXL, IER)
   CALL FSUNMASSSPBCGSSETPRECTYPE(PRETYPE, IER)
   CALL FSUNMASSSPBCGSSETMAXL(MAXL, IER)
   CALL FSUNMASSSPFGMRSETGSTYPE(GSTYPE, IER)
   CALL FSUNMASSSPFGMRSETPRECTYPE(PRETYPE, IER)
   CALL FSUNMASSSPGMRSETGSTYPE(GSTYPE, IER)
   CALL FSUNMASSSPGMRSETPRECTYPE(PRETYPE, IER)
   CALL FSUNMASSSPTFQMRSETPRECTYPE(PRETYPE, IER)
   CALL FSUNMASSSPTFQMRSETMAXL(MAXL, IER)

where again the call sequences are described in the appropriate
sections of the Chapter :ref:`SUNLinSol`.




.. _FInterface.SUNNonlinSol:

SUNNONLINSOL module initialization
--------------------------------------

If using a non-default nonlinear solver method, the user must make a call
of the form

.. code::

   CALL FSUNNEWTONINIT(4, IER)
   CALL FSUNFIXEDPOINTINIT(4, M, IER)

in which the specific arguments are as described in the
appropriate section of the Chapter :ref:`SUNNonlinSol`.

Once one of these has been initialized, its solver parameters may be
modified using a call to the functions

.. code::

   CALL FSUNNEWTONSETMAXITERS(4, MAXITERS, IER)
   CALL FSUNFIXEDPOINTSETMAXITERS(4, MAXITERS, IER)

where again the call sequences are described in the appropriate
sections of the Chapter :ref:`SUNNonlinSol`.




.. _FInterface.Problem:

Problem specification
--------------------------------------

To set various problem and solution parameters and allocate internal
memory, the user must call :f:func:`FARKMALLOC()`.


.. f:subroutine:: FARKMALLOC(T0, Y0, IMEX, IATOL, RTOL, ATOL, IOUT, ROUT, IPAR, RPAR, IER)

   Initializes the Fortran interface to the ARKStep solver, providing
   interfaces to the C routines :c:func:`ARKStepCreate()` and
   :c:func:`ARKStepSetUserData()`, as well as one of :c:func:`ARKStepSStolerances()` or
   :c:func:`ARKStepSVtolerances()`.

   **Arguments:**
      * *T0* (``realtype``, input) -- initial value of :math:`t`.
      * *Y0* (``realtype``, input) -- array of initial conditions.
      * *IMEX* (``int``, input) -- flag denoting basic integration
	method: 0 = implicit, 1 = explicit, 2 = ImEx.
      * *IATOL* (``int``, input) -- type for absolute tolerance input
	*ATOL*: 1 = scalar, 2 = array, 3 = user-supplied function; the
	user must subsequently call :f:func:`FARKEWTSET()` and supply
	a routine :f:func:`FARKEWT()` to compute the error weight vector.
      * *RTOL* (``realtype``, input) -- scalar relative tolerance.
      * *ATOL* (``realtype``, input) -- scalar or array absolute tolerance.
      * *IOUT* (``long int``, input/output) -- array of length 29 for integer optional outputs.
      * *ROUT* (``realtype``, input/output) -- array of length 6 for real optional outputs.
      * *IPAR* (``long int``, input/output) -- array of user integer data, which will be passed
        unmodified to all user-provided routines.
      * *RPAR* (``realtype``, input/output) -- array with user real data, which will be passed
        unmodified to all user-provided routines.
      * *IER* (``int``, output) -- return flag (0 success, :math:`\ne 0` failure).

   **Notes:** Modifications to the user data arrays *IPAR* and *RPAR*
   inside a user-provided routine will be propagated to all
   subsequent calls to such routines. The optional outputs
   associated with the main ARKStep integrator are listed in
   :ref:`FInterface.IOUTTable` and :ref:`FInterface.ROUTTable`, in
   the section :ref:`FInterface.OptionalOutputs`.


As an alternative to providing tolerances in the call to
:f:func:`FARKMALLOC()`, the user may provide a routine to compute the
error weights used in the WRMS norm evaluations.  If supplied, it must
have the following form:

.. f:subroutine:: FARKEWT(Y, EWT, IPAR, RPAR, IER)

   It must set the positive components of the error weight
   vector *EWT* for the calculation of the WRMS norm of *Y*.

   **Arguments:**
      * *Y* (``realtype``, input) -- array containing state variables.
      * *EWT* (``realtype``, output) -- array containing the error weight vector.
      * *IPAR* (``long int``, input) -- array containing the integer user data that was passed
        to :f:func:`FARKMALLOC()`.
      * *RPAR* (``realtype``, input) -- array containing the real user data that was passed to
        :f:func:`FARKMALLOC()`.
      * *IER* (``int``, output) -- return flag (0 success, :math:`\ne 0` failure).


If the :f:func:`FARKEWT()` routine is provided, then, following the
call to :f:func:`FARKMALLOC()`, the user must call the function
:f:func:`FARKEWTSET()`.

.. f:subroutine:: FARKEWTSET(FLAG, IER)

   Informs FARKODE to use the user-supplied :f:func:`FARKEWT()` function.

   **Arguments:**
      * *FLAG* (``int``, input) -- flag, use "1" to denoting to use :f:func:`FARKEWT()`.
      * *IER* (``int``, output) -- return flag (0 success, :math:`\ne 0` failure).




.. _FInterface.OptionalInputs:

Setting optional inputs
--------------------------------------

Unlike ARKStep's C interface, that provides separate functions for
setting each optional input, FARKODE uses only three functions, that
accept keywords to specify which optional input should be set to the
provided value.  These routines are :f:func:`FARKSETIIN()`,
:f:func:`FARKSETRIN()`, and :f:func:`FARKSETVIN()` and are further
described below.


.. f:subroutine:: FARKSETIIN(KEY, IVAL, IER)

   Specification routine to pass optional integer inputs
   to the :f:func:`FARKODE()` solver.

   **Arguments:**
      * *KEY* (quoted string, input) -- which optional input
        is set (see :ref:`FInterface.IINOptionTable`).
      * *IVAL* (``long int``, input) -- the integer input value to be used.
      * *IER* (``int``, output) -- return flag (0 success, :math:`\ne 0` failure).



.. _FInterface.IINOptionTable:

Table: Keys for setting FARKODE integer optional inputs
--------------------------------------------------------

.. cssclass:: table-bordered

=======================  =========================================
Key                      ARKStep routine
=======================  =========================================
``ORDER``                :c:func:`ARKStepSetOrder()`
``DENSE_ORDER``          :c:func:`ARKStepSetDenseOrder()`
``LINEAR``               :c:func:`ARKStepSetLinear()`
``NONLINEAR``            :c:func:`ARKStepSetNonlinear()`
``EXPLICIT``             :c:func:`ARKStepSetExplicit()`
``IMPLICIT``             :c:func:`ARKStepSetImplicit()`
``IMEX``                 :c:func:`ARKStepSetImEx()`
``IRK_TABLE_NUM``        :c:func:`ARKStepSetTableNum()`
``ERK_TABLE_NUM``        :c:func:`ARKStepSetTableNum()`
``ARK_TABLE_NUM`` *(a)*  :c:func:`ARKStepSetTableNum()`
``MAX_NSTEPS``           :c:func:`ARKStepSetMaxNumSteps()`
``HNIL_WARNS``           :c:func:`ARKStepSetMaxHnilWarns()`
``PREDICT_METHOD``       :c:func:`ARKStepSetPredictorMethod()`
``MAX_ERRFAIL``          :c:func:`ARKStepSetMaxErrTestFails()`
``MAX_CONVFAIL``         :c:func:`ARKStepSetMaxConvFails()`
``MAX_NITERS``           :c:func:`ARKStepSetMaxNonlinIters()`
``ADAPT_SMALL_NEF``      :c:func:`ARKStepSetSmallNumEFails()`
``LSETUP_MSBP``          :c:func:`ARKStepSetMaxStepsBetweenLSet()`
``MAX_CONSTR_FAIL``      :c:func:`ARKStepSetMaxNumConstrFails()`
=======================  =========================================

*(a)* When setting ``ARK_TABLE_NUM``, pass in *IVAL* as an array of
length 2, specifying the IRK table number first, then the ERK table
number.  The integer specifiers for each table may be found in the
section :ref:`Constants`, or in the ARKode header files
``arkode_butcher_dirk.h`` and ``arkode_butcher_erk.h``.


.. f:subroutine:: FARKSETRIN(KEY, RVAL, IER)

   Specification routine to pass optional real inputs
   to the :f:func:`FARKODE()` solver.

   **Arguments:**
      * *KEY* (quoted string, input) -- which optional input
        is set (see :ref:`FInterface.RINOptionTable`).
      * *RVAL* (``realtype``, input) -- the real input value to be used.
      * *IER* (``int``, output) -- return flag (0 success, :math:`\ne 0` failure).



.. _FInterface.RINOptionTable:

Table: Keys for setting FARKODE real optional inputs
-----------------------------------------------------

.. cssclass:: table-bordered

=================  =========================================
Key                ARKStep routine
=================  =========================================
``INIT_STEP``      :c:func:`ARKStepSetInitStep()`
``MAX_STEP``       :c:func:`ARKStepSetMaxStep()`
``MIN_STEP``       :c:func:`ARKStepSetMinStep()`
``STOP_TIME``      :c:func:`ARKStepSetStopTime()`
``NLCONV_COEF``    :c:func:`ARKStepSetNonlinConvCoef()`
``ADAPT_CFL``      :c:func:`ARKStepSetCFLFraction()`
``ADAPT_SAFETY``   :c:func:`ARKStepSetSafetyFactor()`
``ADAPT_BIAS``     :c:func:`ARKStepSetErrorBias()`
``ADAPT_GROWTH``   :c:func:`ARKStepSetMaxGrowth()`
``ADAPT_ETAMX1``   :c:func:`ARKStepSetMaxFirstGrowth()`
``ADAPT_BOUNDS``   :c:func:`ARKStepSetFixedStepBounds()`
``ADAPT_ETAMXF``   :c:func:`ARKStepSetMaxEFailGrowth()`
``ADAPT_ETACF``    :c:func:`ARKStepSetMaxCFailGrowth()`
``NONLIN_CRDOWN``  :c:func:`ARKStepSetNonlinCRDown()`
``NONLIN_RDIV``    :c:func:`ARKStepSetNonlinRDiv()`
``LSETUP_DGMAX``   :c:func:`ARKStepSetDeltaGammaMax()`
``FIXED_STEP``     :c:func:`ARKStepSetFixedStep()`
=================  =========================================


.. f:subroutine:: FARKSETVIN(KEY, VVAL, IER)

   Specification routine to pass optional vector inputs
   to the :f:func:`FARKODE()` solver.

   **Arguments:**
      * *KEY* (quoted string, input) -- which optional input
        is set (see :ref:`FInterface.VINOptionTable`).
      * *VVAL* (``realtype*``, input) -- the input vector of real values to be used.
      * *IER* (``int``, output) -- return flag (0 success, :math:`\ne 0` failure).



.. _FInterface.VINOptionTable:

Table: Keys for setting FARKODE vector optional inputs
-------------------------------------------------------

.. cssclass:: table-bordered

=================  =========================================
Key                ARKStep routine
=================  =========================================
``CONSTR_VEC``     :c:func:`ARKStepSetConstraints()`
=================  =========================================




If a user wishes to reset all of the options to their default values,
they may call the routine :f:func:`FARKSETDEFAULTS()`.

.. f:subroutine:: FARKSETDEFAULTS(IER)

   Specification routine to reset all FARKODE optional
   inputs to their default values.

   **Arguments:**
      * *IER* (``int``, output) -- return flag (0 success, :math:`\ne 0` failure).




Optional advanced FARKODE inputs
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

FARKODE supplies additional routines to specify optional advanced
inputs to the :c:func:`ARKStepEvolve()` solver.  These are summarized below,
and the user is referred to their C routine counterparts for more
complete information.



.. f:subroutine:: FARKSETERKTABLE(S, Q, P, C, A, B, BEMBED, IER)

   Interface to the routine :c:func:`ARKStepSetTables()`.

   **Arguments:**
      * *S* (``int``, input) -- number of stages in the table.
      * *Q* (``int``, input) -- global order of accuracy of the method.
      * *P* (``int``, input) -- global order of accuracy of the embedding.
      * *C* (``realtype``, input) -- array of length *S* containing the stage times.
      * *A* (``realtype``, input) -- array of length *S*S* containing the ERK coefficients
        (stored in row-major, "C", order).
      * *B* (``realtype``, input) -- array of length *S* containing the solution coefficients.
      * *BEMBED* (``realtype``, input) -- array of length *S* containing the embedding
        coefficients.
      * *IER* (``int``, output) -- return flag (0 success, :math:`\ne 0` failure).


.. f:subroutine:: FARKSETIRKTABLE(S, Q, P, C, A, B, BEMBED, IER)

   Interface to the routine :c:func:`ARKStepSetTables()`.

   **Arguments:**
      * *S* (``int``, input) -- number of stages in the table.
      * *Q* (``int``, input) -- global order of accuracy of the method.
      * *P* (``int``, input) -- global order of accuracy of the embedding.
      * *C* (``realtype``, input) -- array of length *S* containing the stage times.
      * *A* (``realtype``, input) -- array of length *S*S* containing the IRK coefficients
        (stored in row-major, "C", order).
      * *B* (``realtype``, input) -- array of length *S* containing the solution coefficients.
      * *BEMBED* (``realtype``, input) -- array of length *S* containing the embedding
        coefficients.
      * *IER* (``int``, output) -- return flag (0 success, :math:`\ne 0` failure).


.. f:subroutine:: FARKSETARKTABLES(S, Q, P, CI, CE, AI, AE, BI, BE, B2I, B2E, IER)

   Interface to the routine :c:func:`ARKStepSetTables()`.

   **Arguments:**
      * *S* (``int``, input) -- number of stages in the table.
      * *Q* (``int``, input) -- global order of accuracy of the method.
      * *P* (``int``, input) -- global order of accuracy of the embedding.
      * *CI* (``realtype``, input) -- array of length *S* containing
	the implicit stage times.
      * *CE* (``realtype``, input) -- array of length *S* containing
	the explicit stage times.
      * *AI* (``realtype``, input) -- array of length *S*S* containing the IRK coefficients
        (stored in row-major, "C", order).
      * *AE* (``realtype``, input) -- array of length *S*S* containing the ERK coefficients
        (stored in row-major, "C", order).
      * *BI* (``realtype``, input) -- array of length *S* containing
	the implicit solution coefficients.
      * *BE* (``realtype``, input) -- array of length *S* containing
	the explicit solution coefficients.
      * *B2I* (``realtype``, input) -- array of length *S* containing
	the implicit embedding coefficients.
      * *B2E* (``realtype``, input) -- array of length *S* containing
	the explicit embedding coefficients.
      * *IER* (``int``, output) -- return flag (0 success, :math:`\ne 0` failure).


.. f:subroutine:: FARKSETRESTOLERANCE(IATOL, ATOL, IER)

   Interface to the routines :c:func:`ARKStepResStolerance()` and :c:func:`ARKStepResVtolerance()`.

   **Arguments:**
      * *IATOL* (``int``, input) -- type for absolute residual tolerance input
	*ATOL*: 1 = scalar, 2 = array.
      * *ATOL* (``realtype``, input) -- scalar or array absolute residual tolerance.
      * *IER* (``int``, output) -- return flag (0 success, :math:`\ne 0` failure).



Additionally, a user may set the accuracy-based step size adaptivity
strategy (and it's associated parameters) through a call to
:f:func:`FARKSETADAPTIVITYMETHOD()`, as described below.

.. f:subroutine:: FARKSETADAPTIVITYMETHOD(IMETHOD, IDEFAULT, IPQ, PARAMS, IER)

   Specification routine to set the step size adaptivity strategy and
   parameters within the :f:func:`FARKODE()` solver.  Interfaces with
   the C routine :c:func:`ARKStepSetAdaptivityMethod()`.

   **Arguments:**
      * *IMETHOD* (``int``, input) -- choice of adaptivity method.
      * *IDEFAULT* (``int``, input) -- flag denoting whether to use
	default parameters (1) or that customized parameters will be
	supplied (1).
      * *IPQ* (``int``, input) -- flag denoting whether to use
	the embedding order of accuracy (0) or the method order of
	accuracy (1) within step adaptivity algorithm.
      * *PARAMS* (``realtype``, input) -- array of 3 parameters to be
	used within the adaptivity strategy.
      * *IER* (``int``, output) -- return flag (0 success, :math:`\ne 0` failure).


Lastly, the user may provide functions to aid/replace those within
ARKStep for handling adaptive error control and explicit stability.
The former of these is designed for advanced users who wish to
investigate custom step adaptivity approaches as opposed to using any
of those built-in to ARKStep.  In ARKStep's C/C++ interface, this would be
provided by a function of type :c:func:`ARKAdaptFn()`; in the Fortran
interface this is provided through the user-supplied function:

.. f:subroutine:: FARKADAPT(Y, T, H1, H2, H3, E1, E2, E3, Q, P, HNEW, IPAR, RPAR, IER)

   It must set the new step size *HNEW* based on the three previous
   steps (*H1*, *H2*, *H3*) and the three previous error estimates
   (*E1*, *E2*, *E3*).

   **Arguments:**
      * *Y* (``realtype``, input) -- array containing state variables.
      * *T* (``realtype``, input) -- current value of the independent variable.
      * *H1* (``realtype``, input) -- current step size.
      * *H2* (``realtype``, input) -- previous step size.
      * *H3* (``realtype``, input) -- previous-previous step size.
      * *E1* (``realtype``, input) -- estimated temporal error in current step.
      * *E2* (``realtype``, input) -- estimated temporal error in previous step.
      * *E3* (``realtype``, input) -- estimated temporal error in previous-previous step.
      * *Q* (``int``, input) -- global order of accuracy for RK method.
      * *P* (``int``, input) -- global order of accuracy for RK embedded method.
      * *HNEW* (``realtype``, output) -- array containing the error weight vector.
      * *IPAR* (``long int``, input) -- array containing the integer
	user data that was passed to :f:func:`FARKMALLOC()`.
      * *RPAR* (``realtype``, input) -- array containing the real user
	data that was passed to :f:func:`FARKMALLOC()`.
      * *IER* (``int``, output) -- return flag (0 success, :math:`\ne 0` failure).


This routine is enabled by a call to the activation routine:

.. f:subroutine:: FARKADAPTSET(FLAG, IER)

   Informs FARKODE to use the user-supplied :f:func:`FARKADAPT()` function.

   **Arguments:**
      * *FLAG* (``int``, input) -- flag, use "1" to denoting to use
	:f:func:`FARKADAPT()`, or use "0" to denote a return to the
        default adaptivity strategy.
      * *IER* (``int``, output) -- return flag (0 success, :math:`\ne
	0` failure).

   Note: The call to :f:func:`FARKADAPTSET()` must occur *after* the call
   to :f:func:`FARKMALLOC()`.


Similarly, if either an explicit or mixed implicit-explicit
integration method is to be employed, the user may specify a function
to provide the maximum explicitly-stable step for their problem.
Again, in the C/C++ interface this would be a function of type
:c:func:`ARKExpStabFn()`, while in ARKStep's Fortran interface this
must be given through the user-supplied function:

.. f:subroutine:: FARKEXPSTAB(Y, T, HSTAB, IPAR, RPAR, IER)

   It must set the maximum explicitly-stable step size, *HSTAB*, based
   on the current solution, *Y*.

   **Arguments:**
      * *Y* (``realtype``, input) -- array containing state variables.
      * *T* (``realtype``, input) -- current value of the independent variable.
      * *HSTAB* (``realtype``, output) -- maximum explicitly-stable step size.
      * *IPAR* (``long int``, input) -- array containing the integer user data that was passed
        to :f:func:`FARKMALLOC()`.
      * *RPAR* (``realtype``, input) -- array containing the real user data that was passed to
        :f:func:`FARKMALLOC()`.
      * *IER* (``int``, output) -- return flag (0 success, :math:`\ne 0` failure).


This routine is enabled by a call to the activation routine:

.. f:subroutine:: FARKEXPSTABSET(FLAG, IER)

   Informs FARKODE to use the user-supplied :f:func:`FARKEXPSTAB()` function.

   **Arguments:**
      * *FLAG* (``int``, input) -- flag, use "1" to denoting to use
	:f:func:`FARKEXPSTAB()`, or use "0" to denote a return to the
        default error-based stability strategy.
      * *IER* (``int``, output) -- return flag (0 success, :math:`\ne
	0` failure).

   Note: The call to :f:func:`FARKEXPSTABSET()` must occur *after* the call
   to :f:func:`FARKMALLOC()`.




.. _FInterface.NonlinearSolver:

Nonlinear solver module specification
----------------------------------------------

To use a non-default nonlinear solver algorithm, then after it has
been initialized in step :ref:`FInterface.SUNNonlinSol` above, the
user of FARKODE must attach it to ARKSTEP by calling the
:f:func:`FARKNLSINIT()` routine:


.. f:subroutine:: FARKNLSINIT(IER)

   Interfaces with the :c:func:`ARKStepSetNonlinearSolver()` function to
   specify use of a non-default nonlinear solver module.

   **Arguments:**
      * *IER* (``int``, output) -- return flag (0 if success, -1 if a memory allocation
        error occurred, -2 for an illegal input).



.. _FInterface.LinearSolver:

System linear solver interface specification
----------------------------------------------

To attach the linear solver (and optionally the matrix) object(s)
initialized in steps :ref:`FInterface.SUNMatrix` and
:ref:`FInterface.SUNLinSol` above, the user of FARKODE must
initialize the linear solver interface.  To attach any SUNLINSOL
object (and optional SUNMATRIX object) to ARKStep, following calls to
initialize the SUNLINSOL (and SUNMATRIX) object(s) in steps
:ref:`FInterface.SUNMatrix` and :ref:`FInterface.SUNLinSol` above, the
user must call the :f:func:`FARKLSINIT()` routine:


.. f:subroutine:: FARKLSINIT(IER)

   Interfaces with the :c:func:`ARKStepSetLinearSolver()` function to
   attach a linear solver object (and optionally a matrix object) to ARKStep.

   **Arguments:**
      * *IER* (``int``, output) -- return flag (0 if success, -1 if a memory allocation
        error occurred, -2 for an illegal input).


.. _FInterface.Direct:

Matrix-based linear solvers
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

As an option when using ARKSTEP with either the SUNLINSOL_DENSE or
SUNLINSOL_LAPACKDENSE linear solver modules, the user may supply a
routine that computes a dense approximation of the system Jacobian
:math:`J = \frac{\partial f^I}{\partial y}`.  If supplied, it must
have the following form:

.. f:subroutine:: FARKDJAC(NEQ, T, Y, FY, DJAC, H, IPAR, RPAR, WK1, WK2, WK3, IER)

   Interface to provide a user-supplied dense Jacobian approximation
   function (of type :c:func:`ARKLsJacFn()`), to be used by the
   SUNLINSOL_DENSE or SUNLINSOL_LAPACKDENSE solver modules.

   **Arguments:**
      * *NEQ* (``long int``, input) -- size of the ODE system.
      * *T* (``realtype``, input) -- current value of the independent variable.
      * *Y* (``realtype``, input) -- array containing values of the dependent state variables.
      * *FY* (``realtype``, input) -- array containing values of the dependent state derivatives.
      * *DJAC* (``realtype`` of size (NEQ,NEQ), output) -- 2D array containing the Jacobian entries.
      * *H* (``realtype``, input) -- current step size.
      * *IPAR* (``long int``, input) -- array containing integer user data that was passed to
        :f:func:`FARKMALLOC()`.
      * *RPAR* (``realtype``, input) -- array containing real user data that was passed to
        :f:func:`FARKMALLOC()`.
      * *WK1*, *WK2*, *WK3*  (``realtype``, input) -- array containing temporary workspace
        of same size as *Y*.
      * *IER* (``int``, output) -- return flag (0 if success, >0 if a recoverable error
        occurred, <0 if an unrecoverable error occurred).

   **Notes:** Typically this routine will use only *NEQ*, *T*, *Y*, and
   *DJAC*. It must compute the Jacobian and store it column-wise in *DJAC*.


If the above routine uses difference quotient approximations, it may
need to access the error weight array *EWT* in the calculation of
suitable increments. The array *EWT* can be obtained by calling
:f:func:`FARKGETERRWEIGHTS()` using one of the work arrays as
temporary storage for *EWT*. It may also need the unit roundoff, which
can be obtained as the optional output *ROUT(6)*, passed from the
calling program to this routine using either *RPAR* or a common block.

If the :f:func:`FARKDJAC()` routine is provided, then, following the
call to :f:func:`FARKLSINIT()`, the user must call the routine
:f:func:`FARKDENSESETJAC()`:


.. f:subroutine:: FARKDENSESETJAC(FLAG, IER)

   Interface to the :c:func:`ARKStepSetJacFn()` function, specifying
   to use the user-supplied routine :f:func:`FARKDJAC()` for the
   Jacobian approximation.

   **Arguments:**
      * *FLAG* (``int``, input) -- any nonzero value specifies to use
	:f:func:`FARKDJAC()`.
      * *IER* (``int``, output) -- return flag (0 if success,
	:math:`\ne 0` if an error occurred).


As an option when using ARKStep with either the SUNLINSOL_BAND or
SUNLINSOL_LAPACKBAND linear solver modules, the user may supply a
routine that computes a banded approximation of the linear system
Jacobian :math:`J = \frac{\partial f^I}{\partial y}`. If supplied, it
must have the following form:

.. f:subroutine:: FARKBJAC(NEQ, MU, ML, MDIM, T, Y, FY, BJAC, H, IPAR, RPAR, WK1, WK2, WK3, IER)

   Interface to provide a user-supplied band Jacobian approximation
   function (of type :c:func:`ARKLsJacFn()`), to be used by the
   SUNLINSOL_BAND or SUNLINSOL_LAPACKBAND solver modules.

   **Arguments:**
      * *NEQ* (``long int``, input) -- size of the ODE system.
      * *MU*   (``long int``, input) -- upper half-bandwidth.
      * *ML*   (``long int``, input) -- lower half-bandwidth.
      * *MDIM* (``long int``, input) -- leading dimension of *BJAC* array.
      * *T*    (``realtype``, input) -- current value of the independent variable.
      * *Y*    (``realtype``, input) -- array containing dependent state variables.
      * *FY*   (``realtype``, input) -- array containing dependent state derivatives.
      * *BJAC* (``realtype`` of size *(MDIM,NEQ)*, output) -- 2D array
	containing the Jacobian entries.
      * *H*    (``realtype``, input) -- current step size.
      * *IPAR* (``long int``, input) -- array containing integer user data that was passed to
        :f:func:`FARKMALLOC()`.
      * *RPAR* (``realtype``, input) -- array containing real user data that was passed to
        :f:func:`FARKMALLOC()`.
      * *WK1*, *WK2*, *WK3*  (``realtype``, input) -- array containing temporary workspace
        of same size as *Y*.
      * *IER* (``int``, output) -- return flag (0 if success, >0 if a recoverable error
        occurred, <0 if an unrecoverable error occurred).

   **Notes:**
   Typically this routine will use only *NEQ*, *MU*, *ML*, *T*, *Y*, and
   *BJAC*. It must load the *MDIM* by *N* array *BJAC* with the Jacobian
   matrix at the current :math:`(t,y)` in band form.  Store in
   *BJAC(k,j)* the Jacobian element :math:`J_{i,j}` with
   *k = i - j + MU + 1* (or *k = 1, ..., ML+MU+1*) and *j = 1, ..., N*.


If the above routine uses difference quotient approximations, it may
need to use the error weight array *EWT* in the calculation of
suitable increments. The array *EWT* can be obtained by calling
:f:func:`FARKGETERRWEIGHTS()` using one of the work
arrays as temporary storage for *EWT*. It may also need the unit
roundoff, which can be obtained as the optional output *ROUT(6)*,
passed from the calling program to this routine using either *RPAR*
or a common block.

If the :f:func:`FARKBJAC()` routine is provided, then, following the
call to :f:func:`FARKLSINIT()`, the user must call the routine
:f:func:`FARKBANDSETJAC()`.


.. f:subroutine:: FARKBANDSETJAC(FLAG, IER)

   Interface to the :c:func:`ARKStepSetJacFn()` function, specifying
   to use the user-supplied routine :f:func:`FARKBJAC()` for the
   Jacobian approximation.

   **Arguments:**
      * *FLAG* (``int``, input) -- any nonzero value specifies to use
        :f:func:`FARKBJAC()`.
      * *IER* (``int``, output) -- return flag (0 if success,
	:math:`\ne 0` if an error occurred).


When using ARKStep with either the SUNLINSOL_KLU or
SUNLINSOL_SUPERLUMT sparse direct linear solver modules, the user must
supply a routine that computes a sparse approximation of the system
Jacobian :math:`J = \frac{\partial f^I}{\partial y}`.  Both the KLU
and SuperLU_MT solvers allow specification of :math:`J` in either
compressed-sparse-column (CSC) format or compressed-sparse-row (CSR)
format.  The sparse Jacobian approximation function must have the
following form:


.. f:subroutine:: FARKSPJAC(T, Y, FY, N, NNZ, JDATA, JINDEXVALS, JINDEXPTRS, H, IPAR, RPAR, WK1, WK2, WK3, IER)

   Interface to provide a user-supplied sparse Jacobian approximation
   function (of type :c:func:`ARKLsJacFn()`), to be used by the
   SUNLINSOL_KLU or SUNLINSOL_SUPERLUMT solver modules.

   **Arguments:**
      * *T* (``realtype``, input) -- current value of the independent variable.
      * *Y* (``realtype``, input) -- array containing values of the dependent state variables.
      * *FY* (``realtype``, input) -- array containing values of the dependent state derivatives.
      * *N* (``sunindextype``, input) -- number of matrix rows and columns in Jacobian.
      * *NNZ* (``sunindextype``, input) -- allocated length of nonzero storage in Jacobian.
      * *JDATA* (``realtype`` of size NNZ, output) -- nonzero values in Jacobian.
      * *JINDEXVALS* (``sunindextype`` of size NNZ, output) -- row *[CSR: column]* indices for each
	nonzero Jacobian entry.
      * *JINDEXPTRS* (``sunindextype`` of size N+1, output) -- indices of where
	each column's *[CSR: row's]* nonzeros begin in data array; last entry points
	just past end of data values.
      * *H* (``realtype``, input) -- current step size.
      * *IPAR* (``long int``, input) -- array containing integer user data that was passed to
        :f:func:`FARKMALLOC()`.
      * *RPAR* (``realtype``, input) -- array containing real user data that was passed to
        :f:func:`FARKMALLOC()`.
      * *WK1*, *WK2*, *WK3*  (``realtype``, input) -- array containing temporary workspace
        of same size as *Y*.
      * *IER* (``int``, output) -- return flag (0 if success, >0 if a recoverable error
        occurred, <0 if an unrecoverable error occurred).

   **Notes:** due to the internal storage format of the
   SUNMATRIX_SPARSE module, the matrix-specific integer parameters and
   arrays are all of type ``sunindextype`` -- the index precision
   (32-bit vs 64-bit signed integers) specified during the SUNDIALS
   build.  It is assumed that the user's Fortran codes are constructed
   to have matching type to how SUNDIALS was installed.

If the above routine uses difference quotient approximations to
compute the nonzero entries, it may need to access the error weight
array *EWT* in the calculation of suitable increments. The array *EWT*
can be obtained by calling :f:func:`FARKGETERRWEIGHTS()` using one of
the work arrays as temporary storage for *EWT*.  It may also need the
unit roundoff, which can be obtained as the optional output *ROUT(6)*,
passed from the calling program to this routine using either *RPAR* or
a common block.

When supplying the :f:func:`FARKSPJAC()` routine, following the call
to :f:func:`FARKLSINIT()`, the user must call the routine
:f:func:`FARKSPARSESETJAC()`.


.. f:subroutine:: FARKSPARSESETJAC(IER)

   Interface to the :c:func:`ARKStepSetJacFn()` function,
   specifying that the user-supplied routine :f:func:`FARKSPJAC()` has
   been provided for the Jacobian approximation.

   **Arguments:**
      * *IER* (``int``, output) -- return flag (0 if success,
	:math:`\ne 0` if an error occurred).



.. _FInterface.Iterative:

Iterative linear solvers
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

As described in the section :ref:`Mathematics.Error.Linear`, a user
may adjust the linear solver tolerance scaling factor
:math:`\epsilon_L`.  Fortran users may adjust this value by calling
the function :f:func:`FARKLSSETEPSLIN()`:

.. f:subroutine:: FARKLSSETEPSLIN(EPLIFAC, IER)

   Interface to the function :c:func:`ARKStepSetEpsLin()` to
   specify the linear solver tolerance scale factor :math:`\epsilon_L`
   for the Newton system linear solver.

   This routine must be called *after* :f:func:`FARKLSINIT()`.

   **Arguments:**
      * *EPLIFAC* (``realtype``, input) -- value to use for
        :math:`\epsilon_L`.  Passing a value of 0 indicates to use the
        default value (0.05).
      * *IER*  (``int``, output) -- return flag  (0 if success,
	:math:`\ne 0` if an error).


Optional user-supplied routines :f:func:`FARKJTSETUP()` and
:f:func:`FARKJTIMES()` may be provided to compute the product
of the system Jacobian :math:`J = \frac{\partial f^I}{\partial y}` and
a given vector :math:`v`.  If these are supplied, then following the
call to :f:func:`FARKLSINIT()`, the user must call the
:f:func:`FARKLSSETJAC()` routine with *FLAG* :math:`\ne 0`:

.. f:subroutine:: FARKLSSETJAC(FLAG, IER)

   Interface to the function :c:func:`ARKStepSetJacTimes()` to
   specify use of the user-supplied Jacobian-times-vector setup and
   product functions, :f:func:`FARKJTSETUP()` and
   :f:func:`FARKJTIMES()`, respectively.

   This routine must be called *after* :f:func:`FARKLSINIT()`.

   **Arguments:**
      * *FLAG* (``int``, input) -- flag denoting use of user-supplied
        Jacobian-times-vector routines.  A nonzero value specifies to
        use these the user-supplied routines, a zero value specifies
        not to use these.
      * *IER*  (``int``, output) -- return flag  (0 if success,
	:math:`\ne 0` if an error).


Similarly, optional user-supplied routines :f:func:`FARKPSET()` and
:f:func:`FARKPSOL()` may be provided to perform preconditioning of the
iterative linear solver (note: the SUNLINSOL module must have been
configured with preconditioning enabled).  If these routines are
supplied, then following the call to :f:func:`FARKLSINIT()` the
user must call the routine :f:func:`FARKLSSETPREC()` with *FLAG*
:math:`\ne 0`:

.. f:subroutine:: FARKLSSETPREC(FLAG, IER)

   Interface to the function :c:func:`ARKStepSetPreconditioner()` to
   specify use of the user-supplied preconditioner setup and solve
   functions, :f:func:`FARKPSET()` and :f:func:`FARKPSOL()`, respectively.

   This routine must be called *after* :f:func:`FARKLSINIT()`.

   **Arguments:**
      * *FLAG* (``int``, input) -- flag denoting use of user-supplied
        preconditioning routines.  A nonzero value specifies to
        use these the user-supplied routines, a zero value specifies
        not to use these.
      * *IER*  (``int``, output) -- return flag  (0 if success,
	:math:`\ne 0` if an error).


With treatment of the linear systems by any of the Krylov iterative
solvers, there are four optional user-supplied routines --
:f:func:`FARKJTSETUP()`, :f:func:`FARKJTIMES()`, :f:func:`FARKPSET()`
and :f:func:`FARKPSOL()`. The specifications of these functions are
given below.


As an option when using iterative linear solvers, the user
may supply a routine that computes the product of the system Jacobian
:math:`J = \frac{\partial f^I}{\partial y}` and a given vector
:math:`v`.  If supplied, it must have the following form:



.. f:subroutine:: FARKJTIMES(V, FJV, T, Y, FY, H, IPAR, RPAR, WORK, IER)

   Interface to provide a user-supplied Jacobian-times-vector product
   approximation function (corresponding to a C interface routine of
   type :c:func:`ARKLsJacTimesVecFn()`), to be used by one of the
   Krylov iterative linear solvers.

   **Arguments:**
      * *V*    (``realtype``, input) -- array containing the vector to multiply.
      * *FJV*  (``realtype``, output) -- array containing resulting product vector.
      * *T*    (``realtype``, input) -- current value of the independent variable.
      * *Y*    (``realtype``, input) -- array containing dependent state variables.
      * *FY*   (``realtype``, input) -- array containing dependent state derivatives.
      * *H*    (``realtype``, input) -- current step size.
      * *IPAR* (``long int``, input) -- array containing integer user data that was passed to
        :f:func:`FARKMALLOC()`.
      * *RPAR* (``realtype``, input) -- array containing real user data that was passed to
        :f:func:`FARKMALLOC()`.
      * *WORK* (``realtype``, input) -- array containing temporary workspace of same size as
        *Y*.
      * *IER*  (``int``, output) -- return flag  (0 if success, :math:`\ne 0` if an error).

   **Notes:**
   Typically this routine will use only *T*, *Y*, *V*, and
   *FJV*.  It must compute the product vector :math:`Jv`, where
   :math:`v` is given in *V*, and the product is stored in *FJV*.


If the user's Jacobian-times-vector product routine requires that any
Jacobian related data be evaluated or preprocessed, then the following
routine can be used for the evaluation and preprocessing of this data:



.. f:subroutine:: FARKJTSETUP(T, Y, FY, H, IPAR, RPAR, IER)

   Interface to setup data for use in a user-supplied
   Jacobian-times-vector product approximation function (corresponding
   to a C interface routine of type
   :c:func:`ARKLJacTimesSetupFn()`).

   **Arguments:**
      * *T*    (``realtype``, input) -- current value of the independent variable.
      * *Y*    (``realtype``, input) -- array containing dependent state variables.
      * *FY*   (``realtype``, input) -- array containing dependent state derivatives.
      * *H*    (``realtype``, input) -- current step size.
      * *IPAR* (``long int``, input) -- array containing integer user data that was passed to
        :f:func:`FARKMALLOC()`.
      * *RPAR* (``realtype``, input) -- array containing real user data that was passed to
        :f:func:`FARKMALLOC()`.
      * *IER*  (``int``, output) -- return flag  (0 if success, :math:`\ne 0` if an error).

   **Notes:**
   Typically this routine will use only *T* and *Y*, and store
   the results in either the arrays *IPAR* and *RPAR*, or in a Fortran
   module or common block.


If preconditioning is to be included, the following routine must be
supplied, for solution of the preconditioner linear system:


.. f:subroutine:: FARKPSOL(T,Y,FY,R,Z,GAMMA,DELTA,LR,IPAR,RPAR,VT,IER)

   User-supplied preconditioner solve routine (of type
   :c:func:`ARKLsPrecSolveFn()`).

   **Arguments:**
      * *T* (``realtype``, input) -- current value of the independent variable.
      * *Y* (``realtype``, input) -- current dependent state variable array.
      * *FY* (``realtype``, input) -- current dependent state variable derivative array.
      * *R* (``realtype``, input) -- right-hand side array.
      * *Z* (``realtype``, output) -- solution array.
      * *GAMMA* (``realtype``, input) -- Jacobian scaling factor.
      * *DELTA* (``realtype``, input) -- desired residual tolerance.
      * *LR* (``int``, input) -- flag denoting to solve the right or left preconditioner
        system: 1 = left preconditioner, 2 = right preconditioner.
      * *IPAR* (``long int``, input/output) -- array containing integer user data that was passed to
        :f:func:`FARKMALLOC()`.
      * *RPAR* (``realtype``, input/output) -- array containing real user data that was passed to
        :f:func:`FARKMALLOC()`.
      * *IER*  (``int``, output) -- return flag  (0 if success, >0 if a recoverable
        failure, <0 if a non-recoverable failure).

   **Notes:**
   Typically this routine will use only *T*, *Y*, *GAMMA*, *R*,
   *LR*, and *Z*.  It must solve the preconditioner linear system
   :math:`Pz = r`.  The preconditioner (or the product of the left and
   right preconditioners if both are nontrivial) should be an
   approximation to the matrix  :math:`M - \gamma J`, where
   :math:`M` is the system mass matrix, :math:`\gamma` is the input
   GAMMA, and :math:`J = \frac{\partial f^I}{\partial y}`.


If the user's preconditioner requires that any Jacobian related data be evaluated
or preprocessed, then the following routine can be used for the evaluation and
preprocessing of the preconditioner:

.. f:subroutine:: FARKPSET(T,Y,FY,JOK,JCUR,GAMMA,H,IPAR,RPAR,IER)

   User-supplied preconditioner setup routine (of type
   :c:func:`ARKLsPrecSetupFn()`).

   **Arguments:**
      * *T* (``realtype``, input) -- current value of the independent variable.
      * *Y* (``realtype``, input) -- current dependent state variable array.
      * *FY* (``realtype``, input) -- current dependent state variable derivative array.
      * *JOK* (``int``, input) -- flag indicating whether Jacobian-related data needs to be
        recomputed: 0 = recompute, 1 = reuse with the current value of *GAMMA*.
      * *JCUR* (``realtype``, output) -- return flag to denote if
	Jacobian data was recomputed (1=yes, 0=no).
      * *GAMMA* (``realtype``, input) -- Jacobian scaling factor.
      * *H* (``realtype``, input) -- current step size.
      * *IPAR* (``long int``, input/output) -- array containing integer user data that was passed to
        :f:func:`FARKMALLOC()`.
      * *RPAR* (``realtype``, input/output) -- array containing real user data that was passed to
        :f:func:`FARKMALLOC()`.
      * *IER*  (``int``, output) -- return flag  (0 if success, >0 if a recoverable
        failure, <0 if a non-recoverable failure).

   **Notes:**
   This routine must set up the preconditioner :math:`P` to be used in
   the subsequent call to :f:func:`FARKPSOL()`.  The preconditioner (or
   the product of the left and right preconditioners if using both)
   should be an approximation to the matrix  :math:`M - \gamma J`,
   where :math:`M` is the system mass matrix, :math:`\gamma` is the
   input *GAMMA*, and :math:`J = \frac{\partial f^I}{\partial y}`.


Notes:

(a) If the user's :f:func:`FARKJTSETUP()`, :f:func:`FARKJTIMES()` or
    :f:func:`FARKPSET()` routines use difference quotient
    approximations, they may need to use the error weight array *EWT*
    and/or the unit roundoff, in the calculation of suitable
    increments.  Also, if :f:func:`FARKPSOL()` uses an iterative
    method in its solution, the residual vector :math:`\rho = r - Pz`
    of the system should be made less than :math:`\delta =` *DELTA* in
    the weighted l2 norm, i.e.

    .. math::
       \left(\sum_i \left(\rho_i\, EWT_i\right)^2 \right)^{1/2} < \delta.

(b) If needed in :f:func:`FARKJTSETUP()` :f:func:`FARKJTIMES()`,
    :f:func:`FARKPSOL()`, or :f:func:`FARKPSET()`, the error weight
    array *EWT* can be obtained by calling
    :f:func:`FARKGETERRWEIGHTS()` using a user-allocated array as
    temporary storage for *EWT*.

(c) If needed in :f:func:`FARKJTSETUP()` :f:func:`FARKJTIMES()`,
    :f:func:`FARKPSOL()`, or :f:func:`FARKPSET()`, the unit roundoff
    can be obtained as the optional output *ROUT(6)* (available after
    the call to :f:func:`FARKMALLOC()`) and can be passed using either
    the *RPAR* user data array or a common block.



.. _FInterface.MassLinearSolver:

Mass matrix linear solver interface specification
----------------------------------------------------

To attach the mass matrix linear solver (and optionally the mass
matrix) object(s) initialized in steps :ref:`FInterface.SUNMatrix` and
:ref:`FInterface.SUNLinSol` above, the user of FARKODE must
initialize the mass-matrix linear solver interface.  To attach any
SUNLINSOL object (and optional SUNMATRIX object) to the mass-matrix
solver interface, following calls to initialize the SUNLINSOL (and
SUNMATRIX) object(s) in steps :ref:`FInterface.SUNMatrix` and
:ref:`FInterface.SUNLinSol` above, the user must call the
:f:func:`FARKLSMASSINIT()` routine:


.. f:subroutine:: FARKLSMASSINIT(TIME_DEP, IER)

   Interfaces with the :c:func:`ARKStepSetMassLinearSolver()` function to
   attach a linear solver object (and optionally a matrix object) to
   ARKStep's mass-matrix linear solver interface.

   **Arguments:**
      * *TIME_DEP* (``int``, input) -- flag indicating whether the
        mass matrix is time-dependent (1) or not (0).
        *Currently, only values of "0" are supported*
      * *IER* (``int``, output) -- return flag (0 if success, -1 if a memory allocation
        error occurred, -2 for an illegal input).




Matrix-based mass matrix linear solvers
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

When using the mass-matrix linear solver interface with the
SUNLINSOL_DENSE or SUNLINSOL_LAPACKDENSE mass matrix linear solver
modules, the user must supply a routine that computes the dense mass
matrix :math:`M`.  This routine must have the following form:


.. f:subroutine:: FARKDMASS(NEQ, T, DMASS, IPAR, RPAR, WK1, WK2, WK3, IER)

   Interface to provide a user-supplied dense mass matrix computation
   function (of type :c:func:`ARKLsMassFn()`), to be used by the
   SUNLINSOL_DENSE or SUNLINSOL_LAPACKDENSE solver modules.

   **Arguments:**
      * *NEQ* (``long int``, input) -- size of the ODE system.
      * *T* (``realtype``, input) -- current value of the independent variable.
      * *DMASS* (``realtype`` of size (NEQ,NEQ), output) -- 2D array
	containing the mass matrix entries.
      * *IPAR* (``long int``, input) -- array containing integer user data that was passed to
        :f:func:`FARKMALLOC()`.
      * *RPAR* (``realtype``, input) -- array containing real user data that was passed to
        :f:func:`FARKMALLOC()`.
      * *WK1*, *WK2*, *WK3*  (``realtype``, input) -- array containing temporary workspace
        of same size as *Y*.
      * *IER* (``int``, output) -- return flag (0 if success, >0 if a recoverable error
        occurred, <0 if an unrecoverable error occurred).

   **Notes:** Typically this routine will use only *NEQ*, *T*, and
   *DMASS*. It must compute the mass matrix and store it column-wise in *DMASS*.

To indicate that the :f:func:`FARKDMASS()` routine has been provided, then,
following the call to :f:func:`FARKLSMASSINIT()`, the user must call
the routine :f:func:`FARKDENSESETMASS()`:


.. f:subroutine:: FARKDENSESETMASS(IER)

   Interface to the :c:func:`ARKStepSetMassFn()` function,
   specifying to use the user-supplied routine :f:func:`FARKDMASS()`
   for the mass matrix calculation.

   **Arguments:**
      * *IER* (``int``, output) -- return flag (0 if success,
	:math:`\ne 0` if an error occurred).


When using the mass-matrix linear solver interface with the
SUNLINSOL_BAND or SUNLINSOL_LAPACKBAND mass matrix linear solver
modules, the user must supply a routine that computes the banded mass
matrix :math:`M`.  This routine must have the following form:

.. f:subroutine:: FARKBMASS(NEQ, MU, ML, MDIM, T, BMASS, IPAR, RPAR, WK1, WK2, WK3, IER)

   Interface to provide a user-supplied band mass matrix calculation
   function (of type :c:func:`ARKLsMassFn()`), to be used by the
   SUNLINSOL_BAND or SUNLINSOL_LAPACKBAND solver modules.

   **Arguments:**
      * *NEQ* (``long int``, input) -- size of the ODE system.
      * *MU*   (``long int``, input) -- upper half-bandwidth.
      * *ML*   (``long int``, input) -- lower half-bandwidth.
      * *MDIM* (``long int``, input) -- leading dimension of *BMASS* array.
      * *T*    (``realtype``, input) -- current value of the independent variable.
      * *BMASS* (``realtype`` of size *(MDIM,NEQ)*, output) -- 2D array
	containing the mass matrix entries.
      * *IPAR* (``long int``, input) -- array containing integer user data that was passed to
        :f:func:`FARKMALLOC()`.
      * *RPAR* (``realtype``, input) -- array containing real user data that was passed to
        :f:func:`FARKMALLOC()`.
      * *WK1*, *WK2*, *WK3*  (``realtype``, input) -- array containing temporary workspace
        of same size as *Y*.
      * *IER* (``int``, output) -- return flag (0 if success, >0 if a recoverable error
        occurred, <0 if an unrecoverable error occurred).

   **Notes:**
   Typically this routine will use only *NEQ*, *MU*, *ML*, *T*, and
   *BMASS*. It must load the *MDIM* by *N* array *BMASS* with the mass
   matrix at the current :math:`(t)` in band form.  Store in
   *BMASS(k,j)* the mass matrix element :math:`M_{i,j}` with
   *k = i - j + MU + 1* (or *k = 1, ..., ML+MU+1*) and *j = 1, ..., N*.


To indicate that the :f:func:`FARKBMASS()` routine has been provided, then,
following the call to :f:func:`FARKLSMASSINIT()`, the user must call the routine
:f:func:`FARKBANDSETMASS()`:

.. f:subroutine:: FARKBANDSETMASS(IER)

   Interface to the :c:func:`ARKStepSetMassFn()` function, specifying
   to use the user-supplied routine :f:func:`FARKBMASS()` for the mass
   matrix calculation.

   **Arguments:**
      * *IER* (``int``, output) -- return flag (0 if success,
	:math:`\ne 0` if an error occurred).


When using the mass-matrix linear solver interface with the
SUNLINSOL_KLU or SUNLINSOL_SUPERLUMT mass matrix linear solver
modules, the user must supply a routine that computes the sparse mass
matrix :math:`M`. Both the KLU and SuperLU_MT solver interfaces
support the compressed-sparse-column (CSC) and compressed-sparse-row
(CSR) matrix formats.  The desired format must have been specified to
the :f:func:`FSUNSPARSEMASSMATINIT()` function when initializing the
sparse mass matrix.  The user-provided routine to compute :math:`M`
must have the following form:


.. f:subroutine:: FARKSPMASS(T, N, NNZ, MDATA, MINDEXVALS, MINDEXPTRS, IPAR, RPAR, WK1, WK2, WK3, IER)

   Interface to provide a user-supplied sparse mass matrix approximation
   function (of type :c:func:`ARKLsMassFn()`), to be used by the
   SUNLINSOL_KLU or SUNLINSOL_SUPERLUMT solver modules.

   **Arguments:**
      * *T* (``realtype``, input) -- current value of the independent variable.
      * *N* (``sunindextype``, input) -- number of mass matrix rows and columns.
      * *NNZ* (``sunindextype``, input) -- allocated length of nonzero storage
	in mass matrix.
      * *MDATA* (``realtype`` of size NNZ, output) -- nonzero values
	in mass matrix.
      * *MINDEXVALS* (``sunindextype`` of size NNZ, output) -- row *[CSR: column]* indices for each
	nonzero mass matrix entry.
      * *MINDEXPTRS* (``sunindextype`` of size N+1, output) -- indices of where
	each column's *[CSR: row's]* nonzeros begin in data array; last entry points
	just past end of data values.
      * *IPAR* (``long int``, input) -- array containing integer user data that was passed to
        :f:func:`FARKMALLOC()`.
      * *RPAR* (``realtype``, input) -- array containing real user data that was passed to
        :f:func:`FARKMALLOC()`.
      * *WK1*, *WK2*, *WK3*  (``realtype``, input) -- array containing temporary workspace
        of same size as *Y*.
      * *IER* (``int``, output) -- return flag (0 if success, >0 if a recoverable error
        occurred, <0 if an unrecoverable error occurred).

   **Notes:** due to the internal storage format of the
   SUNMATRIX_SPARSE module, the matrix-specific integer parameters and
   arrays are all of type ``sunindextype`` -- the index precision
   (32-bit vs 64-bit signed integers) specified during the SUNDIALS
   build.  It is assumed that the user's Fortran codes are constructed
   to have matching type to how SUNDIALS was installed.


To indicate that the :f:func:`FARKSPMASS()` routine has been provided, then,
following the call to :f:func:`FARKLSMASSINIT()`, the user must call the routine
:f:func:`FARKSPARSESETMASS()`:


.. f:subroutine:: FARKSPARSESETMASS(IER)

   Interface to the :c:func:`ARKStepSetMassFn()` function,
   specifying that the user-supplied routine :f:func:`FARKSPMASS()` has
   been provided for the mass matrix calculation.

   **Arguments:**
      * *IER* (``int``, output) -- return flag (0 if success,
	:math:`\ne 0` if an error occurred).



Iterative mass matrix linear solvers
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

As described in the section :ref:`Mathematics.Error.Linear`, a user
may adjust the linear solver tolerance scaling factor
:math:`\epsilon_L`.  Fortran users may adjust this value for the mass
matrix linear solver by calling the function
:f:func:`FARKLSSETMASSEPSLIN()`:

.. f:subroutine:: FARKLSSETMASSEPSLIN(EPLIFAC, IER)

   Interface to the function :c:func:`ARKStepSetMassEpsLin()` to
   specify the linear solver tolerance scale factor :math:`\epsilon_L`
   for the mass matrix linear solver.

   This routine must be called *after* :f:func:`FARKLSMASSINIT()`.

   **Arguments:**
      * *EPLIFAC* (``realtype``, input) -- value to use for
        :math:`\epsilon_L`.  Passing a value of 0 indicates to use the
        default value (0.05).
      * *IER*  (``int``, output) -- return flag  (0 if success,
	:math:`\ne 0` if an error).


With treatment of the mass matrix linear systems by any of the Krylov
iterative solvers, there are two required user-supplied routines,
:f:func:`FARKMTSETUP()` and :f:func:`FARKMTIMES()`, and there are two
optional user-supplied routines, :f:func:`FARKMASSPSET()` and
:f:func:`FARKMASSPSOL()`. The specifications of these functions are given below.

The required routines when using a Krylov iterative mass matrix linear
solver perform setup and computation of the product of the system mass
matrix :math:`M` and a given vector :math:`v`.  The product routine
must have the following form:


.. f:subroutine:: FARKMTIMES(V, MV, T, IPAR, RPAR, IER)

   Interface to a user-supplied mass-matrix-times-vector product
   approximation function (corresponding to a C interface routine of
   type :c:func:`ARKLsMassTimesVecFn()`), to be used by one of the
   Krylov iterative linear solvers.

   **Arguments:**
      * *V*    (``realtype``, input) -- array containing the vector to multiply.
      * *MV*   (``realtype``, output) -- array containing resulting product vector.
      * *T*    (``realtype``, input) -- current value of the independent variable.
      * *IPAR* (``long int``, input) -- array containing integer user data that was passed to
        :f:func:`FARKMALLOC()`.
      * *RPAR* (``realtype``, input) -- array containing real user data that was passed to
        :f:func:`FARKMALLOC()`.
      * *IER*  (``int``, output) -- return flag  (0 if success, :math:`\ne 0` if an error).

   **Notes:**
   Typically this routine will use only *T*, *V*, and
   *MV*.  It must compute the product vector :math:`Mv`, where
   :math:`v` is given in *V*, and the product is stored in *MV*.


If the user's mass-matrix-times-vector product routine requires that
any mass matrix data be evaluated or preprocessed, then the following
routine can be used for the evaluation and preprocessing of this data:


.. f:subroutine:: FARKMTSETUP(T, IPAR, RPAR, IER)

   Interface to a user-supplied mass-matrix-times-vector setup
   function (corresponding to a C interface routine of type
   :c:func:`ARKLsMassTimesSetupFn()`).

   **Arguments:**
      * *T*    (``realtype``, input) -- current value of the independent variable.
      * *IPAR* (``long int``, input) -- array containing integer user data that was passed to
        :f:func:`FARKMALLOC()`.
      * *RPAR* (``realtype``, input) -- array containing real user data that was passed to
        :f:func:`FARKMALLOC()`.
      * *IER*  (``int``, output) -- return flag  (0 if success, :math:`\ne 0` if an error).

   **Notes:**
   Typically this routine will use only *T*, and store
   the results in either the arrays *IPAR* and *RPAR*, or in a Fortran
   module or common block.  If no mass matrix setup is needed, this
   routine should just set *IER* to 0 and return.


To indicate that these routines have been supplied by the user, then,
following the call to :f:func:`FARKLSMASSINIT()`, the user must
call the routine :f:func:`FARKLSSETMASS()`:

.. f:subroutine:: FARKLSSETMASS(IER)

   Interface to the function :c:func:`ARKStepSetMassTimes()` to
   specify use of the user-supplied mass-matrix-times-vector setup and
   product functions :f:func:`FARKMTSETUP()` and :f:func:`FARKMTIMES()`.

   This routine must be called *after* :f:func:`FARKLSMASSINIT()`.

   **Arguments:**
      * *IER*  (``int``, output) -- return flag  (0 if success,
	:math:`\ne 0` if an error).


Two optional user-supplied preconditioning routines may be supplied to
help accelerate convergence of the Krylov mass matrix linear solver.
If preconditioning was selected when enabling the Krylov solver
(i.e. the solver was set up with *IPRETYPE* :math:`\ne 0`), then the
user must also call the routine :f:func:`FARKLSSETMASSPREC()` with
*FLAG* :math:`\ne 0`:


.. f:subroutine:: FARKLSSETMASSPREC(FLAG, IER)

   Interface to the function :c:func:`ARKStepSetMassPreconditioner()` to
   specify use of the user-supplied preconditioner setup and solve
   functions, :f:func:`FARKMASSPSET()` and :f:func:`FARKMASSPSOL()`,
   respectively.

   This routine must be called *after* :f:func:`FARKLSMASSINIT()`.

   **Arguments:**
      * *FLAG* (``int``, input) -- flag denoting use of user-supplied
        preconditioning routines.
      * *IER*  (``int``, output) -- return flag  (0 if success,
	:math:`\ne 0` if an error).


In addition, the user must provide the following two routines to
implement the preconditioner setup and solve functions to be used
within the solve.


.. f:subroutine:: FARKMASSPSET(T,IPAR,RPAR,IER)

   User-supplied preconditioner setup routine (of type
   :c:func:`ARKLsMassPrecSetupFn()`).

   **Arguments:**
      * *T* (``realtype``, input) -- current value of the independent
	variable.
      * *IPAR* (``long int``, input/output) -- array containing
	integer user data that was passed to :f:func:`FARKMALLOC()`.
      * *RPAR* (``realtype``, input/output) -- array containing real
	user data that was passed to :f:func:`FARKMALLOC()`.
      * *IER*  (``int``, output) -- return flag  (0 if success, >0 if
	a recoverable failure, <0 if a non-recoverable failure).

   **Notes:**
   This routine must set up the preconditioner :math:`P` to be used in
   the subsequent call to :f:func:`FARKMASSPSOL()`.  The
   preconditioner (or the product of the left and right
   preconditioners if using both) should be an approximation to the
   system mass matrix, :math:`M`.


.. f:subroutine:: FARKMASSPSOL(T,R,Z,DELTA,LR,IPAR,RPAR,IER)

   User-supplied preconditioner solve routine (of type
   :c:func:`ARKLsMassPrecSolveFn()`).

   **Arguments:**
      * *T* (``realtype``, input) -- current value of the independent
	variable.
      * *R* (``realtype``, input) -- right-hand side array.
      * *Z* (``realtype``, output) -- solution array.
      * *DELTA* (``realtype``, input) -- desired residual tolerance.
      * *LR* (``int``, input) -- flag denoting to solve the right or
	left preconditioner system: 1 = left preconditioner, 2 = right
	preconditioner.
      * *IPAR* (``long int``, input/output) -- array containing
	integer user data that was passed to :f:func:`FARKMALLOC()`.
      * *RPAR* (``realtype``, input/output) -- array containing real
	user data that was passed to :f:func:`FARKMALLOC()`.
      * *IER*  (``int``, output) -- return flag  (0 if success, >0 if
	a recoverable failure, <0 if a non-recoverable failure).

   **Notes:**
   Typically this routine will use only *T*, *R*, *LR*, and *Z*.  It
   must solve the preconditioner linear system :math:`Pz = r`.  The
   preconditioner (or the product of the left and right
   preconditioners if both are nontrivial) should be an approximation
   to the system mass matrix :math:`M`.


Notes:

(a) If the user's :f:func:`FARKMASSPSOL()` uses an iterative method in
    its solution, the residual vector :math:`\rho = r - Pz` of the
    system should be made less than :math:`\delta =` *DELTA* in the
    weighted l2 norm, i.e.

    .. math::
       \left(\sum_i \left(\rho_i\, EWT_i\right)^2 \right)^{1/2} < \delta.

(b) If needed in :f:func:`FARKMTIMES()`, :f:func:`FARKMTSETUP()`,
    :f:func:`FARKMASSPSOL()`, or :f:func:`FARKMASSPSET()`, the error
    weight array *EWT* can be obtained by calling
    :f:func:`FARKGETERRWEIGHTS()` using a user-allocated array as
    temporary storage for *EWT*.

(c) If needed in :f:func:`FARKMTIMES()`, :f:func:`FARKMTSETUP()`,
    :f:func:`FARKMASSPSOL()`, or :f:func:`FARKMASSPSET()`, the unit
    roundoff can be obtained as the optional output *ROUT(6)*
    (available after the call to :f:func:`FARKMALLOC()`) and can be
    passed using either the *RPAR* user data array or a common block.





.. _FInterface.Solution:

Problem solution
-----------------------

Carrying out the integration is accomplished by making calls to
:f:func:`FARKODE()`.


.. f:subroutine:: FARKODE(TOUT, T, Y, ITASK, IER)

   Fortran interface to the C routine :c:func:`ARKStepEvolve()`
   for performing the solve, along with many of the ARK*Get*
   routines for reporting on solver statistics.

   **Arguments:**
      * *TOUT* (``realtype``, input) -- next value of :math:`t` at
	which a solution is desired.

      * *T* (``realtype``, output) -- value of independent variable
	that corresponds to the output *Y*

      * *Y* (``realtype``, output) -- array containing dependent state
	variables on output.

      * *ITASK* (``int``, input) -- task indicator :

        * 1 = normal mode (overshoot *TOUT* and interpolate)

        * 2 = one-step mode (return after each internal step taken)

        * 3 = normal 'tstop' mode (like 1, but integration never
          proceeds past *TSTOP*, which must be specified through a
          preceding call to :f:func:`FARKSETRIN()` using the key
          *STOP_TIME*)

        * 4 = one step 'tstop' mode (like 2, but integration never
	  goes past *TSTOP*).

      * *IER* (int, output) -- completion flag:

	* 0 = success,

	* 1 = tstop return,

	* 2 = root return,

	* values -1, ..., -10 are failure modes (see :c:func:`ARKStepEvolve()` and
          :ref:`Constants`).

   **Notes:**
   The current values of the optional outputs are immediately
   available in *IOUT* and *ROUT* upon return from this function (see
   :ref:`FInterface.IOUTTable` and :ref:`FInterface.ROUTTable`).

   A full description of error flags and output behavior of the solver
   (values filled in for *T* and *Y*) is provided in the description
   of :c:func:`ARKStepEvolve()`.




.. _FInterface.AdditionalOutput:

Additional solution output
---------------------------------------

After a successful return from :f:func:`FARKODE()`, the routine
:f:func:`FARKDKY()` may be used to obtain a derivative of the solution,
of order up to 3, at any :math:`t` within the last step taken.


.. f:subroutine:: FARKDKY(T, K, DKY, IER)

   Fortran interface to the C routine :f:func:`ARKDKY()` for
   interpolating output of the solution or its derivatives at any
   point within the last step taken.

   **Arguments:**
      * *T* (``realtype``, input) -- time at which solution derivative
	is desired, within the interval :math:`[t_n-h,t_n]`.
      * *K* (``int``, input) -- derivative order :math:`(0 \le k \le 3)`.
      * *DKY* (``realtype``, output) -- array containing the computed
	*K*-th derivative of :math:`y`.
      * *IER* (``int``, output) -- return flag (0 if success, <0 if an
	illegal argument).



.. _FInterface.ReInit:

Problem reinitialization
---------------------------------------

To re-initialize the ARKStep solver for the solution of a new
problem of the same size as one already solved, the user must call
:f:func:`FARKREINIT()`:


.. f:subroutine:: FARKREINIT(T0, Y0, IMEX, IATOL, RTOL, ATOL, IER)

   Re-initializes the Fortran interface to the ARKStep solver.

   **Arguments:**  The arguments have the same names and meanings as those of
   :f:func:`FARKMALLOC()`.

   **Notes:**
   This routine performs no memory allocation, instead using the
   existing memory created by the previous :f:func:`FARKMALLOC()`
   call.  The call to specify the linear system solution method may
   or may not be needed.


Following a call to :f:func:`FARKREINIT()` if the choice of linear
solver is being changed then a user must make a call to create the
alternate SUNLINSOL module and then attach it to ARKStep, as shown
above.  If only linear solver parameters are being modified, then
these calls may be made without re-attaching to ARKStep.



.. _FInterface.Resize:

Resizing the ODE system
-----------------------------

For simulations involving changes to the number of equations and
unknowns in the ODE system (e.g. when solving a spatially-adaptive
PDE), the :f:func:`FARKODE()` integrator may be "resized" between
integration steps, through calls to the :f:func:`FARKRESIZE()`
function, that interfaces with the C routine :c:func:`ARKStepResize()`.
This function modifies ARKStep's internal memory structures to use the
new problem size, without destruction of the temporal adaptivity
heuristics.  It is assumed that the dynamical time scales before and
after the vector resize will be comparable, so that all time-stepping
heuristics prior to calling :c:func:`FARKRESIZE` remain valid after
the call.  If instead the dynamics should be re-calibrated, the
FARKODE memory structure should be deleted with a call to
:f:func:`FARKFREE()`, and re-created with a call to
:f:func:`FARKMALLOC()`.


.. f:subroutine:: FARKRESIZE(T0, Y0, HSCALE, ITOL, RTOL, ATOL, IER)

   Re-initializes the Fortran interface to the ARKStep solver for a
   differently-sized ODE system.

   **Arguments:**
      * *T0* (``realtype``, input) -- initial value of the independent
	variable :math:`t`.

      * *Y0* (``realtype``, input) -- array of dependent-variable
	initial conditions.

      * *HSCALE* (``realtype``, input) -- desired step size scale factor:

        * 1.0 is the default,

        * any value <= 0.0 results in the default.

      * *ITOL* (``int``, input) -- flag denoting that a new relative
	tolerance and vector of absolute tolerances are supplied in
	the *RTOL* and *ATOL* arguments:

        * 0 = retain the current scalar-valued relative and absolute
	  tolerances, or the user-supplied error weight function,
	  :f:func:`FARKEWT()`.

        * 1 = *RTOL* contains the new scalar-valued relative tolerance
          and *ATOL* contains a new array of absolute tolerances.

      * *RTOL* (``realtype``, input) -- scalar relative tolerance.

      * *ATOL* (``realtype``, input) -- array of absolute tolerances.

      * *IER* (``int``, output) -- return flag (0 success, :math:`\ne 0` failure).

   **Notes:**
   This routine performs the opposite set of of operations as
   :f:func:`FARKREINIT()`: it does not reinitialize any of the
   time-step heuristics, but it does perform memory reallocation.


Following a call to :f:func:`FARKRESIZE()`, the internal data
structures for all linear solver and matrix objects will be the
incorrect size.  Hence, calls must be made to re-create the linear
system solver, mass matrix solver, linear system matrix, and mass
matrix, followed by calls to attach the updated objects to ARKStep.

If any user-supplied linear solver helper routines were used (Jacobian
evaluation, Jacobian-vector product, mass matrix evaluation,
mass-matrix-vector product, preconditioning, etc.), then the
relevant "set" routines to specify their usage must be called again
**following** the re-specification of the linear solver module(s).





.. _FInterface.Deallocation:

Memory deallocation
---------------------------------------

To free the internal memory created by :f:func:`FARKMALLOC()`,
:f:func:`FARKLSINIT()`, :f:func:`FARKLSMASSINIT()`, and the SUNMATRIX,
SUNLINSOL and SUNNONLINSOL objects, the user may call
:f:func:`FARKFREE()`, as follows:


.. f:subroutine:: FARKFREE()

   Frees the internal memory created by :f:func:`FARKMALLOC()`.

   **Arguments:** None.
