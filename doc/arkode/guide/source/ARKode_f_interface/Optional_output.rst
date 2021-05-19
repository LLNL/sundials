..
   Programmer(s): Daniel R. Reynolds @ SMU
   ----------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2002-2021, Lawrence Livermore National Security
   and Southern Methodist University.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   ----------------------------------------------------------------

:tocdepth: 3


.. _FInterface.OptionalOutputs:

FARKODE optional output
==============================

We note that the optional inputs to FARKODE have already been
described in the section :ref:`FInterface.OptionalInputs`.



IOUT and ROUT arrays
----------------------------

In the Fortran interface, the optional outputs from the
:f:func:`FARKODE()` solver are accessed not through individual
functions, but rather through a pair of user-allocated arrays, *IOUT*
(having ``long int`` type) of dimension at least 36, and *ROUT*
(having ``realtype`` type) of dimension at least 6.  These arrays must
be allocated by the user program that calls :f:func:`FARKODE()`, that
passes them through the Fortran interface as arguments to
:f:func:`FARKMALLOC()`.  Following this call, :f:func:`FARKODE()` will
modify the entries of these arrays to contain all optional output
values provided to a Fortran user.

In the following tables, :ref:`FInterface.IOUTTable` and
:ref:`FInterface.ROUTTable`, we list the entries in these
arrays by index, naming them according to their role with the main
ARKStep solver, and list the relevant ARKStep C/C++ function that is
actually called to extract the output value.  Similarly, optional
integer output values that are specific to the ARKLS linear solver
interface are listed in :ref:`FInterface.LsIOUTTable`.

For more details on the optional inputs and outputs to ARKStep, see
the sections :ref:`ARKStep_CInterface.OptionalInputs` and
:ref:`ARKStep_CInterface.OptionalOutputs`.



.. _FInterface.IOUTTable:

Table: Optional FARKODE integer outputs
----------------------------------------

.. cssclass:: table-bordered

==============  ===============  =========================================================
*IOUT* Index    Optional output  ARKStep function
==============  ===============  =========================================================
1               LENRW            :c:func:`ARKStepGetWorkSpace()`
2               LENIW            :c:func:`ARKStepGetWorkSpace()`
3               NST              :c:func:`ARKStepGetNumSteps()`
4               NST_STB          :c:func:`ARKStepGetNumExpSteps()`
5               NST_ACC          :c:func:`ARKStepGetNumAccSteps()`
6               NST_ATT          :c:func:`ARKStepGetNumStepAttempts()`
7               NFE              :c:func:`ARKStepGetNumRhsEvals()` (num :math:`f^E` calls)
8               NFI              :c:func:`ARKStepGetNumRhsEvals()` (num :math:`f^I` calls)
9               NSETUPS          :c:func:`ARKStepGetNumLinSolvSetups()`
10              NETF             :c:func:`ARKStepGetNumErrTestFails()`
11              NNI              :c:func:`ARKStepGetNumNonlinSolvIters()`
12              NCFN             :c:func:`ARKStepGetNumNonlinSolvConvFails()`
13              NGE              :c:func:`ARKStepGetNumGEvals()`
==============  ===============  =========================================================



.. _FInterface.ROUTTable:

Table: Optional FARKODE real outputs
-------------------------------------

.. cssclass:: table-bordered

==============  ===============  =======================================================================
*ROUT* Index    Optional output  ARKStep function
==============  ===============  =======================================================================
1               H0U              :c:func:`ARKStepGetActualInitStep()`
2               HU               :c:func:`ARKStepGetLastStep()`
3               HCUR             :c:func:`ARKStepGetCurrentStep()`
4               TCUR             :c:func:`ARKStepGetCurrentTime()`
5               TOLSF            :c:func:`ARKStepGetTolScaleFactor()`
6               UROUND           ``UNIT_ROUNDOFF`` (see the section :ref:`ARKStep_CInterface.DataTypes`)
==============  ===============  =======================================================================



.. _FInterface.LsIOUTTable:

Table: Optional ARKLS interface outputs
----------------------------------------

.. cssclass:: table-bordered

==============  ===============  ===================================================
*IOUT* Index    Optional output  ARKStep function
==============  ===============  ===================================================
14              LENRWLS          :c:func:`ARKLsGetWorkSpace()`
15              LENIWLS          :c:func:`ARKLsGetWorkSpace()`
16              LSTF             :c:func:`ARKLsGetLastFlag()`
17              NFELS            :c:func:`ARKLsGetNumRhsEvals()`
18              NJE              :c:func:`ARKLsGetNumJacEvals()`
19              NJTS             :c:func:`ARKLsGetNumJTSetupEvals()`
20              NJTV             :c:func:`ARKLsGetNumJtimesEvals()`
21              NPE              :c:func:`ARKLsGetNumPrecEvals()`
22              NPS              :c:func:`ARKLsGetNumPrecSolves()`
23              NLI              :c:func:`ARKLsGetNumLinIters()`
24              NCFL             :c:func:`ARKLsGetNumConvFails()`
==============  ===============  ===================================================



.. _FInterface.LsMassIOUTTable:

Table: Optional ARKLS mass interface outputs
---------------------------------------------

.. cssclass:: table-bordered

==============  ===============  ===================================================
*IOUT* Index    Optional output  ARKStep function
==============  ===============  ===================================================
25              LENRWMS          :c:func:`ARKLsGetMassWorkSpace()`
26              LENIWMS          :c:func:`ARKLsGetMassWorkSpace()`
27              LSTMF            :c:func:`ARKLsGetLastMassFlag()`
28              NMSET            :c:func:`ARKLsGetNumMassSetups()`
29              NMSOL            :c:func:`ARKLsGetNumMassSolves()`
30              NMTSET           :c:func:`ARKLsGetNumMTSetups()`
31              NMMUL            :c:func:`ARKLsGetNumMassMult()`
32              NMPE             :c:func:`ARKLsGetNumMassPrecEvals()`
33              NMPS             :c:func:`ARKLsGetNumMassPrecSolves()`
34              NMLI             :c:func:`ARKLsGetNumMassIters()`
35              NMCFL            :c:func:`ARKLsGetNumMassConvFails()`
==============  ===============  ===================================================



.. _FInterface.ConstrIOUTTable:

Table: Optional ARKode constraints outputs
-------------------------------------------

.. cssclass:: table-bordered

==============  ===============  ===================================================
*IOUT* Index    Optional output  ARKStep function
==============  ===============  ===================================================
36              CONSTRFAILS      :c:func:`ARKStepGetNumConstrFails()`
==============  ===============  ===================================================



Additional optional output routines
---------------------------------------------

In addition to the optional inputs communicated through FARKSET*
calls and the optional outputs extracted from *IOUT* and *ROUT*,
the following user-callable routines are available.


To obtain the error weight array *EWT*, containing the
multiplicative error weights used in the WRMS norms, the user may call
the routine :f:func:`FARKGETERRWEIGHTS()` as follows:


.. f:subroutine:: FARKGETERRWEIGHTS(EWT, IER)

   Retrieves the current error weight vector (interfaces
   with :c:func:`ARKStepGetErrWeights()`).

   **Arguments:**
      * *EWT* (``realtype``, output) -- array containing the error
	weight vector.
      * *IER*  (``int``, output) -- return flag  (0 if success,
	:math:`\ne 0` if an error).

   **Notes:**
   The array *EWT* must have already been allocated by the user, of
   the same size as the solution array *Y*.



Similarly, to obtain the estimated local truncation errors, following
a successful call to :f:func:`FARKODE()`, the user may call the
routine :f:func:`FARKGETESTLOCALERR()` as follows:


.. f:subroutine:: FARKGETESTLOCALERR(ELE, IER)

   Retrieves the current local truncation error estimate
   vector (interfaces with :c:func:`ARKStepGetEstLocalErrors()`).

   **Arguments:**
      * *ELE* (``realtype``, output) -- array with the estimated local
	truncation error vector.
      * *IER*  (``int``, output) -- return flag  (0 if success,
	:math:`\ne 0` if an error).

   **Notes:**
   The array *ELE* must have already been allocated by the user, of
   the same size as the solution array *Y*.
