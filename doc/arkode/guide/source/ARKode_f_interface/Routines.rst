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


.. _FInterface.Routines:

FARKODE routines
===========================

In this section, we list the full set of user-callable functions
comprising the FARKODE solver interface.  For each function, we list
the corresponding ARKStep functions, to provide a mapping between the
two solver interfaces.  Further documentation on each FARKODE function
is provided in the following sections, :ref:`FInterface.Usage`,
:ref:`FInterface.OptionalOutputs`, :ref:`FInterface.Rootfinding` and
:ref:`FInterface.Preconditioning`.  Additionally, all Fortran and C
functions below are hyperlinked to their definitions in the
documentation, for simplified access.



Interface to the NVECTOR modules
----------------------------------

* :f:func:`FNVINITS()` (defined by NVECTOR_SERIAL) interfaces to
  :c:func:`N_VNewEmpty_Serial()`.

* :f:func:`FNVINITP()` (defined by NVECTOR_PARALLEL) interfaces to
  :c:func:`N_VNewEmpty_Parallel()`.

* :f:func:`FNVINITOMP()` (defined by NVECTOR_OPENMP) interfaces to
  :c:func:`N_VNewEmpty_OpenMP()`.

* :f:func:`FNVINITPTS()` (defined by NVECTOR_PTHREADS) interfaces to
  :c:func:`N_VNewEmpty_Pthreads()`.

* :f:func:`FNVINITPH()` (defined by NVECTOR_PARHYP) interfaces to
  :c:func:`N_VNewEmpty_ParHyp()`.


Interface to the SUNMATRIX modules
---------------------------------------

* :f:func:`FSUNBANDMATINIT()` (defined by SUNMATRIX_BAND) interfaces
  to :c:func:`SUNBandMatrix()`.

* :f:func:`FSUNDENSEMATINIT()` (defined by SUNMATRIX_DENSE) interfaces
  to :c:func:`SUNDenseMatrix()`.

* :f:func:`FSUNSPARSEMATINIT()` (defined by SUNMATRIX_SPARSE) interfaces
  to :c:func:`SUNSparseMatrix()`.

  
Interface to the SUNLINSOL modules
------------------------------------------

* :f:func:`FSUNBANDLINSOLINIT()` (defined by SUNLINSOL_BAND)
  interfaces to :c:func:`SUNLinSol_Band()`.

* :f:func:`FSUNDENSELINSOLINIT()` (defined by SUNLINSOL_DENSE)
  interfaces to :c:func:`SUNLinSol_Dense()`.

* :f:func:`FSUNKLUINIT()` (defined by SUNLINSOL_KLU)
  interfaces to :c:func:`SUNLinSol_KLU()`.

* :f:func:`FSUNKLUREINIT()` (defined by SUNLINSOL_KLU)
  interfaces to :c:func:`SUNLinSol_KLUReinit()`.

* :f:func:`FSUNLAPACKBANDINIT()` (defined by SUNLINSOL_LAPACKBAND)
  interfaces to :c:func:`SUNLinSol_LapackBand()`.

* :f:func:`FSUNLAPACKDENSEINIT()` (defined by SUNLINSOL_LAPACKDENSE)
  interfaces to :c:func:`SUNLinSol_LapackDense()`.

* :f:func:`FSUNPCGINIT()` (defined by SUNLINSOL_PCG)
  interfaces to :c:func:`SUNLinSol_PCG()`.

* :f:func:`FSUNSPBCGSINIT()` (defined by SUNLINSOL_SPBCGS)
  interfaces to :c:func:`SUNLinSol_SPBCGS()`.

* :f:func:`FSUNSPFGMRINIT()` (defined by SUNLINSOL_SPFGMR)
  interfaces to :c:func:`SUNLinSol_SPFGMR()`.

* :f:func:`FSUNSPGMRINIT()` (defined by SUNLINSOL_SPGMR)
  interfaces to :c:func:`SUNLinSol_SPGMR()`.

* :f:func:`FSUNSPTFQMRINIT()` (defined by SUNLINSOL_SPTFQMR)
  interfaces to :c:func:`SUNLinSol_SPTFQMR()`.

* :f:func:`FSUNSUPERLUMTINIT()` (defined by SUNLINSOL_SUPERLUMT)
  interfaces to :c:func:`SUNLinSol_SuperLUMT()`.



  
Interface to the SUNNONLINSOL modules
------------------------------------------

* :f:func:`FSUNNEWTONINIT()` (defined by SUNNONLINSOL_NEWTON)
  interfaces to :c:func:`SUNNonlinSol_Newton()`.

* :f:func:`FSUNNEWTONSETMAXITERS()` (defined by SUNNONLINSOL_NEWTON)
  interfaces to :c:func:`SUNNonlinSolSetMaxIters()` for a
  SUNNONLINSOL_NEWTON object.

* :f:func:`FSUNFIXEDPOINTINIT()` (defined by SUNNONLINSOL_FIXEDPOINT)
  interfaces to :c:func:`SUNNonlinSol_Newton()`.

* :f:func:`FSUNFIXEDPOINTSETMAXITERS()` (defined by SUNNONLINSOL_FIXEDPOINT)
  interfaces to :c:func:`SUNNonlinSolSetMaxIters()` for a
  SUNNONLINSOL_FIXEDPOINT object.



  
Interface to the main ARKODE module
--------------------------------------

* :f:func:`FARKMALLOC()` interfaces to :c:func:`ARKStepCreate()` and
  :c:func:`ARKStepSetUserData()`, as well as one of :c:func:`ARKStepSStolerances()` or :c:func:`ARKStepSVtolerances()`.

* :f:func:`FARKREINIT()` interfaces to :c:func:`ARKStepReInit()`.

* :f:func:`FARKRESIZE()` interfaces to :c:func:`ARKStepResize()`.

* :f:func:`FARKSETIIN()` and :f:func:`FARKSETRIN()` interface to the
  ARKStepSet* and ARKStepSet* functions (see :ref:`ARKStep_CInterface.OptionalInputs`).

* :f:func:`FARKEWTSET()` interfaces to :c:func:`ARKStepWFtolerances()`.

* :f:func:`FARKADAPTSET()` interfaces to :c:func:`ARKStepSetAdaptivityFn()`.

* :f:func:`FARKEXPSTABSET()` interfaces to :c:func:`ARKStepSetStabilityFn()`.

* :f:func:`FARKSETERKTABLE()` interfaces to :c:func:`ARKStepSetTables()`.

* :f:func:`FARKSETIRKTABLE()` interfaces to :c:func:`ARKStepSetTables()`.

* :f:func:`FARKSETARKTABLES()` interfaces to :c:func:`ARKStepSetTables()`.

* :f:func:`FARKSETRESTOLERANCE()` interfaces to either
  :c:func:`ARKStepResStolerance()` and :c:func:`ARKStepResVtolerance()`

..
   * :f:func:`FARKSETDIAGNOSTICS()` interfaces to :c:func:`ARKStepSetDiagnostics()`.

* :f:func:`FARKODE()` interfaces to :c:func:`ARKStepEvolve()`, the
  ARKStepGet* functions (see :ref:`ARKStep_CInterface.OptionalOutputs`),
  and to the optional output functions for the selected linear
  solver module (see :ref:`ARKStep_CInterface.OptionalOutputs`).

* :f:func:`FARKDKY()` interfaces to the interpolated output function
  :c:func:`ARKStepGetDky()`.

* :f:func:`FARKGETERRWEIGHTS()` interfaces to
  :c:func:`ARKStepGetErrWeights()`.

* :f:func:`FARKGETESTLOCALERR()` interfaces to
  :c:func:`ARKStepGetEstLocalErrors()`.

* :f:func:`FARKFREE()` interfaces to :c:func:`ARKStepFree()`.



Interface to the system nonlinear solver interface
----------------------------------------------------

* :f:func:`FARKNLSINIT()` interfaces to :c:func:`ARKStepSetNonlinearSolver()`.


     
Interface to the system linear solver interfaces
--------------------------------------------------

* :f:func:`FARKLSINIT()` interfaces to :c:func:`ARKStepSetLinearSolver()`.

* :f:func:`FARKDENSESETJAC()` interfaces to :c:func:`ARKStepSetJacFn()`.

* :f:func:`FARKBANDSETJAC()` interfaces to :c:func:`ARKStepSetJacFn()`.

* :f:func:`FARKSPARSESETJAC()` interfaces to :c:func:`ARKStepSetJacFn()`.

* :f:func:`FARKLSSETEPSLIN()` interfaces to :c:func:`ARKStepSetEpsLin()`.

* :f:func:`FARKLSSETJAC()` interfaces to :c:func:`ARKStepSetJacTimes()`.

* :f:func:`FARKLSSETPREC()` interfaces to :c:func:`ARKStepSetPreconditioner()`.



Interface to the mass matrix linear solver interfaces
-------------------------------------------------------

* :f:func:`FARKLSMASSINIT()` interfaces to :c:func:`ARKStepSetMassLinearSolver()`.

* :f:func:`FARKDENSESETMASS()` interfaces to :c:func:`ARKStepSetMassFn()`.

* :f:func:`FARKBANDSETMASS()` interfaces to :c:func:`ARKStepSetMassFn()`.

* :f:func:`FARKSPARSESETMASS()` interfaces to :c:func:`ARKStepSetMassFn()`.

* :f:func:`FARKLSSETMASSEPSLIN()` interfaces to :c:func:`ARKStepSetMassEpsLin()`.

* :f:func:`FARKLSSETMASS()` interfaces to :c:func:`ARKStepSetMassTimes()`.

* :f:func:`FARKLSSETMASSPREC()` interfaces to :c:func:`ARKStepSetMassPreconditioner()`.



.. _FInterface.UserSupplied:


User-supplied routines
---------------------------------------

As with the native C interface, the FARKODE solver interface requires
user-supplied functions to specify the ODE problem to be solved.  In
contrast to the case of direct use of ARKStep, and of most Fortran ODE
solvers, the names of all user-supplied routines here are fixed, in
order to maximize portability for the resulting mixed-language program.
As a result, whether using a purely implicit, purely explicit, or
mixed implicit-explicit solver, routines for both :math:`f^E(t,y)` and
:math:`f^I(t,y)` must be provided by the user (though either of which
may do nothing):

.. cssclass:: table-bordered

+---------------------------+-----------------------------------+
| FARKODE routine           | ARKStep interface                 |
| (FORTRAN, user-supplied)  | function type                     |
+===========================+===================================+
| :f:func:`FARKIFUN()`      | :c:func:`ARKRhsFn()`              |
+---------------------------+-----------------------------------+
| :f:func:`FARKEFUN()`      | :c:func:`ARKRhsFn()`              |
+---------------------------+-----------------------------------+

In addition, as with the native C interface a user may provide
additional routines to assist in the solution process.  Each of the
following user-supplied routines is activated by calling the specified
"activation" routine, with the exception of :f:func:`FARKSPJAC()`
which is required whenever a sparse matrix solver is used:

.. cssclass:: table-bordered

+--------------------------+-----------------------------------+-------------------------------+
| FARKODE routine          | ARKStep interface                 | FARKODE "activation" routine  |
| (FORTRAN, user-supplied) | function type                     |                               |
+==========================+===================================+===============================+
| :f:func:`FARKDJAC()`     | :c:func:`ARKLsJacFn()`            | :f:func:`FARKDENSESETJAC()`   |
+--------------------------+-----------------------------------+-------------------------------+
| :f:func:`FARKBJAC()`     | :c:func:`ARKLsJacFn()`            | :f:func:`FARKBANDSETJAC()`    |
+--------------------------+-----------------------------------+-------------------------------+
| :f:func:`FARKSPJAC()`    | :c:func:`ARKLsJacFn()`            | :f:func:`FARKSPARSESETJAC()`  |
+--------------------------+-----------------------------------+-------------------------------+
| :f:func:`FARKDMASS()`    | :c:func:`ARKLsMassFn()`           | :f:func:`FARKDENSESETMASS()`  |
+--------------------------+-----------------------------------+-------------------------------+
| :f:func:`FARKBMASS()`    | :c:func:`ARKLsMassFn()`           | :f:func:`FARKBANDSETMASS()`   |
+--------------------------+-----------------------------------+-------------------------------+
| :f:func:`FARKSPMASS()`   | :c:func:`ARKLsMassFn()`           | :f:func:`FARKSPARSESETMASS()` |
+--------------------------+-----------------------------------+-------------------------------+
| :f:func:`FARKPSET()`     | :c:func:`ARKLsPrecSetupFn()`      | :f:func:`FARKLSSETPREC()`     |
+--------------------------+-----------------------------------+-------------------------------+
| :f:func:`FARKPSOL()`     | :c:func:`ARKLsPrecSolveFn()`      | :f:func:`FARKLSSETPREC()`     |
+--------------------------+-----------------------------------+-------------------------------+
| :f:func:`FARKJTSETUP()`  | :c:func:`ARKLsJacTimesSetupFn()`  | :f:func:`FARKLSSETJAC()`      |
+--------------------------+-----------------------------------+-------------------------------+
| :f:func:`FARKJTIMES()`   | :c:func:`ARKLsJacTimesVecFn()`    | :f:func:`FARKLSSETJAC()`      |
+--------------------------+-----------------------------------+-------------------------------+
| :f:func:`FARKMASSPSET()` | :c:func:`ARKLsMassPrecSetupFn()`  | :f:func:`FARKLSSETMASSPREC()` |
+--------------------------+-----------------------------------+-------------------------------+
| :f:func:`FARKMASSPSOL()` | :c:func:`ARKLsMassPrecSolveFn()`  | :f:func:`FARKLSSETMASSPREC()` |
+--------------------------+-----------------------------------+-------------------------------+
| :f:func:`FARKMTSETUP()`  | :c:func:`ARKLsMassTimesSetupFn()` | :f:func:`FARKLSSETMASS()`     |
+--------------------------+-----------------------------------+-------------------------------+
| :f:func:`FARKMTIMES()`   | :c:func:`ARKLsMassTimesVecFn()`   | :f:func:`FARKLSSETMASS()`     |
+--------------------------+-----------------------------------+-------------------------------+
| :f:func:`FARKEWT()`      | :c:func:`ARKEwtFn()`              | :f:func:`FARKEWTSET()`        |
+--------------------------+-----------------------------------+-------------------------------+
| :f:func:`FARKADAPT()`    | :c:func:`ARKAdaptFn()`            | :f:func:`FARKADAPTSET()`      |
+--------------------------+-----------------------------------+-------------------------------+
| :f:func:`FARKEXPSTAB()`  | :c:func:`ARKExpStabFn()`          | :f:func:`FARKEXPSTABSET()`    |
+--------------------------+-----------------------------------+-------------------------------+
