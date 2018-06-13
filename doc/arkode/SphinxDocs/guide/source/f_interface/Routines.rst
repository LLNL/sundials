..
   Programmer(s): Daniel R. Reynolds @ SMU
   ----------------------------------------------------------------
   Copyright (c) 2013, Southern Methodist University.
   All rights reserved.
   For details, see the LICENSE file.
   ----------------------------------------------------------------

:tocdepth: 3


.. _FInterface.Routines:

FARKODE routines
===========================

In this section, we list the full set of user-callable functions
comprising the FARKODE solver interface.  For each function, we list
the corresponding ARKode functions, to provide a mapping between the
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
  interfaces to :c:func:`SUNBandLinearSolver()`.

* :f:func:`FSUNDENSELINSOLINIT()` (defined by SUNLINSOL_DENSE)
  interfaces to :c:func:`SUNDenseLinearSolver()`.

* :f:func:`FSUNKLUINIT()` (defined by SUNLINSOL_KLU)
  interfaces to :c:func:`SUNKLU()`.

* :f:func:`FSUNKLUREINIT()` (defined by SUNLINSOL_KLU)
  interfaces to :c:func:`SUNKLUReinit()`.

* :f:func:`FSUNLAPACKBANDINIT()` (defined by SUNLINSOL_LAPACKBAND)
  interfaces to :c:func:`SUNLapackBand()`.

* :f:func:`FSUNLAPACKDENSEINIT()` (defined by SUNLINSOL_LAPACKDENSE)
  interfaces to :c:func:`SUNLapackDense()`.

* :f:func:`FSUNPCGINIT()` (defined by SUNLINSOL_PCG)
  interfaces to :c:func:`SUNPCG()`.

* :f:func:`FSUNSPBCGSINIT()` (defined by SUNLINSOL_SPBCGS)
  interfaces to :c:func:`SUNSPBCGS()`.

* :f:func:`FSUNSPFGMRINIT()` (defined by SUNLINSOL_SPFGMR)
  interfaces to :c:func:`SUNSPFGMR()`.

* :f:func:`FSUNSPGMRINIT()` (defined by SUNLINSOL_SPGMR)
  interfaces to :c:func:`SUNSPGMR()`.

* :f:func:`FSUNSPTFQMRINIT()` (defined by SUNLINSOL_SPTFQMR)
  interfaces to :c:func:`SUNSPTFQMR()`.

* :f:func:`FSUNSUPERLUMTINIT()` (defined by SUNLINSOL_SUPERLUMT)
  interfaces to :c:func:`SUNSuperLUMT()`.

     

Interface to the main ARKODE module
--------------------------------------

* :f:func:`FARKMALLOC()` interfaces to :c:func:`ARKodeCreate()`,
  :c:func:`ARKodeSetUserData()`, and :c:func:`ARKStepCreate()`, as well
  as one of :c:func:`ARKodeSStolerances()` or :c:func:`ARKodeSVtolerances()`.

* :f:func:`FARKREINIT()` interfaces to :c:func:`ARKStepReInit()`.

* :f:func:`FARKRESIZE()` interfaces to :c:func:`ARKodeResize()`.

* :f:func:`FARKSETIIN()` and :f:func:`FARKSETRIN()` interface to the
  ARKodeSet* and ARKStepSet* functions (see :ref:`CInterface.OptionalInputs`).

* :f:func:`FARKEWTSET()` interfaces to :c:func:`ARKodeWFtolerances()`.

* :f:func:`FARKADAPTSET()` interfaces to :c:func:`ARKStepSetAdaptivityFn()`.

* :f:func:`FARKEXPSTABSET()` interfaces to :c:func:`ARKStepSetStabilityFn()`.

* :f:func:`FARKSETERKTABLE()` interfaces to :c:func:`ARKStepSetERKTable()`.

* :f:func:`FARKSETIRKTABLE()` interfaces to :c:func:`ARKStepSetIRKTable()`.

* :f:func:`FARKSETARKTABLES()` interfaces to :c:func:`ARKStepSetARKTables()`.

* :f:func:`FARKSETRESTOLERANCE()` interfaces to either
  :c:func:`ARKodeResStolerance()` and :c:func:`ARKodeResVtolerance()`

..
   * :f:func:`FARKSETDIAGNOSTICS()` interfaces to :c:func:`ARKodeSetDiagnostics()`.

* :f:func:`FARKODE()` interfaces to :c:func:`ARKode()`, the
  ARKodeGet* functions (see :ref:`CInterface.OptionalOutputs`), 
  and to the optional output functions for the selected linear
  solver module (see :ref:`CInterface.OptionalOutputs`). 

* :f:func:`FARKDKY()` interfaces to the interpolated output function
  :c:func:`ARKodeGetDky()`.

* :f:func:`FARKGETERRWEIGHTS()` interfaces to
  :c:func:`ARKodeGetErrWeights()`.

* :f:func:`FARKGETESTLOCALERR()` interfaces to
  :c:func:`ARKStepGetEstLocalErrors()`.

* :f:func:`FARKFREE()` interfaces to :c:func:`ARKodeFree()`.



Interface to the system linear solver interfaces
--------------------------------------------------

* :f:func:`FARKDLSINIT()` interfaces to :c:func:`ARKDlsSetLinearSolver()`.

* :f:func:`FARKDENSESETJAC()` interfaces to :c:func:`ARKDlsSetJacFn()`.

* :f:func:`FARKBANDSETJAC()` interfaces to :c:func:`ARKDlsSetJacFn()`.

* :f:func:`FARKSPARSESETJAC()` interfaces to :c:func:`ARKDlsSetJacFn()`.

* :f:func:`FARKSPILSINIT()` interfaces to :c:func:`ARKSpilsSetLinearSolver()`

* :f:func:`FARKSPILSSETEPSLIN()` interfaces to :c:func:`ARKSpilsSetEpsLin()`.

* :f:func:`FARKSPILSSETJAC()` interfaces to :c:func:`ARKSpilsSetJacTimes()`.

* :f:func:`FARKSPILSSETPREC()` interfaces to :c:func:`ARKSpilsSetPreconditioner()`.



Interface to the mass matrix linear solver interfaces
-------------------------------------------------------

* :f:func:`FARKDLSMASSINIT()` interfaces to :c:func:`ARKDlsSetMassLinearSolver()`.

* :f:func:`FARKDENSESETMASS()` interfaces to :c:func:`ARKDlsSetMassFn()`. 

* :f:func:`FARKBANDSETMASS()` interfaces to :c:func:`ARKDlsSetMassFn()`. 

* :f:func:`FARKSPARSESETMASS()` interfaces to :c:func:`ARKDlsSetMassFn()`. 

* :f:func:`FARKSPILSMASSINIT()` interfaces to :c:func:`ARKSpilsSetMassLinearSolver()`.

* :f:func:`FARKSPILSSETMASSEPSLIN()` interfaces to :c:func:`ARKSpilsSetMassEpsLin()`.

* :f:func:`FARKSPILSSETMASS()` interfaces to :c:func:`ARKSpilsSetMassTimes()`. 

* :f:func:`FARKSPILSSETMASSPREC()` interfaces to
  :c:func:`ARKSpilsSetMassPreconditioner()`. 



.. _FInterface.UserSupplied:


User-supplied routines
---------------------------------------

As with the native C interface, the FARKode solver interface requires
user-supplied functions to specify the ODE problem to be solved.  In
contrast to the case of direct use of ARKode, and of most Fortran ODE
solvers, the names of all user-supplied routines here are fixed, in
order to maximize portability for the resulting mixed-language program. 
As a result, whether using a purely implicit, purely explicit, or
mixed implicit-explicit solver, routines for both :math:`f_E(t,y)` and
:math:`f_I(t,y)` must be provided by the user (though either of which
may do nothing): 

.. cssclass:: table-bordered

+---------------------------+-----------------------------------+
| FARKODE routine           | ARKode interface                  |
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

+--------------------------+--------------------------------------+----------------------------------+
| FARKODE routine          | ARKode interface                     | FARKODE "activation" routine     |
| (FORTRAN, user-supplied) | function type                        |                                  |
+==========================+======================================+==================================+
| :f:func:`FARKDJAC()`     | :c:func:`ARKDlsJacFn()`              | :f:func:`FARKDENSESETJAC()`      |
+--------------------------+--------------------------------------+----------------------------------+
| :f:func:`FARKBJAC()`     | :c:func:`ARKDlsJacFn()`              | :f:func:`FARKBANDSETJAC()`       |
+--------------------------+--------------------------------------+----------------------------------+
| :f:func:`FARKSPJAC()`    | :c:func:`ARKDlsJacFn()`              | :f:func:`FARKSPARSESETJAC()`     |
+--------------------------+--------------------------------------+----------------------------------+
| :f:func:`FARKDMASS()`    | :c:func:`ARKDlsMassFn()`             | :f:func:`FARKDENSESETMASS()`     |
+--------------------------+--------------------------------------+----------------------------------+
| :f:func:`FARKBMASS()`    | :c:func:`ARKDlsMassFn()`             | :f:func:`FARKBANDSETMASS()`      |
+--------------------------+--------------------------------------+----------------------------------+
| :f:func:`FARKSPMASS()`   | :c:func:`ARKDlsMassFn()`             | :f:func:`FARKSPARSESETMASS()`    |
+--------------------------+--------------------------------------+----------------------------------+
| :f:func:`FARKPSET()`     | :c:func:`ARKSpilsPrecSetupFn()`      | :f:func:`FARKSPILSSETPREC()`     |
+--------------------------+--------------------------------------+----------------------------------+
| :f:func:`FARKPSOL()`     | :c:func:`ARKSpilsPrecSolveFn()`      | :f:func:`FARKSPILSSETPREC()`     |
+--------------------------+--------------------------------------+----------------------------------+
| :f:func:`FARKJTSETUP()`  | :c:func:`ARKSpilsJacTimesSetupFn()`  | :f:func:`FARKSPILSSETJAC()`      |
+--------------------------+--------------------------------------+----------------------------------+
| :f:func:`FARKJTIMES()`   | :c:func:`ARKSpilsJacTimesVecFn()`    | :f:func:`FARKSPILSSETJAC()`      |
+--------------------------+--------------------------------------+----------------------------------+
| :f:func:`FARKMASSPSET()` | :c:func:`ARKSpilsMassPrecSetupFn()`  | :f:func:`FARKSPILSSETMASSPREC()` |
+--------------------------+--------------------------------------+----------------------------------+
| :f:func:`FARKMASSPSOL()` | :c:func:`ARKSpilsMassPrecSolveFn()`  | :f:func:`FARKSPILSSETMASSPREC()` |
+--------------------------+--------------------------------------+----------------------------------+
| :f:func:`FARKMTSETUP()`  | :c:func:`ARKSpilsMassTimesSetupFn()` | :f:func:`FARKSPILSSETMASS()`     |
+--------------------------+--------------------------------------+----------------------------------+
| :f:func:`FARKMTIMES()`   | :c:func:`ARKSpilsMassTimesVecFn()`   | :f:func:`FARKSPILSSETMASS()`     |
+--------------------------+--------------------------------------+----------------------------------+
| :f:func:`FARKEWT()`      | :c:func:`ARKEwtFn()`                 | :f:func:`FARKEWTSET()`           |
+--------------------------+--------------------------------------+----------------------------------+
| :f:func:`FARKADAPT()`    | :c:func:`ARKAdaptFn()`               | :f:func:`FARKADAPTSET()`         |
+--------------------------+--------------------------------------+----------------------------------+
| :f:func:`FARKEXPSTAB()`  | :c:func:`ARKExpStabFn()`             | :f:func:`FARKEXPSTABSET()`       |
+--------------------------+--------------------------------------+----------------------------------+
