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


.. _SUNLinSol_SPBCGS:

The SUNLinSol_SPBCGS Module
======================================

The SPBCGS (Scaled, Preconditioned, Bi-Conjugate Gradient,
Stabilized [V1992]_) implementation of the ``SUNLinearSolver`` module
provided with SUNDIALS, SUNLinSol_SPBCGS, is an iterative linear
solver that is designed to be compatible with any ``N_Vector``
implementation (serial, threaded, parallel, and user-supplied) that
supports a minimal subset of operations (:c:func:`N_VClone()`,
:c:func:`N_VDotProd()`, :c:func:`N_VScale()`,
:c:func:`N_VLinearSum()`, :c:func:`N_VProd()`, :c:func:`N_VDiv()`, and
:c:func:`N_VDestroy()`).  Unlike the SPGMR and SPFGMR algorithms,
SPBCGS requires a fixed amount of memory that does not increase with
the number of allowed iterations.


.. _SUNLinSol_SPBCGS.Usage:

SUNLinSol_SPBCGS Usage
------------------------

The header file to be included when using this module
is ``sunlinsol/sunlinsol_spbcgs.h``.  The SUNLinSol_SPBCGS module
is accessible from all SUNDIALS solvers *without*
linking to the ``libsundials_sunlinsolspbcgs`` module library.


The module SUNLinSol_SPBCGS provides the following
user-callable routines:


.. c:function:: SUNLinearSolver SUNLinSol_SPBCGS(N_Vector y, int pretype, int maxl)

   This constructor function creates and allocates memory for a SPBCGS
   ``SUNLinearSolver``.  Its arguments are an ``N_Vector``, the desired
   type of preconditioning, and the number of linear iterations to allow.

   This routine will perform consistency checks to ensure that it is
   called with a consistent ``N_Vector`` implementation (i.e. that it
   supplies the requisite vector operations).  If ``y`` is
   incompatible, then this routine will return ``NULL``.

   A ``maxl`` argument that is :math:`\le0` will result in the default
   value (5).

   Allowable inputs for ``pretype`` are ``PREC_NONE`` (0),
   ``PREC_LEFT`` (1), ``PREC_RIGHT`` (2) and ``PREC_BOTH`` (3);
   any other integer input will result in the default (no
   preconditioning).  We note that some SUNDIALS solvers are designed
   to only work with left preconditioning (IDA and IDAS) and others
   with only right preconditioning (KINSOL). While it is possible to
   configure a SUNLinSol_SPBCGS object to use any of the
   preconditioning options with these solvers, this use mode is not
   supported and may result in inferior performance.


.. c:function:: int SUNLinSol_SPBCGSSetPrecType(SUNLinearSolver S, int pretype)

   This function updates the type of preconditioning to use.  Supported
   values are ``PREC_NONE`` (0), ``PREC_LEFT`` (1),
   ``PREC_RIGHT`` (2), and ``PREC_BOTH`` (3).

   This routine will return with one of the error codes
   ``SUNLS_ILL_INPUT`` (illegal ``pretype``), ``SUNLS_MEM_NULL``
   (``S`` is ``NULL``), or ``SUNLS_SUCCESS``.


.. c:function:: int SUNLinSol_SPBCGSSetMaxl(SUNLinearSolver S, int maxl)

   This function updates the number of linear solver iterations to allow.

   A ``maxl`` argument that is :math:`\le0` will result in the default
   value (5).

   This routine will return with one of the error codes
   ``SUNLS_MEM_NULL`` (``S`` is ``NULL``) or ``SUNLS_SUCCESS``.


.. c:function:: int SUNLinSolSetInfoFile_SPBCGS(SUNLinearSolver LS, FILE* info_file)

   The function :c:func:`SUNLinSolSetInfoFile_SPBCGS()` sets the
   output file where all informative (non-error) messages should be directed.

   **Arguments:**
      * *LS* -- a SUNLinSol object
      * *info_file* -- pointer to output file (``stdout`` by default);
         a ``NULL`` input will disable output

   **Return value:**
      * *SUNLS_SUCCESS* if successful
      * *SUNLS_MEM_NULL* if the SUNLinearSolver memory was ``NULL``
      * *SUNLS_ILL_INPUT* if SUNDIALS was not built with monitoring enabled

   **Notes:**
   This function is intended for users that wish to monitor the linear
   solver progress. By default, the file pointer is set to ``stdout``.

   **SUNDIALS must be built with the CMake option
   ``SUNDIALS_BUILD_WITH_MONITORING``, to utilize this function.**
   See section :ref:`Installation.CMake.Options` for more information.


.. c:function:: int SUNLinSolSetPrintLevel_SPBCGS(SUNLinearSolver LS, int print_level)

   The function :c:func:`SUNLinSolSetPrintLevel_SPBCGS()` specifies the
   level of verbosity of the output.

   **Arguments:**
      * *LS* -- a SUNLinSol object
      * *print_level* -- flag indicating level of verbosity;
        must be one of:

         * 0, no information is printed (default)
         * 1, for each linear iteration the residual norm is printed

   **Return value:**
      * *SUNLS_SUCCESS* if successful
      * *SUNLS_MEM_NULL* if the SUNLinearSolver memory was ``NULL``
      * *SUNLS_ILL_INPUT* if SUNDIALS was not built with monitoring enabled, or
        if the print level value was invalid

   **Notes:**
   This function is intended for users that wish to monitor the linear
   solver progress. By default, the print level is 0.

   **SUNDIALS must be built with the CMake option
   ``SUNDIALS_BUILD_WITH_MONITORING``, to utilize this function.**
   See section :ref:`Installation.CMake.Options` for more information.


For backwards compatibility, we also provide the wrapper functions,
each with identical input and output arguments to the routines that
they wrap:

.. c:function:: SUNLinearSolver SUNSPBCGS(N_Vector y, int pretype, int maxl)

   Wrapper function for :c:func:`SUNLinSol_SPBCGS()`

.. c:function:: int SUNSPBCGSSetPrecType(SUNLinearSolver S, int pretype)

   Wrapper function for :c:func:`SUNLinSol_SPBCGSSetPrecType()`

.. c:function:: int SUNSPBCGSSetMaxl(SUNLinearSolver S, int maxl)

   Wrapper function for :c:func:`SUNLinSol_SPBCGSSetMaxl()`




For solvers that include a Fortran interface module, the
SUNLinSol_SPBCGS module also includes the Fortran-callable
function :f:func:`FSUNSPBCGSInit()` to initialize this
SUNLinSol_SPBCGS module for a given SUNDIALS solver.

.. f:subroutine:: FSUNSPBCGSInit(CODE, PRETYPE, MAXL, IER)

   Initializes a SPBCGS ``SUNLinearSolver`` structure for
   use in a SUNDIALS package.

   This routine must be called *after* the ``N_Vector`` object has
   been initialized.

   **Arguments:**
      * *CODE* (``int``, input) -- flag denoting the SUNDIALS solver
        this matrix will be used for: CVODE=1, IDA=2, KINSOL=3, ARKode=4.
      * *PRETYPE* (``int``, input) -- flag denoting type of
        preconditioning to use: none=0, left=1, right=2, both=3.
      * *MAXL* (``int``, input) -- number of SPBCGS iterations to allow.
      * *IER* (``int``, output) -- return flag (0 success, -1 for failure).

Additionally, when using ARKode with a non-identity mass matrix, the
Fortran-callable function :f:func:`FSUNMassSPBCGSInit()` initializes
this SUNLinSol_SPBCGS module for solving mass matrix linear systems.

.. f:subroutine:: FSUNMassSPBCGSInit(PRETYPE, MAXL, IER)

   Initializes a SPBCGS ``SUNLinearSolver`` structure for use in
   solving mass matrix systems in ARKode.

   This routine must be called *after* the ``N_Vector`` object has
   been initialized.

   **Arguments:**
      * *PRETYPE* (``int``, input) -- flag denoting type of
        preconditioning to use: none=0, left=1, right=2, both=3.
      * *MAXL* (``int``, input) -- number of SPBCGS iterations to allow.
      * *IER* (``int``, output) -- return flag (0 success, -1 for failure).

The :c:func:`SUNLinSol_SPBCGSSetPrecType()` and :c:func:`SUNLinSol_SPBCGSSetMaxl()`
routines also support Fortran interfaces for the system and mass
matrix solvers:

.. f:subroutine:: FSUNSPBCGSSetPrecType(CODE, PRETYPE, IER)

   Fortran interface to :c:func:`SUNLinSol_SPBCGSSetPrecType()` for system
   linear solvers.

   This routine must be called *after* :f:func:`FSUNSPBCGSInit()` has
   been called.

   **Arguments:** all should have type ``int``, and have meanings
   identical to those listed above.


.. f:subroutine:: FSUNMassSPBCGSSetPrecType(PRETYPE, IER)

   Fortran interface to :c:func:`SUNLinSol_SPBCGSSetPrecType()` for mass matrix
   linear solvers in ARKode.

   This routine must be called *after* :f:func:`FSUNMassSPBCGSInit()` has
   been called.

   **Arguments:** all should have type ``int``, and have meanings
   identical to those listed above.


.. f:subroutine:: FSUNSPBCGSSetMaxl(CODE, MAXL, IER)

   Fortran interface to :c:func:`SUNLinSol_SPBCGSSetMaxl()` for system
   linear solvers.

   This routine must be called *after* :f:func:`FSUNSPBCGSInit()` has
   been called.

   **Arguments:** all should have type ``int``, and have meanings
   identical to those listed above.


.. f:subroutine:: FSUNMassSPBCGSSetMaxl(MAXL, IER)

   Fortran interface to :c:func:`SUNLinSol_SPBCGSSetMaxl()` for mass matrix
   linear solvers in ARKode.

   This routine must be called *after* :f:func:`FSUNMassSPBCGSInit()` has
   been called.

   **Arguments:** all should have type ``int``, and have meanings
   identical to those listed above.




.. _SUNLinSol_SPBCGS.Description:

SUNLinSol_SPBCGS Description
-------------------------------

The SUNLinSol_SPBCGS module defines the *content* field of a
``SUNLinearSolver`` to be the following structure:

.. code-block:: c

   struct _SUNLinearSolverContent_SPBCGS {
     int maxl;
     int pretype;
     int numiters;
     realtype resnorm;
     int last_flag;
     ATimesFn ATimes;
     void* ATData;
     PSetupFn Psetup;
     PSolveFn Psolve;
     void* PData;
     N_Vector s1;
     N_Vector s2;
     N_Vector r;
     N_Vector r_star;
     N_Vector p;
     N_Vector q;
     N_Vector u;
     N_Vector Ap;
     N_Vector vtemp;
     int      print_level;
     FILE*    info_file;
   };

These entries of the *content* field contain the following
information:

* ``maxl`` - number of SPBCGS iterations to allow (default is 5),

* ``pretype`` - flag for type of preconditioning to employ
  (default is none),

* ``numiters`` - number of iterations from the most-recent solve,

* ``resnorm`` - final linear residual norm from the most-recent
  solve,

* ``last_flag`` - last error return flag from an internal
  function,

* ``ATimes`` - function pointer to perform :math:`Av` product,

* ``ATData`` - pointer to structure for ``ATimes``,

* ``Psetup`` - function pointer to preconditioner setup routine,

* ``Psolve`` - function pointer to preconditioner solve routine,

* ``PData`` - pointer to structure for ``Psetup`` and ``Psolve``,

* ``s1, s2`` - vector pointers for supplied scaling matrices
  (default is ``NULL``),

* ``r`` - a ``N_Vector`` which holds the current scaled,
  preconditioned linear system residual,

* ``r_star`` - a ``N_Vector`` which holds the initial scaled,
  preconditioned linear system residual,

* ``p, q, u, Ap, vtemp`` - ``N_Vector`` used for workspace by the
  SPBCGS algorithm.

* ``print_level`` - controls the amount of information to be printed to the info file

* ``info_file``   - the file where all informative (non-error) messages will be directed


This solver is constructed to perform the following operations:

* During construction all ``N_Vector`` solver data is allocated, with
  vectors cloned from a template ``N_Vector`` that is input, and
  default solver parameters are set.

* User-facing "set" routines may be called to modify default
  solver parameters.

* Additional "set" routines are called by the SUNDIALS solver
  that interfaces with SUNLinSol_SPBCGS to supply the ``ATimes``,
  ``PSetup``, and ``Psolve`` function pointers and ``s1`` and ``s2``
  scaling vectors.

* In the "initialize" call, the solver parameters are checked
  for validity.

* In the "setup" call, any non-``NULL`` ``PSetup`` function is
  called.  Typically, this is provided by the SUNDIALS solver itself,
  that translates between the generic ``PSetup`` function and the
  solver-specific routine (solver-supplied or user-supplied).

* In the "solve" call the SPBCGS iteration is performed.  This
  will include scaling and preconditioning if those options have been
  supplied.

The SUNLinSol_SPBCGS module defines implementations of all
"iterative" linear solver operations listed in the section
:ref:`SUNLinSol.API`:

* ``SUNLinSolGetType_SPBCGS``

* ``SUNLinSolInitialize_SPBCGS``

* ``SUNLinSolSetATimes_SPBCGS``

* ``SUNLinSolSetPreconditioner_SPBCGS``

* ``SUNLinSolSetScalingVectors_SPBCGS``

* ``SUNLinSolSetup_SPBCGS``

* ``SUNLinSolSolve_SPBCGS``

* ``SUNLinSolNumIters_SPBCGS``

* ``SUNLinSolResNorm_SPBCGS``

* ``SUNLinSolResid_SPBCGS``

* ``SUNLinSolLastFlag_SPBCGS``

* ``SUNLinSolSpace_SPBCGS``

* ``SUNLinSolFree_SPBCGS``
