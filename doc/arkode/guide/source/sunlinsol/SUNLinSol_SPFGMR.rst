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


.. _SUNLinSol_SPFGMR:

The SUNLinSol_SPFGMR Module
======================================

The SPFGMR (Scaled, Preconditioned, Flexible, Generalized Minimum
Residual [S1993]_) implementation of the ``SUNLinearSolver`` module
provided with SUNDIALS, SUNLinSol_SPFGMR, is an iterative linear
solver that is designed to be compatible with any ``N_Vector``
implementation (serial, threaded, parallel, and user-supplied) that
supports a minimal subset of operations (:c:func:`N_VClone()`,
:c:func:`N_VDotProd()`, :c:func:`N_VScale()`,
:c:func:`N_VLinearSum()`, :c:func:`N_VProd()`, :c:func:`N_VConst()`,
:c:func:`N_VDiv()`, and :c:func:`N_VDestroy()`).  Unlike the other
Krylov iterative linear solvers supplied with SUNDIALS, FGMRES is
specifically designed to work with a changing preconditioner
(e.g. from an iterative method).


.. _SUNLinSol_SPFGMR.Usage:

SUNLinSol_SPFGMR Usage
-------------------------

The header file to be included when using this module
is ``sunlinsol/sunlinsol_spfgmr.h``.  The SUNLinSol_SPFGMR module is
accessible from all SUNDIALS solvers *without*
linking to the ``libsundials_sunlinsolspfgmr`` module library.


The module SUNLinSol_SPFGMR provides the following
user-callable routines:


.. c:function:: SUNLinearSolver SUNLinSol_SPFGMR(N_Vector y, int pretype, int maxl)

   This constructor function creates and allocates memory for a SPFGMR
   ``SUNLinearSolver``.  Its arguments are an ``N_Vector``, a flag
   indicating to use preconditioning, and the number of Krylov basis
   vectors to use.

   This routine will perform consistency checks to ensure that it is
   called with a consistent ``N_Vector`` implementation (i.e. that it
   supplies the requisite vector operations).  If ``y`` is
   incompatible, then this routine will return ``NULL``.

   A ``maxl`` argument that is :math:`\le0` will result in the default
   value (5).

   Since the FGMRES algorithm is designed to only support right
   preconditioning, then any of the ``pretype``
   inputs ``PREC_LEFT`` (1), ``PREC_RIGHT`` (2), or ``PREC_BOTH``
   (3) will result in use of ``PREC_RIGHT``;  any other integer input
   will result in the default (no preconditioning).
   We note that some SUNDIALS solvers are designed
   to only work with left preconditioning (IDA and IDAS). While it is
   possible to use a right-preconditioned SUNLinSol_SPFGMR object for
   these packages, this use mode is not supported and may result in
   inferior performance.

.. c:function:: int SUNLinSol_SPFGMRSetPrecType(SUNLinearSolver S, int pretype)

   This function updates the flag indicating use of preconditioning.
   Since the FGMRES algorithm is designed to only support right
   preconditioning, then any of the ``pretype``
   inputs ``PREC_LEFT`` (1), ``PREC_RIGHT`` (2), or ``PREC_BOTH``
   (3) will result in use of ``PREC_RIGHT``;  any other integer input
   will result in the default (no preconditioning).

   This routine will return with one of the error codes
   ``SUNLS_MEM_NULL`` (``S`` is ``NULL``) or ``SUNLS_SUCCESS``.


.. c:function:: int SUNLinSol_SPFGMRSetGSType(SUNLinearSolver S, int gstype)

   This function sets the type of Gram-Schmidt orthogonalization to
   use.  Supported values are ``MODIFIED_GS`` (1) and
   ``CLASSICAL_GS`` (2).  Any other integer input will result in a
   failure, returning error code ``SUNLS_ILL_INPUT``.

   This routine will return with one of the error codes
   ``SUNLS_ILL_INPUT`` (illegal ``gstype``), ``SUNLS_MEM_NULL``
   (``S`` is ``NULL``), or ``SUNLS_SUCCESS``.


.. c:function:: int SUNLinSol_SPFGMRSetMaxRestarts(SUNLinearSolver S, int maxrs)

   This function sets the number of FGMRES restarts to
   allow.  A negative input will result in the default of 0.

   This routine will return with one of the error codes
   ``SUNLS_MEM_NULL`` (``S`` is ``NULL``) or ``SUNLS_SUCCESS``.


.. c:function:: int SUNLinSolSetInfoFile_SPFGMR(SUNLinearSolver LS, FILE* info_file)

   The function :c:func:`SUNLinSolSetInfoFile_SPFGMR()` sets the
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


.. c:function:: int SUNLinSolSetPrintLevel_SPFGMR(SUNLinearSolver LS, int print_level)

   The function :c:func:`SUNLinSolSetPrintLevel_SPFGMR()` specifies the
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

.. c:function:: SUNLinearSolver SUNSPFGMR(N_Vector y, int pretype, int maxl)

   Wrapper function for :c:func:`SUNLinSol_SPFGMR()`

.. c:function:: int SUNSPFGMRSetPrecType(SUNLinearSolver S, int pretype)

   Wrapper function for :c:func:`SUNLinSol_SPFGMRSetPrecType()`

.. c:function:: int SUNSPFGMRSetGSType(SUNLinearSolver S, int gstype)

   Wrapper function for :c:func:`SUNLinSol_SPFGMRSetGSType()`

.. c:function:: int SUNSPFGMRSetMaxRestarts(SUNLinearSolver S, int maxrs)

   Wrapper function for :c:func:`SUNLinSol_SPFGMRSetMaxRestarts()`


For solvers that include a Fortran interface module, the
SUNLinSol_SPFGMR module also includes the Fortran-callable
function :f:func:`FSUNSPFGMRInit()` to initialize this
SUNLinSol_SPFGMR module for a given SUNDIALS solver.

.. f:subroutine:: FSUNSPFGMRInit(CODE, PRETYPE, MAXL, IER)

   Initializes a SPFGMR ``SUNLinearSolver`` structure for
   use in a SUNDIALS package.

   This routine must be called *after* the ``N_Vector`` object has
   been initialized.

   **Arguments:**
      * *CODE* (``int``, input) -- flag denoting the SUNDIALS solver
        this matrix will be used for: CVODE=1, IDA=2, KINSOL=3, ARKode=4.
      * *PRETYPE* (``int``, input) -- flag denoting whether to use
        preconditioning: no=0, yes=1.
      * *MAXL* (``int``, input) -- number of FGMRES basis vectors to use.
      * *IER* (``int``, output) -- return flag (0 success, -1 for failure).

Additionally, when using ARKode with a non-identity mass matrix, the
Fortran-callable function :f:func:`FSUNMassSPFGMRInit()` initializes
this SUNLinSol_SPFGMR module for solving mass matrix linear systems.

.. f:subroutine:: FSUNMassSPFGMRInit(PRETYPE, MAXL, IER)

   Initializes a SPFGMR ``SUNLinearSolver`` structure for use in
   solving mass matrix systems in ARKode.

   This routine must be called *after* the ``N_Vector`` object has
   been initialized.

   **Arguments:**
      * *PRETYPE* (``int``, input) -- flag denoting whether to use
        preconditioning: no=0, yes=1.
      * *MAXL* (``int``, input) -- number of FGMRES basis vectors to use.
      * *IER* (``int``, output) -- return flag (0 success, -1 for failure).


The :c:func:`SUNLinSol_SPFGMRSetGSType()`, :c:func:`SUNLinSol_SPFGMRSetPrecType()`
and :c:func:`SUNLinSol_SPFGMRSetMaxRestarts()` routines also support Fortran
interfaces for the system and mass matrix solvers:

.. f:subroutine:: FSUNSPFGMRSetGSType(CODE, GSTYPE, IER)

   Fortran interface to :c:func:`SUNLinSol_SPFGMRSetGSType()` for system
   linear solvers.

   This routine must be called *after* :f:func:`FSUNSPFGMRInit()` has
   been called.

   **Arguments:** all should have type ``int``, and have meanings
   identical to those listed above.

.. f:subroutine:: FSUNMassSPFGMRSetGSType(GSTYPE, IER)

   Fortran interface to :c:func:`SUNLinSol_SPFGMRSetGSType()` for mass matrix
   linear solvers in ARKode.

   This routine must be called *after* :f:func:`FSUNMassSPFGMRInit()` has
   been called.

   **Arguments:** all should have type ``int``, and have meanings
   identical to those listed above.

.. f:subroutine:: FSUNSPFGMRSetPrecType(CODE, PRETYPE, IER)

   Fortran interface to :c:func:`SUNLinSol_SPFGMRSetPrecType()` for system
   linear solvers.

   This routine must be called *after* :f:func:`FSUNSPFGMRInit()` has
   been called.

   **Arguments:** all should have type ``int``, and have meanings
   identical to those listed above.

.. f:subroutine:: FSUNMassSPFGMRSetPrecType(PRETYPE, IER)

   Fortran interface to :c:func:`SUNLinSol_SPFGMRSetPrecType()` for mass matrix
   linear solvers in ARKode.

   This routine must be called *after* :f:func:`FSUNMassSPFGMRInit()` has
   been called.

   **Arguments:** all should have type ``int``, and have meanings
   identical to those listed above.


.. f:subroutine:: FSUNSPFGMRSetMaxRS(CODE, MAXRS, IER)

   Fortran interface to :c:func:`SUNLinSol_SPFGMRSetMaxRS()` for system
   linear solvers.

   This routine must be called *after* :f:func:`FSUNSPFGMRInit()` has
   been called.

   **Arguments:** all should have type ``int``, and have meanings
   identical to those listed above.

.. f:subroutine:: FSUNMassSPFGMRSetMaxRS(MAXRS, IER)

   Fortran interface to :c:func:`SUNLinSol_SPFGMRSetMaxRS()` for mass matrix
   linear solvers in ARKode.

   This routine must be called *after* :f:func:`FSUNMassSPFGMRInit()` has
   been called.

   **Arguments:** all should have type ``int``, and have meanings
   identical to those listed above.





.. _SUNLinSol_SPFGMR.Description:

SUNLinSol_SPFGMR Description
---------------------------------


The SUNLinSol_SPFGMR module defines the *content* field of a
``SUNLinearSolver`` to be the following structure:

.. code-block:: c

   struct _SUNLinearSolverContent_SPFGMR {
     int maxl;
     int pretype;
     int gstype;
     int max_restarts;
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
     N_Vector *V;
     N_Vector *Z;
     realtype **Hes;
     realtype *givens;
     N_Vector xcor;
     realtype *yg;
     N_Vector vtemp;
     int      print_level;
     FILE*    info_file;
   };

These entries of the *content* field contain the following
information:

* ``maxl`` - number of FGMRES basis vectors to use (default is 5),

* ``pretype`` - flag for use of preconditioning (default is none),

* ``gstype`` - flag for type of Gram-Schmidt orthogonalization
  (default is modified Gram-Schmidt),

* ``max_restarts`` - number of FGMRES restarts to allow
  (default is 0),

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

* ``V`` - the array of Krylov basis vectors
  :math:`v_1, \ldots, v_{\text{maxl}+1}`, stored in
  ``V[0], ..., V[maxl]``. Each :math:`v_i` is a vector of type ``N_Vector``,

* ``Z`` - the array of preconditioned Krylov basis vectors
  :math:`z_1, \ldots, z_{\text{maxl}+1}`, stored in
  ``Z[0], ..., Z[maxl]``. Each :math:`z_i` is a vector of type ``N_Vector``,

* ``Hes`` - the :math:`(\text{maxl}+1)\times\text{maxl}`
  Hessenberg matrix. It is stored row-wise so that the (i,j)th
  element is given by ``Hes[i][j]``,

* ``givens`` - a length :math:`2\,\text{maxl}` array which represents
  the Givens rotation matrices that arise in the FGMRES
  algorithm. These matrices are :math:`F_0, F_1, \ldots, F_j`, where

  .. math::

     F_i = \begin{bmatrix}
        1 &        &   &     &      &   &        &   \\
          & \ddots &   &     &      &   &        &   \\
          &        & 1 &     &      &   &        &   \\
          &        &   & c_i & -s_i &   &        &   \\
          &        &   & s_i &  c_i &   &        &   \\
          &        &   &     &      & 1 &        &   \\
          &        &   &     &      &   & \ddots &   \\
          &        &   &     &      &   &        & 1\end{bmatrix},

  are represented in the ``givens`` vector as
  ``givens[0]`` :math:`= c_0`,
  ``givens[1]`` :math:`= s_0`,
  ``givens[2]`` :math:`= c_1`,
  ``givens[3]`` :math:`= s_1`, :math:`\ldots`,
  ``givens[2j]`` :math:`= c_j`,
  ``givens[2j+1]`` :math:`= s_j`,

* ``xcor`` - a vector which holds the scaled, preconditioned
  correction to the initial guess,

* ``yg`` - a length :math:`(\text{maxl}+1)` array of ``realtype``
  values used to hold "short" vectors (e.g. :math:`y` and :math:`g`),

* ``vtemp`` - temporary vector storage.

* ``print_level`` - controls the amount of information to be printed to the info file

* ``info_file``   - the file where all informative (non-error) messages will be directed


This solver is constructed to perform the following operations:

* During construction, the ``xcor`` and ``vtemp`` arrays are cloned
  from a template ``N_Vector`` that is input, and default solver
  parameters are set.

* User-facing "set" routines may be called to modify default
  solver parameters.

* Additional "set" routines are called by the SUNDIALS solver
  that interfaces with SUNLinSol_SPFGMR to supply the
  ``ATimes``, ``PSetup``, and ``Psolve`` function pointers and
  ``s1`` and ``s2`` scaling vectors.

* In the "initialize" call, the remaining solver data is
  allocated (``V``, ``Hes``, ``givens``, and ``yg`` )

* In the "setup" call, any non-``NULL`` ``PSetup`` function is called.
  Typically, this is provided by the SUNDIALS solver itself, that
  translates between the generic ``PSetup`` function and the
  solver-specific routine (solver-supplied or user-supplied).

* In the "solve" call, the FGMRES iteration is performed.  This
  will include scaling, preconditioning, and restarts if those options
  have been supplied.

The SUNLinSol_SPFGMR module defines implementations of all
"iterative" linear solver operations listed in the section
:ref:`SUNLinSol.API`:

* ``SUNLinSolGetType_SPFGMR``

* ``SUNLinSolInitialize_SPFGMR``

* ``SUNLinSolSetATimes_SPFGMR``

* ``SUNLinSolSetPreconditioner_SPFGMR``

* ``SUNLinSolSetScalingVectors_SPFGMR``

* ``SUNLinSolSetup_SPFGMR``

* ``SUNLinSolSolve_SPFGMR``

* ``SUNLinSolNumIters_SPFGMR``

* ``SUNLinSolResNorm_SPFGMR``

* ``SUNLinSolResid_SPFGMR``

* ``SUNLinSolLastFlag_SPFGMR``

* ``SUNLinSolSpace_SPFGMR``

* ``SUNLinSolFree_SPFGMR``
