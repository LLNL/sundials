..
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

.. _SUNonlinSol_PetscSNES:

================================================
The SUNNonlinearSolver_PetscSNES implementation
================================================

This section describes the SUNNonlinSol interface to the PETSc SNES nonlinear
solver(s).To enable the SUNonlinSol_PetscSNES module, SUNDIALS must be
configured to use PETSc. Instructions on how to do thus are given in Chapter
:ref:`Installation.CMake.ExternalLibraries.PETSc`. To access the
SUNNonlinSol_PetscSNES module, include the header file
``sunnonlinsol/sunnonlinsol_petscsnes.h``. The library to link to is
``libsundials_sunnonlinsolpetsc.lib`` where ``.lib`` is typically ``.so`` for
shared libaries and ``.a`` for static libraries. Users of the
``SUNNonlinearSolver_PetscSNES`` should also see the section
:ref:`NVectors.NVPETSc` which discusses the NVECTOR interface to the PETSc Vec
API.

.. _SUNNonlinSolPetscSNES.Description:

SUNNonlinearSolver_PetscSNES description
----------------------------------------

The ``SUNNonlinearSolver_PetscSNES`` implementation allows users to utilize a
PETSc SNES nonlinear solver to solve the nonlinear systems that arise in the
SUNDIALS integrators. Since SNES uses the KSP linear solver interface underneath
it, the ``SUNNonlinearSolver_PetscSNES`` implementation does not interface with
SUNDIALS linear solvers. Instead, users should set nonlinear solver options,
linear solver options, and preconditioner options through the PETSc SNES, KSP,
and PC APIs.

*Important usage notes for the ``SUNNonlinearSolver_PetscSNES`` implementation
are provided below:*

* The ``SUNNonlinearSolver_PetscSNES`` implementation handles calling
  ``SNESSetFunction`` at construction. The actual residual function :math:`F(y)`
  is set by the SUNDIALS integrator when the ``SUNNonlinearSolver_PetscSNES``
  object is attached to it. Therefore, a user should not call ``SNESSetFunction``
  on a ``SNES`` object that is being used with ``SUNNonlinearSolver_PetscSNES``.
  For these reasons, it is recommended, although not always necessary, that the
  user calls ``SUNNonlinSol_PetscSNES`` with the new ``SNES`` object immediately
  after calling ``SNESCreate``.

* The number of nonlinear iterations is tracked by SUNDIALS separately from the
  count kept by SNES. As such, the function ``SUNNonlinSolGetNumIters`` reports
  the cumulative number of iterations across the lifetime of the
  ``SUNNonlinearSolver`` object.

* Some "converged" and "diverged" convergence reasons returned by SNES are
  treated as recoverable convergence failures by SUNDIALS. Therefore, the count of
  convergence failures returned by ``SUNNonlinSolGetNumConvFails`` will reflect
  the number of recoverable convergence failures as determined by SUNDIALS, and
  may differ from the count returned by ``SNESGetNonlinearStepFailures``.

* The ``SUNNonlinearSolver_PetscSNES`` module is not currently compatible with
  the CVODES or IDAS staggered or simultaneous sensitivity strategies.

.. _SUNNonlinSolPetscSNES.functions:

SUNNonlinearSolver_PetscSNES functions
--------------------------------------

The ``SUNNonlinearSolver_PetscSNES`` module provides the following constructor
for creating a ``SUNNonlinearSolver`` object.

.. c:function:: SUNNonlinearSolver SUNNonlinSol_PetscSNES(N_Vector y, SNES snes)

  The function ``SUNNonlinSol_PetscSNES`` creates a ``SUNNonlinearSolver``
  object that wraps a PETSc ``SNES`` object for use with SUNDIALS. This
  will call ``SNESSetFunction`` on the provided ``SNES`` object.

  **Arguments:**
    * *snes* -- a PETSc ``SNES`` object
    * *y* -- a ``N_Vector`` object of type NVECTOR_PETSC that is used as a template
      for the residual vector

  **Return value:** a SUNNonlinSol object if the constructor exits
   successfully, otherwise it will be ``NULL``.

  *This function calls ``SNESSetFunction`` and will overwrite whatever
  function was previously set. Users should not call ``SNESSetFunction``
  on the ``SNES`` object provided to the constructor.*

The ``SUNNonlinSol_PetscSNES`` module implements all of the functions defined in
sections :ref:`SUNNonlinSol.CoreFn` through :ref:`SUNNonlinSol.GetFn` except for
``SUNNonlinSolSetup``, ``SUNNonlinSolSetLSetupFn``, ``SUNNonlinSolSetLSolveFn``,
``SUNNonlinSolSetConvTestFn``, and ``SUNNonlinSolSetMaxIters``.

The ``SUNNonlinSol_PetscSNES`` functions have the same names as those defined by
the generic ``SUNNonlinearSolver`` API with ``_PetscSNES`` appended to the
function name. Unless using the ``SUNNonlinSol_PetscSNES`` module as a
standalone nonlinear solver the generic functions defined in sections
:ref:`SUNNonlinSol.CoreFn` through :ref:`SUNNonlinSol.GetFn` should be called in
favor of the ``SUNNonlinSol_PetscSNES`` specific implementations.

The ``SUNNonlinSol_PetscSNES`` module also defines the following additional
user-callable functions.

.. c:function:: int SUNNonlinSolGetSNES_PetscSNES(SUNNonlinearSolver NLS, SNES* snes)

  The function ``SUNNonlinSolGetSNES_PetscSNES`` gets the ``SNES`` object that
  was wrapped.

  **Arguments:**
    * *NLS* -- a ``SUNNonlinearSolver`` object
    * *snes* -- a pointer to a PETSc ``SNES`` object that will be set upon return

  **Return value:** The return value (of type ``int``) should be zero
  for a successful call, and a negative value for a failure.

.. c:function:: int SUNNonlinSolGetPetscError_PetscSNES(SUNNonlinearSolver NLS, PestcErrorCode* error)

  The function ``SUNNonlinSolGetPetscError_PetscSNES`` gets the last error code
  returned by the last internal call to a PETSc API function.

  **Arguments:**
    * *NLS* -- a ``SUNNonlinearSolver`` object
    * *error* -- a pointer to a PETSc error integer that will be set upon return

  **Return value:** The return value (of type ``int``) should be zero
  for a successful call, and a negative value for a failure.

.. c:function:: int SUNNonlinSolGetSysFn_PetscSNES(SUNNonlinearSolver NLS, SUNNonlinSolSysFn* SysFn)

  The function ``SUNNonlinSolGetSysFn_PetscSNES`` returns the residual
  function that defines the nonlinear system.

  **Arguments:**
    * *NLS* -- a ``SUNNonlinearSolver`` object
    * *SysFn* -- the function defining the nonlinear system

  **Return value:** The return value (of type ``int``) should be zero
  for a successful call, and a negative value for a failure.

.. _SUNNonlinSolPetscSNES.Content:

SUNNonlinearSolver_PetscSNES content
------------------------------------

The *content* field of the SUNNonlinSol_Newton module is the following
structure.

.. code-block:: c

  struct _SUNNonlinearSolverContent_PetscSNES {
    int sysfn_last_err;
    PetscErrorCode petsc_last_err;
    long int nconvfails;
    long int nni;
    void *imem;
    SNES snes;
    Vec r;
    N_Vector y, f;
    SUNNonlinSolSysFn Sys;
  };

These entries of the *content* field contain the following information:

* ``sysfn_last_err``  -- last error returned by the system defining function,
* ``petsc_last_err``  -- last error returned by PETSc
* ``nconvfails``      -- number of nonlinear converge failures (recoverable or not),
* ``nni``             -- number of nonlinear iterations,
* ``imem``            -- SUNDIALS integrator memory,
* ``snes``            -- PETSc ``SNES`` object,
* ``r``               -- the nonlinear residual,
* ``y``               -- wrapper for PETSc vectors used in the system function,
* ``f``               -- wrapper for PETSc vectors used in the system function,
* ``Sys``             -- nonlinear system definining function.
