..
   Programmer(s): Daniel R. Reynolds @ SMU
   ----------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2002-2019, Lawrence Livermore National Security
   and Southern Methodist University.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   ----------------------------------------------------------------

:tocdepth: 3

.. _Introduction:

Introduction
============

The ARKode infrastructure provides adaptive-step time integration
modules for stiff, nonstiff and mixed stiff/nonstiff systems of
ordinary differential equations (ODEs).  ARKode itself is structured
to support a wide range of one-step (but multi-stage) methods,
allowing for rapid development of parallel implementations of
state-of-the-art time integration methods.  At present, ARKode is
packaged with two time-stepping modules, *ARKStep* and *ERKStep*.


*ARKStep* supports ODE systems posed in split, linearly-implicit form,

.. math::
   M \dot{y} = f^E(t,y) + f^I(t,y),  \qquad y(t_0) = y_0,
   :label: ODE_split_linearly_implicit

where :math:`t` is the independent variable, :math:`y` is the set of
dependent variables (in :math:`\mathbb{R}^N`), :math:`M` is a
user-specified, nonsingular operator from :math:`\mathbb{R}^N` to
:math:`\mathbb{R}^N`, and the right-hand side function is partitioned
into up to two components:

- :math:`f^E(t,y)` contains the "nonstiff" time scale components to be
  integrated explicitly, and
- :math:`f^I(t,y)`  contains the "stiff" time scale components to be
  integrated implicitly.

Either of these operators may be disabled, allowing for fully
explicit, fully implicit, or combination implicit-explicit (ImEx) time
integration.

The algorithms used in ARKStep are adaptive- and fixed-step additive
Runge Kutta methods. Such methods are defined through combining two
complementary Runge-Kutta methods: one explicit (ERK) and the other
diagonally implicit (DIRK).  Through appropriately partitioning the
ODE right-hand side into explicit and implicit components
:eq:`ODE_split_linearly_implicit`, such methods have the potential to
enable accurate and efficient time integration of stiff, nonstiff, and
mixed stiff/nonstiff systems of ordinary differential equations.  A
key feature allowing for high efficiency of these methods is that only
the components in :math:`f^I(t,y)` must be solved implicitly, allowing
for splittings tuned for use with optimal implicit solver algorithms.

This framework allows for significant freedom over the constitutive
methods used for each component, and ARKode is packaged with a wide
array of built-in methods for use.  These built-in Butcher tables
include adaptive explicit methods of orders 2-8, adaptive implicit
methods of orders 2-5, and adaptive ImEx methods of orders 3-5.


*ERKStep* focuses specifically on problems posed in explicit form,

.. math::
   \dot{y} = f(t,y),  \qquad y(t_0) = y_0.
   :label: ODE_explicit

allowing for increased computational efficiency and memory savings.
The algorithms used in ERKStep are adaptive- and fixed-step explicit
Runge Kutta methods.   As with ARKStep, the ERKStep module is packaged
with adaptive explicit methods of orders 2-8.


For problems that include nonzero implicit term :math:`f^I(t,y)`, the
resulting implicit system (assumed nonlinear, unless specified
otherwise) is solved approximately at each integration step, using a
modified Newton method, inexact Newton method, or an
accelerated fixed-point solver.  For the Newton-based methods and the
serial or threaded NVECTOR modules in SUNDIALS, ARKode may use a
variety of linear solvers provided with SUNDIALS, including both
direct (dense, band, or sparse) and preconditioned Krylov iterative
(GMRES [SS1986]_, BiCGStab [V1992]_, TFQMR [F1993]_, FGMRES [S1993]_,
or PCG [HS1952]_) linear solvers.  When used with the MPI-based
parallel, PETSc, *hypre*, CUDA, and Raja NVECTOR modules, or a
user-provided vector data structure, only the Krylov solvers are
available, although a user may supply their own linear solver for any
data structures if desired.  For the serial or threaded vector
structures, we provide a banded preconditioner module called ARKBANDPRE
that may be used with the Krylov solvers, while for the MPI-based
parallel vector structure there is a preconditioner module called
ARKBBDPRE which provides a band-block-diagonal preconditioner.
Additionally, a user may supply more optimal, problem-specific
preconditioner routines.




Changes from previous versions
--------------------------------

Changes in v4.0.0-dev.0
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

An additional NVECTOR implementation, NVECTOR_MANYVECTOR, was
created to support flexible partitioning of solution data among
different processing elements (e.g., CPU + GPU) or for multi-physics
problems that couple distinct MPI-based simulations together (see the
Section :ref:`NVectors.ManyVector` for more details).  This
implementation is accompanied by additions to user documentation and
SUNDIALS examples.

Eleven new optional vector operations have been added to the NVECTOR API to
support the new NVECTOR_MANYVECTOR implementation (see the Chapter :ref:`NVectors`
for more details). Two of the operations, :c:func:`N_VGetCommunicator()` and
:c:func:`N_VGetLength()`, must be implemented by subvectors that are combined to
create an NVECTOR_MANYVECTOR, but are not used outside of this context. The
remaining nine operations are optional local reduction operations intended to
eliminate unnecessary latency when performing vector reduction operations
(norms, etc.) on distributed memory systems. The optional local reduction vector
operations are
:c:func:`N_VDotProdLocal()`,
:c:func:`N_VMaxNormLocal()`,
:c:func:`N_VMinLocal()`,
:c:func:`N_VL1NormLocal()`,
:c:func:`N_VWSqrSumLocal()`,
:c:func:`N_VWSqrSumMaskLocal()`,
:c:func:`N_VInvTestLocal()`,
:c:func:`N_VConstrMaskLocal()`, and
:c:func:`N_VMinQuotientLocal()`.
If an NVECTOR implementation defines any of the local operations as ``NULL``,
then the NVECTOR_MANYVECTOR will call standard NVECTOR operations to complete
the computation.

A new SUNMatrix and SUNLinearSolver implementation was added to facilitate 
the use of the SuperLU_DIST library with SUNDIALS.

A new operation, ``SUNMatMatvecSetup``, was added to the SUNMatrix API. Users
who have implemented custom ``SUNMatrix`` modules will need to at least update
their code to set the corresponding ops structure member, ``matvecsetup``, to ``NULL``.

The generic SUNMatrix API now defines error codes to be returned by SUNMatrix operations.
Operations which return an integer flag indiciating success/failure may return different
values than previously.

The MRIStep module has been updated to support explicit, implicit, or IMEX
methods as the fast integrator using the ARKStep module. As a result some
function signatures have been changed including :c:func:`MRIStepCreate` which
now takes an ARKStep memory structure for the fast integration as an input.

Changes in v3.1.0
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

An additional NVECTOR implementation was added for the
Tpetra vector from the Trilinos library to facilitate interoperability
between SUNDIALS and Trilinos. This implementation is accompanied by
additions to user documentation and SUNDIALS examples.

A bug was fixed where a nonlinear solver object could be freed twice in some use
cases.

The ``EXAMPLES_ENABLE_RAJA`` CMake option has been removed. The option ``EXAMPLES_ENABLE_CUDA``
enables all examples that use CUDA including the RAJA examples with a CUDA back end
(if the RAJA NVECTOR is enabled).

The implementation header file `arkode_impl.h` is no longer installed. This means users
who are directly manipulating the ``ARKodeMem`` structure will need to update their code
to use ARKode's public API.

Python is no longer required to run ``make test`` and ``make test_install``.

Fixed a bug in ``ARKodeButcherTable_Write`` when printing a Butcher table
without an embedding.

Changes in v3.0.2
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Added information on how to contribute to SUNDIALS and a contributing agreement.

Changes in v3.0.1
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A bug in ARKode where single precision builds would fail to compile has been fixed.


Changes in v3.0.0
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The ARKode library has been entirely rewritten to support a modular
approach to one-step methods, which should allow rapid research and
development of novel integration methods without affecting existing
solver functionality.  To support this, the existing ARK-based methods
have been encapsulated inside the new ``ARKStep`` time-stepping
module. Two new time-stepping modules have been added:

* The ``ERKStep`` module provides an optimized implementation for explicit
  Runge-Kutta methods with reduced storage and number of calls to the ODE
  right-hand side function.

* The ``MRIStep`` module implements two-rate explicit-explicit multirate
  infinitesimal step methods utilizing different step sizes for slow
  and fast processes in an additive splitting.

This restructure has resulted in numerous small changes to the user
interface, particularly the suite of "Set" routines for user-provided
solver parameters and and "Get" routines to access solver statistics,
that are now prefixed with the name of time-stepping module (e.g., ``ARKStep``
or ``ERKStep``) instead of ``ARKode``.  Aside from affecting the names of these
routines, user-level changes have been kept to a minimum.  However, we recommend
that users consult both this documentation and the ARKode example programs for
further details on the updated infrastructure.

As part of the ARKode restructuring an :c:type:`ARKodeButcherTable` structure
has been added for storing Butcher tables. Functions for creating new Butcher
tables and checking their analytic order are provided along with other utility
routines. For more details see :ref:`ARKodeButcherTable`.

Two changes were made in the initial step size algorithm:

* Fixed an efficiency bug where an extra call to the right hand side function was made.

* Changed the behavior of the algorithm if the max-iterations case is hit.
  Before the algorithm would exit with the step size calculated on the
  penultimate iteration. Now it will exit with the step size calculated
  on the final iteration.

ARKode's dense output infrastructure has been improved to support
higher-degree Hermite polynomial interpolants (up to degree 5) over
the last successful time step.

ARKode's previous direct and iterative linear solver interfaces, ARKDLS and
ARKSPILS, have been merged into a single unified linear solver interface, ARKLS,
to support any valid SUNLINSOL module. This includes ``DIRECT`` and
``ITERATIVE`` types as well as the new ``MATRIX_ITERATIVE`` type. Details
regarding how ARKLS utilizes linear solvers of each type as well as discussion
regarding intended use cases for user-supplied SUNLinSol implementations are
included in the chapter :ref:`SUNLinSol`. All ARKode examples programs and the
standalone linear solver examples have been updated to use the unified linear
solver interface.

The user interface for the new ARKLS module is very similar to the previous
ARKDLS and ARKSPILS interfaces. Additionally, we note that Fortran users will
need to enlarge their ``iout`` array of optional integer outputs, and update the
indices that they query for certain linear-solver-related statistics.

The names of all constructor routines for SUNDIALS-provided SUNLinSol
implementations have been updated to follow the naming convention
``SUNLinSol_*`` where ``*`` is the name of the linear solver. The new names are
``SUNLinSol_Band``, ``SUNLinSol_Dense``, ``SUNLinSol_KLU``,
``SUNLinSol_LapackBand``, ``SUNLinSol_LapackDense``, ``SUNLinSol_PCG``,
``SUNLinSol_SPBCGS``, ``SUNLinSol_SPFGMR``, ``SUNLinSol_SPGMR``,
``SUNLinSol_SPTFQMR``, and ``SUNLinSol_SuperLUMT``.  Solver-specific "set"
routine names have been similarly standardized.  To minimize challenges in user
migration to the new names, the previous routine names may still be used; these
will be deprecated in future releases, so we recommend that users migrate to the
new names soon. All ARKode example programs and the standalone linear solver
examples have been updated to use the new naming convention.

The ``SUNBandMatrix`` constructor has been simplified to remove the
storage upper bandwidth argument.

SUNDIALS integrators have been updated to utilize generic nonlinear solver
modules defined through the SUNNONLINSOL API. This API will ease the addition of
new nonlinear solver options and allow for external or user-supplied nonlinear
solvers. The SUNNONLINSOL API and SUNDIALS provided modules are described in
:ref:`SUNNonlinSol` and follow the same object oriented design and
implementation used by the NVector, SUNMatrix, and SUNLinSol modules. Currently
two SUNNONLINSOL implementations are provided, SUNNonlinSol_Newton and
SUNNonlinSol_FixedPoint. These replicate the previous integrator specific
implementations of a Newton iteration and an accelerated fixed-point iteration,
respectively. Example programs using each of these nonlinear solver modules in a
standalone manner have been added and all ARKode example programs have been
updated to use generic SUNNonlinSol modules.

As with previous versions, ARKode will use the Newton solver (now
provided by SUNNonlinSol_Newton) by default.  Use of the
:c:func:`ARKStepSetLinear()` routine (previously named
``ARKodeSetLinear``) will indicate that the problem is
linearly-implicit, using only a single Newton iteration per implicit
stage.  Users wishing to switch to the accelerated fixed-point solver
are now required to create a SUNNonlinSol_FixedPoint object and attach
that to ARKode, instead of calling the previous
``ARKodeSetFixedPoint`` routine.  See the documentation sections
:ref:`ARKStep_CInterface.Skeleton`,
:ref:`ARKStep_CInterface.NonlinearSolvers`, and
:ref:`SUNNonlinSol_FixedPoint` for further details, or the serial C
example program ``ark_brusselator_fp.c`` for an example.

Three fused vector operations and seven vector array operations have been added
to the NVECTOR API. These *optional* operations are disabled by default and may
be activated by calling vector specific routines after creating an NVector (see
:ref:`NVectors.Description` for more details). The new operations are intended
to increase data reuse in vector operations, reduce parallel communication on
distributed memory systems, and lower the number of kernel launches on systems
with accelerators. The fused operations are ``N_VLinearCombination``,
``N_VScaleAddMulti``, and ``N_VDotProdMulti``, and the vector array operations
are ``N_VLinearCombinationVectorArray``, ``N_VScaleVectorArray``,
``N_VConstVectorArray``, ``N_VWrmsNormVectorArray``,
``N_VWrmsNormMaskVectorArray``, ``N_VScaleAddMultiVectorArray``, and
``N_VLinearCombinationVectorArray``. If an NVector implementation defines any of
these operations as ``NULL``, then standard NVector operations will
automatically be called as necessary to complete the computation.

Multiple changes to the CUDA NVECTOR were made:

  * Changed the ``N_VMake_Cuda`` function to take a host data pointer and a device
    data pointer instead of an ``N_VectorContent_Cuda`` object.

  * Changed ``N_VGetLength_Cuda`` to return the global vector length instead of
    the local vector length.

  * Added ``N_VGetLocalLength_Cuda`` to return the local vector length.

  * Added ``N_VGetMPIComm_Cuda`` to return the MPI communicator used.

  * Removed the accessor functions in the namespace ``suncudavec``.

  * Added the ability to set the ``cudaStream_t`` used for execution of the CUDA
    NVECTOR kernels. See the function ``N_VSetCudaStreams_Cuda``.

  * Added ``N_VNewManaged_Cuda``, ``N_VMakeManaged_Cuda``, and ``N_VIsManagedMemory_Cuda``
    functions to accommodate using managed memory with the CUDA NVECTOR.

Multiple changes to the RAJA NVECTOR were made:

  * Changed ``N_VGetLength_Raja`` to return the global vector length instead of
    the local vector length.

  * Added ``N_VGetLocalLength_Raja`` to return the local vector length.

  * Added ``N_VGetMPIComm_Raja`` to return the MPI communicator used.

  * Removed the accessor functions in the namespace ``sunrajavec``.

A new NVECTOR implementation for leveraging OpenMP 4.5+ device offloading has
been added, NVECTOR_OpenMPDEV. See :ref:`NVectors.OpenMPDEV` for more details.


Changes in v2.2.1
^^^^^^^^^^^^^^^^^^^^^^^

Fixed a bug in the CUDA NVECTOR where the ``N_VInvTest`` operation could
write beyond the allocated vector data.

Fixed library installation path for multiarch systems. This fix changes the default
library installation path to ``CMAKE_INSTALL_PREFIX/CMAKE_INSTALL_LIBDIR``
from ``CMAKE_INSTALL_PREFIX/lib``. ``CMAKE_INSTALL_LIBDIR`` is automatically
set, but is available as a CMAKE option that can modified.


Changes in v2.2.0
^^^^^^^^^^^^^^^^^^^^^^^

Fixed a problem with setting ``sunindextype`` which would occur with
some compilers (e.g. armclang) that did not define ``__STDC_VERSION__``.

Added hybrid MPI/CUDA and MPI/RAJA vectors to allow use of more than
one MPI rank when using a GPU system.  The vectors assume one GPU
device per MPI rank.

Changed the name of the RAJA NVECTOR library to
``libsundials_nveccudaraja.lib`` from
``libsundials_nvecraja.lib`` to better reflect that we only support CUDA
as a backend for RAJA currently.

Several changes were made to the build system:

* CMake 3.1.3 is now the minimum required CMake version.

* Deprecate the behavior of the ``SUNDIALS_INDEX_TYPE`` CMake option and
  added the ``SUNDIALS_INDEX_SIZE`` CMake option to select the ``sunindextype``
  integer size.

* The native CMake FindMPI module is now used to locate an MPI
  installation.

* If MPI is enabled and MPI compiler wrappers are not set, the build system
  will check if ``CMAKE_<language>_COMPILER`` can compile MPI programs before
  trying to locate and use an MPI installation.

* The previous options for setting MPI compiler wrappers and the executable
  for running MPI programs have been have been depreated. The new options that
  align with those used in native CMake FindMPI module are
  ``MPI_C_COMPILER``, ``MPI_CXX_COMPILER``, ``MPI_Fortran_COMPILER``,
  and ``MPIEXEC_EXECUTABLE``.

* When a Fortran name-mangling scheme is needed (e.g., ``LAPACK_ENABLE``
  is ``ON``) the build system will infer the scheme from the Fortran
  compiler. If a Fortran compiler is not available or the inferred or default
  scheme needs to be overridden, the advanced options
  ``SUNDIALS_F77_FUNC_CASE`` and ``SUNDIALS_F77_FUNC_UNDERSCORES`` can
  be used to manually set the name-mangling scheme and bypass trying to infer
  the scheme.

* Parts of the main CMakeLists.txt file were moved to new files in the
  ``src`` and ``example`` directories to make the CMake configuration file
  structure more modular.



Changes in v2.1.2
^^^^^^^^^^^^^^^^^^^^^^^

Updated the minimum required version of CMake to 2.8.12 and enabled
using rpath by default to locate shared libraries on OSX.

Fixed Windows specific problem where sunindextype was not correctly
defined when using 64-bit integers for the SUNDIALS index type. On Windows
sunindextype is now defined as the MSVC basic type ``__int64``.

Added sparse SUNMatrix "Reallocate" routine to allow specification of
the nonzero storage.

Updated the KLU SUNLinearSolver module to set constants for the two
reinitialization types, and fixed a bug in the full reinitialization
approach where the sparse SUNMatrix pointer would go out of scope on
some architectures.

Updated the "ScaleAdd" and "ScaleAddI" implementations in the
sparse SUNMatrix module to more optimally handle the case where the
target matrix contained sufficient storage for the sum, but had the
wrong sparsity pattern.  The sum now occurs in-place, by performing
the sum backwards in the existing storage.  However, it is still more
efficient if the user-supplied Jacobian routine allocates storage for
the sum :math:`I+\gamma J` or :math:`M+\gamma J` manually (with zero
entries if needed).

Changed LICENSE install path to ``instdir/include/sundials``.



Changes in v2.1.1
^^^^^^^^^^^^^^^^^^^

Fixed a potential memory leak in the SPGMR and SPFGMR linear solvers:
if "Initialize" was called multiple times then the solver memory was
reallocated (without being freed).

Fixed a minor bug in the ARKReInit routine, where a flag was
incorrectly set to indicate that the problem had been resized (instead
of just re-initialized).

Fixed C++11 compiler errors/warnings about incompatible use of string
literals.

Updated KLU SUNLinearSolver module to use a ``typedef`` for the
precision-specific solve function to be used (to avoid compiler
warnings).

Added missing typecasts for some ``(void*)`` pointers (again, to avoid
compiler warnings).

Bugfix in ``sunmatrix_sparse.c`` where we had used ``int`` instead of
``sunindextype`` in one location.

Added missing ``#include <stdio.h>`` in NVECTOR and SUNMATRIX header files.

Added missing prototype for ``ARKSpilsGetNumMTSetups``.

Fixed an indexing bug in the CUDA NVECTOR implementation of
``N_VWrmsNormMask`` and revised the RAJA NVECTOR implementation of
``N_VWrmsNormMask`` to work with mask arrays using values other than
zero or one. Replaced ``double`` with ``realtype`` in the RAJA vector
test functions.

Fixed compilation issue with GCC 7.3.0 and Fortran programs that do
not require a SUNMatrix or SUNLinearSolver module (e.g. iterative
linear solvers, explicit methods, fixed point solver, etc.).


Changes in v2.1.0
^^^^^^^^^^^^^^^^^^^

Added NVECTOR print functions that write vector data to a specified
file (e.g. ``N_VPrintFile_Serial``).

Added ``make test`` and ``make test_install`` options to the build
system for testing SUNDIALS after building with ``make`` and
installing with ``make install`` respectively.


Changes in v2.0.0
^^^^^^^^^^^^^^^^^^^

All interfaces to matrix structures and linear solvers have been
reworked, and all example programs have been updated.  The goal of the
redesign of these interfaces was to provide more encapsulation and
ease in interfacing custom linear solvers and interoperability with
linear solver libraries.

Specific changes include:

* Added generic SUNMATRIX module with three provided implementations:
  dense, banded and sparse.  These replicate previous SUNDIALS Dls and
  Sls matrix structures in a single object-oriented API.

* Added example problems demonstrating use of generic SUNMATRIX modules.

* Added generic SUNLINEARSOLVER module with eleven provided
  implementations: dense, banded, LAPACK dense, LAPACK band, KLU,
  SuperLU_MT, SPGMR, SPBCGS, SPTFQMR, SPFGMR, PCG.  These replicate
  previous SUNDIALS generic linear solvers in a single object-oriented
  API.

* Added example problems demonstrating use of generic SUNLINEARSOLVER modules.

* Expanded package-provided direct linear solver (Dls) interfaces and
  scaled, preconditioned, iterative linear solver (Spils) interfaces
  to utilize generic SUNMATRIX and SUNLINEARSOLVER objects.

* Removed package-specific, linear solver-specific, solver modules
  (e.g. CVDENSE, KINBAND, IDAKLU, ARKSPGMR) since their functionality
  is entirely replicated by the generic Dls/Spils interfaces and
  SUNLINEARSOLVER/SUNMATRIX modules.  The exception is CVDIAG, a
  diagonal approximate Jacobian solver available to CVODE and CVODES.

* Converted all SUNDIALS example problems to utilize new generic
  SUNMATRIX and SUNLINEARSOLVER objects, along with updated Dls and
  Spils linear solver interfaces.

* Added Spils interface routines to ARKode, CVODE, CVODES, IDA and
  IDAS to allow specification of a user-provided "JTSetup" routine.
  This change supports users who wish to set up data structures for
  the user-provided Jacobian-times-vector ("JTimes") routine, and
  where the cost of one JTSetup setup per Newton iteration can be
  amortized between multiple JTimes calls.

Two additional NVECTOR implementations were added -- one for CUDA and
one for RAJA vectors.  These vectors are supplied to provide very
basic support for running on GPU architectures.  Users are advised
that these vectors both move all data to the GPU device upon
construction, and speedup will only be realized if the user also
conducts the right-hand-side function evaluation on the device. In
addition, these vectors assume the problem fits on one GPU. Further
information about RAJA, users are referred to the web site,
`https://software.llnl.gov/RAJA/ <https://software.llnl.gov/RAJA/>`_.
These additions are accompanied by additions to various interface
functions and to user documentation.

All indices for data structures were updated to a new ``sunindextype``
that can be configured to be a 32- or 64-bit integer data index type.
``sunindextype`` is defined to be ``int32_t`` or ``int64_t`` when
portable types are supported, otherwise it is defined as ``int`` or
``long int``. The Fortran interfaces continue to use ``long int`` for
indices, except for their sparse matrix interface that now uses the
new ``sunindextype``.  This new flexible capability for index types
includes interfaces to PETSc, *hypre*, SuperLU_MT, and KLU with either
32-bit or 64-bit capabilities depending how the user configures
SUNDIALS.

To avoid potential namespace conflicts, the macros defining
``booleantype`` values ``TRUE`` and ``FALSE`` have been changed to
``SUNTRUE`` and ``SUNFALSE`` respectively.

Temporary vectors were removed from preconditioner setup and solve
routines for all packages.  It is assumed that all necessary data
for user-provided preconditioner operations will be allocated and
stored in user-provided data structures.

The file ``include/sundials_fconfig.h`` was added.  This file contains
SUNDIALS type information for use in Fortran programs.

Added functions SUNDIALSGetVersion and SUNDIALSGetVersionNumber to get
SUNDIALS release version information at runtime.

The build system was expanded to support many of the xSDK-compliant keys.
The xSDK is a movement in scientific software to provide a foundation for the
rapid and efficient production of high-quality,
sustainable extreme-scale scientific applications.  More information can
be found at, `https://xsdk.info <https://xsdk.info>`_.

In addition, numerous changes were made to the build system.
These include the addition of separate ``BLAS_ENABLE`` and ``BLAS_LIBRARIES``
CMake variables, additional error checking during CMake configuration,
minor bug fixes, and renaming CMake options to enable/disable examples
for greater clarity and an added option to enable/disable Fortran 77 examples.
These changes included changing ``ENABLE_EXAMPLES`` to ``ENABLE_EXAMPLES_C``,
changing ``CXX_ENABLE`` to ``EXAMPLES_ENABLE_CXX``, changing ``F90_ENABLE`` to
``EXAMPLES_ENABLE_F90``, and adding an ``EXAMPLES_ENABLE_F77`` option.

Corrections and additions were made to the examples, to
installation-related files, and to the user documentation.




Changes in v1.1.0
^^^^^^^^^^^^^^^^^^^

We have included numerous bugfixes and enhancements since the
v1.0.2 release.

The bugfixes include:

* For each linear solver, the various solver performance counters are
  now initialized to 0 in both the solver specification function and
  in the solver's ``linit`` function.  This ensures that these solver
  counters are initialized upon linear solver instantiation as well as
  at the beginning of the problem solution.

* The choice of the method vs embedding the Billington and TRBDF2
  explicit Runge-Kutta methods were swapped, since in those the
  lower-order coefficients result in an A-stable method, while the
  higher-order coefficients do not.  This change results in
  significantly improved robustness when using those methods.

* A bug was fixed for the situation where a user supplies a vector of
  absolute tolerances, and also uses the vector Resize() functionality.

* A bug was fixed wherein a user-supplied Butcher table without an
  embedding is supplied, and the user is running with either fixed
  time steps (or they do adaptivity manually); previously this had
  resulted in an error since the embedding order was below 1.

* Numerous aspects of the documentation were fixed and/or clarified.


The feature changes/enhancements include:

* Two additional NVECTOR implementations were added -- one for Hypre
  (parallel) ParVector vectors, and one for PETSc vectors.  These
  additions are accompanied by additions to various interface
  functions and to user documentation.

* Each NVECTOR module now includes a function, ``N_VGetVectorID``,
  that returns the NVECTOR module name.

* A memory leak was fixed in the banded preconditioner and
  banded-block-diagonal preconditioner interfaces.  In addition,
  updates were done to return integers from linear solver and
  preconditioner 'free' routines.

* The Krylov linear solver Bi-CGstab was enhanced by removing a
  redundant dot product.  Various additions and corrections were made
  to the interfaces to the sparse solvers KLU and SuperLU_MT,
  including support for CSR format when using KLU.

* The ARKode implicit predictor algorithms were updated: methods 2 and
  3 were improved slightly, a new predictor approach was added, and
  the default choice was modified.

* The underlying sparse matrix structure was enhanced to allow both
  CSR and CSC matrices, with CSR supported by the KLU linear solver
  interface.  ARKode interfaces to the KLU solver from both C and
  Fortran were updated to enable selection of sparse matrix type, and a
  Fortran-90 CSR example program was added.

* The missing :c:func:`ARKSpilsGetNumMtimesEvals()` function was added
  -- this had been included in the previous documentation but had not
  been implemented.

* The handling of integer codes for specifying built-in ARKode Butcher
  tables was enhanced.  While a global numbering system is still used,
  methods now have #defined names to simplify the user interface and to
  streamline incorporation of new Butcher tables into ARKode.

* The maximum number of Butcher table stages was increased from 8 to
  15 to accommodate very high order methods, and an 8th-order adaptive
  ERK method was added.

* Support was added for the explicit and implicit methods in an
  additive Runge-Kutta method to utilize different stage times,
  solution and embedding coefficients, to support new SSP-ARK
  methods.

* The FARKODE interface was extended to include a routine to set
  scalar/array-valued residual tolerances, to support Fortran
  applications with non-identity mass-matrices.







Reading this User Guide
----------------------------

This user guide is a combination of general usage instructions and
specific example programs.  We expect that some readers will want to
concentrate on the general instructions, while others will refer
mostly to the examples, and the organization is intended to
accommodate both styles.

The structure of this document is as follows:

* In the next section we provide a thorough presentation of the
  underlying :ref:`mathematics <Mathematics>` used within the ARKode
  family of solvers.

* We follow this with an overview of how the source code for ARKode is
  :ref:`organized <Organization>`.

* The largest section follows, providing a full account of the
  ARKStep module user interface, including a description of all
  user-accessible functions and outlines for usage in serial and
  parallel applications. Since ARKode is written in C, we first
  present a section on :ref:`using ARKStep for C and C++ applications
  <ARKStep_CInterface>`, followed with a separate section on
  :ref:`using ARKode within Fortran applications <FortranInterface>`.

* The much smaller section describing the ERKStep time-stepping
  module, :ref:`using ERKStep for C and C++ applications
  <ERKStep_CInterface>`, follows.

* Subsequent sections discuss shared features between ARKode
  and the rest of the SUNDIALS library:
  :ref:`vector data structures <NVectors>`,
  :ref:`matrix data structures <SUNMatrix>`,
  :ref:`linear solver data structures <SUNLinSol>`, and the
  :ref:`installation procedure <Installation>`.

* The final sections catalog the full set of :ref:`ARKode constants
  <Constants>`, that are used for both input specifications and return
  codes, and the full set of :ref:`Butcher tables <Butcher>` that are
  packaged with ARKode.



SUNDIALS Release License
----------------------------

All SUNDIALS packages are released open source, under the BSD 3-Clause
license. The only requirements of the license are preservation of copyright and
a standard disclaimer of liability. The full text of the license and an
additional notice are provided below and may also be found in the LICENSE and
NOTICE files provided with all SUNDIALS packages.

**PLEASE NOTE**  If you are using SUNDIALS with any third party
libraries linked in (e.g., LAPACK, KLU, SuperLU_MT, PETSc, or
*hypre*), be sure to review the respective license of the package as
that license may have more restrictive terms than the SUNDIALS
license.  For example, if someone builds SUNDIALS with a statically
linked KLU, the build is subject to terms of the more-restrictive LGPL
license (which is what KLU is released with) and *not* the SUNDIALS
BSD license anymore.


BSD 3-Clause License
^^^^^^^^^^^^^^^^^^^^

Copyright (c) 2002-2019, Lawrence Livermore National Security and Southern
Methodist University.

All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.

* Neither the name of the copyright holder nor the names of its
  contributors may be used to endorse or promote products derived from
  this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS ''AS IS''
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

Additional Notice
^^^^^^^^^^^^^^^^^

This work was produced under the auspices of the U.S. Department of
Energy by Lawrence Livermore National Laboratory under Contract
DE-AC52-07NA27344.

This work was prepared as an account of work sponsored by an agency of
the United States Government. Neither the United States Government nor
Lawrence Livermore National Security, LLC, nor any of their employees
makes any warranty, expressed or implied, or assumes any legal liability
or responsibility for the accuracy, completeness, or usefulness of any
information, apparatus, product, or process disclosed, or represents that
its use would not infringe privately owned rights.

Reference herein to any specific commercial product, process, or service
by trade name, trademark, manufacturer, or otherwise does not necessarily
constitute or imply its endorsement, recommendation, or favoring by the
United States Government or Lawrence Livermore National Security, LLC.

The views and opinions of authors expressed herein do not necessarily
state or reflect those of the United States Government or Lawrence
Livermore National Security, LLC, and shall not be used for advertising
or product endorsement purposes.

SUNDIALS Release Numbers
^^^^^^^^^^^^^^^^^^^^^^^^
LLNL-CODE-667205  (ARKODE)

UCRL-CODE-155951  (CVODE)

UCRL-CODE-155950  (CVODES)

UCRL-CODE-155952  (IDA)

UCRL-CODE-237203  (IDAS)

LLNL-CODE-665877  (KINSOL)
