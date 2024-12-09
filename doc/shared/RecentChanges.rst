**Major Features**

Added a time-stepping module to ARKODE for low storage Runge--Kutta methods, 
:ref:`LSRKStep <ARKODE.Usage.LSRKStep>`.  This currently supports five explicit low-storage 
methods: the second-order Runge--Kutta--Chebyshev and Runge--Kutta--Legendre methods, 
and the second- through fourth-order optimal strong stability preserving Runge--Kutta methods.  
All methods include embeddings for temporal adaptivity.

Added an operator splitting module,
:ref:`SplittingStep <ARKODE.Usage.SplittingStep>`, and forcing method module,
:ref:`ForcingStep <ARKODE.Usage.ForcingStep>`, to ARKODE. These modules support
a broad range of operator-split time integration methods for multiphysics
applications.

**New Features and Enhancements**

Added the :c:func:`ARKodeSetStepDirection` and :c:func:`ARKodeGetStepDirection`
functions to change and query the direction of integration.

Added the :c:type:`SUNStepper` base class to represent a generic solution
procedure for IVPs. This is used by the
:ref:`SplittingStep <ARKODE.Usage.SplittingStep>` and
:ref:`ForcingStep <ARKODE.Usage.ForcingStep>` modules of ARKODE. A SUNStepper
can be created from an ARKODE memory block with the new function
:c:func:`ARKodeCreateSUNStepper`. To enable interoperability with
:c:type:`MRIStepInnerStepper`, the function
:c:func:`MRIStepInnerStepper_CreateFromSUNStepper` was added.

The following DIRK schemes now have coefficients accurate to quad precision:

* ``ARKODE_BILLINGTON_3_3_2``

* ``ARKODE_KVAERNO_4_2_3``

* ``ARKODE_CASH_5_2_4``

* ``ARKODE_CASH_5_3_4``

* ``ARKODE_KVAERNO_5_3_4``

* ``ARKODE_KVAERNO_7_4_5``

The default value of :cmakeop:`CMAKE_CUDA_ARCHITECTURES` is no longer set to
``70`` and is now determined automatically by CMake. The previous default was
only valid for Volta GPUs while the automatically selected value will vary
across compilers and compiler versions. As such, users are encouraged to
override this value with the architecture for their system.

The Trilinos Tpetra NVector interface has been updated to utilize CMake
imported targets added in Trilinos 14 to improve support for different Kokkos
backends with Trilinos. As such, Trilinos 14 or newer is required and the
``Trilinos_INTERFACE_*`` CMake options have been removed.

Example programs using *hypre* have been updated to support v2.20 and newer.

The build system has been updated to utilize the CMake LAPACK imported target
which should ease building SUNDIALS with LAPACK libraries that require setting
specific linker flags e.g., MKL.

Added support for multirate time step adaptivity controllers, based on the
recently introduced :c:type:`SUNAdaptController` base class, to ARKODE's MRIStep module.
As a part of this, we added embeddings for existing MRI-GARK methods, as well as
support for embedded MERK and IMEX-MRI-SR methods.  Added new default MRI methods
for temporally adaptive versus fixed-step runs.  Added the function
:c:func:`MRIStepGetNumInnerStepperFails` to retrieve the number of recoverable
failures reported by the MRIStepInnerStepper.

Added functionality to ARKODE to accumulate a temporal error
estimate over multiple time steps.  See the routines
:c:func:`ARKodeSetAccumulatedErrorType`, :c:func:`ARKodeResetAccumulatedError`,
and :c:func:`ARKodeGetAccumulatedError` for details.

Added a utility routine to wrap any valid ARKODE integrator for use as an MRIStep
inner stepper object, :c:func:`ARKodeCreateMRIStepInnerStepper`.

**Bug Fixes**

Fixed a bug where :c:func:`CVodeSetProjFailEta` would ignore the `eta`
parameter.

Fixed a bug in the SPTFQMR linear solver where recoverable preconditioner errors
were reported as unrecoverable.

Fixed a `bug <https://github.com/LLNL/sundials/issues/581>`__ in the sparse
matrix implementation of :c:func:`SUNMatScaleAddI` which caused out of bounds
writes unless ``indexvals`` were in ascending order for each row/column.

Fixed :c:func:`ARKodeResize` not using the default ``hscale`` when an argument
of ``0`` was provided.

Fixed the loading of ARKStep's default first order explicit method.

Fixed a bug in ARKODE when enabling rootfinding with fixed step sizes and the
initial value of the rootfinding function is zero. In this case, uninitialized
right-hand side data was used to compute a state value near the initial
condition to determine if any rootfinding functions are initially active.

Fixed a CMake bug regarding usage of missing "print_warning" macro
that was only triggered when the deprecated ``CUDA_ARCH`` option was used.

Fixed a memory leak that could occur if :c:func:`ARKodeSetDefaults` is called
repeatedly.

Fixed compilation errors when building the Trilinos Teptra NVector with CUDA
support.

Fixed loading the default IMEX-MRI method if :c:func:`ARKodeSetOrder` is used to
specify a third or fourth order method. Previously, the default second order method
was loaded in both cases.

Fixed a bug in MRIStep where the data supplied to the Hermite interpolation module did
not include contributions from the fast right-hand side function. With this fix, users
will see one additional fast right-hand side function evaluation per slow step with the
Hermite interpolation option.

Fixed potential memory leaks and out of bounds array accesses that could occur
in the ARKODE Lagrange interpolation module when changing the method order or
polynomial degree after re-initializing an integrator.

Fixed a CMake configuration issue related to aliasing an ``ALIAS`` target when
using ``ENABLE_KLU=ON`` in combination with a static-only build of SuiteSparse.

Fixed a CMake issue which caused third-party CMake variables to be unset.
Users may see more options in the CMake GUI now as a result of the fix.
See details in GitHub Issue `#538 <https://github.com/LLNL/sundials/issues/538>`__.

**Deprecation Notices**

Deprecated the ARKStep-specific utility routine for wrapping an ARKStep instance
as an MRIStep inner stepper object, :c:func:`ARKStepCreateMRIStepInnerStepper`. Use
:c:func:`ARKodeCreateMRIStepInnerStepper` instead.

The ARKODE stepper specific functions to retrieve the number of right-hand side
function evaluations have been deprecated. Use :c:func:`ARKodeGetNumRhsEvals`
instead.
