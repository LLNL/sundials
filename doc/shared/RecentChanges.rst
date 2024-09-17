**Major Features**

**New Features and Enhancements**

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
recently introduced `SUNAdaptController` base class, to ARKODE's MRIStep module.
Added new default MRI methods for temporally adaptive versus fixed-step runs.

Added functionality to ARKODE to accumulate a temporal error
estimate over multiple time steps.  See the routines
:c:func:`ARKodeSetAccumulatedErrorType`, :c:func:`ARKodeResetAccumulatedError`,
and :c:func:`ARKodeGetAccumulatedError` for details.

**Bug Fixes**

Fixed c:func:`ARKodeResize` not using the default ``hscale`` when an argument of
``0`` was provided.

Fixed the loading of ARKStep's default first order explicit method.

Fixed a CMake bug regarding usage of missing "print_warning" macro
that was only triggered when the deprecated ``CUDA_ARCH`` option was used.

Fixed a memory leak that could occur if :c:func:`ARKodeSetDefaults` is called
repeatedly.

Fixed compilation errors when building the Trilinos Teptra NVector with CUDA
support.

Fixed loading the default IMEX-MRI method if :c:func:`ARKodeSetOrder` is used to
specify a third or fourth order method. Previously, the default second order method
was loaded in both cases.


**Deprecation Notices**
