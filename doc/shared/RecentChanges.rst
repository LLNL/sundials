**Major Features**

**New Features and Enhancements**

The default value of :cmakeop:`CMAKE_CUDA_ARCHITECTURES` is no longer set to
``70`` and is now determined automatically by CMake. The previous default was
only valid for Volta GPUs while the automatically selected value will vary
across compilers and compiler versions. As such, users are encouraged to
override this value with the architecture for their system.

Added support for multirate time step adaptivity controllers, based on the
recently introduced `SUNAdaptController` base class, to ARKODE's MRIStep module.
Added new default MRI methods for temporally adaptive versus fixed-step runs.

Added functionality to ARKODE to accumulate a temporal error
estimate over multiple time steps.  See the routines
:c:func:`ARKodeSetAccumulatedErrorType`, :c:func:`ARKodeResetAccumulatedError`,
and :c:func:`ARKodeGetAccumulatedError` for details.

**Bug Fixes**

Fixed the loading of ARKStep's default first order explicit method.

Fixed a CMake bug regarding usage of missing "print_warning" macro
that was only triggered when the deprecated ``CUDA_ARCH`` option was used.

Fixed loading the default IMEX-MRI method if :c:func:`ARKodeSetOrder` is used to
specify a third or fourth order method. Previously, the default second order method
was loaded in both cases.

**Deprecation Notices**
