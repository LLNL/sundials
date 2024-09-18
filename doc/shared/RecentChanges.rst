**Major Features**

Added an operator splitting module,
:ref:`SplittingStep <ARKODE.Usage.SplittingStep>`, and forcing method module,
:ref:`ForcingStep <ARKODE.Usage.ForcingStep>`, to ARKODE. These modules support
a broad range of operator-split time integration methods for multiphysics
applications.

**New Features and Enhancements**

Added the :c:func:`ARKodeSetStepDirection` and :c:func:`ARKodeGetStepDirection`
functions to change and query the direction of integration.

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

**Deprecation Notices**
