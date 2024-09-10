**Major Features**

**New Features and Enhancements**

The default value of :cmakeop:`CMAKE_CUDA_ARCHITECTURES` is no longer set to
``70`` and is now determined automatically by CMake. The previous default was
only valid for Volta GPUs while the automatically selected value will vary
across compilers and compiler versions. As such, users are encouraged to
override this value with the architecture for their system.

**Bug Fixes**

Fixed the loading of ARKStep's default first order explicit method.

Fixed loading the default IMEX-MRI method if :c:func:`ARKodeSetOrder` is used to
specify a third or fourth order method. Previously, the default second order
method was loaded in both cases.

Fixed a CMake bug regarding usage of missing "print_warning" macro
that was only triggered when the deprecated ``CUDA_ARCH`` option was used.

**Deprecation Notices**
