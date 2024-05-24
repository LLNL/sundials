**Major Features**

Created shared user interface functions for ARKODE to allow more uniform control
over time-stepping algorithms, improved extensibility, and simplified code
maintenance. The corresponding stepper-specific user-callable functions are now
deprecated and will be removed in a future major release.

Added CMake infrastructure that enables externally maintained addons/plugins to
be *optionally* built with SUNDIALS. See :ref:`Contributing` for details.

**New Features and Enhancements**

Added support for Kokkos Kernels v4.

Added the following Runge-Kutta Butcher tables

* ``ARKODE_FORWARD_EULER_1_1``
* ``ARKODE_RALSTON_EULER_2_1_2``
* ``ARKODE_EXPLICIT_MIDPOINT_EULER_2_1_2``
* ``ARKODE_BACKWARD_EULER_1_1``
* ``ARKODE_IMPLICIT_MIDPOINT_1_2``
* ``ARKODE_IMPLICIT_TRAPEZOIDAL_2_2``

Added the following MRI coupling tables

* ``ARKODE_MRI_GARK_FORWARD_EULER``
* ``ARKODE_MRI_GARK_RALSTON2``
* ``ARKODE_MRI_GARK_RALSTON3``
* ``ARKODE_MRI_GARK_BACKWARD_EULER``
* ``ARKODE_MRI_GARK_IMPLICIT_MIDPOINT``
* ``ARKODE_IMEX_MRI_GARK_EULER``
* ``ARKODE_IMEX_MRI_GARK_TRAPEZOIDAL``
* ``ARKODE_IMEX_MRI_GARK_MIDPOINT``

Added :c:func:`ARKodeButcherTable_ERKIDToName` and
:c:func:`ARKodeButcherTable_DIRKIDToName` to convert a Butcher table ID to a
string representation.

Users may now disable interpolated output in ARKODE by passing
``ARK_INTERP_NONE`` to :c:func:`ARKodeSetInterpolantType`. When interpolation is
disabled, rootfinding is not supported, implicit methods must use the trivial
predictor (the default option), and interpolation at stop times cannot be used
(interpolating at stop times is disabled by default). With interpolation
disabled, calling :c:func:`ARKodeEvolve` in ``ARK_NORMAL`` mode will return at
or past the requested output time (setting a stop time may still be used to halt
the integrator at a specific time). Disabling interpolation will reduce the
memory footprint of an integrator by two or more state vectors (depending on the
interpolant type and degree) which can be beneficial when interpolation is not
needed e.g., when integrating to a final time without output in between or using
an explicit fast time scale integrator with an MRI method.

Added "Resize" capability to ARKODE's SPRKStep time-stepping module.

**Bug Fixes**

Updated the CMake variable ``HIP_PLATFORM`` default to ``amd`` as the previous
default, ``hcc``, is no longer recognized in ROCm 5.7.0 or newer. The new
default is also valid in older version of ROCm (at least back to version 4.3.1).

Renamed the DPCPP value for the :cmakeop:`SUNDIALS_GINKGO_BACKENDS` CMake option
to ``SYCL`` to match Ginkgo's updated naming convention.

Changed the CMake version compatibility mode for SUNDIALS to ``AnyNewerVersion``
instead of ``SameMajorVersion``. This fixes the issue seen `here
<https://github.com/AMReX-Codes/amrex/pull/3835>`_.

Fixed a CMake bug that caused an MPI linking error for our C++ examples in some
instances. Fixes `GitHub Issue #464
<https://github.com/LLNL/sundials/issues/464>`_.

Fixed the runtime library installation path for windows systems. This fix
changes the default library installation path from
``CMAKE_INSTALL_PREFIX/CMAKE_INSTALL_LIBDIR`` to
``CMAKE_INSTALL_PREFIX/CMAKE_INSTALL_BINDIR``.

Fixed conflicting ``.lib`` files between shared and static libs when using
``MSVC`` on Windows

Fixed invalid ``SUNDIALS_EXPORT`` generated macro when building both shared and
static libs.

Fixed a bug in some Fortran examples where ``c_null_ptr`` was passed as an
argument to a function pointer instead of ``c_null_funptr``. This caused
compilation issues with the Cray Fortran compiler.

Fixed a bug in the HIP execution policies where ``WARP_SIZE`` would not be set
with ROCm 6.0.0 or newer.

Fixed a bug that caused error messages to be cut off in some cases. Fixes
`GitHub Issue #461 <https://github.com/LLNL/sundials/issues/461>`_.

Fixed a memory leak when an error handler was added to a
:c:type:`SUNContext`. Fixes `GitHub Issue #466
<https://github.com/LLNL/sundials/issues/466>`_.

Fixed a bug where :c:func:`MRIStepEvolve` would not handle a recoverable error
produced from evolving the inner stepper.

Added missing ``SetRootDirection`` and ``SetNoInactiveRootWarn`` functions to
ARKODE's SPRKStep time-stepping module.

Fixed a bug in :c:func:`ARKodeSPRKTable_Create` where the coefficient arrays
were not allocated.

Fix bug on LLP64 platforms (like Windows 64-bit) where ``KLU_INDEXTYPE`` could be
32 bits wide even if ``SUNDIALS_INT64_T`` is defined.

Check if size of ``SuiteSparse_long`` is 8 if the size of ``sunindextype`` is 8
when using KLU.

**Deprecation Notices**

Numerous ARKODE stepper-specific functions are now deprecated in favor of
ARKODE-wide functions.

Deprecated the `ARKStepSetOptimalParams` function. Since this function does not have an
ARKODE-wide equivalent, instructions have been added to the user guide for how
to retain the current functionality using other user-callable functions.