**New Features**

Added CMake infrastructure that enables externally maintained addons/plugins
to be *optionally* built with SUNDIALS. See :ref:`Contributing` for details.

Added Multirate time step adaptivity controllers, based on the recently introduced
`SUNAdaptController` base class, to ARKODE's MRIStep module.

Added functionality to ARKStep and ERKStep to accumulate a temporal error
estimate over multiple time steps.  See the routines :c:func:`ARKStepSetAccumulatedErrorType`,
:c:func:`ARKStepResetAccumulatedError`, :c:func:`ARKStepGetAccumulatedError`,
:c:func:`ERKStepSetAccumulatedErrorType`, :c:func:`ERKStepResetAccumulatedError`,
and :c:func:`ERKStepGetAccumulatedError` for details.

**Bug Fixes**

Updated the CMake variable ``HIP_PLATFORM`` default to ``amd`` as the previous
default, ``hcc``, is no longer recognized in ROCm 5.7.0 or newer. The new
default is also valid in older version of ROCm (at least back to version 4.3.1).

Fixed a bug in the HIP execution policies where ``WARP_SIZE`` would not be set
with ROCm 6.0.0 or newer.

Changed the CMake version compatibility mode for SUNDIALS to ``AnyNewerVersion``
instead of ``SameMajorVersion``. This fixes the issue seen
`here <https://github.com/AMReX-Codes/amrex/pull/3835>`_.

Fixed a bug in some Fortran examples where ``c_null_ptr`` was passed as an argument
to a function pointer instead of ``c_null_funptr``. This caused compilation issues
with the Cray Fortran compiler.

Fixed a bug where :c:func:`MRIStepEvolve` would not handle a recoverable error
produced from evolving the inner stepper.

Added support for Kokkos Kernels v4.

Fixed a bug that caused error messages to be cut off in some cases. Fixes `GitHub Issue #461 <https://github.com/LLNL/sundials/issues/461>`_.
