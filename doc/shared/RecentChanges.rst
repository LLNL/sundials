**New Features**

Created shared user interface for ARKODE user-callable routines, to allow more
uniform control over time-stepping algorithms, improved extensibility, and
simplified code maintenance.  Marked the corresponding stepper-specific
user-callable routines as deprecated; these will be removed in a future major
release.

Added "Resize" capability, as well as missing ``SetRootDirection`` and
``SetNoInactiveRootWarn`` functions, to ARKODE's SPRKStep time-stepping module.

Deprecated ``ARKStepSetOptimalParams`` function; added instructions to user guide
for users who wish to retain the current functionality.

Added CMake infrastructure that enables externally maintained addons/plugins
to be *optionally* built with SUNDIALS. See :ref:`Contributing` for details.

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

Fixed a memory leak when an error handler was added to a :c:type:`SUNContext`. Fixes `GitHub Issue #466 <https://github.com/LLNL/sundials/issues/466>`_.

Fixed a CMake bug that caused an MPI linking error for our C++ examples in some instances. Fixes `GitHub Issue #464 <https://github.com/LLNL/sundials/issues/464>`_.

Fixed a bug in :c:func:`ARKodeSPRKTable_Create` where the coefficient arrays
where not allocated.

Fix bug on LLP64 platforms (like Windows 64-bit) where ``KLU_INDEXTYPE`` could be
32 bits wide even if ``SUNDIALS_INT64_T`` is defined.

Check if size of ``SuiteSparse_long`` is 8 if the size of ``sunindextype`` is 8
when using KLU.
