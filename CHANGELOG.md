# SUNDIALS Changelog

## Changes to SUNDIALS in release X.X.X

Updated the CMake variable `HIP_PLATFORM` default to `amd` as the previous
default, `hcc`, is no longer recognized in ROCm 5.7.0 or newer. The new default
is also valid in older version of ROCm (at least back to version 4.3.1).

Fixed a bug in the HIP execution policies where `WARP_SIZE` would not be set
with ROCm 6.0.0 or newer.

## Changes to SUNDIALS in release v7.0.0

### Major Feature

SUNDIALS now has more robust and uniform error handling. Non-release builds will
be built with additional error checking by default. See the
[Error Checking](https://sundials.readthedocs.io/en/latest/sundials/Errors_link.html)
section in the user guide for details.

### Breaking Changes

#### Minimum C Standard

SUNDIALS now requires using a compiler that supports a subset of the C99
standard. Note with the Microsoft C/C++ compiler the subset of C99 features
utilized by SUNDIALS are available starting with [Visual Studio 2015](https://learn.microsoft.com/en-us/cpp/overview/visual-cpp-language-conformance?view=msvc-170#c-standard-library-features-1).

#### Deprecated Types and Functions Removed

The previously deprecated types `realtype` and `booleantype` were removed from
`sundials_types.h` and replaced with `sunrealtype` and `sunbooleantype`. The
deprecated names for these types can be used by including the header file
`sundials_types_deprecated.h` but will be fully removed in the next major
release. Functions, types and header files that were previously deprecated have
also been removed.

#### Error Handling Changes

With the addition of the new error handling capability, the `*SetErrHandlerFn`
and `*SetErrFile` functions in CVODE(S), IDA(S), ARKODE, and KINSOL have been
removed. Users of these functions can use the functions
`SUNContext_PushErrHandler`, and `SUNLogger_SetErrorFilename` instead. For
further details see the
[Error Checking](https://sundials.readthedocs.io/en/latest/sundials/Errors_link.html)
and
[Logging](https://sundials.readthedocs.io/en/latest/sundials/Logging_link.html)
sections in the documentation.

In addition the following names/symbols were replaced by ``SUN_ERR_*`` codes:

| Removed                        | Replaced with ``SUNErrCode``      |
|:-------------------------------|:----------------------------------|
| `SUNLS_SUCCESS`                | `SUN_SUCCESS`                     |
| `SUNLS_UNRECOV_FAILURE`        | no replacement (value was unused) |
| `SUNLS_MEM_NULL`               | `SUN_ERR_ARG_CORRUPT`             |
| `SUNLS_ILL_INPUT`              | `SUN_ERR_ARG_*`                   |
| `SUNLS_MEM_FAIL`               | `SUN_ERR_MEM_FAIL`                |
| `SUNLS_PACKAGE_FAIL_UNREC`     | `SUN_ERR_EXT_FAIL`                |
| `SUNLS_VECTOROP_ERR`           | `SUN_ERR_OP_FAIL`                 |
| `SUN_NLS_SUCCESS`              | `SUN_SUCCESS`                     |
| `SUN_NLS_MEM_NULL`             | `SUN_ERR_ARG_CORRUPT`             |
| `SUN_NLS_MEM_FAIL`             | `SUN_ERR_MEM_FAIL`                |
| `SUN_NLS_ILL_INPUT`            | `SUN_ERR_ARG_*`                   |
| `SUN_NLS_VECTOROP_ERR`         | `SUN_ERR_OP_FAIL`                 |
| `SUN_NLS_EXT_FAIL`             | `SUN_ERR_EXT_FAIL`                |
| `SUNMAT_SUCCESS`               | `SUN_SUCCESS`                     |
| `SUNMAT_ILL_INPUT`             | `SUN_ERR_ARG_*`                   |
| `SUNMAT_MEM_FAIL`              | `SUN_ERR_MEM_FAIL`                |
| `SUNMAT_OPERATION_FAIL`        | `SUN_ERR_OP_FAIL`                 |
| `SUNMAT_MATVEC_SETUP_REQUIRED` | `SUN_ERR_OP_FAIL`                 |

The following functions have had their signature updated to ensure they can
leverage the new SUNDIALS error handling capabilities.

```c
// From sundials_futils.h
SUNDIALSFileOpen
SUNDIALSFileClose

// From sundials_memory.h
SUNMemorNewEmpty
SUNMemoryHelper_Alias
SUNMemoryHelper_Wrap

// From sundials_nvector.h
N_VNewVectorArray
```

#### SUNComm Type Added

We have replaced the use of a type-erased (i.e., `void*`) pointer to a
communicator in place of `MPI_Comm` throughout the SUNDIALS API with a
`SUNComm`, which is just a typedef to an `int` in builds without MPI
and a typedef to a `MPI_Comm` in builds with MPI. As a result:

- All users will need to update their codes because the call to
  `SUNContext_Create` now takes a `SUNComm` instead
  of type-erased pointer to a communicator. For non-MPI codes,
  pass `SUN_COMM_NULL` to the `comm` argument instead of
  `NULL`. For MPI codes, pass the `MPI_Comm` directly.

- The same change must be made for calls to
  `SUNLogger_Create` or `SUNProfiler_Create`.

- Some users will need to update their calls to `N_VGetCommunicator`, and
  update any custom `N_Vector` implementations that provide
  `N_VGetCommunicator`, since it now returns a `SUNComm`.

The change away from type-erased pointers for `SUNComm` fixes problems like the
one described in [GitHub Issue #275](https://github.com/LLNL/sundials/issues/275).

The SUNLogger is now always MPI-aware if MPI is enabled in SUNDIALS and the
`SUNDIALS_LOGGING_ENABLE_MPI` CMake option and macro definition were removed
accordingly.

#### SUNDIALS Core Library

Users now need to link to `sundials_core` in addition to the libraries already
linked to. This will be picked up automatically in projects that use the
SUNDIALS CMake target. The library `sundials_generic` has been superseded by
`sundials_core` and is no longer available. This fixes some duplicate symbol
errors on Windows when linking to multiple SUNDIALS libraries.

#### Fortran Interface Modules Streamlined

We have streamlined the Fortran modules that need to be included by users by combining
the SUNDIALS core into one Fortran module, `fsundials_core_mod`. Modules for
implementations of the core APIs still exist (e.g., for the Dense linear solver there
is `fsunlinsol_dense_mod`) as do the modules for the SUNDIALS packages (e.g., `fcvode_mod`).
The following modules are the ones that have been consolidated into `fsundials_core_mod`:

```
fsundials_adaptcontroller_mod
fsundials_context_mod
fsundials_futils_mod
fsundials_linearsolver_mod
fsundials_logger_mod
fsundials_matrix_mod
fsundials_nonlinearsolver_mod
fsundials_nvector_mod
fsundials_profiler_mod
fsundials_types_mod
```

### Deprecation notice

The functions in `sundials_math.h` will be deprecated in the next release.

```c
  sunrealtype SUNRpowerI(sunrealtype base, int exponent);
  sunrealtype SUNRpowerR(sunrealtype base, sunrealtype exponent);
  sunbooleantype SUNRCompare(sunrealtype a, sunrealtype b);
  sunbooleantype SUNRCompareTol(sunrealtype a, sunrealtype b, sunrealtype tol);
  sunrealtype SUNStrToReal(const char* str);
```

Additionally, the following header files (and everything in them) will be
deprecated -- users who rely on these are recommended to transition to the
corresponding `SUNMatrix` and `SUNLinearSolver` modules:

```c
sundials_direct.h
sundials_dense.h
sundials_band.h
```

### Minor Changes

Fixed [#329](https://github.com/LLNL/sundials/issues/329) so that C++20 aggregate initialization can be used.

Fixed integer overflow in the internal SUNDIALS hashmap. This resolves
[#409](https://github.com/LLNL/sundials/issues/409) and
[#249](https://github.com/LLNL/sundials/issues/249).

The `CMAKE_BUILD_TYPE` defaults to `RelWithDebInfo` mode now i.e., SUNDIALS
will be built with optimizations and debugging symbols enabled by default.
Previously the build type was unset by default so no optimization or debugging
flags were set.

The advanced CMake options to override the inferred LAPACK name-mangling scheme
have been updated from `SUNDIALS_F77_FUNC_CASE` and
`SUNDIALS_F77_FUNC_UNDERSCORES` to `SUNDIALS_LAPACK_CASE` and
`SUNDIALS_LAPACK_UNDERSCORES`, respectively.

Converted most previous Fortran 77 and 90 examples to use SUNDIALS' Fortran 2003
interface.

## Changes to SUNDIALS in release 6.7.0

Added the `SUNAdaptController` base class, ported ARKODE's internal
implementations of time step controllers into implementations of this class,
and updated ARKODE to use these objects instead of its own implementations.
Added `ARKStepSetAdaptController` and `ERKStepSetAdaptController` routines
so that users can modify controller parameters, or even provide custom
implementations.

Added the routines `ARKStepSetAdaptivityAdjustment` and
`ERKStepSetAdaptivityAdjustment`, that allow users to adjust the
value for the method order supplied to the temporal adaptivity controllers.
The ARKODE default for this adjustment has been -1 since its initial
release, but for some applications a value of 0 is more appropriate.
Users who notice that their simulations encounter a large number of
temporal error test failures may want to experiment with adjusting this value.

Added the third order ERK method `ARKODE_SHU_OSHER_3_2_3`, the fourth order
ERK method `ARKODE_SOFRONIOU_SPALETTA_5_3_4`, the sixth order ERK method
`ARKODE_VERNER_9_5_6`, the seventh order ERK method `ARKODE_VERNER_10_6_7`,
the eighth order ERK method `ARKODE_VERNER_13_7_8`, and the ninth order ERK
method `ARKODE_VERNER_16_8_9`.

ARKStep, ERKStep, MRIStep, and SPRKStep were updated to remove a potentially
unnecessary right-hand side evaluation at the end of an integration. ARKStep was
additionally updated to remove extra right-hand side evaluations when using an
explicit method or an implicit method with an explicit first stage.

Improved computational complexity of `SUNMatScaleAddI_Sparse` from `O(M*N)` to
`O(NNZ)`.

Added Fortran support for the LAPACK dense `SUNLinearSolver` implementation.

Fixed a regression introduced by the stop time bug fix in v6.6.1 where ARKODE,
CVODE, CVODES, IDA, and IDAS would return at the stop time rather than the
requested output time if the stop time was reached in the same step in which the
output time was passed.

Fixed a bug in ERKStep where methods with `c[s-1] = 1` but `a[s-1,j] != b[j]`
were incorrectly treated as having the first same as last (FSAL) property.

Fixed a bug in ARKODE where `ARKStepSetInterpolateStopTime` would return an
interpolated solution at the stop time in some cases when interpolation was
disabled.

Fixed a bug in `ARKStepSetTableNum` wherein it did not recognize
`ARKODE_ARK2_ERK_3_1_2` and `ARKODE_ARK2_DIRK_3_1_2` as a valid additive
Runge--Kutta Butcher table pair.

Fixed a bug in `MRIStepCoupling_Write` where explicit coupling tables were not
written to the output file pointer.

The `MRIStepInnerStepper` class in MRIStep was updated to make supplying an
`MRIStepInnerFullRhsFn` optional.

Fixed scaling bug in `SUNMatScaleAddI_Sparse` for non-square matrices.

Changed the `SUNProfiler` so that it does not rely on `MPI_WTime` in any case.
This fixes [GitHub Issue #312](https://github.com/LLNL/sundials/issues/312).

Fixed missing soversions in some `SUNLinearSolver` and `SUNNonlinearSolver`
CMake targets.

Renamed some internal types in CVODES and IDAS to allow both packages to be
built together in the same binary.

### ARKODE Changes in v5.7.0

Added the :c:type:`SUNAdaptController` base class, ported ARKODE's internal
implementations of time step controllers into implementations of this class,
and updated ARKODE to use these objects instead of its own implementations.  Added
:c:func:`ARKStepSetAdaptController` and :c:func:`ERKStepSetAdaptController`
routines so that users can modify controller parameters, or even provide custom
implementations.

Added the routines :c:func:`ARKStepSetAdaptivityAdjustment` and
:c:func:`ERKStepSetAdaptivityAdjustment`, that allow users to adjust the
value for the method order supplied to the temporal adaptivity controllers.
The ARKODE default for this adjustment has been :math:`-1` since its initial
release, but for some applications a value of :math:`0` is more appropriate.
Users who notice that their simulations encounter a large number of
temporal error test failures may want to experiment with adjusting this value.

Added the third order ERK method ``ARKODE_SHU_OSHER_3_2_3``, the fourth order
ERK method ``ARKODE_SOFRONIOU_SPALETTA_5_3_4``, the sixth order ERK method
``ARKODE_VERNER_9_5_6``, the seventh order ERK method ``ARKODE_VERNER_10_6_7``,
the eighth order ERK method ``ARKODE_VERNER_13_7_8``, and the ninth order ERK
method ``ARKODE_VERNER_16_8_9``.

ARKStep, ERKStep, MRIStep, and SPRKStep were updated to remove a potentially
unnecessary right-hand side evaluation at the end of an integration. ARKStep was
additionally updated to remove extra right-hand side evaluations when using an
explicit method or an implicit method with an explicit first stage.

Improved computational complexity of ``SUNMatScaleAddI_Sparse`` from ``O(M*N)``
to ``O(NNZ)``.

Added Fortran support for the LAPACK  dense ``SUNLinearSolver`` implementation.

Fixed a regression introduced by the stop time bug fix in v6.6.1 where ARKODE
steppers would return at the stop time rather than the requested output time if
the stop time was reached in the same step in which the output time was passed.

Fixed a bug in ERKStep where methods with :math:`c_s = 1` but
:math:`a_{s,j} \neq b_j` were incorrectly treated as having the first same as
last (FSAL) property.

Fixed a bug in ARKODE where :c:func:`ARKStepSetInterpolateStopTime` would return
an interpolated solution at the stop time in some cases when interpolation was
disabled.

Fixed a bug in :c:func:`ARKStepSetTableNum` wherein it did not recognize
`ARKODE_ARK2_ERK_3_1_2` and `ARKODE_ARK2_DIRK_3_1_2` as a valid additive
Runge--Kutta Butcher table pair.

Fixed a bug in :c:func:`MRIStepCoupling_Write` where explicit coupling tables
were not written to the output file pointer.

The :c:type:`MRIStepInnerStepper` class in MRIStep was updated to make supplying
an :c:func:`MRIStepInnerFullRhsFn` optional.

Fixed scaling bug in ``SUNMatScaleAddI_Sparse`` for non-square matrices.

Changed the ``SUNProfiler`` so that it does not rely on ``MPI_WTime`` in any case.
This fixes `GitHub Issue #312 <https://github.com/LLNL/sundials/issues/312>`_.

Fixed missing soversions in some ``SUNLinearSolver`` and ``SUNNonlinearSolver``
CMake targets.

### CVODE Changes in v6.7.0

Improved computational complexity of ``SUNMatScaleAddI_Sparse`` from ``O(M*N)``
to ``O(NNZ)``.

Added Fortran support for the LAPACK  dense ``SUNLinearSolver`` implementation.

Fixed a regression introduced by the stop time bug fix in v6.6.1 where CVODE
would return at the stop time rather than the requested output time if the stop
time was reached in the same step in which the output time was passed.

Fixed scaling bug in ``SUNMatScaleAddI_Sparse`` for non-square matrices.

Changed the ``SUNProfiler`` so that it does not rely on ``MPI_WTime`` in any case.
This fixes `GitHub Issue #312 <https://github.com/LLNL/sundials/issues/312>`_.

Fixed missing soversions in some ``SUNLinearSolver`` and ``SUNNonlinearSolver``
CMake targets.

### CVODES Changes in v6.7.0

Improved computational complexity of ``SUNMatScaleAddI_Sparse`` from ``O(M*N)``
to ``O(NNZ)``.

Added Fortran support for the LAPACK  dense ``SUNLinearSolver`` implementation.

Fixed a regression introduced by the stop time bug fix in v6.6.1 where CVODE
would return at the stop time rather than the requested output time if the stop
time was reached in the same step in which the output time was passed.

Fixed scaling bug in ``SUNMatScaleAddI_Sparse`` for non-square matrices.

Changed the ``SUNProfiler`` so that it does not rely on ``MPI_WTime`` in any case.
This fixes `GitHub Issue #312 <https://github.com/LLNL/sundials/issues/312>`_.

Fixed missing soversions in some ``SUNLinearSolver`` and ``SUNNonlinearSolver``
CMake targets.

Renamed some internal types in CVODES and IDAS to allow both packages to be
built together in the same binary.

### IDA Changes in v6.7.0

Improved computational complexity of ``SUNMatScaleAddI_Sparse`` from ``O(M*N)``
to ``O(NNZ)``.

Added Fortran support for the LAPACK  dense ``SUNLinearSolver`` implementation.

Fixed a regression introduced by the stop time bug fix in v6.6.1 where CVODE
would return at the stop time rather than the requested output time if the stop
time was reached in the same step in which the output time was passed.

Fixed scaling bug in ``SUNMatScaleAddI_Sparse`` for non-square matrices.

Changed the ``SUNProfiler`` so that it does not rely on ``MPI_WTime`` in any case.
This fixes `GitHub Issue #312 <https://github.com/LLNL/sundials/issues/312>`_.

Fixed missing soversions in some ``SUNLinearSolver`` and ``SUNNonlinearSolver``
CMake targets.

### IDAS Changes in v5.7.0

Improved computational complexity of ``SUNMatScaleAddI_Sparse`` from ``O(M*N)``
to ``O(NNZ)``.

Added Fortran support for the LAPACK  dense ``SUNLinearSolver`` implementation.

Fixed a regression introduced by the stop time bug fix in v6.6.1 where CVODE
would return at the stop time rather than the requested output time if the stop
time was reached in the same step in which the output time was passed.

Fixed scaling bug in ``SUNMatScaleAddI_Sparse`` for non-square matrices.

Changed the ``SUNProfiler`` so that it does not rely on ``MPI_WTime`` in any case.
This fixes `GitHub Issue #312 <https://github.com/LLNL/sundials/issues/312>`_.

Fixed missing soversions in some ``SUNLinearSolver`` and ``SUNNonlinearSolver``
CMake targets.

Renamed some internal types in CVODES and IDAS to allow both packages to be
built together in the same binary.

### KINSOL Changes in v6.7.0

Added Fortran support for the LAPACK  dense ``SUNLinearSolver`` implementation.

Improved computational complexity of ``SUNMatScaleAddI_Sparse`` from ``O(M*N)`` to
``O(NNZ)``.

Fixed scaling bug in ``SUNMatScaleAddI_Sparse`` for non-square matrices.

Fixed missing soversions in some ``SUNLinearSolver`` and ``SUNNonlinearSolver``
CMake targets.



## Changes to SUNDIALS in release 6.6.2

Fixed the build system support for MAGMA when using a NVIDIA HPC SDK installation of CUDA
and fixed the targets used for rocBLAS and rocSPARSE.

### ARKODE Changes in v5.6.2

Fixed the build system support for MAGMA when using a NVIDIA HPC SDK installation of CUDA
and fixed the targets used for rocBLAS and rocSPARSE.

### CVODE Changes in v6.6.2

Fixed the build system support for MAGMA when using a NVIDIA HPC SDK installation of CUDA
and fixed the targets used for rocBLAS and rocSPARSE.

### CVODES Changes in v6.6.2

Fixed the build system support for MAGMA when using a NVIDIA HPC SDK installation of CUDA
and fixed the targets used for rocBLAS and rocSPARSE.

### IDA Changes in v6.6.2

Fixed the build system support for MAGMA when using a NVIDIA HPC SDK installation of CUDA
and fixed the targets used for rocBLAS and rocSPARSE.

### IDAS Changes in v5.6.2

Fixed the build system support for MAGMA when using a NVIDIA HPC SDK installation of CUDA
and fixed the targets used for rocBLAS and rocSPARSE.

### KINSOL Changes in v6.6.2

Fixed the build system support for MAGMA when using a NVIDIA HPC SDK installation of CUDA
and fixed the targets used for rocBLAS and rocSPARSE.



## Changes to SUNDIALS in release 6.6.1

Updated the Tpetra NVector interface to support Trilinos 14.

Fixed a memory leak when destroying a CUDA, HIP, SYCL, or system SUNMemoryHelper
object.

Fixed a bug in ARKODE, CVODE, CVODES, IDA, and IDAS where the stop time may not
be cleared when using normal mode if the requested output time is the same as
the stop time. Additionally, with ARKODE, CVODE, and CVODES this fix removes an
unnecessary interpolation of the solution at the stop time that could occur in
this case.

### ARKODE Changes in v5.6.1

Updated the Tpetra NVector interface to support Trilinos 14.

Fixed a memory leak when destroying a CUDA, HIP, SYCL, or system SUNMemoryHelper
object.

Fixed a bug where the stop time may not be cleared when using normal mode if the
requested output time is the same as the stop time. Additionally, this fix
removes an unnecessary interpolation of the solution at the stop time that could
occur in this case.

### CVODE Changes in v6.6.1

Updated the Tpetra NVector interface to support Trilinos 14.

Fixed a memory leak when destroying a CUDA, HIP, SYCL, or system SUNMemoryHelper
object.

Fixed a bug where the stop time may not be cleared when using normal mode if the
requested output time is the same as the stop time. Additionally, this fix
removes an unnecessary interpolation of the solution at the stop time that could
occur in this case.

### CVODES Changes in v6.6.1

Updated the Tpetra NVector interface to support Trilinos 14.

Fixed a memory leak when destroying a CUDA, HIP, SYCL, or system SUNMemoryHelper
object.

Fixed a bug where the stop time may not be cleared when using normal mode if the
requested output time is the same as the stop time. Additionally, this fix
removes an unnecessary interpolation of the solution at the stop time that could
occur in this case.

### IDA Changes in v6.6.1

Updated the Tpetra NVector interface to support Trilinos 14.

Fixed a memory leak when destroying a CUDA, HIP, SYCL, or system SUNMemoryHelper
object.

Fixed a bug where the stop time may not be cleared when using normal mode if the
requested output time is the same as the stop time.

### IDAS Changes in v5.6.1

Updated the Tpetra NVector interface to support Trilinos 14.

Fixed a memory leak when destroying a CUDA, HIP, SYCL, or system SUNMemoryHelper
object.

Fixed a bug where the stop time may not be cleared when using normal mode if the
requested output time is the same as the stop time.

### KINSOL Changes in v6.6.1

Updated the Tpetra NVector interface to support Trilinos 14.

Fixed a memory leak when destroying a CUDA, HIP, SYCL, or system SUNMemoryHelper
object.

Changed the ``SUNProfiler`` so that it does not rely on ``MPI_WTime`` in any case.
This fixes `GitHub Issue #312 <https://github.com/LLNL/sundials/issues/312>`_.



## Changes to SUNDIALS in release 6.6.0

A new time-stepping module, `SPRKStep`, was added to ARKODE. This time-stepper
provides explicit symplectic partitioned Runge-Kutta methods up to order 10
for separable Hamiltonian systems.

Added support for relaxation Runge-Kutta methods to ERKStep and ARKStep in
ARKODE.

Added the second order IMEX method from Giraldo, Kelly, and Constantinescu 2013
as the default second order IMEX method in ARKStep. The explicit table is given
by `ARKODE_ARK2_ERK_3_1_2` and the implicit table by `ARKODE_ARK2_DIRK_3_1_2`.

Updated CVODE, CVODES and ARKODE default behavior when returning the solution when
the internal time has reached a user-specified stop time.  Previously, the output
solution was interpolated to the value of `tstop`; the default is now to copy the
internal solution vector.  Users who wish to revert to interpolation may call a new
routine `CVodeSetInterpolateStopTime`, `ARKStepSetInterpolateStopTime`,
`ERKStepSetInterpolateStopTime`, or `MRIStepSetInterpolateStopTime`.

A potential bug was fixed when using inequality constraint handling and
calling `ARKStepGetEstLocalErrors` or `ERKStepGetEstLocalErrors` after a failed
step in which an inequality constraint violation occurred. In this case, the
values returned by `ARKStepGetEstLocalErrors` or `ERKStepGetEstLocalErrors` may
have been invalid.

Updated the F2003 utility routines `SUNDIALSFileOpen` and `SUNDIALSFileClose`
to support user specification of `stdout` and `stderr` strings for the output
file names.

### ARKODE Changes in v5.6.0

A new time-stepping module, :ref:`SPRKStep <ARKODE.Mathematics.SPRKStep>`, was
added to ARKODE. This time-stepper provides explicit symplectic partitioned
Runge-Kutta methods up to order 10 for separable Hamiltonian systems.

Added support for relaxation Runge-Kutta methods in ERKStep and ARKStep, see
:numref:`ARKODE.Mathematics.Relaxation`, :numref:`ARKODE.Usage.ERKStep.Relaxation`,
and :numref:`ARKODE.Usage.ARKStep.Relaxation` for more information.

Added the second order IMEX method from :cite:p:`giraldo2013implicit` as the
default second order IMEX method in ARKStep. The explicit table is given by
``ARKODE_ARK2_ERK_3_1_2`` (see :numref:`Butcher.ARK2_ERK`) and the implicit
table by ``ARKODE_ARK2_DIRK_3_1_2`` (see :numref:`Butcher.ARK2_DIRK`).

Updated the default ARKODE behavior when returning the solution when
the internal time has reached a user-specified stop time.  Previously, the output
solution was interpolated to the value of ``tstop``; the default is now to copy the
internal solution vector.  Users who wish to revert to interpolation may call a new
routine :c:func:`ARKStepSetInterpolateStopTime`,
:c:func:`ERKStepSetInterpolateStopTime`, or :c:func:`MRIStepSetInterpolateStopTime`.

A potential bug was fixed when using inequality constraint handling and
calling :c:func:`ARKStepGetEstLocalErrors` or :c:func:`ERKStepGetEstLocalErrors`
after a failed step in which an inequality constraint violation occurred. In
this case, the values returned by :c:func:`ARKStepGetEstLocalErrors` or
:c:func:`ERKStepGetEstLocalErrors` may have been invalid.

Updated the F2003 utility routines :c:func:`SUNDIALSFileOpen` and :c:func:`SUNDIALSFileClose`
to support user specification of ``stdout`` and ``stderr`` strings for the output
file names.

### CVODE Changes in v6.6.0

Updated the default CVODE behavior when returning the solution when
the internal time has reached a user-specified stop time.  Previously, the output
solution was interpolated to the value of ``tstop``; the default is now to copy the
internal solution vector.  Users who wish to revert to interpolation may call the
routine :c:func:`CVodeSetInterpolateStopTime`.

Updated the F2003 utility routines :c:func:`SUNDIALSFileOpen` and :c:func:`SUNDIALSFileClose`
to support user specification of ``stdout`` and ``stderr`` strings for the output
file names.

### CVODES Changes in v6.6.0

Updated the default CVODES behavior when returning the solution when
the internal time has reached a user-specified stop time.  Previously, the output
solution was interpolated to the value of ``tstop``; the default is now to copy the
internal solution vector.  Users who wish to revert to interpolation may call the
routine :c:func:`CVodeSetInterpolateStopTime`.

Updated the F2003 utility routines :c:func:`SUNDIALSFileOpen` and :c:func:`SUNDIALSFileClose`
to support user specification of ``stdout`` and ``stderr`` strings for the output
file names.

### IDA Changes in v6.6.0

Updated the F2003 utility routines :c:func:`SUNDIALSFileOpen` and :c:func:`SUNDIALSFileClose`
to support user specification of ``stdout`` and ``stderr`` strings for the output
file names.

### IDAS Changes in v5.6.0

Updated the F2003 utility routines :c:func:`SUNDIALSFileOpen` and :c:func:`SUNDIALSFileClose`
to support user specification of ``stdout`` and ``stderr`` strings for the output
file names.

### KINSOL Changes in v6.6.0

Updated the F2003 utility routines :c:func:`SUNDIALSFileOpen` and :c:func:`SUNDIALSFileClose`
to support user specification of ``stdout`` and ``stderr`` strings for the output
file names.



## Changes to SUNDIALS in release 6.5.1

Added the functions `ARKStepClearStopTime`, `ERKStepClearStopTime`,
`MRIStepClearStopTime`, `CVodeClearStopTime`, and `IDAClearStopTime` to
disable a previously set stop time.

Fixed build errors when using SuperLU_DIST with ROCM enabled to target AMD GPUs.

Fixed compilation errors in some SYCL examples when using the `icx` compiler.

The default interpolant in ARKODE when using a first order method has been
updated to a linear interpolant to ensure values obtained by the integrator are
returned at the ends of the time interval. To restore the previous behavior of
using a constant interpolant call `ARKStepSetInterpolantDegree`,
`ERKStepSetInterpolantDegree`, or `MRIStepSetInterpolantDegree` and set the
interpolant degree to zero before evolving the problem.

### ARKODE Changes in v5.5.1

Added the functions :c:func:`ARKStepClearStopTime`,
:c:func:`ERKStepClearStopTime`, and :c:func:`MRIStepClearStopTime` to disable a
previously set stop time.

Fixed build errors when using SuperLU_DIST with ROCM enabled to target AMD GPUs.

Fixed compilation errors in some SYCL examples when using the ``icx`` compiler.

The default interpolant in ARKODE when using a first order method has been
updated to a linear interpolant to ensure values obtained by the integrator are
returned at the ends of the time interval. To restore the previous behavior of
using a constant interpolant call :c:func:`ARKStepSetInterpolantDegree`,
:c:func:`ERKStepSetInterpolantDegree`, or :c:func:`MRIStepSetInterpolantDegree`
and set the interpolant degree to zero before evolving the problem.

### CVODE Changes in v6.5.1

Added the function :c:func:`CVodeClearStopTime` to disable a previously set stop
time.

Fixed build errors when using SuperLU_DIST with ROCM enabled to target AMD GPUs.

Fixed compilation errors in some SYCL examples when using the ``icx`` compiler.

### CVODES Changes in v6.5.1

Added the function :c:func:`CVodeClearStopTime` to disable a previously set stop
time.

Fixed build errors when using SuperLU_DIST with ROCM enabled to target AMD GPUs.

Fixed compilation errors in some SYCL examples when using the ``icx`` compiler.

### IDA Changes in v6.5.1

Added the function :c:func:`IDAClearStopTime` to disable a previously set stop
time.

Fixed build errors when using SuperLU_DIST with ROCM enabled to target AMD GPUs.

Fixed compilation errors in some SYCL examples when using the ``icx`` compiler.

### IDAS Changes in v5.5.1

Added the function :c:func:`IDAClearStopTime` to disable a previously set stop
time.

Fixed build errors when using SuperLU_DIST with ROCM enabled to target AMD GPUs.

Fixed compilation errors in some SYCL examples when using the ``icx`` compiler.

### KINSOL Changes in v6.5.1

Fixed build errors when using SuperLU_DIST with ROCM enabled to target AMD GPUs.

Fixed compilation errors in some SYCL examples when using the ``icx`` compiler.





## Changes to SUNDIALS in release 6.5.0

Added the functions `ARKStepGetJac`, `ARKStepGetJacTime`,
`ARKStepGetJacNumSteps`, `MRIStepGetJac`, `MRIStepGetJacTime`,
`MRIStepGetJacNumSteps`, `CVodeGetJac`, `CVodeGetJacTime`,
`CVodeGetJacNumSteps`, `IDAGetJac`, `IDAGetJacCj`, `IDAGetJacTime`,
`IDAGetJacNumSteps`, `KINGetJac`, `KINGetJacNumIters` to assist in
debugging simulations utilizing a matrix-based linear solver.

Added support for the SYCL backend with RAJA 2022.x.y.

Fixed an underflow bug during root finding in ARKODE, CVODE, CVODES, IDA and
IDAS.

Fixed an issue with finding oneMKL when using the `icpx` compiler with the
`-fsycl` flag as the C++ compiler instead of `dpcpp`.

Fixed the shape of the arrays returned by `FN_VGetArrayPointer` functions as well
as the `FSUNDenseMatrix_Data`, `FSUNBandMatrix_Data`, `FSUNSparseMatrix_Data`,
`FSUNSparseMatrix_IndexValues`, and `FSUNSparseMatrix_IndexPointers` functions.
Compiling and running code that uses the SUNDIALS Fortran interfaces with
bounds checking will now work.

Fixed an implicit conversion error in the Butcher table for ESDIRK5(4)7L[2]SA2.

A new capability to keep track of memory allocations made through the `SUNMemoryHelper`
classes has been added. Memory allocation stats can be accessed through the
`SUNMemoryHelper_GetAllocStats` function. See the documentation for
the `SUNMemoryHelper` classes for more details.

Added support for CUDA 12.

### ARKODE Changes in v5.5.0

Added the functions :c:func:`ARKStepGetJac`, :c:func:`ARKStepGetJacTime`,
:c:func:`ARKStepGetJacNumSteps`, :c:func:`MRIStepGetJac`,
:c:func:`MRIStepGetJacTime`, and :c:func:`MRIStepGetJacNumSteps` to assist in
debugging simulations utilizing a matrix-based linear solver.

Added support for the SYCL backend with RAJA 2022.x.y.

Fixed an underflow bug during root finding.

A new capability to keep track of memory allocations made through the ``SUNMemoryHelper``
classes has been added. Memory allocation stats can be accessed through the
:c:func:`SUNMemoryHelper_GetAllocStats` function. See the documentation for
the ``SUNMemoryHelper`` classes for more details.

Added support for CUDA v12.

Fixed an issue with finding oneMKL when using the ``icpx`` compiler with the
``-fsycl`` flag as the C++ compiler instead of ``dpcpp``.

Fixed the shape of the arrays returned by ``FN_VGetArrayPointer`` functions as well
as the ``FSUNDenseMatrix_Data``, ``FSUNBandMatrix_Data``, ``FSUNSparseMatrix_Data``,
``FSUNSparseMatrix_IndexValues``, and ``FSUNSparseMatrix_IndexPointers`` functions.
Compiling and running code that uses the SUNDIALS Fortran interfaces with
bounds checking will now work.

Fixed an implicit conversion error in the Butcher table for ESDIRK5(4)7L[2]SA2.

### CVODE Changes in v6.5.0

Added the functions :c:func:`CVodeGetJac`, :c:func:`CVodeGetJacTime`,
:c:func:`CVodeGetJacNumSteps` to assist in debugging simulations utilizing
a matrix-based linear solver.

Added support for the SYCL backend with RAJA 2022.x.y.

Fixed an underflow bug during root finding.

A new capability to keep track of memory allocations made through the ``SUNMemoryHelper``
classes has been added. Memory allocation stats can be accessed through the
:c:func:`SUNMemoryHelper_GetAllocStats` function. See the documentation for
the ``SUNMemoryHelper`` classes for more details.

Added support for CUDA v12.
Fixed an issue with finding oneMKL when using the ``icpx`` compiler with the
``-fsycl`` flag as the C++ compiler instead of ``dpcpp``.

Fixed the shape of the arrays returned by ``FN_VGetArrayPointer`` functions as well
as the ``FSUNDenseMatrix_Data``, ``FSUNBandMatrix_Data``, ``FSUNSparseMatrix_Data``,
``FSUNSparseMatrix_IndexValues``, and ``FSUNSparseMatrix_IndexPointers`` functions.
Compiling and running code that uses the SUNDIALS Fortran interfaces with
bounds checking will now work.

### CVODES Changes in v6.5.0

Added the functions :c:func:`CVodeGetJac`, :c:func:`CVodeGetJacTime`,
:c:func:`CVodeGetJacNumSteps` to assist in debugging simulations utilizing
a matrix-based linear solver.

Added support for the SYCL backend with RAJA 2022.x.y.

Fixed an underflow bug during root finding.

A new capability to keep track of memory allocations made through the ``SUNMemoryHelper``
classes has been added. Memory allocation stats can be accessed through the
:c:func:`SUNMemoryHelper_GetAllocStats` function. See the documentation for
the ``SUNMemoryHelper`` classes for more details.

Added support for CUDA v12.

Fixed an issue with finding oneMKL when using the ``icpx`` compiler with the
``-fsycl`` flag as the C++ compiler instead of ``dpcpp``.

Fixed the shape of the arrays returned by ``FN_VGetArrayPointer`` functions as well
as the ``FSUNDenseMatrix_Data``, ``FSUNBandMatrix_Data``, ``FSUNSparseMatrix_Data``,
``FSUNSparseMatrix_IndexValues``, and ``FSUNSparseMatrix_IndexPointers`` functions.
Compiling and running code that uses the SUNDIALS Fortran interfaces with
bounds checking will now work.

### IDA Changes in v6.5.0

Added the functions :c:func:`IDAGetJac`, :c:func:`IDAGetJacCj`,
:c:func:`IDAGetJacTime`, :c:func:`IDAGetJacNumSteps` to assist in debugging
simulations utilizing a matrix-based linear solver.

Added support for the SYCL backend with RAJA 2022.x.y.

Fixed an underflow bug during root finding.

A new capability to keep track of memory allocations made through the ``SUNMemoryHelper``
classes has been added. Memory allocation stats can be accessed through the
:c:func:`SUNMemoryHelper_GetAllocStats` function. See the documentation for
the ``SUNMemoryHelper`` classes for more details.

Added support for CUDA v12.

Fixed an issue with finding oneMKL when using the ``icpx`` compiler with the
``-fsycl`` flag as the C++ compiler instead of ``dpcpp``.

Fixed the shape of the arrays returned by ``FN_VGetArrayPointer`` functions as well
as the ``FSUNDenseMatrix_Data``, ``FSUNBandMatrix_Data``, ``FSUNSparseMatrix_Data``,
``FSUNSparseMatrix_IndexValues``, and ``FSUNSparseMatrix_IndexPointers`` functions.
Compiling and running code that uses the SUNDIALS Fortran interfaces with
bounds checking will now work.

### IDAS Changes in v5.5.0

Added the functions :c:func:`IDAGetJac`, :c:func:`IDAGetJacCj`,
:c:func:`IDAGetJacTime`, :c:func:`IDAGetJacNumSteps` to assist in debugging
simulations utilizing a matrix-based linear solver.

Added support for the SYCL backend with RAJA 2022.x.y.

Fixed an underflow bug during root finding.

A new capability to keep track of memory allocations made through the ``SUNMemoryHelper``
classes has been added. Memory allocation stats can be accessed through the
:c:func:`SUNMemoryHelper_GetAllocStats` function. See the documentation for
the ``SUNMemoryHelper`` classes for more details.

Added support for CUDA v12.

Fixed an issue with finding oneMKL when using the ``icpx`` compiler with the
``-fsycl`` flag as the C++ compiler instead of ``dpcpp``.

Fixed the shape of the arrays returned by ``FN_VGetArrayPointer`` functions as well
as the ``FSUNDenseMatrix_Data``, ``FSUNBandMatrix_Data``, ``FSUNSparseMatrix_Data``,
``FSUNSparseMatrix_IndexValues``, and ``FSUNSparseMatrix_IndexPointers`` functions.
Compiling and running code that uses the SUNDIALS Fortran interfaces with
bounds checking will now work.

### KINSOL Changes in v6.5.0

A new capability to keep track of memory allocations made through the ``SUNMemoryHelper``
classes has been added. Memory allocation stats can be accessed through the
:c:func:`SUNMemoryHelper_GetAllocStats` function. See the documentation for
the ``SUNMemoryHelper`` classes for more details.

Added the functions :c:func:`KINGetJac` and :c:func:`KINGetJacNumIters` to
assist in debugging simulations utilizing a matrix-based linear solver.

Added support for the SYCL backend with RAJA 2022.x.y.

Fixed an issue with finding oneMKL when using the ``icpx`` compiler with the
``-fsycl`` flag as the C++ compiler instead of ``dpcpp``.

Fixed the shape of the arrays returned by ``FN_VGetArrayPointer`` functions as well
as the ``FSUNDenseMatrix_Data``, ``FSUNBandMatrix_Data``, ``FSUNSparseMatrix_Data``,
``FSUNSparseMatrix_IndexValues``, and ``FSUNSparseMatrix_IndexPointers`` functions.
Compiling and running code that uses the SUNDIALS Fortran interfaces with
bounds checking will now work.



## Changes to SUNDIALS in release 6.4.1

Fixed a bug with the Kokkos interfaces that would arise when using clang.

Fixed a compilation error with the Intel oneAPI 2022.2 Fortran compiler in the
Fortran 2003 interface test for the serial `N_Vector`.

Fixed a bug in the SUNLINSOL_LAPACKBAND and SUNLINSOL_LAPACKDENSE modules
which would cause the tests to fail on some platforms.

### ARKODE Changes in v5.4.1

Fixed a bug with the Kokkos interfaces that would arise when using clang.

Fixed a compilation error with the Intel oneAPI 2022.2 Fortran compiler in the
Fortran 2003 interface test for the serial ``N_Vector``.

Fixed a bug in the SUNLINSOL_LAPACKBAND and SUNLINSOL_LAPACKDENSE modules
which would cause the tests to fail on some platforms.

### CVODE Changes in v6.4.1

Fixed a bug with the Kokkos interfaces that would arise when using clang.

Fixed a compilation error with the Intel oneAPI 2022.2 Fortran compiler in the
Fortran 2003 interface test for the serial ``N_Vector``.

Fixed a bug in the SUNLINSOL_LAPACKBAND and SUNLINSOL_LAPACKDENSE modules
which would cause the tests to fail on some platforms.

### CVODES Changes in v6.4.1

Fixed a bug with the Kokkos interfaces that would arise when using clang.

Fixed a compilation error with the Intel oneAPI 2022.2 Fortran compiler in the
Fortran 2003 interface test for the serial ``N_Vector``.

Fixed a bug in the SUNLINSOL_LAPACKBAND and SUNLINSOL_LAPACKDENSE modules
which would cause the tests to fail on some platforms.

### IDA  Changes in v6.4.1

Fixed a bug with the Kokkos interfaces that would arise when using clang.

Fixed a compilation error with the Intel oneAPI 2022.2 Fortran compiler in the
Fortran 2003 interface test for the serial ``N_Vector``.

Fixed a bug in the SUNLINSOL_LAPACKBAND and SUNLINSOL_LAPACKDENSE modules
which would cause the tests to fail on some platforms.

### IDAS Changes in v5.4.1

Fixed a bug with the Kokkos interfaces that would arise when using clang.

Fixed a compilation error with the Intel oneAPI 2022.2 Fortran compiler in the
Fortran 2003 interface test for the serial ``N_Vector``.

Fixed a bug in the SUNLINSOL_LAPACKBAND and SUNLINSOL_LAPACKDENSE modules
which would cause the tests to fail on some platforms.

### KINSOL Changes in v6.4.1

Fixed a bug with the Kokkos interfaces that would arise when using clang.

Fixed a compilation error with the Intel oneAPI 2022.2 Fortran compiler in the
Fortran 2003 interface test for the serial ``N_Vector``.

Fixed a bug in the SUNLINSOL_LAPACKBAND and SUNLINSOL_LAPACKDENSE modules
which would cause the tests to fail on some platforms.



## Changes to SUNDIALS in release 6.4.0

CMake 3.18.0 or newer is now required for CUDA support.

A C++14 compliant compiler is now required for C++ based features and examples
e.g., CUDA, HIP, RAJA, Trilinos, SuperLU_DIST, MAGMA, GINKGO, and KOKKOS.

Added support for GPU enabled SuperLU_DIST and SuperLU_DIST v8.x.x. Removed
support for SuperLU_DIST v6.x.x or older. Fix mismatched definition and
declaration bug in SuperLU_DIST matrix constructor.

Added support for the [Ginkgo](https://ginkgo-project.github.io/) linear algebra
library. This support includes new `SUNMatrix` and `SUNLinearSolver`
implementations, see the `SUNMATRIX_GINKGO` and `SUNLINEARSOLVER_GINKGO`
sections in the documentation for more information.

Added new `NVector`, dense `SUNMatrix`, and dense `SUNLinearSolver`
implementations utilizing [Kokkos Ecosystem](https://kokkos.org/) for
performance portability, see the `NVECTOR_KOKKOS`, `SUNMATRIX_KOKKOSDENSE` and
`SUNLINEARSOLVER_KOKKOSDENSE` sections in the documentation for more
information.

Added the functions `ARKStepSetTableName`, `ERKStepSetTableName`,
`MRIStepCoupling_LoadTableByName`, `ARKodeButcherTable_LoadDIRKByName`, and
`ARKodeButcherTable_LoadERKByName` to load a table from a string.

Fixed a bug in the CUDA and HIP vectors where `N_VMaxNorm` would return the
minimum positive floating-point value for the zero vector.

Fixed memory leaks/out of bounds memory accesses in the ARKODE MRIStep module
that could occur when attaching a coupling table after reinitialization with a
different number of stages than originally selected.

Fixed a memory leak in CVODE and CVODES where the projection memory would not be
deallocated when calling `CVodeFree`.

### ARKODE Changes in v5.4.0

CMake 3.18.0 or newer is now required for CUDA support.

A C++14 compliant compiler is now required for C++ based features and examples
e.g., CUDA, HIP, RAJA, Trilinos, SuperLU_DIST, MAGMA, GINKGO, and KOKKOS.

Added support for GPU enabled SuperLU_DIST and SuperLU_DIST v8.x.x. Removed
support for SuperLU_DIST v6.x.x or older. Fix mismatched definition and
declaration bug in SuperLU_DIST matrix constructor.

Added support for the `Ginkgo <https://ginkgo-project.github.io/>`_  linear
algebra library. This support includes new ``SUNMatrix`` and ``SUNLinearSolver``
implementations, see the sections :numref:`SUNMatrix.Ginkgo` and
:numref:`SUNLinSol.Ginkgo`.

Added new ``NVector``, dense ``SUNMatrix``, and dense ``SUNLinearSolver``
implementations utilizing the `Kokkos Ecosystem <https://kokkos.org/>`_ for
performance portability, see sections :numref:`NVectors.Kokkos`,
:numref:`SUNMatrix.Kokkos`, and :numref:`SUNLinSol.Kokkos` for more information.

Added the functions :c:func:`ARKStepSetTableName`,
:c:func:`ERKStepSetTableName`, :c:func:`MRIStepCoupling_LoadTableByName`,
:c:func:`ARKodeButcherTable_LoadDIRKByName`, and
:c:func:`ARKodeButcherTable_LoadERKByName` to load a table from a string.

Fixed a bug in the CUDA and HIP vectors where :c:func:`N_VMaxNorm` would return
the minimum positive floating-point value for the zero vector.

Fixed memory leaks/out of bounds memory accesses in the ARKODE MRIStep module
that could occur when attaching a coupling table after reinitialization with a
different number of stages than originally selected.


### CVODE Changes in v6.4.0

CMake 3.18.0 or newer is now required for CUDA support.

A C++14 compliant compiler is now required for C++ based features and examples
e.g., CUDA, HIP, RAJA, Trilinos, SuperLU_DIST, MAGMA, GINKGO, and KOKKOS.

Added support for GPU enabled SuperLU_DIST and SuperLU_DIST v8.x.x. Removed
support for SuperLU_DIST v6.x.x or older. Fix mismatched definition and
declaration bug in SuperLU_DIST matrix constructor.

Added support for the `Ginkgo <https://ginkgo-project.github.io/>`_  linear
algebra library. This support includes new ``SUNMatrix`` and ``SUNLinearSolver``
implementations, see the sections :numref:`SUNMatrix.Ginkgo` and
:numref:`SUNLinSol.Ginkgo`.

Added new ``NVector``, dense ``SUNMatrix``, and dense ``SUNLinearSolver``
implementations utilizing the `Kokkos Ecosystem <https://kokkos.org/>`_ for
performance portability, see sections :numref:`NVectors.Kokkos`,
:numref:`SUNMatrix.Kokkos`, and :numref:`SUNLinSol.Kokkos` for more information.

Fixed a bug in the CUDA and HIP vectors where :c:func:`N_VMaxNorm` would return
the minimum positive floating-point value for the zero vector.

Fixed a memory leak where the projection memory would not be deallocated when
calling :c:func:`CVodeFree`.

### CVODES Changes in v6.4.0

CMake 3.18.0 or newer is now required for CUDA support.

A C++14 compliant compiler is now required for C++ based features and examples
e.g., CUDA, HIP, RAJA, Trilinos, SuperLU_DIST, MAGMA, GINKGO, and KOKKOS.

Added support for GPU enabled SuperLU_DIST and SuperLU_DIST v8.x.x. Removed
support for SuperLU_DIST v6.x.x or older. Fix mismatched definition and
declaration bug in SuperLU_DIST matrix constructor.

Added support for the `Ginkgo <https://ginkgo-project.github.io/>`_  linear
algebra library. This support includes new ``SUNMatrix`` and ``SUNLinearSolver``
implementations, see the sections :numref:`SUNMatrix.Ginkgo` and
:numref:`SUNLinSol.Ginkgo`.

Added new ``NVector``, dense ``SUNMatrix``, and dense ``SUNLinearSolver``
implementations utilizing the `Kokkos Ecosystem <https://kokkos.org/>`_ for
performance portability, see sections :numref:`NVectors.Kokkos`,
:numref:`SUNMatrix.Kokkos`, and :numref:`SUNLinSol.Kokkos` for more information.

Fixed a bug in the CUDA and HIP vectors where :c:func:`N_VMaxNorm` would return
the minimum positive floating-point value for the zero vector.

Fixed a memory leak where the projection memory would not be deallocated when
calling :c:func:`CVodeFree`.

### IDA Changes in v6.4.0

CMake 3.18.0 or newer is now required for CUDA support.

A C++14 compliant compiler is now required for C++ based features and examples
e.g., CUDA, HIP, RAJA, Trilinos, SuperLU_DIST, MAGMA, GINKGO, and KOKKOS.

Added support for GPU enabled SuperLU_DIST and SuperLU_DIST v8.x.x. Removed
support for SuperLU_DIST v6.x.x or older. Fix mismatched definition and
declaration bug in SuperLU_DIST matrix constructor.

Added support for the `Ginkgo <https://ginkgo-project.github.io/>`_  linear
algebra library. This support includes new ``SUNMatrix`` and ``SUNLinearSolver``
implementations, see the sections :numref:`SUNMatrix.Ginkgo` and
:numref:`SUNLinSol.Ginkgo`.

Added new ``NVector``, dense ``SUNMatrix``, and dense ``SUNLinearSolver``
implementations utilizing the `Kokkos Ecosystem <https://kokkos.org/>`_ for
performance portability, see sections :numref:`NVectors.Kokkos`,
:numref:`SUNMatrix.Kokkos`, and :numref:`SUNLinSol.Kokkos` for more information.

Fixed a bug in the CUDA and HIP vectors where :c:func:`N_VMaxNorm` would return
the minimum positive floating-point value for the zero vector.

### IDAS Changes in v5.4.0

CMake 3.18.0 or newer is now required for CUDA support.

A C++14 compliant compiler is now required for C++ based features and examples
e.g., CUDA, HIP, RAJA, Trilinos, SuperLU_DIST, MAGMA, GINKGO, and KOKKOS.

Added support for GPU enabled SuperLU_DIST and SuperLU_DIST v8.x.x. Removed
support for SuperLU_DIST v6.x.x or older. Fix mismatched definition and
declaration bug in SuperLU_DIST matrix constructor.

Added support for the `Ginkgo <https://ginkgo-project.github.io/>`_  linear
algebra library. This support includes new ``SUNMatrix`` and ``SUNLinearSolver``
implementations, see the sections :numref:`SUNMatrix.Ginkgo` and
:numref:`SUNLinSol.Ginkgo`.

Added new ``NVector``, dense ``SUNMatrix``, and dense ``SUNLinearSolver``
implementations utilizing the `Kokkos Ecosystem <https://kokkos.org/>`_ for
performance portability, see sections :numref:`NVectors.Kokkos`,
:numref:`SUNMatrix.Kokkos`, and :numref:`SUNLinSol.Kokkos` for more information.

Fixed a bug in the CUDA and HIP vectors where :c:func:`N_VMaxNorm` would return
the minimum positive floating-point value for the zero vector.k

### KINSOL Changes in v6.4.0

CMake 3.18.0 or newer is now required for CUDA support.

A C++14 compliant compiler is now required for C++ based features and examples
e.g., CUDA, HIP, RAJA, Trilinos, SuperLU_DIST, MAGMA, GINKGO, and KOKKOS.

Added support for GPU enabled SuperLU_DIST and SuperLU_DIST v8.x.x. Removed
support for SuperLU_DIST v6.x.x or older. Fix mismatched definition and
declaration bug in SuperLU_DIST matrix constructor.

Added support for the `Ginkgo <https://ginkgo-project.github.io/>`_  linear
algebra library. This support includes new ``SUNMatrix`` and ``SUNLinearSolver``
implementations, see the sections :numref:`SUNMatrix.Ginkgo` and
:numref:`SUNLinSol.Ginkgo`.

Added new ``NVector``, dense ``SUNMatrix``, and dense ``SUNLinearSolver``
implementations utilizing the `Kokkos Ecosystem <https://kokkos.org/>`_ for
performance portability, see sections :numref:`NVectors.Kokkos`,
:numref:`SUNMatrix.Kokkos`, and :numref:`SUNLinSol.Kokkos` for more information.

Fixed a bug in the CUDA and HIP vectors where :c:func:`N_VMaxNorm` would return
the minimum positive floating-point value for the zero vector.



## Changes to SUNDIALS in release 6.3.0

Added `GetUserData` functions in each package to retrieve the user data pointer
provided to `SetUserData` functions. See `ARKStepGetUserData`,
`ERKStepGetUserData`, `MRIStepGetUserData`, `CVodeGetUserData`,
`IDAGetUserData`, or `KINGetUserData` for more information.

Fixed a bug in `ERKStepReset`, `ERKStepReInit`, `ARKStepReset`, `ARKStepReInit`,
`MRIStepReset`, and `MRIStepReInit` where a previously-set value of *tstop* (from
a call to `ERKStepSetStopTime`, `ARKStepSetStopTime`, or `MRIStepSetStopTime`,
respectively) would not be cleared.

Updated `MRIStepReset` to call the corresponding `MRIStepInnerResetFn` with the same
(*tR*,*yR*) arguments for the `MRIStepInnerStepper` object that is used to evolve the
MRI "fast" time scale subproblems.

Added a new [example](examples/cvode/serial/cvRocket_dns.c) which
demonstrates using CVODE with a discontinuous right-hand-side function
and rootfinding.

Added a variety of embedded DIRK methods from [Kennedy & Carpenter,
NASA TM-2016-219173, 2016] and [Kennedy & Carpenter, Appl. Numer. Math., 146, 2019] to
ARKODE.

Fixed the unituitive behavior of the `USE_GENERIC_MATH` CMake option which
caused the double precision math functions to be used regardless of the value of
`SUNDIALS_PRECISION`. Now, SUNDIALS will use precision appropriate math
functions when they are available and the user may provide the math library to
link to via the advanced CMake option `SUNDIALS_MATH_LIBRARY`.

Changed `SUNDIALS_LOGGING_ENABLE_MPI` CMake option default to be 'OFF'.

### ARKODE Changes in v5.3.0

Added the functions :c:func:`ARKStepGetUserData`, :c:func:`ERKStepGetUserData`,
and :c:func:`MRIStepGetUserData` to retrieve the user data pointer provided to
:c:func:`ARKStepSetUserData`, :c:func:`ERKStepSetUserData`, and
:c:func:`MRIStepSetUserData`, respectively.

Fixed a bug in :c:func:`ERKStepReset()`, :c:func:`ERKStepReInit()`,
:c:func:`ARKStepReset()`, :c:func:`ARKStepReInit()`, :c:func:`MRIStepReset()`, and
:c:func:`MRIStepReInit()` where a previously-set value of *tstop* (from a call to
:c:func:`ERKStepSetStopTime()`, :c:func:`ARKStepSetStopTime()`, or
:c:func:`MRIStepSetStopTime()`, respectively) would not be cleared.

Updated :c:func:`MRIStepReset()` to call the corresponding
:c:type:`MRIStepInnerResetFn` with the same :math:`(t_R,y_R)` arguments for the
:c:type:`MRIStepInnerStepper` object that is used to evolve the MRI "fast" time
scale subproblems.

Added a variety of embedded DIRK methods from :cite:p:`KenCarp:16` and :cite:p:`KenCarp:19b`.

Fixed the unituitive behavior of the :cmakeop:`USE_GENERIC_MATH` CMake option which
caused the double precision math functions to be used regardless of the value of
:cmakeop:`SUNDIALS_PRECISION`. Now, SUNDIALS will use precision appropriate math
functions when they are available and the user may provide the math library to
link to via the advanced CMake option :cmakeop:`SUNDIALS_MATH_LIBRARY`.

Changed :cmakeop:`SUNDIALS_LOGGING_ENABLE_MPI` CMake option default to be 'OFF'.

### CVODE Changes in v6.3.0

Added the function :c:func:`CVodeGetUserData` to retrieve the user data pointer
provided to :c:func:`CVodeSetUserData`.

Added a new example, ``examples/cvode/serial/cvRocket_dns.c,`` which
demonstrates using CVODE with a discontinuous right-hand-side function
and rootfinding.

Fixed the unituitive behavior of the :cmakeop:`USE_GENERIC_MATH` CMake option which
caused the double precision math functions to be used regardless of the value of
:cmakeop:`SUNDIALS_PRECISION`. Now, SUNDIALS will use precision appropriate math
functions when they are available and the user may provide the math library to
link to via the advanced CMake option :cmakeop:`SUNDIALS_MATH_LIBRARY`.

Changed :cmakeop:`SUNDIALS_LOGGING_ENABLE_MPI` CMake option default to be 'OFF'.

### CVODES Changes in v6.3.0

Added the function :c:func:`CVodeGetUserData` to retrieve the user data pointer
provided to :c:func:`CVodeSetUserData`.

Fixed the unituitive behavior of the :cmakeop:`USE_GENERIC_MATH` CMake option which
caused the double precision math functions to be used regardless of the value of
:cmakeop:`SUNDIALS_PRECISION`. Now, SUNDIALS will use precision appropriate math
functions when they are available and the user may provide the math library to
link to via the advanced CMake option :cmakeop:`SUNDIALS_MATH_LIBRARY`.

Changed :cmakeop:`SUNDIALS_LOGGING_ENABLE_MPI` CMake option default to be 'OFF'.

### IDA Changes in v6.3.0

Added the function :c:func:`IDAGetUserData` to retrieve the user data pointer
provided to :c:func:`IDASetUserData`.

Fixed the unituitive behavior of the :cmakeop:`USE_GENERIC_MATH` CMake option which
caused the double precision math functions to be used regardless of the value of
:cmakeop:`SUNDIALS_PRECISION`. Now, SUNDIALS will use precision appropriate math
functions when they are available and the user may provide the math library to
link to via the advanced CMake option :cmakeop:`SUNDIALS_MATH_LIBRARY`.

Changed :cmakeop:`SUNDIALS_LOGGING_ENABLE_MPI` CMake option default to be 'OFF'.

### IDAS Changes in v5.3.0

Added the function :c:func:`IDAGetUserData` to retrieve the user data pointer
provided to :c:func:`IDASetUserData`.

Fixed the unituitive behavior of the :cmakeop:`USE_GENERIC_MATH` CMake option which
caused the double precision math functions to be used regardless of the value of
:cmakeop:`SUNDIALS_PRECISION`. Now, SUNDIALS will use precision appropriate math
functions when they are available and the user may provide the math library to
link to via the advanced CMake option :cmakeop:`SUNDIALS_MATH_LIBRARY`.

Changed :cmakeop:`SUNDIALS_LOGGING_ENABLE_MPI` CMake option default to be 'OFF'.

### KINSOL Changes in v6.3.0

Added the function :c:func:`KINGetUserData` to retrieve the user data pointer
provided to :c:func:`KINSetUserData`.

Fixed the unituitive behavior of the :cmakeop:`USE_GENERIC_MATH` CMake option which
caused the double precision math functions to be used regardless of the value of
:cmakeop:`SUNDIALS_PRECISION`. Now, SUNDIALS will use precision appropriate math
functions when they are available and the user may provide the math library to
link to via the advanced CMake option :cmakeop:`SUNDIALS_MATH_LIBRARY`.

Changed :cmakeop:`SUNDIALS_LOGGING_ENABLE_MPI` CMake option default to be 'OFF'.



## Changes to SUNDIALS in release 6.2.0

Added the `SUNLogger` API which provides a SUNDIALS-wide
mechanism for logging of errors, warnings, informational output,
and debugging output.

Deprecated the following functions, it is recommended to use the `SUNLogger` API
instead.

* `ARKStepSetDiagnostics`
* `ERKStepSetDiagnostics`
* `MRIStepSetDiagnostics`
* `KINSetInfoFile`
* `SUNNonlinSolSetPrintLevel_Newton`
* `SUNNonlinSolSetInfoFile_Newton`
* `SUNNonlinSolSetPrintLevel_FixedPoint`
* `SUNNonlinSolSetInfoFile_FixedPoint`
* `SUNLinSolSetInfoFile_PCG`
* `SUNLinSolSetPrintLevel_PCG`
* `SUNLinSolSetInfoFile_SPGMR`
* `SUNLinSolSetPrintLevel_SPGMR`
* `SUNLinSolSetInfoFile_SPFGMR`
* `SUNLinSolSetPrintLevel_SPFGMR`
* `SUNLinSolSetInfoFile_SPTFQM`
* `SUNLinSolSetPrintLevel_SPTFQMR`
* `SUNLinSolSetInfoFile_SPBCGS`
* `SUNLinSolSetPrintLevel_SPBCGS`

The `SUNLinSolSetInfoFile_**` and  `SUNNonlinSolSetInfoFile_*` family of
functions are now enabled by setting the CMake option `SUNDIALS_LOGGING_LEVEL`
to a value `>= 3`.

Added the function `SUNProfiler_Reset` to reset the region timings and counters
to zero.

Added the functions `ARKStepPrintAllStats`, `ERKStepPrintAllStats`,
`MRIStepPrintAll`, `CVodePrintAllStats`, `IDAPrintAllStats`, and
`KINPrintAllStats` to output all of the integrator, nonlinear solver, linear
solver, and other statistics in one call. The file `scripts/sundials_csv.py`
contains functions for parsing the comma-separated value output files.

Added functions to CVODE, CVODES, IDA, and IDAS to change the default step size
adaptivity parameters. For more information see the documentation for:

* `CVodeSetEtaFixedStepBounds`
* `CVodeSetEtaMaxFirstStep`
* `CVodeSetEtaMaxEarlyStep`
* `CVodeSetNumStepsEtaMaxEarlyStep`
* `CVodeSetEtaMax`
* `CVodeSetEtaMin`
* `CVodeSetEtaMinErrFail`
* `CVodeSetEtaMaxErrFail`
* `CVodeSetNumFailsEtaMaxErrFail`
* `CVodeSetEtaConvFail`
* `IDASetEtaFixedStepBounds`
* `IDAsetEtaMax`
* `IDASetEtaMin`
* `IDASetEtaLow`
* `IDASetEtaMinErrFail`
* `IDASetEtaConvFail`

Added the functions `CVodeSetDeltaGammaMaxLSetup` and
`CVodeSetDeltaGammaMaxBadJac` in CVODE and CVODES to adjust the `gamma` change
thresholds to require a linear solver setup or Jacobian/precondition update,
respectively.

Added the function `IDASetDetlaCjLSetup` in IDA and IDAS to adjust the parameter
that determines when a change in `c_j` requires calling the linear solver setup
function.

Added the function `MRIStepSetOrder` to select the default MRI method of a given
order.

Added support to CVODES for integrating IVPs with constraints using BDF methods
and projecting the solution onto the constraint manifold with a user defined
projection function. This implementation is accompanied by additions to the
CVODES user documentation and examples.

The behavior of `N_VSetKernelExecPolicy_Sycl` has been updated to be consistent
with the CUDA and HIP vectors. The input execution policies are now cloned and
may be freed after calling `N_VSetKernelExecPolicy_Sycl`. Additionally, `NULL`
inputs are now allowed and, if provided, will reset the vector execution
policies to the defaults.

Fixed the `SUNContext` convenience class for C++ users to disallow copy
construction and allow move construction.

A memory leak in the SYCL vector was fixed where the execution policies were
not freed when the vector was destroyed.

The include guard in `nvector_mpimanyvector.h` has been corrected to enable
using both the ManyVector and MPIManyVector NVector implementations in the same
simulation.

Changed exported SUNDIALS PETSc CMake targets to be INTERFACE IMPORTED instead
of UNKNOWN IMPORTED.

A bug was fixed in the integrator functions to retrieve the number of nonlinear
solver failures. The failure count returned was the number of failed *steps* due
to a nonlinear solver failure i.e., if a nonlinear solve failed with a stale
Jacobian or preconditioner but succeeded after updating the Jacobian or
preconditioner, the initial failure was not included in the nonlinear solver
failure count. The following functions have been updated to return the total
number of nonlinear solver failures:

* `ARKStepGetNumNonlinSolvConvFails`
* `ARKStepGetNonlinSolvStats`
* `MRIStepGetNumNonlinSolvConvFails`
* `MRIStepGetNonlinSolvStats`
* `CVodeGetNumNonlinSolvConvFails`
* `CVodeGetNonlinSolvStats`
* `CVodeGetSensNumNonlinSolvConvFails`
* `CVodeGetSensNonlinSolvStats`
* `CVodeGetStgrSensNumNonlinSolvConvFails`
* `CVodeGetStgrSensNonlinSolvStats`
* `IDAGetNumNonlinSolvConvFails`
* `IDAGetNonlinSolvStats`
* `IDAGetSensNumNonlinSolvConvFails`
* `IDAGetSensNonlinSolvStats`

As such users may see an increase in the number of failures reported from the
above functions. The following functions have been added to retrieve the number
of failed steps due to a nonlinear solver failure i.e., the counts previously
returned by the above functions:

* `ARKStepGetNumStepSolveFails`
* `MRIStepGetNumStepSolveFails`
* `CVodeGetNumStepSolveFails`
* `CVodeGetNumStepSensSolveFails`
* `CVodeGetNumStepStgrSensSolveFails`
* `IDAGetNumStepSolveFails`
* `IDAGetNumStepSensSolveFails`

### ARKODE Changes in v5.2.0

Added the :c:type:`SUNLogger` API which provides a SUNDIALS-wide
mechanism for logging of errors, warnings, informational output,
and debugging output.

Deprecated :c:func:`ARKStepSetDiagnostics`,
:c:func:`MRIStepSetDiagnostics`, :c:func:`ERKStepSetDiagnostics`,
:c:func:`SUNNonlinSolSetPrintLevel_Newton`,
:c:func:`SUNNonlinSolSetInfoFile_Newton`,
:c:func:`SUNNonlinSolSetPrintLevel_FixedPoint`,
:c:func:`SUNNonlinSolSetInfoFile_FixedPoint`,
:c:func:`SUNLinSolSetInfoFile_PCG`, :c:func:`SUNLinSolSetPrintLevel_PCG`,
:c:func:`SUNLinSolSetInfoFile_SPGMR`, :c:func:`SUNLinSolSetPrintLevel_SPGMR`,
:c:func:`SUNLinSolSetInfoFile_SPFGMR`, :c:func:`SUNLinSolSetPrintLevel_SPFGMR`,
:c:func:`SUNLinSolSetInfoFile_SPTFQM`, :c:func:`SUNLinSolSetPrintLevel_SPTFQMR`,
:c:func:`SUNLinSolSetInfoFile_SPBCGS`, :c:func:`SUNLinSolSetPrintLevel_SPBCGS`
it is recommended to use the `SUNLogger` API instead. The ``SUNLinSolSetInfoFile_**``
and ``SUNNonlinSolSetInfoFile_*`` family of functions are now enabled
by setting the CMake option :cmakeop:`SUNDIALS_LOGGING_LEVEL` to a value ``>= 3``.

Added the function :c:func:`SUNProfiler_Reset` to reset the region timings and
counters to zero.

Added the functions :c:func:`ARKStepPrintAllStats`,
:c:func:`ERKStepPrintAllStats`, and :c:func:`MRIStepPrintAll` to output all of
the integrator, nonlinear solver, linear solver, and other statistics in one
call. The file ``scripts/sundials_csv.py`` contains functions for parsing the
comma-separated value output files.

Added the functions :c:func:`ARKStepSetDeduceImplicitRhs` and
:c:func:`MRIStepSetDeduceImplicitRhs` to optionally remove an evaluation of the
implicit right-hand side function after nonlinear solves. See
:numref:`ARKODE.Mathematics.Nonlinear`, for considerations on using this
optimization.

Added the function :c:func:`MRIStepSetOrder` to select the default MRI method of
a given order.

The behavior of :c:func:`N_VSetKernelExecPolicy_Sycl` has been updated to be
consistent with the CUDA and HIP vectors. The input execution policies are now
cloned and may be freed after calling :c:func:`N_VSetKernelExecPolicy_Sycl`.
Additionally, ``NULL`` inputs are now allowed and, if provided, will reset the
vector execution policies to the defaults.

Fixed the :c:type:`SUNContext` convenience class for C++ users to disallow copy
construction and allow move construction.

A memory leak in the SYCL vector was fixed where the execution policies were
not freed when the vector was destroyed.

The include guard in ``nvector_mpimanyvector.h`` has been corrected to enable
using both the ManyVector and MPIManyVector NVector implementations in the same
simulation.

Changed exported SUNDIALS PETSc CMake targets to be INTERFACE IMPORTED instead
of UNKNOWN IMPORTED.

A bug was fixed in the functions
:c:func:`ARKStepGetNumNonlinSolvConvFails`,
:c:func:`ARKStepGetNonlinSolvStats`,
:c:func:`MRIStepGetNumNonlinSolvConvFails`, and
:c:func:`MRIStepGetNonlinSolvStats`
where the number of nonlinear solver failures returned was the number of failed
*steps* due to a nonlinear solver failure i.e., if a nonlinear solve failed with
a stale Jacobian or preconditioner but succeeded after updating the Jacobian or
preconditioner, the initial failure was not included in the nonlinear solver
failure count. These functions have been updated to return the total number of
nonlinear solver failures. As such users may see an increase in the number of
failures reported.

The functions :c:func:`ARKStepGetNumStepSolveFails` and
:c:func:`MRIStepGetNumStepSolveFails` have been added to retrieve the number of
failed steps due to a nonlinear solver failure. The counts returned from these
functions will match those previously returned by
:c:func:`ARKStepGetNumNonlinSolvConvFails`,
:c:func:`ARKStepGetNonlinSolvStats`,
:c:func:`MRIStepGetNumNonlinSolvConvFails`, and
:c:func:`MRIStepGetNonlinSolvStats`.



### CVODE Changes in v6.2.0

Added the :c:type:`SUNLogger` API which provides a SUNDIALS-wide
mechanism for logging of errors, warnings, informational output,
and debugging output.

Deprecated :c:func:`SUNNonlinSolSetPrintLevel_Newton`,
:c:func:`SUNNonlinSolSetInfoFile_Newton`,
:c:func:`SUNNonlinSolSetPrintLevel_FixedPoint`,
:c:func:`SUNNonlinSolSetInfoFile_FixedPoint`,
:c:func:`SUNLinSolSetInfoFile_PCG`, :c:func:`SUNLinSolSetPrintLevel_PCG`,
:c:func:`SUNLinSolSetInfoFile_SPGMR`, :c:func:`SUNLinSolSetPrintLevel_SPGMR`,
:c:func:`SUNLinSolSetInfoFile_SPFGMR`, :c:func:`SUNLinSolSetPrintLevel_SPFGMR`,
:c:func:`SUNLinSolSetInfoFile_SPTFQM`, :c:func:`SUNLinSolSetPrintLevel_SPTFQMR`,
:c:func:`SUNLinSolSetInfoFile_SPBCGS`, :c:func:`SUNLinSolSetPrintLevel_SPBCGS`
it is recommended to use the `SUNLogger` API instead. The ``SUNLinSolSetInfoFile_**``
and ``SUNNonlinSolSetInfoFile_*`` family of functions are now enabled
by setting the CMake option :cmakeop:`SUNDIALS_LOGGING_LEVEL` to a value ``>= 3``.

Added the function :c:func:`SUNProfiler_Reset` to reset the region timings and
counters to zero.

Added the function :c:func:`CVodePrintAllStats` to output all of the integrator,
nonlinear solver, linear solver, and other statistics in one call. The file
``scripts/sundials_csv.py`` contains functions for parsing the comma-separated
value output files.

Added the functions
:c:func:`CVodeSetEtaFixedStepBounds`,
:c:func:`CVodeSetEtaMaxFirstStep`,
:c:func:`CVodeSetEtaMaxEarlyStep`,
:c:func:`CVodeSetNumStepsEtaMaxEarlyStep`,
:c:func:`CVodeSetEtaMax`,
:c:func:`CVodeSetEtaMin`,
:c:func:`CVodeSetEtaMinErrFail`,
:c:func:`CVodeSetEtaMaxErrFail`,
:c:func:`CVodeSetNumFailsEtaMaxErrFail`, and
:c:func:`CVodeSetEtaConvFail` to adjust various parameters controlling changes
in step size.

Added the functions :c:func:`CVodeSetDeltaGammaMaxLSetup` and
:c:func:`CVodeSetDeltaGammaMaxBadJac` to adjust the :math:`\gamma` change
thresholds to require a linear solver setup or Jacobian/precondition update,
respectively.

The behavior of :c:func:`N_VSetKernelExecPolicy_Sycl` has been updated to be
consistent with the CUDA and HIP vectors. The input execution policies are now
cloned and may be freed after calling :c:func:`N_VSetKernelExecPolicy_Sycl`.
Additionally, ``NULL`` inputs are now allowed and, if provided, will reset the
vector execution policies to the defaults.

Fixed the :c:type:`SUNContext` convenience class for C++ users to disallow copy
construction and allow move construction.

A memory leak in the SYCL vector was fixed where the execution policies were
not freed when the vector was destroyed.

The include guard in ``nvector_mpimanyvector.h`` has been corrected to enable
using both the ManyVector and MPIManyVector NVector implementations in the same
simulation.

Changed exported SUNDIALS PETSc CMake targets to be INTERFACE IMPORTED instead
of UNKNOWN IMPORTED.

A bug was fixed in the functions :c:func:`CVodeGetNumNonlinSolvConvFails` and
:c:func:`CVodeGetNonlinSolvStats` where the number of nonlinear solver failures
returned was the number of failed *steps* due to a nonlinear solver failure
i.e., if a nonlinear solve failed with a stale Jacobian or preconditioner but
succeeded after updating the Jacobian or preconditioner, the initial failure was
not included in the nonlinear solver failure count. These functions have been
updated to return the total number of nonlinear solver failures. As such users
may see an increase in the number of failures reported.

The function :c:func:`CVodeGetNumStepSolveFails` has been added to retrieve the
number of failed steps due to a nonlinear solver failure. The count returned by
this function will match those previously returned by
:c:func:`CVodeGetNumNonlinSolvConvFails` and :c:func:`CVodeGetNonlinSolvStats`.

### CVODES Changes in v6.2.0

Added the :c:type:`SUNLogger` API which provides a SUNDIALS-wide
mechanism for logging of errors, warnings, informational output,
and debugging output.

Deprecated :c:func:`SUNNonlinSolSetPrintLevel_Newton`,
:c:func:`SUNNonlinSolSetInfoFile_Newton`,
:c:func:`SUNNonlinSolSetPrintLevel_FixedPoint`,
:c:func:`SUNNonlinSolSetInfoFile_FixedPoint`,
:c:func:`SUNLinSolSetInfoFile_PCG`, :c:func:`SUNLinSolSetPrintLevel_PCG`,
:c:func:`SUNLinSolSetInfoFile_SPGMR`, :c:func:`SUNLinSolSetPrintLevel_SPGMR`,
:c:func:`SUNLinSolSetInfoFile_SPFGMR`, :c:func:`SUNLinSolSetPrintLevel_SPFGMR`,
:c:func:`SUNLinSolSetInfoFile_SPTFQM`, :c:func:`SUNLinSolSetPrintLevel_SPTFQMR`,
:c:func:`SUNLinSolSetInfoFile_SPBCGS`, :c:func:`SUNLinSolSetPrintLevel_SPBCGS`
it is recommended to use the `SUNLogger` API instead. The ``SUNLinSolSetInfoFile_**``
and ``SUNNonlinSolSetInfoFile_*`` family of functions are now enabled
by setting the CMake option :cmakeop:`SUNDIALS_LOGGING_LEVEL` to a value ``>= 3``.

Added the function :c:func:`SUNProfiler_Reset` to reset the region timings and
counters to zero.

Added the function :c:func:`CVodePrintAllStats` to output all of the integrator,
nonlinear solver, linear solver, and other statistics in one call. The file
``scripts/sundials_csv.py`` contains functions for parsing the comma-separated
value output files.

Added support for integrating IVPs with constraints using BDF methods
and projecting the solution onto the constraint manifold with a user
defined projection function. This implementation is accompanied by
additions to user documentation and CVODES examples. See
:c:func:`CVodeSetProjFn` for more information.

Added the functions
:c:func:`CVodeSetEtaFixedStepBounds`,
:c:func:`CVodeSetEtaMaxFirstStep`,
:c:func:`CVodeSetEtaMaxEarlyStep`,
:c:func:`CVodeSetNumStepsEtaMaxEarlyStep`,
:c:func:`CVodeSetEtaMax`,
:c:func:`CVodeSetEtaMin`,
:c:func:`CVodeSetEtaMinErrFail`,
:c:func:`CVodeSetEtaMaxErrFail`,
:c:func:`CVodeSetNumFailsEtaMaxErrFail`, and
:c:func:`CVodeSetEtaConvFail` to adjust various parameters controlling changes
in step size.

Added the functions :c:func:`CVodeSetDeltaGammaMaxLSetup` and
:c:func:`CVodeSetDeltaGammaMaxBadJac` to adjust the :math:`\gamma` change
thresholds to require a linear solver setup or Jacobian/precondition update,
respectively.

The behavior of :c:func:`N_VSetKernelExecPolicy_Sycl` has been updated to be
consistent with the CUDA and HIP vectors. The input execution policies are now
cloned and may be freed after calling :c:func:`N_VSetKernelExecPolicy_Sycl`.
Additionally, ``NULL`` inputs are now allowed and, if provided, will reset the
vector execution policies to the defaults.

Fixed the :c:type:`SUNContext` convenience class for C++ users to disallow copy
construction and allow move construction.

A memory leak in the SYCL vector was fixed where the execution policies were
not freed when the vector was destroyed.

The include guard in ``nvector_mpimanyvector.h`` has been corrected to enable
using both the ManyVector and MPIManyVector NVector implementations in the same
simulation.

Changed exported SUNDIALS PETSc CMake targets to be INTERFACE IMPORTED instead
of UNKNOWN IMPORTED.

A bug was fixed in the functions
:c:func:`CVodeGetNumNonlinSolvConvFails`,
:c:func:`CVodeGetNonlinSolvStats`,
:c:func:`CVodeGetSensNumNonlinSolvConvFails`,
:c:func:`CVodeGetSensNonlinSolvStats`,
:c:func:`CVodeGetStgrSensNumNonlinSolvConvFails`, and
:c:func:`CVodeGetStgrSensNonlinSolvStats`
where the number of nonlinear solver failures returned was the number of failed
*steps* due to a nonlinear solver failure i.e., if a nonlinear solve failed
with a stale Jacobian or preconditioner but succeeded after updating the
Jacobian or preconditioner, the initial failure was not included in the
nonlinear solver failure count. These functions have been updated to return the
total number of nonlinear solver failures. As such users may see an increase in
the number of failures reported.

The functions
:c:func:`CVodeGetNumStepSolveFails`,
:c:func:`CVodeGetNumStepSensSolveFails`, and
:c:func:`CVodeGetNumStepStgrSensSolveFails`
have been added to retrieve the number of failed steps due to a nonlinear solver
failure. The counts returned from these functions will match those previously
returned by
:c:func:`CVodeGetNumNonlinSolvConvFails`,
:c:func:`CVodeGetNonlinSolvStats`,
:c:func:`CVodeGetSensNumNonlinSolvConvFails`,
:c:func:`CVodeGetSensNonlinSolvStats`,
:c:func:`CVodeGetStgrSensNumNonlinSolvConvFails`, and
:c:func:`CVodeGetStgrSensNonlinSolvStats`.

### IDA Changes in v6.2.0

Added the :c:type:`SUNLogger` API which provides a SUNDIALS-wide
mechanism for logging of errors, warnings, informational output,
and debugging output.

Deprecated :c:func:`SUNNonlinSolSetPrintLevel_Newton`,
:c:func:`SUNNonlinSolSetInfoFile_Newton`,
:c:func:`SUNNonlinSolSetPrintLevel_FixedPoint`,
:c:func:`SUNNonlinSolSetInfoFile_FixedPoint`,
:c:func:`SUNLinSolSetInfoFile_PCG`, :c:func:`SUNLinSolSetPrintLevel_PCG`,
:c:func:`SUNLinSolSetInfoFile_SPGMR`, :c:func:`SUNLinSolSetPrintLevel_SPGMR`,
:c:func:`SUNLinSolSetInfoFile_SPFGMR`, :c:func:`SUNLinSolSetPrintLevel_SPFGMR`,
:c:func:`SUNLinSolSetInfoFile_SPTFQM`, :c:func:`SUNLinSolSetPrintLevel_SPTFQMR`,
:c:func:`SUNLinSolSetInfoFile_SPBCGS`, :c:func:`SUNLinSolSetPrintLevel_SPBCGS`
it is recommended to use the `SUNLogger` API instead. The ``SUNLinSolSetInfoFile_**``
and ``SUNNonlinSolSetInfoFile_*`` family of functions are now enabled
by setting the CMake option :cmakeop:`SUNDIALS_LOGGING_LEVEL` to a value ``>= 3``.

Added the function :c:func:`SUNProfiler_Reset` to reset the region timings and
counters to zero.

Added the function :c:func:`IDAPrintAllStats` to output all of the integrator,
nonlinear solver, linear solver, and other statistics in one call. The file
``scripts/sundials_csv.py`` contains functions for parsing the comma-separated
value output files.

Added the function :c:func:`IDASetDetlaCjLSetup` to adjust the parameter that
determines when a change in :math:`c_j` requires calling the linear solver setup
function.

Added the functions :c:func:`IDASetEtaFixedStepBounds`, :c:func:`IDASetEtaMax`,
:c:func:`IDASetEtaMin`, :c:func:`IDASetEtaLow`, :c:func:`IDASetEtaMinErrFail`,
and :c:func:`IDASetEtaConvFail` to adjust various parameters controlling changes
in step size.

Added the function :c:func:`IDASetMinStep` to set a minimum step size.

The behavior of :c:func:`N_VSetKernelExecPolicy_Sycl` has been updated to be
consistent with the CUDA and HIP vectors. The input execution policies are now
cloned and may be freed after calling :c:func:`N_VSetKernelExecPolicy_Sycl`.
Additionally, ``NULL`` inputs are now allowed and, if provided, will reset the
vector execution policies to the defaults.

Fixed the :c:type:`SUNContext` convenience class for C++ users to disallow copy
construction and allow move construction.

A memory leak in the SYCL vector was fixed where the execution policies were
not freed when the vector was destroyed.

The include guard in ``nvector_mpimanyvector.h`` has been corrected to enable
using both the ManyVector and MPIManyVector NVector implementations in the same
simulation.

Changed exported SUNDIALS PETSc CMake targets to be INTERFACE IMPORTED instead
of UNKNOWN IMPORTED.

A bug was fixed in the functions :c:func:`IDAGetNumNonlinSolvConvFails` and
:c:func:`IDAGetNonlinSolvStats` where the number of nonlinear solver failures
returned was the number of failed *steps* due to a nonlinear solver failure
i.e., if a nonlinear solve failed with a stale Jacobian or preconditioner but
succeeded after updating the Jacobian or preconditioner, the initial failure was
not included in the nonlinear solver failure count. These functions have been
updated to return the total number of nonlinear solver failures. As such users
may see an increase in the number of failures reported.

The function :c:func:`IDAGetNumStepSolveFails` has been added to retrieve the
number of failed steps due to a nonlinear solver failure. The count returned by
this function will match those previously returned by
:c:func:`IDAGetNumNonlinSolvConvFails` and :c:func:`IDAGetNonlinSolvStats`.

### IDAS Changes in v5.2.0

Added the :c:type:`SUNLogger` API which provides a SUNDIALS-wide
mechanism for logging of errors, warnings, informational output,
and debugging output.

Deprecated :c:func:`SUNNonlinSolSetPrintLevel_Newton`,
:c:func:`SUNNonlinSolSetInfoFile_Newton`,
:c:func:`SUNNonlinSolSetPrintLevel_FixedPoint`,
:c:func:`SUNNonlinSolSetInfoFile_FixedPoint`,
:c:func:`SUNLinSolSetInfoFile_PCG`, :c:func:`SUNLinSolSetPrintLevel_PCG`,
:c:func:`SUNLinSolSetInfoFile_SPGMR`, :c:func:`SUNLinSolSetPrintLevel_SPGMR`,
:c:func:`SUNLinSolSetInfoFile_SPFGMR`, :c:func:`SUNLinSolSetPrintLevel_SPFGMR`,
:c:func:`SUNLinSolSetInfoFile_SPTFQM`, :c:func:`SUNLinSolSetPrintLevel_SPTFQMR`,
:c:func:`SUNLinSolSetInfoFile_SPBCGS`, :c:func:`SUNLinSolSetPrintLevel_SPBCGS`
it is recommended to use the `SUNLogger` API instead. The ``SUNLinSolSetInfoFile_**``
and ``SUNNonlinSolSetInfoFile_*`` family of functions are now enabled
by setting the CMake option :cmakeop:`SUNDIALS_LOGGING_LEVEL` to a value ``>= 3``.

Added the function :c:func:`SUNProfiler_Reset` to reset the region timings and
counters to zero.

Added the function :c:func:`IDAPrintAllStats` to output all of the integrator,
nonlinear solver, linear solver, and other statistics in one call. The file
``scripts/sundials_csv.py`` contains functions for parsing the comma-separated
value output files.

Added the function :c:func:`IDASetDetlaCjLSetup` to adjust the parameter that
determines when a change in :math:`c_j` requires calling the linear solver setup
function.

Added the functions :c:func:`IDASetEtaFixedStepBounds`, :c:func:`IDASetEtaMax`,
:c:func:`IDASetEtaMin`, :c:func:`IDASetEtaLow`, :c:func:`IDASetEtaMinErrFail`,
and :c:func:`IDASetEtaConvFail` to adjust various parameters controlling changes
in step size.

Added the function :c:func:`IDASetMinStep` to set a minimum step size.

The behavior of :c:func:`N_VSetKernelExecPolicy_Sycl` has been updated to be
consistent with the CUDA and HIP vectors. The input execution policies are now
cloned and may be freed after calling :c:func:`N_VSetKernelExecPolicy_Sycl`.
Additionally, ``NULL`` inputs are now allowed and, if provided, will reset the
vector execution policies to the defaults.

Fixed the :c:type:`SUNContext` convenience class for C++ users to disallow copy
construction and allow move construction.

A memory leak in the SYCL vector was fixed where the execution policies were
not freed when the vector was destroyed.

The include guard in ``nvector_mpimanyvector.h`` has been corrected to enable
using both the ManyVector and MPIManyVector NVector implementations in the same
simulation.

Changed exported SUNDIALS PETSc CMake targets to be INTERFACE IMPORTED instead
of UNKNOWN IMPORTED.

A bug was fixed in the functions
:c:func:`IDAGetNumNonlinSolvConvFails`,
:c:func:`IDAGetNonlinSolvStats`,
:c:func:`IDAGetSensNumNonlinSolvConvFails`, and
:c:func:`IDAGetSensNonlinSolvStats`
where the number of nonlinear solver failures returned was the number of failed
*steps* due to a nonlinear solver failure i.e., if a nonlinear solve failed
with a stale Jacobian or preconditioner but succeeded after updating the
Jacobian or preconditioner, the initial failure was not included in the
nonlinear solver failure count.  These functions have been updated to return the
total number of nonlinear solver failures. As such users may see an increase in
the number of failures reported.

The functions :c:func:`IDAGetNumStepSolveFails` and
:c:func:`IDAGetNumStepSensSolveFails` have been added to retrieve the number of
failed steps due to a nonlinear solver failure. The counts returned from these
functions will match those previously returned by
:c:func:`IDAGetNumNonlinSolvConvFails`,
:c:func:`IDAGetNonlinSolvStats`,
:c:func:`IDAGetSensNumNonlinSolvConvFails`, and
:c:func:`IDAGetSensNonlinSolvStats`.

### KINSOL Changes in v6.2.0

Added the :c:type:`SUNLogger` API which provides a SUNDIALS-wide
mechanism for logging of errors, warnings, informational output,
and debugging output.

Deprecated :c:func:`KINSetInfoFile`, :c:func:`KINSetDebugFile`,
:c:func:`SUNNonlinSolSetPrintLevel_Newton`,
:c:func:`SUNNonlinSolSetInfoFile_Newton`,
:c:func:`SUNNonlinSolSetPrintLevel_FixedPoint`,
:c:func:`SUNNonlinSolSetInfoFile_FixedPoint`,
:c:func:`SUNLinSolSetInfoFile_PCG`, :c:func:`SUNLinSolSetPrintLevel_PCG`,
:c:func:`SUNLinSolSetInfoFile_SPGMR`, :c:func:`SUNLinSolSetPrintLevel_SPGMR`,
:c:func:`SUNLinSolSetInfoFile_SPFGMR`, :c:func:`SUNLinSolSetPrintLevel_SPFGMR`,
:c:func:`SUNLinSolSetInfoFile_SPTFQM`, :c:func:`SUNLinSolSetPrintLevel_SPTFQMR`,
:c:func:`SUNLinSolSetInfoFile_SPBCGS`, :c:func:`SUNLinSolSetPrintLevel_SPBCGS`
it is recommended to use the `SUNLogger` API instead. The ``SUNLinSolSetInfoFile_**``
and ``SUNNonlinSolSetInfoFile_*`` family of functions are now enabled
by setting the CMake option :cmakeop:`SUNDIALS_LOGGING_LEVEL` to a value ``>= 3``.

Added the function :c:func:`SUNProfiler_Reset` to reset the region timings and
counters to zero.

Added the function :c:func:`KINPrintAllStats` to output all of the nonlinear
solver, linear solver, and other statistics in one call.  The file
``scripts/sundials_csv.py`` contains functions for parsing the comma-separated
value output files.

The behavior of :c:func:`N_VSetKernelExecPolicy_Sycl` has been updated to be
consistent with the CUDA and HIP vectors. The input execution policies are now
cloned and may be freed after calling :c:func:`N_VSetKernelExecPolicy_Sycl`.
Additionally, ``NULL`` inputs are now allowed and, if provided, will reset the
vector execution policies to the defaults.

Fixed the :c:type:`SUNContext` convenience class for C++ users to disallow copy
construction and allow move construction.

A memory leak in the SYCL vector was fixed where the execution policies were
not freed when the vector was destroyed.

The include guard in ``nvector_mpimanyvector.h`` has been corrected to enable
using both the ManyVector and MPIManyVector NVector implementations in the same
simulation.

Changed exported SUNDIALS PETSc CMake targets to be INTERFACE IMPORTED instead
of UNKNOWN IMPORTED.



## Changes to SUNDIALS in release 6.1.1

Fixed exported `SUNDIALSConfig.cmake`.

Fixed Fortran interface to `MRIStepInnerStepper` and `MRIStepCoupling`
structures and functions.

Added new Fortran example program,
`examples/arkode/F2003_serial/ark_kpr_mri_f2003.f90` demonstrating MRI
capabilities.

### ARKODE Changes in v5.1.1

Fixed exported ``SUNDIALSConfig.cmake``.

Fixed Fortran interface to :c:type:`MRIStepInnerStepper` and :c:type:`MRIStepCoupling`
structures and functions.

Added new Fortran example program,
``examples/arkode/F2003_serial/ark_kpr_mri_f2003.f90`` demonstrating MRI
capabilities.


### CVODE Changes in v6.1.1

Fixed exported ``SUNDIALSConfig.cmake``.

### CVODES Changes in v6.1.1

Fixed exported ``SUNDIALSConfig.cmake``.

### IDA Changes in v6.1.1

Fixed exported ``SUNDIALSConfig.cmake``.

### IDAS Changes in v5.1.1

Fixed exported ``SUNDIALSConfig.cmake``.

### KINSOL Changes in v6.1.1

Fixed exported ``SUNDIALSConfig.cmake``.


## Changes to SUNDIALS in release 6.1.0

Added new reduction implementations for the CUDA and HIP NVECTORs that use
shared memory (local data storage) instead of atomics. These new implementations
are recommended when the target hardware does not provide atomic support for the
floating point precision that SUNDIALS is being built with. The HIP vector uses
these by default, but the `N_VSetKernelExecPolicy_Cuda` and
`N_VSetKernelExecPolicy_Hip` functions can be used to choose between
different reduction implementations.

`SUNDIALS::<lib>` targets with no static/shared suffix have been added for use
within the build directory (this mirrors the targets exported on installation).

`CMAKE_C_STANDARD` is now set to 99 by default.

Fixed exported `SUNDIALSConfig.cmake` when profiling is enabled without Caliper.

Fixed `sundials_export.h` include in `sundials_config.h`.

Fixed memory leaks in the SUNLINSOL_SUPERLUMT linear solver.

### ARKODE Changes in v5.1.0

Added new reduction implementations for the CUDA and HIP NVECTORs that use
shared memory (local data storage) instead of atomics. These new implementations
are recommended when the target hardware does not provide atomic support for the
floating point precision that SUNDIALS is being built with. The HIP vector uses
these by default, but the :c:func:`N_VSetKernelExecPolicy_Cuda` and
:c:func:`N_VSetKernelExecPolicy_Hip` functions can be used to choose between
different reduction implementations.

``SUNDIALS::<lib>`` targets with no static/shared suffix have been added for use
within the build directory (this mirrors the targets exported on installation).

:cmakeop:`CMAKE_C_STANDARD` is now set to 99 by default.

Fixed exported ``SUNDIALSConfig.cmake`` when profiling is enabled without Caliper.

Fixed ``sundials_export.h`` include in ``sundials_config.h``.

Fixed memory leaks in the SUNLINSOL_SUPERLUMT linear solver.


### CVODE Changes in v6.1.0

Added new reduction implementations for the CUDA and HIP NVECTORs that use
shared memory (local data storage) instead of atomics. These new implementations
are recommended when the target hardware does not provide atomic support for the
floating point precision that SUNDIALS is being built with. The HIP vector uses
these by default, but the :c:func:`N_VSetKernelExecPolicy_Cuda` and
:c:func:`N_VSetKernelExecPolicy_Hip` functions can be used to choose between
different reduction implementations.

``SUNDIALS::<lib>`` targets with no static/shared suffix have been added for use
within the build directory (this mirrors the targets exported on installation).

:cmakeop:`CMAKE_C_STANDARD` is now set to 99 by default.

Fixed exported ``SUNDIALSConfig.cmake`` when profiling is enabled without Caliper.

Fixed ``sundials_export.h`` include in ``sundials_config.h``.

Fixed memory leaks in the SUNLINSOL_SUPERLUMT linear solver.

### CVODES Changes in v6.1.0

Added new reduction implementations for the CUDA and HIP NVECTORs that use
shared memory (local data storage) instead of atomics. These new implementations
are recommended when the target hardware does not provide atomic support for the
floating point precision that SUNDIALS is being built with. The HIP vector uses
these by default, but the :c:func:`N_VSetKernelExecPolicy_Cuda` and
:c:func:`N_VSetKernelExecPolicy_Hip` functions can be used to choose between
different reduction implementations.

``SUNDIALS::<lib>`` targets with no static/shared suffix have been added for use
within the build directory (this mirrors the targets exported on installation).

:cmakeop:`CMAKE_C_STANDARD` is now set to 99 by default.

Fixed exported ``SUNDIALSConfig.cmake`` when profiling is enabled without Caliper.

Fixed ``sundials_export.h`` include in ``sundials_config.h``.

Fixed memory leaks in the SUNLINSOL_SUPERLUMT linear solver.

### IDA Changes in v6.1.0

Added new reduction implementations for the CUDA and HIP NVECTORs that use
shared memory (local data storage) instead of atomics. These new implementations
are recommended when the target hardware does not provide atomic support for the
floating point precision that SUNDIALS is being built with. The HIP vector uses
these by default, but the :c:func:`N_VSetKernelExecPolicy_Cuda` and
:c:func:`N_VSetKernelExecPolicy_Hip` functions can be used to choose between
different reduction implementations.

``SUNDIALS::<lib>`` targets with no static/shared suffix have been added for use
within the build directory (this mirrors the targets exported on installation).

:cmakeop:`CMAKE_C_STANDARD` is now set to 99 by default.

Fixed exported ``SUNDIALSConfig.cmake`` when profiling is enabled without Caliper.

Fixed ``sundials_export.h`` include in ``sundials_config.h``.

Fixed memory leaks in the SUNLINSOL_SUPERLUMT linear solver.

### IDAS Changes in v5.1.0

Added new reduction implementations for the CUDA and HIP NVECTORs that use
shared memory (local data storage) instead of atomics. These new implementations
are recommended when the target hardware does not provide atomic support for the
floating point precision that SUNDIALS is being built with. The HIP vector uses
these by default, but the :c:func:`N_VSetKernelExecPolicy_Cuda` and
:c:func:`N_VSetKernelExecPolicy_Hip` functions can be used to choose between
different reduction implementations.

``SUNDIALS::<lib>`` targets with no static/shared suffix have been added for use
within the build directory (this mirrors the targets exported on installation).

:cmakeop:`CMAKE_C_STANDARD` is now set to 99 by default.

Fixed exported ``SUNDIALSConfig.cmake`` when profiling is enabled without Caliper.

Fixed ``sundials_export.h`` include in ``sundials_config.h``.

Fixed memory leaks in the SUNLINSOL_SUPERLUMT linear solver.

### KINSOL Changes in v6.1.0

Added new reduction implementations for the CUDA and HIP NVECTORs that use
shared memory (local data storage) instead of atomics. These new implementations
are recommended when the target hardware does not provide atomic support for the
floating point precision that SUNDIALS is being built with. The HIP vector uses
these by default, but the :c:func:`N_VSetKernelExecPolicy_Cuda` and
:c:func:`N_VSetKernelExecPolicy_Hip` functions can be used to choose between
different reduction implementations.

``SUNDIALS::<lib>`` targets with no static/shared suffix have been added for use
within the build directory (this mirrors the targets exported on installation).

:cmakeop:`CMAKE_C_STANDARD` is now set to 99 by default.

Fixed exported ``SUNDIALSConfig.cmake`` when profiling is enabled without Caliper.

Fixed ``sundials_export.h`` include in ``sundials_config.h``.

Fixed memory leaks in the SUNLINSOL_SUPERLUMT linear solver.



## Changes to SUNDIALS in release 6.0.0

### SUNContext

SUNDIALS v6.0.0 introduces a new `SUNContext` object on which all other SUNDIALS
objects depend. As such, the constructors for all SUNDIALS packages, vectors,
matrices, linear solvers, nonlinear solvers, and memory helpers have been
updated to accept a context as the last input. Users upgrading to SUNDIALS
v6.0.0 will need to call `SUNContext_Create` to create a context object with
before calling any other SUNDIALS library function, and then provide this object
to other SUNDIALS constructors. The context object has been introduced to allow
SUNDIALS to provide new features, such as the profiling/instrumentation also
introduced in this release, while maintaining thread-safety. See the
documentation section on the `SUNContext` for more details.

A script `upgrade-to-sundials-6-from-5.sh` has been provided with the release
(obtainable from the GitHub release page) to help ease the transition to
SUNDIALS v6.0.0. The script will add a `SUNCTX_PLACEHOLDER` argument to all of
the calls to SUNDIALS constructors that now require a `SUNContext` object. It
can also update deprecated SUNDIALS constants/types to the new names. It can be
run like this:

```
> ./upgrade-to-sundials-6-from-5.sh <files to update>
```

### SUNProfiler

A capability to profile/instrument SUNDIALS library code has been added. This
can be enabled with the CMake option `SUNDIALS_BUILD_WITH_PROFILING`. A built-in
profiler will be used by default, but the
[Caliper](https://github.com/LLNL/Caliper) library can also be used instead with
the CMake option `ENABLE_CALIPER`. See the documentation section on profiling
for more details.  **WARNING**: Profiling will impact performance, and should be
enabled judiciously.

### SUNMemoryHelper

The `SUNMemoryHelper` functions `Alloc`, `Dealloc`, and `Copy` have been updated
to accept an opaque handle as the last input. At a minimum, existing
`SUNMemoryHelper` implementations will need to update these functions to accept
the additional argument. Typically, this handle is the execution stream (e.g., a
CUDA/HIP stream or SYCL queue) for the operation. The CUDA, HIP, and SYCL
`SUNMemoryHelper` implementations have been updated accordingly. Additionally,
the constructor for the SYCL implementation has been updated to remove the SYCL
queue as an input.

### NVector

Two new optional vector operations, `N_VDotProdMultiLocal` and
`N_VDotProdMultiAllReduce`, have been added to support low-synchronization
methods for Anderson acceleration.

The CUDA, HIP, and SYCL execution policies have been moved from the `sundials`
namespace to the `sundials::cuda`, `sundials::hip`, and `sundials::sycl`
namespaces respectively. Accordingly, the prefixes "Cuda", "Hip", and "Sycl"
have been removed from the execution policy classes and methods.

The `Sundials` namespace used by the Trilinos Tpetra NVector has been replaced
with the `sundials::trilinos::nvector_tpetra` namespace.

The serial, PThreads, PETSc, *hypre*, Parallel, OpenMP_DEV, and OpenMP vector
functions `N_VCloneVectorArray_*` and `N_VDestroyVectorArray_*` have been
deprecated. The generic `N_VCloneVectorArray` and `N_VDestroyVectorArray`
functions should be used instead.

The previously deprecated constructor `N_VMakeWithManagedAllocator_Cuda` and
the function `N_VSetCudaStream_Cuda` have been removed and replaced with
`N_VNewWithMemHelp_Cuda` and `N_VSetKerrnelExecPolicy_Cuda` respectively.

The previously deprecated macros `PVEC_REAL_MPI_TYPE` and
`PVEC_INTEGER_MPI_TYPE` have been removed and replaced with
`MPI_SUNREALTYPE` and `MPI_SUNINDEXTYPE` respectively.

### SUNLinearSolver

The following previously deprecated functions have been removed

| Removed                   | Replaced with                    |
|:--------------------------|:---------------------------------|
| `SUNBandLinearSolver`     | `SUNLinSol_Band`                 |
| `SUNDenseLinearSolver`    | `SUNLinSol_Dense`                |
| `SUNKLU`                  | `SUNLinSol_KLU`                  |
| `SUNKLUReInit`            | `SUNLinSol_KLUReInit`            |
| `SUNKLUSetOrdering`       | `SUNLinSol_KLUSetOrdering`       |
| `SUNLapackBand`           | `SUNLinSol_LapackBand`           |
| `SUNLapackDense`          | `SUNLinSol_LapackDense`          |
| `SUNPCG`                  | `SUNLinSol_PCG`                  |
| `SUNPCGSetPrecType`       | `SUNLinSol_PCGSetPrecType`       |
| `SUNPCGSetMaxl`           | `SUNLinSol_PCGSetMaxl`           |
| `SUNSPBCGS`               | `SUNLinSol_SPBCGS`               |
| `SUNSPBCGSSetPrecType`    | `SUNLinSol_SPBCGSSetPrecType`    |
| `SUNSPBCGSSetMaxl`        | `SUNLinSol_SPBCGSSetMaxl`        |
| `SUNSPFGMR`               | `SUNLinSol_SPFGMR`               |
| `SUNSPFGMRSetPrecType`    | `SUNLinSol_SPFGMRSetPrecType`    |
| `SUNSPFGMRSetGSType`      | `SUNLinSol_SPFGMRSetGSType`      |
| `SUNSPFGMRSetMaxRestarts` | `SUNLinSol_SPFGMRSetMaxRestarts` |
| `SUNSPGMR`                | `SUNLinSol_SPGMR`                |
| `SUNSPGMRSetPrecType`     | `SUNLinSol_SPGMRSetPrecType`     |
| `SUNSPGMRSetGSType`       | `SUNLinSol_SPGMRSetGSType`       |
| `SUNSPGMRSetMaxRestarts`  | `SUNLinSol_SPGMRSetMaxRestarts`  |
| `SUNSPTFQMR`              | `SUNLinSol_SPTFQMR`              |
| `SUNSPTFQMRSetPrecType`   | `SUNLinSol_SPTFQMRSetPrecType`   |
| `SUNSPTFQMRSetMaxl`       | `SUNLinSol_SPTFQMRSetMaxl`       |
| `SUNSuperLUMT`            | `SUNLinSol_SuperLUMT`            |
| `SUNSuperLUMTSetOrdering` | `SUNLinSol_SuperLUMTSetOrdering` |

### Fortran Interfaces

The ARKODE, CVODE, IDA, and KINSOL Fortran 77 interfaces have been removed. See
the "SUNDIALS Fortran Interface" section in the user guides and the F2003
example programs for more details using the SUNDIALS Fortran 2003 module
interfaces.

### ARKODE

The ARKODE MRIStep module has been extended to support implicit-explicit (IMEX)
multirate infinitesimal generalized additive Runge-Kutta (MRI-GARK) methods. As
such, `MRIStepCreate` has been updated to include arguments for the slow
explicit and slow implicit ODE right-hand side functions. `MRIStepCreate` has
also been updated to require attaching an `MRIStepInnerStepper` for evolving the
fast time scale. `MRIStepReInit` has been similarly updated to take explicit
and implicit right-hand side functions as input. Codes using explicit or
implicit MRI methods will need to update `MRIStepCreate` and `MRIStepReInit`
calls to pass `NULL` for either the explicit or implicit right-hand side
function as appropriate. If ARKStep is used as the fast time scale integrator,
codes will need to call `ARKStepCreateMRIStepInnerStepper` to wrap the ARKStep
memory as an `MRIStepInnerStepper` object. Additionally, `MRIStepGetNumRhsEvals`
has been updated to return the number of slow implicit and explicit function
evaluations. The coupling table structure `MRIStepCouplingMem` and the
functions `MRIStepCoupling_Alloc` and `MRIStepCoupling_Create` have also
been updated to support IMEX-MRI-GARK methods.

The deprecated functions `MRIStepGetCurrentButcherTables` and
`MRIStepWriteButcher` and the utility functions `MRIStepSetTable` and
`MRIStepSetTableNum` have been removed. Users wishing to create an MRI-GARK
method from a Butcher table should use `MRIStepCoupling_MIStoMRI` to create
the corresponding MRI coupling table and attach it with `MRIStepSetCoupling`.

The implementation of solve-decoupled implicit MRI-GARK methods has been updated
to remove extraneous slow implicit function calls and reduce the memory
requirements.

Deprecated ARKODE nonlinear solver predictors: specification of the ARKStep
"bootstrap" or "minimum correction" predictors (options 4 and 5 from
`ARKStepSetPredictorMethod`), or MRIStep "bootstrap" predictor (option 4 from
`MRIStepSetPredictorMethod`), will output a deprecation warning message.
These options will be removed in a future release.

The previously deprecated functions `ARKStepSetMaxStepsBetweenLSet` and
`ARKStepSetMaxStepsBetweenJac` have been removed and replaced with
`ARKStepSetLSetupFrequency` and `ARKStepSetMaxStepsBetweenJac` respectively.

### CVODE

The previously deprecated function `CVodeSetMaxStepsBetweenJac` has been removed
and replaced with `CVodeSetJacEvalFrequency`.

### CVODES

Added a new function `CVodeGetLinSolveStats` to get the CVODES linear solver
statistics as a group.

Added a new function, `CVodeSetMonitorFn`, that takes a user-function
to be called by CVODES after every `nst` successfully completed time-steps.
This is intended to provide a way of monitoring the CVODES statistics
throughout the simulation.

The previously deprecated function `CVodeSetMaxStepsBetweenJac` has been removed
and replaced with `CVodeSetJacEvalFrequency`.

### KINSOL

New orthogonalization methods were added for use within Anderson acceleration
in KINSOL. See the "Anderson Acceleration QR Factorization" subsection within
the mathematical considerations chapter of the user guide and the
`KINSetOrthAA` function documentation for more details.

### Deprecations

In addition to the deprecations noted elsewhere, many constants, types, and
functions have been renamed so that they are properly namespaced. The old names
have been deprecated and will be removed in SUNDIALS v7.0.0.

The following constants, macros, and  typedefs are now deprecated:

| Deprecated Name            | New Name                          |
|:---------------------------|:----------------------------------|
| `realtype`                 | `sunrealtype`                     |
| `booleantype`              | `sunbooleantype`                  |
| `RCONST`                   | `SUN_RCONST`                      |
| `BIG_REAL`                 | `SUN_BIG_REAL`                    |
| `SMALL_REAL`               | `SUN_SMALL_REAL`                  |
| `UNIT_ROUNDOFF`            | `SUN_UNIT_ROUNDOFF`               |
| `PREC_NONE`                | `SUN_PREC_NONE`                   |
| `PREC_LEFT`                | `SUN_PREC_LEFT`                   |
| `PREC_RIGHT`               | `SUN_PREC_RIGHT`                  |
| `PREC_BOTH`                | `SUN_PREC_BOTH`                   |
| `MODIFIED_GS`              | `SUN_MODIFIED_GS`                 |
| `CLASSICAL_GS`             | `SUN_CLASSICAL_GS`                |
| `ATimesFn`                 | `SUNATimesFn`                     |
| `PSetupFn`                 | `SUNPSetupFn`                     |
| `PSolveFn`                 | `SUNPSolveFn`                     |
| `DlsMat`                   | `SUNDlsMat`                       |
| `DENSE_COL`                | `SUNDLS_DENSE_COL`                |
| `DENSE_ELEM`               | `SUNDLS_DENSE_ELEM`               |
| `BAND_COL`                 | `SUNDLS_BAND_COL`                 |
| `BAND_COL_ELEM`            | `SUNDLS_BAND_COL_ELEM`            |
| `BAND_ELEM`                | `SUNDLS_BAND_ELEM`                |
| `SDIRK_2_1_2`              | `ARKODE_SDIRK_2_1_2`              |
| `BILLINGTON_3_3_2`         | `ARKODE_BILLINGTON_3_3_2`         |
| `TRBDF2_3_3_2`             | `ARKODE_TRBDF2_3_3_2`             |
| `KVAERNO_4_2_3`            | `ARKODE_KVAERNO_4_2_3`            |
| `ARK324L2SA_DIRK_4_2_3`    | `ARKODE_ARK324L2SA_DIRK_4_2_3`    |
| `CASH_5_2_4`               | `ARKODE_CASH_5_2_4`               |
| `CASH_5_3_4`               | `ARKODE_CASH_5_3_4`               |
| `SDIRK_5_3_4`              | `ARKODE_SDIRK_5_3_4`              |
| `KVAERNO_5_3_4`            | `ARKODE_KVAERNO_5_3_4`            |
| `ARK436L2SA_DIRK_6_3_4`    | `ARKODE_ARK436L2SA_DIRK_6_3_4`    |
| `KVAERNO_7_4_5`            | `ARKODE_KVAERNO_7_4_5`            |
| `ARK548L2SA_DIRK_8_4_5`    | `ARKODE_ARK548L2SA_DIRK_8_4_5`    |
| `ARK437L2SA_DIRK_7_3_4`    | `ARKODE_ARK437L2SA_DIRK_7_3_4`    |
| `ARK548L2SAb_DIRK_8_4_5`   | `ARKODE_ARK548L2SAb_DIRK_8_4_5`   |
| `MIN_DIRK_NUM`             | `ARKODE_MIN_DIRK_NUM`             |
| `MAX_DIRK_NUM`             | `ARKODE_MAX_DIRK_NUM`             |
| `MIS_KW3`                  | `ARKODE_MIS_KW3`                  |
| `MRI_GARK_ERK33a`          | `ARKODE_MRI_GARK_ERK33a`          |
| `MRI_GARK_ERK45a`          | `ARKODE_MRI_GARK_ERK45a`          |
| `MRI_GARK_IRK21a`          | `ARKODE_MRI_GARK_IRK21a`          |
| `MRI_GARK_ESDIRK34a`       | `ARKODE_MRI_GARK_ESDIRK34a`       |
| `MRI_GARK_ESDIRK46a`       | `ARKODE_MRI_GARK_ESDIRK46a`       |
| `IMEX_MRI_GARK3a`          | `ARKODE_IMEX_MRI_GARK3a`          |
| `IMEX_MRI_GARK3b`          | `ARKODE_IMEX_MRI_GARK3b`          |
| `IMEX_MRI_GARK4`           | `ARKODE_IMEX_MRI_GARK4`           |
| `MIN_MRI_NUM`              | `ARKODE_MIN_MRI_NUM`              |
| `MAX_MRI_NUM`              | `ARKODE_MAX_MRI_NUM`              |
| `DEFAULT_MRI_TABLE_3`      | `MRISTEP_DEFAULT_TABLE_3`         |
| `DEFAULT_EXPL_MRI_TABLE_3` | `MRISTEP_DEFAULT_EXPL_TABLE_3`    |
| `DEFAULT_EXPL_MRI_TABLE_4` | `MRISTEP_DEFAULT_EXPL_TABLE_4`    |
| `DEFAULT_IMPL_SD_TABLE_2`  | `MRISTEP_DEFAULT_IMPL_SD_TABLE_2` |
| `DEFAULT_IMPL_SD_TABLE_3`  | `MRISTEP_DEFAULT_IMPL_SD_TABLE_3` |
| `DEFAULT_IMPL_SD_TABLE_4`  | `MRISTEP_DEFAULT_IMPL_SD_TABLE_4` |
| `DEFAULT_IMEX_SD_TABLE_3`  | `MRISTEP_DEFAULT_IMEX_SD_TABLE_3` |
| `DEFAULT_IMEX_SD_TABLE_4`  | `MRISTEP_DEFAULT_IMEX_SD_TABLE_4` |
| `HEUN_EULER_2_1_2`         | `ARKODE_HEUN_EULER_2_1_2`         |
| `BOGACKI_SHAMPINE_4_2_3`   | `ARKODE_BOGACKI_SHAMPINE_4_2_3`   |
| `ARK324L2SA_ERK_4_2_3`     | `ARKODE_ARK324L2SA_ERK_4_2_3`     |
| `ZONNEVELD_5_3_4`          | `ARKODE_ZONNEVELD_5_3_4`          |
| `ARK436L2SA_ERK_6_3_4`     | `ARKODE_ARK436L2SA_ERK_6_3_4`     |
| `SAYFY_ABURUB_6_3_4`       | `ARKODE_SAYFY_ABURUB_6_3_4`       |
| `CASH_KARP_6_4_5`          | `ARKODE_CASH_KARP_6_4_5`          |
| `FEHLBERG_6_4_5`           | `ARKODE_FEHLBERG_6_4_5`           |
| `DORMAND_PRINCE_7_4_5`     | `ARKODE_DORMAND_PRINCE_7_4_5`     |
| `ARK548L2SA_ERK_8_4_5`     | `ARKODE_ARK548L2SA_ERK_8_4_5`     |
| `VERNER_8_5_6`             | `ARKODE_VERNER_8_5_6`             |
| `FEHLBERG_13_7_8`          | `ARKODE_FEHLBERG_13_7_8`          |
| `KNOTH_WOLKE_3_3`          | `ARKODE_KNOTH_WOLKE_3_3`          |
| `ARK437L2SA_ERK_7_3_4`     | `ARKODE_ARK437L2SA_ERK_7_3_4`     |
| `ARK548L2SAb_ERK_8_4_5`    | `ARKODE_ARK548L2SAb_ERK_8_4_5`    |
| `MIN_ERK_NUM`              | `ARKODE_MIN_ERK_NUM`              |
| `MAX_ERK_NUM`              | `ARKODE_MAX_ERK_NUM`              |
| `DEFAULT_ERK_2`            | `ARKSTEP_DEFAULT_ERK_2`           |
| `DEFAULT_ERK_3`            | `ARKSTEP_DEFAULT_ERK_3`           |
| `DEFAULT_ERK_4`            | `ARKSTEP_DEFAULT_ERK_4`           |
| `DEFAULT_ERK_5`            | `ARKSTEP_DEFAULT_ERK_5`           |
| `DEFAULT_ERK_6`            | `ARKSTEP_DEFAULT_ERK_6`           |
| `DEFAULT_ERK_8`            | `ARKSTEP_DEFAULT_ERK_8`           |
| `DEFAULT_DIRK_2`           | `ARKSTEP_DEFAULT_DIRK_2`          |
| `DEFAULT_DIRK_3`           | `ARKSTEP_DEFAULT_DIRK_3`          |
| `DEFAULT_DIRK_4`           | `ARKSTEP_DEFAULT_DIRK_4`          |
| `DEFAULT_DIRK_5`           | `ARKSTEP_DEFAULT_DIRK_5`          |
| `DEFAULT_ARK_ETABLE_3`     | `ARKSTEP_DEFAULT_ARK_ETABLE_3`    |
| `DEFAULT_ARK_ETABLE_4`     | `ARKSTEP_DEFAULT_ARK_ETABLE_4`    |
| `DEFAULT_ARK_ETABLE_5`     | `ARKSTEP_DEFAULT_ARK_ETABLE_4`    |
| `DEFAULT_ARK_ITABLE_3`     | `ARKSTEP_DEFAULT_ARK_ITABLE_3`    |
| `DEFAULT_ARK_ITABLE_4`     | `ARKSTEP_DEFAULT_ARK_ITABLE_4`    |
| `DEFAULT_ARK_ITABLE_5`     | `ARKSTEP_DEFAULT_ARK_ITABLE_5`    |
| `DEFAULT_ERK_2`            | `ERKSTEP_DEFAULT_2`               |
| `DEFAULT_ERK_3`            | `ERKSTEP_DEFAULT_3`               |
| `DEFAULT_ERK_4`            | `ERKSTEP_DEFAULT_4`               |
| `DEFAULT_ERK_5`            | `ERKSTEP_DEFAULT_5`               |
| `DEFAULT_ERK_6`            | `ERKSTEP_DEFAULT_6`               |
| `DEFAULT_ERK_8`            | `ERKSTEP_DEFAULT_8`               |

In addition, the following functions are now deprecated (compile-time warnings
will be thrown if supported by the compiler):

| Deprecated Name               | New Name                     |
|:------------------------------|:-----------------------------|
| `CVSpilsSetLinearSolver`      | `CVodeSetLinearSolver`       |
| `CVSpilsSetEpsLin`            | `CVodeSetEpsLin`             |
| `CVSpilsSetPreconditioner`    | `CVodeSetPreconditioner`     |
| `CVSpilsSetJacTimes`          | `CVodeSetJacTimes`           |
| `CVSpilsGetWorkSpace`         | `CVodeGetLinWorkSpace`       |
| `CVSpilsGetNumPrecEvals`      | `CVodeGetNumPrecEvals`       |
| `CVSpilsGetNumPrecSolves`     | `CVodeGetNumPrecSolves`      |
| `CVSpilsGetNumLinIters`       | `CVodeGetNumLinIters`        |
| `CVSpilsGetNumConvFails`      | `CVodeGetNumConvFails`       |
| `CVSpilsGetNumJTSetupEvals`   | `CVodeGetNumJTSetupEvals`    |
| `CVSpilsGetNumJtimesEvals`    | `CVodeGetNumJtimesEvals`     |
| `CVSpilsGetNumRhsEvals`       | `CVodeGetNumLinRhsEvals`     |
| `CVSpilsGetLastFlag`          | `CVodeGetLastLinFlag`        |
| `CVSpilsGetReturnFlagName`    | `CVodeGetLinReturnFlagName`  |
| `CVSpilsSetLinearSolverB`     | `CVodeSetLinearSolverB`      |
| `CVSpilsSetEpsLinB`           | `CVodeSetEpsLinB`            |
| `CVSpilsSetPreconditionerB`   | `CVodeSetPreconditionerB`    |
| `CVSpilsSetPreconditionerBS`  | `CVodeSetPreconditionerBS`   |
| `CVSpilsSetJacTimesB`         | `CVodeSetJacTimesB`          |
| `CVSpilsSetJacTimesBS`        | `CVodeSetJacTimesBS`         |
| `CVDlsSetLinearSolver`        | `CVodeSetLinearSolver`       |
| `CVDlsSetJacFn`               | `CVodeSetJacFn`              |
| `CVDlsGetWorkSpace`           | `CVodeGetLinWorkSpace`       |
| `CVDlsGetNumJacEvals`         | `CVodeGetNumJacEvals`        |
| `CVDlsGetNumRhsEvals`         | `CVodeGetNumLinRhsEvals`     |
| `CVDlsGetLastFlag`            | `CVodeGetLastLinFlag`        |
| `CVDlsGetReturnFlagName`      | `CVodeGetLinReturnFlagName`  |
| `CVDlsSetLinearSolverB`       | `CVodeSetLinearSolverB`      |
| `CVDlsSetJacFnB`              | `CVodeSetJacFnB`             |
| `CVDlsSetJacFnBS`             | `CVodeSetJacFnBS`            |
| `CVDlsSetLinearSolver`        | `CVodeSetLinearSolver`       |
| `CVDlsSetJacFn`               | `CVodeSetJacFn`              |
| `CVDlsGetWorkSpace`           | `CVodeGetLinWorkSpace`       |
| `CVDlsGetNumJacEvals`         | `CVodeGetNumJacEvals`        |
| `CVDlsGetNumRhsEvals`         | `CVodeGetNumLinRhsEvals`     |
| `CVDlsGetLastFlag`            | `CVodeGetLastLinFlag`        |
| `CVDlsGetReturnFlagName`      | `CVodeGetLinReturnFlagName`  |
| `KINDlsSetLinearSolver`       | `KINSetLinearSolver`         |
| `KINDlsSetJacFn`              | `KINSetJacFn`                |
| `KINDlsGetWorkSpace`          | `KINGetLinWorkSpace`         |
| `KINDlsGetNumJacEvals`        | `KINGetNumJacEvals`          |
| `KINDlsGetNumFuncEvals`       | `KINGetNumLinFuncEvals`      |
| `KINDlsGetLastFlag`           | `KINGetLastLinFlag`          |
| `KINDlsGetReturnFlagName`     | `KINGetLinReturnFlagName`    |
| `KINSpilsSetLinearSolver`     | `KINSetLinearSolver`         |
| `KINSpilsSetPreconditioner`   | `KINSetPreconditioner`       |
| `KINSpilsSetJacTimesVecFn`    | `KINSetJacTimesVecFn`        |
| `KINSpilsGetWorkSpace`        | `KINGetLinWorkSpace`         |
| `KINSpilsGetNumPrecEvals`     | `KINGetNumPrecEvals`         |
| `KINSpilsGetNumPrecSolves`    | `KINGetNumPrecSolves`        |
| `KINSpilsGetNumLinIters`      | `KINGetNumLinIters`          |
| `KINSpilsGetNumConvFails`     | `KINGetNumLinConvFails`      |
| `KINSpilsGetNumJtimesEvals`   | `KINGetNumJtimesEvals`       |
| `KINSpilsGetNumFuncEvals`     | `KINGetNumLinFuncEvals`      |
| `KINSpilsGetLastFlag`         | `KINGetLastLinFlag`          |
| `KINSpilsGetReturnFlagName`   | `KINGetLinReturnFlagName`    |
| `IDASpilsSetLinearSolver`     | `IDASetLinearSolver`         |
| `IDASpilsSetPreconditioner`   | `IDASetPreconditioner`       |
| `IDASpilsSetJacTimes`         | `IDASetJacTimes`             |
| `IDASpilsSetEpsLin`           | `IDASetEpsLin`               |
| `IDASpilsSetIncrementFactor`  | `IDASetIncrementFactor`      |
| `IDASpilsGetWorkSpace`        | `IDAGetLinWorkSpace`         |
| `IDASpilsGetNumPrecEvals`     | `IDAGetNumPrecEvals`         |
| `IDASpilsGetNumPrecSolves`    | `IDAGetNumPrecSolves`        |
| `IDASpilsGetNumLinIters`      | `IDAGetNumLinIters`          |
| `IDASpilsGetNumConvFails`     | `IDAGetNumLinConvFails`      |
| `IDASpilsGetNumJTSetupEvals`  | `IDAGetNumJTSetupEvals`      |
| `IDASpilsGetNumJtimesEvals`   | `IDAGetNumJtimesEvals`       |
| `IDASpilsGetNumResEvals`      | `IDAGetNumLinResEvals`       |
| `IDASpilsGetLastFlag`         | `IDAGetLastLinFlag`          |
| `IDASpilsGetReturnFlagName`   | `IDAGetLinReturnFlagName`    |
| `IDASpilsSetLinearSolverB`    | `IDASetLinearSolverB`        |
| `IDASpilsSetEpsLinB`          | `IDASetEpsLinB`              |
| `IDASpilsSetIncrementFactorB` | `IDASetIncrementFactorB`     |
| `IDASpilsSetPreconditionerB`  | `IDASetPreconditionerB`      |
| `IDASpilsSetPreconditionerBS` | `IDASetPreconditionerBS`     |
| `IDASpilsSetJacTimesB`        | `IDASetJacTimesB`            |
| `IDASpilsSetJacTimesBS`       | `IDASetJacTimesBS`           |
| `IDADlsSetLinearSolver`       | `IDASetLinearSolver`         |
| `IDADlsSetJacFn`              | `IDASetJacFn`                |
| `IDADlsGetWorkSpace`          | `IDAGetLinWorkSpace`         |
| `IDADlsGetNumJacEvals`        | `IDAGetNumJacEvals`          |
| `IDADlsGetNumResEvals`        | `IDAGetNumLinResEvals`       |
| `IDADlsGetLastFlag`           | `IDAGetLastLinFlag`          |
| `IDADlsGetReturnFlagName`     | `IDAGetLinReturnFlagName`    |
| `IDADlsSetLinearSolverB`      | `IDASetLinearSolverB`        |
| `IDADlsSetJacFnB`             | `IDASetJacFnB`               |
| `IDADlsSetJacFnBS`            | `IDASetJacFnBS`              |
| `DenseGETRF`                  | `SUNDlsMat_DenseGETRF`       |
| `DenseGETRS`                  | `SUNDlsMat_DenseGETRS`       |
| `denseGETRF`                  | `SUNDlsMat_denseGETRF`       |
| `denseGETRS`                  | `SUNDlsMat_denseGETRS`       |
| `DensePOTRF`                  | `SUNDlsMat_DensePOTRF`       |
| `DensePOTRS`                  | `SUNDlsMat_DensePOTRS`       |
| `densePOTRF`                  | `SUNDlsMat_densePOTRF`       |
| `densePOTRS`                  | `SUNDlsMat_densePOTRS`       |
| `DenseGEQRF`                  | `SUNDlsMat_DenseGEQRF`       |
| `DenseORMQR`                  | `SUNDlsMat_DenseORMQR`       |
| `denseGEQRF`                  | `SUNDlsMat_denseGEQRF`       |
| `denseORMQR`                  | `SUNDlsMat_denseORMQR`       |
| `DenseCopy`                   | `SUNDlsMat_DenseCopy`        |
| `denseCopy`                   | `SUNDlsMat_denseCopy`        |
| `DenseScale`                  | `SUNDlsMat_DenseScale`       |
| `denseScale`                  | `SUNDlsMat_denseScale`       |
| `denseAddIdentity`            | `SUNDlsMat_denseAddIdentity` |
| `DenseMatvec`                 | `SUNDlsMat_DenseMatvec`      |
| `denseMatvec`                 | `SUNDlsMat_denseMatvec`      |
| `BandGBTRF`                   | `SUNDlsMat_BandGBTRF`        |
| `bandGBTRF`                   | `SUNDlsMat_bandGBTRF`        |
| `BandGBTRS`                   | `SUNDlsMat_BandGBTRS`        |
| `bandGBTRS`                   | `SUNDlsMat_bandGBTRS`        |
| `BandCopy`                    | `SUNDlsMat_BandCopy`         |
| `bandCopy`                    | `SUNDlsMat_bandCopy`         |
| `BandScale`                   | `SUNDlsMat_BandScale`        |
| `bandScale`                   | `SUNDlsMat_bandScale`        |
| `bandAddIdentity`             | `SUNDlsMat_bandAddIdentity`  |
| `BandMatvec`                  | `SUNDlsMat_BandMatvec`       |
| `bandMatvec`                  | `SUNDlsMat_bandMatvec`       |
| `ModifiedGS`                  | `SUNModifiedGS`              |
| `ClassicalGS`                 | `SUNClassicalGS`             |
| `QRfact`                      | `SUNQRFact`                  |
| `QRsol`                       | `SUNQRsol`                   |
| `DlsMat_NewDenseMat`          | `SUNDlsMat_NewDenseMat`      |
| `DlsMat_NewBandMat`           | `SUNDlsMat_NewBandMat`       |
| `DestroyMat`                  | `SUNDlsMat_DestroyMat`       |
| `NewIntArray`                 | `SUNDlsMat_NewIntArray`      |
| `NewIndexArray`               | `SUNDlsMat_NewIndexArray`    |
| `NewRealArray`                | `SUNDlsMat_NewRealArray`     |
| `DestroyArray`                | `SUNDlsMat_DestroyArray`     |
| `AddIdentity`                 | `SUNDlsMat_AddIdentity`      |
| `SetToZero`                   | `SUNDlsMat_SetToZero`        |
| `PrintMat`                    | `SUNDlsMat_PrintMat`         |
| `newDenseMat`                 | `SUNDlsMat_newDenseMat`      |
| `newBandMat`                  | `SUNDlsMat_newBandMat`       |
| `destroyMat`                  | `SUNDlsMat_destroyMat`       |
| `newIntArray`                 | `SUNDlsMat_newIntArray`      |
| `newIndexArray`               | `SUNDlsMat_newIndexArray`    |
| `newRealArray`                | `SUNDlsMat_newRealArray`     |
| `destroyArray`                | `SUNDlsMat_destroyArray`     |

In addition, the entire `sundials_lapack.h` header file is now deprecated for
removal in SUNDIALS v7.0.0. Note, this header file is not needed to use the
SUNDIALS LAPACK linear solvers.

### ARKODE Changes in v5.0.0

**SUNContext**

SUNDIALS v6.0.0 introduces a new :c:type:`SUNContext` object on which all other
SUNDIALS objects depend. As such, the constructors for all SUNDIALS packages,
vectors, matrices, linear solvers, nonlinear solvers, and memory helpers have
been updated to accept a context as the last input. Users upgrading to SUNDIALS
v6.0.0 will need to call :c:func:`SUNContext_Create` to create a context object
with before calling any other SUNDIALS library function, and then provide this
object to other SUNDIALS constructors. The context object has been introduced to
allow SUNDIALS to provide new features, such as the profiling/instrumentation
also introduced in this release, while maintaining thread-safety. See the
documentation section on the :c:type:`SUNContext` for more details.

A script ``upgrade-to-sundials-6-from-5.sh`` has been provided with the release
(obtainable from the GitHub release page) to help ease the transition to
SUNDIALS v6.0.0. The script will add a ``SUNCTX_PLACEHOLDER`` argument to all of
the calls to SUNDIALS constructors that now require a ``SUNContext`` object. It
can also update deprecated SUNDIALS constants/types to the new names. It can be
run like this:

.. code-block::

   > ./upgrade-to-sundials-6-from-5.sh <files to update>

**SUNProfiler**

A capability to profile/instrument SUNDIALS library code has been added. This
can be enabled with the CMake option :cmakeop:`SUNDIALS_BUILD_WITH_PROFILING`. A
built-in profiler will be used by default, but the `Caliper
<https://github.com/LLNL/Caliper>`_ library can also be used instead with the
CMake option :cmakeop:`ENABLE_CALIPER`. See the documentation section on
profiling for more details.  **WARNING**: Profiling will impact performance, and
should be enabled judiciously.

**SUNMemoryHelper**

The :c:type:`SUNMemoryHelper` functions :c:func:`SUNMemoryHelper_Alloc`,
:c:func:`SUNMemoryHelper_Dealloc`, and :c:func:`SUNMemoryHelper_Copy` have been
updated to accept an opaque handle as the last input. At a minimum, user-defined
:c:type:`SUNMemoryHelper` implementations will need to update these functions to
accept the additional argument. Typically, this handle is the execution stream
(e.g., a CUDA/HIP stream or SYCL queue) for the operation. The :ref:`CUDA
<SUNMemory.CUDA>`, :ref:`HIP <SUNMemory.HIP>`, and :ref:`SYCL <SUNMemory.SYCL>`
implementations have been updated accordingly. Additionally, the constructor
:c:func:`SUNMemoryHelper_Sycl` has been updated to remove the SYCL queue as an
input.

**NVector**

Two new optional vector operations, :c:func:`N_VDotProdMultiLocal` and
:c:func:`N_VDotProdMultiAllReduce`, have been added to support
low-synchronization methods for Anderson acceleration.

The CUDA, HIP, and SYCL execution policies have been moved from the ``sundials``
namespace to the ``sundials::cuda``, ``sundials::hip``, and ``sundials::sycl``
namespaces respectively. Accordingly, the prefixes "Cuda", "Hip", and "Sycl"
have been removed from the execution policy classes and methods.

The ``Sundials`` namespace used by the Trilinos Tpetra NVector has been replaced
with the ``sundials::trilinos::nvector_tpetra`` namespace.

The serial, PThreads, PETSc, *hypre*, Parallel, OpenMP_DEV, and OpenMP vector
functions ``N_VCloneVectorArray_*`` and ``N_VDestroyVectorArray_*`` have been
deprecated. The generic :c:func:`N_VCloneVectorArray` and
:c:func:`N_VDestroyVectorArray` functions should be used instead.

The previously deprecated constructor ``N_VMakeWithManagedAllocator_Cuda`` and
the function ``N_VSetCudaStream_Cuda`` have been removed and replaced with
:c:func:`N_VNewWithMemHelp_Cuda` and :c:func:`N_VSetKerrnelExecPolicy_Cuda`
respectively.

The previously deprecated macros ``PVEC_REAL_MPI_TYPE`` and
``PVEC_INTEGER_MPI_TYPE`` have been removed and replaced with
``MPI_SUNREALTYPE`` and ``MPI_SUNINDEXTYPE`` respectively.

**SUNLinearSolver**

The following previously deprecated functions have been removed:

+-----------------------------+------------------------------------------+
| Removed                     | Replacement                              |
+=============================+==========================================+
| ``SUNBandLinearSolver``     | :c:func:`SUNLinSol_Band`                 |
+-----------------------------+------------------------------------------+
| ``SUNDenseLinearSolver``    | :c:func:`SUNLinSol_Dense`                |
+-----------------------------+------------------------------------------+
| ``SUNKLU``                  | :c:func:`SUNLinSol_KLU`                  |
+-----------------------------+------------------------------------------+
| ``SUNKLUReInit``            | :c:func:`SUNLinSol_KLUReInit`            |
+-----------------------------+------------------------------------------+
| ``SUNKLUSetOrdering``       | :c:func:`SUNLinSol_KLUSetOrdering`       |
+-----------------------------+------------------------------------------+
| ``SUNLapackBand``           | :c:func:`SUNLinSol_LapackBand`           |
+-----------------------------+------------------------------------------+
| ``SUNLapackDense``          | :c:func:`SUNLinSol_LapackDense`          |
+-----------------------------+------------------------------------------+
| ``SUNPCG``                  | :c:func:`SUNLinSol_PCG`                  |
+-----------------------------+------------------------------------------+
| ``SUNPCGSetPrecType``       | :c:func:`SUNLinSol_PCGSetPrecType`       |
+-----------------------------+------------------------------------------+
| ``SUNPCGSetMaxl``           | :c:func:`SUNLinSol_PCGSetMaxl`           |
+-----------------------------+------------------------------------------+
| ``SUNSPBCGS``               | :c:func:`SUNLinSol_SPBCGS`               |
+-----------------------------+------------------------------------------+
| ``SUNSPBCGSSetPrecType``    | :c:func:`SUNLinSol_SPBCGSSetPrecType`    |
+-----------------------------+------------------------------------------+
| ``SUNSPBCGSSetMaxl``        | :c:func:`SUNLinSol_SPBCGSSetMaxl`        |
+-----------------------------+------------------------------------------+
| ``SUNSPFGMR``               | :c:func:`SUNLinSol_SPFGMR`               |
+-----------------------------+------------------------------------------+
| ``SUNSPFGMRSetPrecType``    | :c:func:`SUNLinSol_SPFGMRSetPrecType`    |
+-----------------------------+------------------------------------------+
| ``SUNSPFGMRSetGSType``      | :c:func:`SUNLinSol_SPFGMRSetGSType`      |
+-----------------------------+------------------------------------------+
| ``SUNSPFGMRSetMaxRestarts`` | :c:func:`SUNLinSol_SPFGMRSetMaxRestarts` |
+-----------------------------+------------------------------------------+
| ``SUNSPGMR``                | :c:func:`SUNLinSol_SPGMR`                |
+-----------------------------+------------------------------------------+
| ``SUNSPGMRSetPrecType``     | :c:func:`SUNLinSol_SPGMRSetPrecType`     |
+-----------------------------+------------------------------------------+
| ``SUNSPGMRSetGSType``       | :c:func:`SUNLinSol_SPGMRSetGSType`       |
+-----------------------------+------------------------------------------+
| ``SUNSPGMRSetMaxRestarts``  | :c:func:`SUNLinSol_SPGMRSetMaxRestarts`  |
+-----------------------------+------------------------------------------+
| ``SUNSPTFQMR``              | :c:func:`SUNLinSol_SPTFQMR`              |
+-----------------------------+------------------------------------------+
| ``SUNSPTFQMRSetPrecType``   | :c:func:`SUNLinSol_SPTFQMRSetPrecType`   |
+-----------------------------+------------------------------------------+
| ``SUNSPTFQMRSetMaxl``       | :c:func:`SUNLinSol_SPTFQMRSetMaxl`       |
+-----------------------------+------------------------------------------+
| ``SUNSuperLUMT``            | :c:func:`SUNLinSol_SuperLUMT`            |
+-----------------------------+------------------------------------------+
| ``SUNSuperLUMTSetOrdering`` | :c:func:`SUNLinSol_SuperLUMTSetOrdering` |
+-----------------------------+------------------------------------------+

**ARKODE**

The MRIStep module has been extended to support implicit-explicit (ImEx)
multirate infinitesimal generalized additive Runge--Kutta (MRI-GARK) methods. As
such, :c:func:`MRIStepCreate` has been updated to include arguments for the slow
explicit and slow implicit ODE right-hand side functions.
:c:func:`MRIStepCreate` has also been updated to require attaching an
MRIStepInnerStepper for evolving the fast time scale. :c:func:`MRIStepReInit`
has been similarly updated to take explicit and implicit right-hand side
functions as input. Codes using explicit or implicit MRI methods will need to
update :c:func:`MRIStepCreate` and :c:func:`MRIStepReInit` calls to pass
``NULL`` for either the explicit or implicit right-hand side function as
appropriate. If ARKStep is used as the fast time scale integrator, codes will
need to call :c:func:`ARKStepCreateMRIStepInnerStepper` to wrap the ARKStep
memory as an MRIStepInnerStepper object. Additionally,
:c:func:`MRIStepGetNumRhsEvals` has been updated to return the number of slow
implicit and explicit function evaluations. The coupling table structure
:c:type:`MRIStepCouplingMem` and the functions :c:func:`MRIStepCoupling_Alloc`
and :c:func:`MRIStepCoupling_Create` have also been updated to support
IMEX-MRI-GARK methods.

The deprecated functions ``MRIStepGetCurrentButcherTables`` and
``MRIStepWriteButcher`` and the utility functions ``MRIStepSetTable`` and
``MRIStepSetTableNum`` have been removed. Users wishing to create an MRI-GARK
method from a Butcher table should use :c:func:`MRIStepCoupling_MIStoMRI` to
create the corresponding MRI coupling table and attach it with
:c:func:`MRIStepSetCoupling`.

The implementation of solve-decoupled implicit MRI-GARK methods has been updated
to remove extraneous slow implicit function calls and reduce the memory
requirements.

The previously deprecated functions ``ARKStepSetMaxStepsBetweenLSet`` and
``ARKStepSetMaxStepsBetweenJac`` have been removed and replaced with
:c:func:`ARKStepSetLSetupFrequency` and :c:func:`ARKStepSetMaxStepsBetweenJac`
respectively.

The ARKODE Fortran 77 interface has been removed. See :numref:`SUNDIALS.Fortran`
and the F2003 example programs for more details using the SUNDIALS Fortran 2003
module interfaces.

**Deprecations**

In addition to the deprecations noted elsewhere, many constants, types, and
functions have been renamed so that they are properly namespaced. The old names
have been deprecated and will be removed in SUNDIALS v7.0.0.

The following constants, macros, and typedefs are now deprecated:

+------------------------------+-------------------------------------+
| Deprecated Name              | New Name                            |
+==============================+=====================================+
| ``realtype``                 | ``sunrealtype``                     |
+------------------------------+-------------------------------------+
| ``booleantype``              | ``sunbooleantype``                  |
+------------------------------+-------------------------------------+
| ``RCONST``                   | ``SUN_RCONST``                      |
+------------------------------+-------------------------------------+
| ``BIG_REAL``                 | ``SUN_BIG_REAL``                    |
+------------------------------+-------------------------------------+
| ``SMALL_REAL``               | ``SUN_SMALL_REAL``                  |
+------------------------------+-------------------------------------+
| ``UNIT_ROUNDOFF``            | ``SUN_UNIT_ROUNDOFF``               |
+------------------------------+-------------------------------------+
| ``PREC_NONE``                | ``SUN_PREC_NONE``                   |
+------------------------------+-------------------------------------+
| ``PREC_LEFT``                | ``SUN_PREC_LEFT``                   |
+------------------------------+-------------------------------------+
| ``PREC_RIGHT``               | ``SUN_PREC_RIGHT``                  |
+------------------------------+-------------------------------------+
| ``PREC_BOTH``                | ``SUN_PREC_BOTH``                   |
+------------------------------+-------------------------------------+
| ``MODIFIED_GS``              | ``SUN_MODIFIED_GS``                 |
+------------------------------+-------------------------------------+
| ``CLASSICAL_GS``             | ``SUN_CLASSICAL_GS``                |
+------------------------------+-------------------------------------+
| ``ATimesFn``                 | ``SUNATimesFn``                     |
+------------------------------+-------------------------------------+
| ``PSetupFn``                 | ``SUNPSetupFn``                     |
+------------------------------+-------------------------------------+
| ``PSolveFn``                 | ``SUNPSolveFn``                     |
+------------------------------+-------------------------------------+
| ``DlsMat``                   | ``SUNDlsMat``                       |
+------------------------------+-------------------------------------+
| ``DENSE_COL``                | ``SUNDLS_DENSE_COL``                |
+------------------------------+-------------------------------------+
| ``DENSE_ELEM``               | ``SUNDLS_DENSE_ELEM``               |
+------------------------------+-------------------------------------+
| ``BAND_COL``                 | ``SUNDLS_BAND_COL``                 |
+------------------------------+-------------------------------------+
| ``BAND_COL_ELEM``            | ``SUNDLS_BAND_COL_ELEM``            |
+------------------------------+-------------------------------------+
| ``BAND_ELEM``                | ``SUNDLS_BAND_ELEM``                |
+------------------------------+-------------------------------------+
| ``SDIRK_2_1_2``              | ``ARKODE_SDIRK_2_1_2``              |
+------------------------------+-------------------------------------+
| ``BILLINGTON_3_3_2``         | ``ARKODE_BILLINGTON_3_3_2``         |
+------------------------------+-------------------------------------+
| ``TRBDF2_3_3_2``             | ``ARKODE_TRBDF2_3_3_2``             |
+------------------------------+-------------------------------------+
| ``KVAERNO_4_2_3``            | ``ARKODE_KVAERNO_4_2_3``            |
+------------------------------+-------------------------------------+
| ``ARK324L2SA_DIRK_4_2_3``    | ``ARKODE_ARK324L2SA_DIRK_4_2_3``    |
+------------------------------+-------------------------------------+
| ``CASH_5_2_4``               | ``ARKODE_CASH_5_2_4``               |
+------------------------------+-------------------------------------+
| ``CASH_5_3_4``               | ``ARKODE_CASH_5_3_4``               |
+------------------------------+-------------------------------------+
| ``SDIRK_5_3_4``              | ``ARKODE_SDIRK_5_3_4``              |
+------------------------------+-------------------------------------+
| ``KVAERNO_5_3_4``            | ``ARKODE_KVAERNO_5_3_4``            |
+------------------------------+-------------------------------------+
| ``ARK436L2SA_DIRK_6_3_4``    | ``ARKODE_ARK436L2SA_DIRK_6_3_4``    |
+------------------------------+-------------------------------------+
| ``KVAERNO_7_4_5``            | ``ARKODE_KVAERNO_7_4_5``            |
+------------------------------+-------------------------------------+
| ``ARK548L2SA_DIRK_8_4_5``    | ``ARKODE_ARK548L2SA_DIRK_8_4_5``    |
+------------------------------+-------------------------------------+
| ``ARK437L2SA_DIRK_7_3_4``    | ``ARKODE_ARK437L2SA_DIRK_7_3_4``    |
+------------------------------+-------------------------------------+
| ``ARK548L2SAb_DIRK_8_4_5``   | ``ARKODE_ARK548L2SAb_DIRK_8_4_5``   |
+------------------------------+-------------------------------------+
| ``MIN_DIRK_NUM``             | ``ARKODE_MIN_DIRK_NUM``             |
+------------------------------+-------------------------------------+
| ``MAX_DIRK_NUM``             | ``ARKODE_MAX_DIRK_NUM``             |
+------------------------------+-------------------------------------+
| ``MIS_KW3``                  | ``ARKODE_MIS_KW3``                  |
+------------------------------+-------------------------------------+
| ``MRI_GARK_ERK33a``          | ``ARKODE_MRI_GARK_ERK33a``          |
+------------------------------+-------------------------------------+
| ``MRI_GARK_ERK45a``          | ``ARKODE_MRI_GARK_ERK45a``          |
+------------------------------+-------------------------------------+
| ``MRI_GARK_IRK21a``          | ``ARKODE_MRI_GARK_IRK21a``          |
+------------------------------+-------------------------------------+
| ``MRI_GARK_ESDIRK34a``       | ``ARKODE_MRI_GARK_ESDIRK34a``       |
+------------------------------+-------------------------------------+
| ``MRI_GARK_ESDIRK46a``       | ``ARKODE_MRI_GARK_ESDIRK46a``       |
+------------------------------+-------------------------------------+
| ``IMEX_MRI_GARK3a``          | ``ARKODE_IMEX_MRI_GARK3a``          |
+------------------------------+-------------------------------------+
| ``IMEX_MRI_GARK3b``          | ``ARKODE_IMEX_MRI_GARK3b``          |
+------------------------------+-------------------------------------+
| ``IMEX_MRI_GARK4``           | ``ARKODE_IMEX_MRI_GARK4``           |
+------------------------------+-------------------------------------+
| ``MIN_MRI_NUM``              | ``ARKODE_MIN_MRI_NUM``              |
+------------------------------+-------------------------------------+
| ``MAX_MRI_NUM``              | ``ARKODE_MAX_MRI_NUM``              |
+------------------------------+-------------------------------------+
| ``DEFAULT_MRI_TABLE_3``      | ``MRISTEP_DEFAULT_TABLE_3``         |
+------------------------------+-------------------------------------+
| ``DEFAULT_EXPL_MRI_TABLE_3`` | ``MRISTEP_DEFAULT_EXPL_TABLE_3``    |
+------------------------------+-------------------------------------+
| ``DEFAULT_EXPL_MRI_TABLE_4`` | ``MRISTEP_DEFAULT_EXPL_TABLE_4``    |
+------------------------------+-------------------------------------+
| ``DEFAULT_IMPL_SD_TABLE_2``  | ``MRISTEP_DEFAULT_IMPL_SD_TABLE_2`` |
+------------------------------+-------------------------------------+
| ``DEFAULT_IMPL_SD_TABLE_3``  | ``MRISTEP_DEFAULT_IMPL_SD_TABLE_3`` |
+------------------------------+-------------------------------------+
| ``DEFAULT_IMPL_SD_TABLE_4``  | ``MRISTEP_DEFAULT_IMPL_SD_TABLE_4`` |
+------------------------------+-------------------------------------+
| ``DEFAULT_IMEX_SD_TABLE_3``  | ``MRISTEP_DEFAULT_IMEX_SD_TABLE_3`` |
+------------------------------+-------------------------------------+
| ``DEFAULT_IMEX_SD_TABLE_4``  | ``MRISTEP_DEFAULT_IMEX_SD_TABLE_4`` |
+------------------------------+-------------------------------------+
| ``HEUN_EULER_2_1_2``         | ``ARKODE_HEUN_EULER_2_1_2``         |
+------------------------------+-------------------------------------+
| ``BOGACKI_SHAMPINE_4_2_3``   | ``ARKODE_BOGACKI_SHAMPINE_4_2_3``   |
+------------------------------+-------------------------------------+
| ``ARK324L2SA_ERK_4_2_3``     | ``ARKODE_ARK324L2SA_ERK_4_2_3``     |
+------------------------------+-------------------------------------+
| ``ZONNEVELD_5_3_4``          | ``ARKODE_ZONNEVELD_5_3_4``          |
+------------------------------+-------------------------------------+
| ``ARK436L2SA_ERK_6_3_4``     | ``ARKODE_ARK436L2SA_ERK_6_3_4``     |
+------------------------------+-------------------------------------+
| ``SAYFY_ABURUB_6_3_4``       | ``ARKODE_SAYFY_ABURUB_6_3_4``       |
+------------------------------+-------------------------------------+
| ``CASH_KARP_6_4_5``          | ``ARKODE_CASH_KARP_6_4_5``          |
+------------------------------+-------------------------------------+
| ``FEHLBERG_6_4_5``           | ``ARKODE_FEHLBERG_6_4_5``           |
+------------------------------+-------------------------------------+
| ``DORMAND_PRINCE_7_4_5``     | ``ARKODE_DORMAND_PRINCE_7_4_5``     |
+------------------------------+-------------------------------------+
| ``ARK548L2SA_ERK_8_4_5``     | ``ARKODE_ARK548L2SA_ERK_8_4_5``     |
+------------------------------+-------------------------------------+
| ``VERNER_8_5_6``             | ``ARKODE_VERNER_8_5_6``             |
+------------------------------+-------------------------------------+
| ``FEHLBERG_13_7_8``          | ``ARKODE_FEHLBERG_13_7_8``          |
+------------------------------+-------------------------------------+
| ``KNOTH_WOLKE_3_3``          | ``ARKODE_KNOTH_WOLKE_3_3``          |
+------------------------------+-------------------------------------+
| ``ARK437L2SA_ERK_7_3_4``     | ``ARKODE_ARK437L2SA_ERK_7_3_4``     |
+------------------------------+-------------------------------------+
| ``ARK548L2SAb_ERK_8_4_5``    | ``ARKODE_ARK548L2SAb_ERK_8_4_5``    |
+------------------------------+-------------------------------------+
| ``MIN_ERK_NUM``              | ``ARKODE_MIN_ERK_NUM``              |
+------------------------------+-------------------------------------+
| ``MAX_ERK_NUM``              | ``ARKODE_MAX_ERK_NUM``              |
+------------------------------+-------------------------------------+
| ``DEFAULT_ERK_2``            | ``ARKSTEP_DEFAULT_ERK_2``           |
+------------------------------+-------------------------------------+
| ``DEFAULT_ERK_3``            | ``ARKSTEP_DEFAULT_ERK_3``           |
+------------------------------+-------------------------------------+
| ``DEFAULT_ERK_4``            | ``ARKSTEP_DEFAULT_ERK_4``           |
+------------------------------+-------------------------------------+
| ``DEFAULT_ERK_5``            | ``ARKSTEP_DEFAULT_ERK_5``           |
+------------------------------+-------------------------------------+
| ``DEFAULT_ERK_6``            | ``ARKSTEP_DEFAULT_ERK_6``           |
+------------------------------+-------------------------------------+
| ``DEFAULT_ERK_8``            | ``ARKSTEP_DEFAULT_ERK_8``           |
+------------------------------+-------------------------------------+
| ``DEFAULT_DIRK_2``           | ``ARKSTEP_DEFAULT_DIRK_2``          |
+------------------------------+-------------------------------------+
| ``DEFAULT_DIRK_3``           | ``ARKSTEP_DEFAULT_DIRK_3``          |
+------------------------------+-------------------------------------+
| ``DEFAULT_DIRK_4``           | ``ARKSTEP_DEFAULT_DIRK_4``          |
+------------------------------+-------------------------------------+
| ``DEFAULT_DIRK_5``           | ``ARKSTEP_DEFAULT_DIRK_5``          |
+------------------------------+-------------------------------------+
| ``DEFAULT_ARK_ETABLE_3``     | ``ARKSTEP_DEFAULT_ARK_ETABLE_3``    |
+------------------------------+-------------------------------------+
| ``DEFAULT_ARK_ETABLE_4``     | ``ARKSTEP_DEFAULT_ARK_ETABLE_4``    |
+------------------------------+-------------------------------------+
| ``DEFAULT_ARK_ETABLE_5``     | ``ARKSTEP_DEFAULT_ARK_ETABLE_4``    |
+------------------------------+-------------------------------------+
| ``DEFAULT_ARK_ITABLE_3``     | ``ARKSTEP_DEFAULT_ARK_ITABLE_3``    |
+------------------------------+-------------------------------------+
| ``DEFAULT_ARK_ITABLE_4``     | ``ARKSTEP_DEFAULT_ARK_ITABLE_4``    |
+------------------------------+-------------------------------------+
| ``DEFAULT_ARK_ITABLE_5``     | ``ARKSTEP_DEFAULT_ARK_ITABLE_5``    |
+------------------------------+-------------------------------------+
| ``DEFAULT_ERK_2``            | ``ERKSTEP_DEFAULT_2``               |
+------------------------------+-------------------------------------+
| ``DEFAULT_ERK_3``            | ``ERKSTEP_DEFAULT_3``               |
+------------------------------+-------------------------------------+
| ``DEFAULT_ERK_4``            | ``ERKSTEP_DEFAULT_4``               |
+------------------------------+-------------------------------------+
| ``DEFAULT_ERK_5``            | ``ERKSTEP_DEFAULT_5``               |
+------------------------------+-------------------------------------+
| ``DEFAULT_ERK_6``            | ``ERKSTEP_DEFAULT_6``               |
+------------------------------+-------------------------------------+
| ``DEFAULT_ERK_8``            | ``ERKSTEP_DEFAULT_8``               |
+------------------------------+-------------------------------------+

In addition, the following functions are now deprecated (compile-time warnings
will be thrown if supported by the compiler):

+---------------------------------+--------------------------------+
| Deprecated Name                 | New Name                       |
+=================================+================================+
| ``DenseGETRF``                  | ``SUNDlsMat_DenseGETRF``       |
+---------------------------------+--------------------------------+
| ``DenseGETRS``                  | ``SUNDlsMat_DenseGETRS``       |
+---------------------------------+--------------------------------+
| ``denseGETRF``                  | ``SUNDlsMat_denseGETRF``       |
+---------------------------------+--------------------------------+
| ``denseGETRS``                  | ``SUNDlsMat_denseGETRS``       |
+---------------------------------+--------------------------------+
| ``DensePOTRF``                  | ``SUNDlsMat_DensePOTRF``       |
+---------------------------------+--------------------------------+
| ``DensePOTRS``                  | ``SUNDlsMat_DensePOTRS``       |
+---------------------------------+--------------------------------+
| ``densePOTRF``                  | ``SUNDlsMat_densePOTRF``       |
+---------------------------------+--------------------------------+
| ``densePOTRS``                  | ``SUNDlsMat_densePOTRS``       |
+---------------------------------+--------------------------------+
| ``DenseGEQRF``                  | ``SUNDlsMat_DenseGEQRF``       |
+---------------------------------+--------------------------------+
| ``DenseORMQR``                  | ``SUNDlsMat_DenseORMQR``       |
+---------------------------------+--------------------------------+
| ``denseGEQRF``                  | ``SUNDlsMat_denseGEQRF``       |
+---------------------------------+--------------------------------+
| ``denseORMQR``                  | ``SUNDlsMat_denseORMQR``       |
+---------------------------------+--------------------------------+
| ``DenseCopy``                   | ``SUNDlsMat_DenseCopy``        |
+---------------------------------+--------------------------------+
| ``denseCopy``                   | ``SUNDlsMat_denseCopy``        |
+---------------------------------+--------------------------------+
| ``DenseScale``                  | ``SUNDlsMat_DenseScale``       |
+---------------------------------+--------------------------------+
| ``denseScale``                  | ``SUNDlsMat_denseScale``       |
+---------------------------------+--------------------------------+
| ``denseAddIdentity``            | ``SUNDlsMat_denseAddIdentity`` |
+---------------------------------+--------------------------------+
| ``DenseMatvec``                 | ``SUNDlsMat_DenseMatvec``      |
+---------------------------------+--------------------------------+
| ``denseMatvec``                 | ``SUNDlsMat_denseMatvec``      |
+---------------------------------+--------------------------------+
| ``BandGBTRF``                   | ``SUNDlsMat_BandGBTRF``        |
+---------------------------------+--------------------------------+
| ``bandGBTRF``                   | ``SUNDlsMat_bandGBTRF``        |
+---------------------------------+--------------------------------+
| ``BandGBTRS``                   | ``SUNDlsMat_BandGBTRS``        |
+---------------------------------+--------------------------------+
| ``bandGBTRS``                   | ``SUNDlsMat_bandGBTRS``        |
+---------------------------------+--------------------------------+
| ``BandCopy``                    | ``SUNDlsMat_BandCopy``         |
+---------------------------------+--------------------------------+
| ``bandCopy``                    | ``SUNDlsMat_bandCopy``         |
+---------------------------------+--------------------------------+
| ``BandScale``                   | ``SUNDlsMat_BandScale``        |
+---------------------------------+--------------------------------+
| ``bandScale``                   | ``SUNDlsMat_bandScale``        |
+---------------------------------+--------------------------------+
| ``bandAddIdentity``             | ``SUNDlsMat_bandAddIdentity``  |
+---------------------------------+--------------------------------+
| ``BandMatvec``                  | ``SUNDlsMat_BandMatvec``       |
+---------------------------------+--------------------------------+
| ``bandMatvec``                  | ``SUNDlsMat_bandMatvec``       |
+---------------------------------+--------------------------------+
| ``ModifiedGS``                  | ``SUNModifiedGS``              |
+---------------------------------+--------------------------------+
| ``ClassicalGS``                 | ``SUNClassicalGS``             |
+---------------------------------+--------------------------------+
| ``QRfact``                      | ``SUNQRFact``                  |
+---------------------------------+--------------------------------+
| ``QRsol``                       | ``SUNQRsol``                   |
+---------------------------------+--------------------------------+
| ``DlsMat_NewDenseMat``          | ``SUNDlsMat_NewDenseMat``      |
+---------------------------------+--------------------------------+
| ``DlsMat_NewBandMat``           | ``SUNDlsMat_NewBandMat``       |
+---------------------------------+--------------------------------+
| ``DestroyMat``                  | ``SUNDlsMat_DestroyMat``       |
+---------------------------------+--------------------------------+
| ``NewIntArray``                 | ``SUNDlsMat_NewIntArray``      |
+---------------------------------+--------------------------------+
| ``NewIndexArray``               | ``SUNDlsMat_NewIndexArray``    |
+---------------------------------+--------------------------------+
| ``NewRealArray``                | ``SUNDlsMat_NewRealArray``     |
+---------------------------------+--------------------------------+
| ``DestroyArray``                | ``SUNDlsMat_DestroyArray``     |
+---------------------------------+--------------------------------+
| ``AddIdentity``                 | ``SUNDlsMat_AddIdentity``      |
+---------------------------------+--------------------------------+
| ``SetToZero``                   | ``SUNDlsMat_SetToZero``        |
+---------------------------------+--------------------------------+
| ``PrintMat``                    | ``SUNDlsMat_PrintMat``         |
+---------------------------------+--------------------------------+
| ``newDenseMat``                 | ``SUNDlsMat_newDenseMat``      |
+---------------------------------+--------------------------------+
| ``newBandMat``                  | ``SUNDlsMat_newBandMat``       |
+---------------------------------+--------------------------------+
| ``destroyMat``                  | ``SUNDlsMat_destroyMat``       |
+---------------------------------+--------------------------------+
| ``newIntArray``                 | ``SUNDlsMat_newIntArray``      |
+---------------------------------+--------------------------------+
| ``newIndexArray``               | ``SUNDlsMat_newIndexArray``    |
+---------------------------------+--------------------------------+
| ``newRealArray``                | ``SUNDlsMat_newRealArray``     |
+---------------------------------+--------------------------------+
| ``destroyArray``                | ``SUNDlsMat_destroyArray``     |
+---------------------------------+--------------------------------+

In addition, the entire ``sundials_lapack.h`` header file is now deprecated for
removal in SUNDIALS v7.0.0. Note, this header file is not needed to use the
SUNDIALS LAPACK linear solvers.


### CVODE Changes in v6.0.0

**SUNContext**

SUNDIALS v6.0.0 introduces a new :c:type:`SUNContext` object on which all other
SUNDIALS objects depend. As such, the constructors for all SUNDIALS packages,
vectors, matrices, linear solvers, nonlinear solvers, and memory helpers have
been updated to accept a context as the last input. Users upgrading to SUNDIALS
v6.0.0 will need to call :c:func:`SUNContext_Create` to create a context object
with before calling any other SUNDIALS library function, and then provide this
object to other SUNDIALS constructors. The context object has been introduced to
allow SUNDIALS to provide new features, such as the profiling/instrumentation
also introduced in this release, while maintaining thread-safety. See the
documentation section on the :c:type:`SUNContext` for more details.

A script ``upgrade-to-sundials-6-from-5.sh`` has been provided with the release
(obtainable from the GitHub release page) to help ease the transition to
SUNDIALS v6.0.0. The script will add a ``SUNCTX_PLACEHOLDER`` argument to all of
the calls to SUNDIALS constructors that now require a ``SUNContext`` object. It
can also update deprecated SUNDIALS constants/types to the new names. It can be
run like this:

.. code-block::

   > ./upgrade-to-sundials-6-from-5.sh <files to update>

**SUNProfiler**

A capability to profile/instrument SUNDIALS library code has been added. This
can be enabled with the CMake option :cmakeop:`SUNDIALS_BUILD_WITH_PROFILING`. A
built-in profiler will be used by default, but the `Caliper
<https://github.com/LLNL/Caliper>`_ library can also be used instead with the
CMake option :cmakeop:`ENABLE_CALIPER`. See the documentation section on
profiling for more details.  **WARNING**: Profiling will impact performance, and
should be enabled judiciously.

**SUNMemoryHelper**

The :c:type:`SUNMemoryHelper` functions :c:func:`SUNMemoryHelper_Alloc`,
:c:func:`SUNMemoryHelper_Dealloc`, and :c:func:`SUNMemoryHelper_Copy` have been
updated to accept an opaque handle as the last input. At a minimum, user-defined
:c:type:`SUNMemoryHelper` implementations will need to update these functions to
accept the additional argument. Typically, this handle is the execution stream
(e.g., a CUDA/HIP stream or SYCL queue) for the operation. The :ref:`CUDA
<SUNMemory.CUDA>`, :ref:`HIP <SUNMemory.HIP>`, and :ref:`SYCL <SUNMemory.SYCL>`
implementations have been updated accordingly. Additionally, the constructor
:c:func:`SUNMemoryHelper_Sycl` has been updated to remove the SYCL queue as an
input.

**NVector**

Two new optional vector operations, :c:func:`N_VDotProdMultiLocal` and
:c:func:`N_VDotProdMultiAllReduce`, have been added to support
low-synchronization methods for Anderson acceleration.

The CUDA, HIP, and SYCL execution policies have been moved from the ``sundials``
namespace to the ``sundials::cuda``, ``sundials::hip``, and ``sundials::sycl``
namespaces respectively. Accordingly, the prefixes "Cuda", "Hip", and "Sycl"
have been removed from the execution policy classes and methods.

The ``Sundials`` namespace used by the Trilinos Tpetra NVector has been replaced
with the ``sundials::trilinos::nvector_tpetra`` namespace.

The serial, PThreads, PETSc, *hypre*, Parallel, OpenMP_DEV, and OpenMP vector
functions ``N_VCloneVectorArray_*`` and ``N_VDestroyVectorArray_*`` have been
deprecated. The generic :c:func:`N_VCloneVectorArray` and
:c:func:`N_VDestroyVectorArray` functions should be used instead.

The previously deprecated constructor ``N_VMakeWithManagedAllocator_Cuda`` and
the function ``N_VSetCudaStream_Cuda`` have been removed and replaced with
:c:func:`N_VNewWithMemHelp_Cuda` and :c:func:`N_VSetKerrnelExecPolicy_Cuda`
respectively.

The previously deprecated macros ``PVEC_REAL_MPI_TYPE`` and
``PVEC_INTEGER_MPI_TYPE`` have been removed and replaced with
``MPI_SUNREALTYPE`` and ``MPI_SUNINDEXTYPE`` respectively.

**SUNLinearSolver**

The following previously deprecated functions have been removed:

+-----------------------------+------------------------------------------+
| Removed                     | Replacement                              |
+=============================+==========================================+
| ``SUNBandLinearSolver``     | :c:func:`SUNLinSol_Band`                 |
+-----------------------------+------------------------------------------+
| ``SUNDenseLinearSolver``    | :c:func:`SUNLinSol_Dense`                |
+-----------------------------+------------------------------------------+
| ``SUNKLU``                  | :c:func:`SUNLinSol_KLU`                  |
+-----------------------------+------------------------------------------+
| ``SUNKLUReInit``            | :c:func:`SUNLinSol_KLUReInit`            |
+-----------------------------+------------------------------------------+
| ``SUNKLUSetOrdering``       | :c:func:`SUNLinSol_KLUSetOrdering`       |
+-----------------------------+------------------------------------------+
| ``SUNLapackBand``           | :c:func:`SUNLinSol_LapackBand`           |
+-----------------------------+------------------------------------------+
| ``SUNLapackDense``          | :c:func:`SUNLinSol_LapackDense`          |
+-----------------------------+------------------------------------------+
| ``SUNPCG``                  | :c:func:`SUNLinSol_PCG`                  |
+-----------------------------+------------------------------------------+
| ``SUNPCGSetPrecType``       | :c:func:`SUNLinSol_PCGSetPrecType`       |
+-----------------------------+------------------------------------------+
| ``SUNPCGSetMaxl``           | :c:func:`SUNLinSol_PCGSetMaxl`           |
+-----------------------------+------------------------------------------+
| ``SUNSPBCGS``               | :c:func:`SUNLinSol_SPBCGS`               |
+-----------------------------+------------------------------------------+
| ``SUNSPBCGSSetPrecType``    | :c:func:`SUNLinSol_SPBCGSSetPrecType`    |
+-----------------------------+------------------------------------------+
| ``SUNSPBCGSSetMaxl``        | :c:func:`SUNLinSol_SPBCGSSetMaxl`        |
+-----------------------------+------------------------------------------+
| ``SUNSPFGMR``               | :c:func:`SUNLinSol_SPFGMR`               |
+-----------------------------+------------------------------------------+
| ``SUNSPFGMRSetPrecType``    | :c:func:`SUNLinSol_SPFGMRSetPrecType`    |
+-----------------------------+------------------------------------------+
| ``SUNSPFGMRSetGSType``      | :c:func:`SUNLinSol_SPFGMRSetGSType`      |
+-----------------------------+------------------------------------------+
| ``SUNSPFGMRSetMaxRestarts`` | :c:func:`SUNLinSol_SPFGMRSetMaxRestarts` |
+-----------------------------+------------------------------------------+
| ``SUNSPGMR``                | :c:func:`SUNLinSol_SPGMR`                |
+-----------------------------+------------------------------------------+
| ``SUNSPGMRSetPrecType``     | :c:func:`SUNLinSol_SPGMRSetPrecType`     |
+-----------------------------+------------------------------------------+
| ``SUNSPGMRSetGSType``       | :c:func:`SUNLinSol_SPGMRSetGSType`       |
+-----------------------------+------------------------------------------+
| ``SUNSPGMRSetMaxRestarts``  | :c:func:`SUNLinSol_SPGMRSetMaxRestarts`  |
+-----------------------------+------------------------------------------+
| ``SUNSPTFQMR``              | :c:func:`SUNLinSol_SPTFQMR`              |
+-----------------------------+------------------------------------------+
| ``SUNSPTFQMRSetPrecType``   | :c:func:`SUNLinSol_SPTFQMRSetPrecType`   |
+-----------------------------+------------------------------------------+
| ``SUNSPTFQMRSetMaxl``       | :c:func:`SUNLinSol_SPTFQMRSetMaxl`       |
+-----------------------------+------------------------------------------+
| ``SUNSuperLUMT``            | :c:func:`SUNLinSol_SuperLUMT`            |
+-----------------------------+------------------------------------------+
| ``SUNSuperLUMTSetOrdering`` | :c:func:`SUNLinSol_SuperLUMTSetOrdering` |
+-----------------------------+------------------------------------------+

**CVODE**

The previously deprecated function ``CVodeSetMaxStepsBetweenJac`` has been
removed and replaced with :c:func:`CVodeSetJacEvalFrequency`.

The CVODE Fortran 77 interface has been removed. See :numref:`SUNDIALS.Fortran`
and the F2003 example programs for more details using the SUNDIALS Fortran 2003
module interfaces.

**Deprecations**

In addition to the deprecations noted elsewhere, many constants, types, and
functions have been renamed so that they are properly namespaced. The old names
have been deprecated and will be removed in SUNDIALS v7.0.0.

The following constants, macros, and typedefs are now deprecated:

+------------------------------+-------------------------------------+
| Deprecated Name              | New Name                            |
+==============================+=====================================+
| ``realtype``                 | ``sunrealtype``                     |
+------------------------------+-------------------------------------+
| ``booleantype``              | ``sunbooleantype``                  |
+------------------------------+-------------------------------------+
| ``RCONST``                   | ``SUN_RCONST``                      |
+------------------------------+-------------------------------------+
| ``BIG_REAL``                 | ``SUN_BIG_REAL``                    |
+------------------------------+-------------------------------------+
| ``SMALL_REAL``               | ``SUN_SMALL_REAL``                  |
+------------------------------+-------------------------------------+
| ``UNIT_ROUNDOFF``            | ``SUN_UNIT_ROUNDOFF``               |
+------------------------------+-------------------------------------+
| ``PREC_NONE``                | ``SUN_PREC_NONE``                   |
+------------------------------+-------------------------------------+
| ``PREC_LEFT``                | ``SUN_PREC_LEFT``                   |
+------------------------------+-------------------------------------+
| ``PREC_RIGHT``               | ``SUN_PREC_RIGHT``                  |
+------------------------------+-------------------------------------+
| ``PREC_BOTH``                | ``SUN_PREC_BOTH``                   |
+------------------------------+-------------------------------------+
| ``MODIFIED_GS``              | ``SUN_MODIFIED_GS``                 |
+------------------------------+-------------------------------------+
| ``CLASSICAL_GS``             | ``SUN_CLASSICAL_GS``                |
+------------------------------+-------------------------------------+
| ``ATimesFn``                 | ``SUNATimesFn``                     |
+------------------------------+-------------------------------------+
| ``PSetupFn``                 | ``SUNPSetupFn``                     |
+------------------------------+-------------------------------------+
| ``PSolveFn``                 | ``SUNPSolveFn``                     |
+------------------------------+-------------------------------------+
| ``DlsMat``                   | ``SUNDlsMat``                       |
+------------------------------+-------------------------------------+
| ``DENSE_COL``                | ``SUNDLS_DENSE_COL``                |
+------------------------------+-------------------------------------+
| ``DENSE_ELEM``               | ``SUNDLS_DENSE_ELEM``               |
+------------------------------+-------------------------------------+
| ``BAND_COL``                 | ``SUNDLS_BAND_COL``                 |
+------------------------------+-------------------------------------+
| ``BAND_COL_ELEM``            | ``SUNDLS_BAND_COL_ELEM``            |
+------------------------------+-------------------------------------+
| ``BAND_ELEM``                | ``SUNDLS_BAND_ELEM``                |
+------------------------------+-------------------------------------+

In addition, the following functions are now deprecated (compile-time warnings
will be thrown if supported by the compiler):

+---------------------------------+--------------------------------+
| Deprecated Name                 | New Name                       |
+=================================+================================+
| ``CVSpilsSetLinearSolver``      | ``CVodeSetLinearSolver``       |
+---------------------------------+--------------------------------+
| ``CVSpilsSetEpsLin``            | ``CVodeSetEpsLin``             |
+---------------------------------+--------------------------------+
| ``CVSpilsSetPreconditioner``    | ``CVodeSetPreconditioner``     |
+---------------------------------+--------------------------------+
| ``CVSpilsSetJacTimes``          | ``CVodeSetJacTimes``           |
+---------------------------------+--------------------------------+
| ``CVSpilsGetWorkSpace``         | ``CVodeGetLinWorkSpace``       |
+---------------------------------+--------------------------------+
| ``CVSpilsGetNumPrecEvals``      | ``CVodeGetNumPrecEvals``       |
+---------------------------------+--------------------------------+
| ``CVSpilsGetNumPrecSolves``     | ``CVodeGetNumPrecSolves``      |
+---------------------------------+--------------------------------+
| ``CVSpilsGetNumLinIters``       | ``CVodeGetNumLinIters``        |
+---------------------------------+--------------------------------+
| ``CVSpilsGetNumConvFails``      | ``CVodeGetNumConvFails``       |
+---------------------------------+--------------------------------+
| ``CVSpilsGetNumJTSetupEvals``   | ``CVodeGetNumJTSetupEvals``    |
+---------------------------------+--------------------------------+
| ``CVSpilsGetNumJtimesEvals``    | ``CVodeGetNumJtimesEvals``     |
+---------------------------------+--------------------------------+
| ``CVSpilsGetNumRhsEvals``       | ``CVodeGetNumLinRhsEvals``     |
+---------------------------------+--------------------------------+
| ``CVSpilsGetLastFlag``          | ``CVodeGetLastLinFlag``        |
+---------------------------------+--------------------------------+
| ``CVSpilsGetReturnFlagName``    | ``CVodeGetLinReturnFlagName``  |
+---------------------------------+--------------------------------+
| ``CVDlsSetLinearSolver``        | ``CVodeSetLinearSolver``       |
+---------------------------------+--------------------------------+
| ``CVDlsSetJacFn``               | ``CVodeSetJacFn``              |
+---------------------------------+--------------------------------+
| ``CVDlsGetWorkSpace``           | ``CVodeGetLinWorkSpace``       |
+---------------------------------+--------------------------------+
| ``CVDlsGetNumJacEvals``         | ``CVodeGetNumJacEvals``        |
+---------------------------------+--------------------------------+
| ``CVDlsGetNumRhsEvals``         | ``CVodeGetNumLinRhsEvals``     |
+---------------------------------+--------------------------------+
| ``CVDlsGetLastFlag``            | ``CVodeGetLastLinFlag``        |
+---------------------------------+--------------------------------+
| ``CVDlsGetReturnFlagName``      | ``CVodeGetLinReturnFlagName``  |
+---------------------------------+--------------------------------+
| ``DenseGETRF``                  | ``SUNDlsMat_DenseGETRF``       |
+---------------------------------+--------------------------------+
| ``DenseGETRS``                  | ``SUNDlsMat_DenseGETRS``       |
+---------------------------------+--------------------------------+
| ``denseGETRF``                  | ``SUNDlsMat_denseGETRF``       |
+---------------------------------+--------------------------------+
| ``denseGETRS``                  | ``SUNDlsMat_denseGETRS``       |
+---------------------------------+--------------------------------+
| ``DensePOTRF``                  | ``SUNDlsMat_DensePOTRF``       |
+---------------------------------+--------------------------------+
| ``DensePOTRS``                  | ``SUNDlsMat_DensePOTRS``       |
+---------------------------------+--------------------------------+
| ``densePOTRF``                  | ``SUNDlsMat_densePOTRF``       |
+---------------------------------+--------------------------------+
| ``densePOTRS``                  | ``SUNDlsMat_densePOTRS``       |
+---------------------------------+--------------------------------+
| ``DenseGEQRF``                  | ``SUNDlsMat_DenseGEQRF``       |
+---------------------------------+--------------------------------+
| ``DenseORMQR``                  | ``SUNDlsMat_DenseORMQR``       |
+---------------------------------+--------------------------------+
| ``denseGEQRF``                  | ``SUNDlsMat_denseGEQRF``       |
+---------------------------------+--------------------------------+
| ``denseORMQR``                  | ``SUNDlsMat_denseORMQR``       |
+---------------------------------+--------------------------------+
| ``DenseCopy``                   | ``SUNDlsMat_DenseCopy``        |
+---------------------------------+--------------------------------+
| ``denseCopy``                   | ``SUNDlsMat_denseCopy``        |
+---------------------------------+--------------------------------+
| ``DenseScale``                  | ``SUNDlsMat_DenseScale``       |
+---------------------------------+--------------------------------+
| ``denseScale``                  | ``SUNDlsMat_denseScale``       |
+---------------------------------+--------------------------------+
| ``denseAddIdentity``            | ``SUNDlsMat_denseAddIdentity`` |
+---------------------------------+--------------------------------+
| ``DenseMatvec``                 | ``SUNDlsMat_DenseMatvec``      |
+---------------------------------+--------------------------------+
| ``denseMatvec``                 | ``SUNDlsMat_denseMatvec``      |
+---------------------------------+--------------------------------+
| ``BandGBTRF``                   | ``SUNDlsMat_BandGBTRF``        |
+---------------------------------+--------------------------------+
| ``bandGBTRF``                   | ``SUNDlsMat_bandGBTRF``        |
+---------------------------------+--------------------------------+
| ``BandGBTRS``                   | ``SUNDlsMat_BandGBTRS``        |
+---------------------------------+--------------------------------+
| ``bandGBTRS``                   | ``SUNDlsMat_bandGBTRS``        |
+---------------------------------+--------------------------------+
| ``BandCopy``                    | ``SUNDlsMat_BandCopy``         |
+---------------------------------+--------------------------------+
| ``bandCopy``                    | ``SUNDlsMat_bandCopy``         |
+---------------------------------+--------------------------------+
| ``BandScale``                   | ``SUNDlsMat_BandScale``        |
+---------------------------------+--------------------------------+
| ``bandScale``                   | ``SUNDlsMat_bandScale``        |
+---------------------------------+--------------------------------+
| ``bandAddIdentity``             | ``SUNDlsMat_bandAddIdentity``  |
+---------------------------------+--------------------------------+
| ``BandMatvec``                  | ``SUNDlsMat_BandMatvec``       |
+---------------------------------+--------------------------------+
| ``bandMatvec``                  | ``SUNDlsMat_bandMatvec``       |
+---------------------------------+--------------------------------+
| ``ModifiedGS``                  | ``SUNModifiedGS``              |
+---------------------------------+--------------------------------+
| ``ClassicalGS``                 | ``SUNClassicalGS``             |
+---------------------------------+--------------------------------+
| ``QRfact``                      | ``SUNQRFact``                  |
+---------------------------------+--------------------------------+
| ``QRsol``                       | ``SUNQRsol``                   |
+---------------------------------+--------------------------------+
| ``DlsMat_NewDenseMat``          | ``SUNDlsMat_NewDenseMat``      |
+---------------------------------+--------------------------------+
| ``DlsMat_NewBandMat``           | ``SUNDlsMat_NewBandMat``       |
+---------------------------------+--------------------------------+
| ``DestroyMat``                  | ``SUNDlsMat_DestroyMat``       |
+---------------------------------+--------------------------------+
| ``NewIntArray``                 | ``SUNDlsMat_NewIntArray``      |
+---------------------------------+--------------------------------+
| ``NewIndexArray``               | ``SUNDlsMat_NewIndexArray``    |
+---------------------------------+--------------------------------+
| ``NewRealArray``                | ``SUNDlsMat_NewRealArray``     |
+---------------------------------+--------------------------------+
| ``DestroyArray``                | ``SUNDlsMat_DestroyArray``     |
+---------------------------------+--------------------------------+
| ``AddIdentity``                 | ``SUNDlsMat_AddIdentity``      |
+---------------------------------+--------------------------------+
| ``SetToZero``                   | ``SUNDlsMat_SetToZero``        |
+---------------------------------+--------------------------------+
| ``PrintMat``                    | ``SUNDlsMat_PrintMat``         |
+---------------------------------+--------------------------------+
| ``newDenseMat``                 | ``SUNDlsMat_newDenseMat``      |
+---------------------------------+--------------------------------+
| ``newBandMat``                  | ``SUNDlsMat_newBandMat``       |
+---------------------------------+--------------------------------+
| ``destroyMat``                  | ``SUNDlsMat_destroyMat``       |
+---------------------------------+--------------------------------+
| ``newIntArray``                 | ``SUNDlsMat_newIntArray``      |
+---------------------------------+--------------------------------+
| ``newIndexArray``               | ``SUNDlsMat_newIndexArray``    |
+---------------------------------+--------------------------------+
| ``newRealArray``                | ``SUNDlsMat_newRealArray``     |
+---------------------------------+--------------------------------+
| ``destroyArray``                | ``SUNDlsMat_destroyArray``     |
+---------------------------------+--------------------------------+

In addition, the entire ``sundials_lapack.h`` header file is now deprecated for
removal in SUNDIALS v7.0.0. Note, this header file is not needed to use the
SUNDIALS LAPACK linear solvers.

### CVODES Changes in v6.0.0

**SUNContext**

SUNDIALS v6.0.0 introduces a new :c:type:`SUNContext` object on which all other
SUNDIALS objects depend. As such, the constructors for all SUNDIALS packages,
vectors, matrices, linear solvers, nonlinear solvers, and memory helpers have
been updated to accept a context as the last input. Users upgrading to SUNDIALS
v6.0.0 will need to call :c:func:`SUNContext_Create` to create a context object
with before calling any other SUNDIALS library function, and then provide this
object to other SUNDIALS constructors. The context object has been introduced to
allow SUNDIALS to provide new features, such as the profiling/instrumentation
also introduced in this release, while maintaining thread-safety. See the
documentation section on the :c:type:`SUNContext` for more details.

A script ``upgrade-to-sundials-6-from-5.sh`` has been provided with the release
(obtainable from the GitHub release page) to help ease the transition to
SUNDIALS v6.0.0. The script will add a ``SUNCTX_PLACEHOLDER`` argument to all of
the calls to SUNDIALS constructors that now require a ``SUNContext`` object. It
can also update deprecated SUNDIALS constants/types to the new names. It can be
run like this:

.. code-block::

   > ./upgrade-to-sundials-6-from-5.sh <files to update>

**SUNProfiler**

A capability to profile/instrument SUNDIALS library code has been added. This
can be enabled with the CMake option :cmakeop:`SUNDIALS_BUILD_WITH_PROFILING`. A
built-in profiler will be used by default, but the `Caliper
<https://github.com/LLNL/Caliper>`_ library can also be used instead with the
CMake option :cmakeop:`ENABLE_CALIPER`. See the documentation section on
profiling for more details.  **WARNING**: Profiling will impact performance, and
should be enabled judiciously.

**SUNMemoryHelper**

The :c:type:`SUNMemoryHelper` functions :c:func:`SUNMemoryHelper_Alloc`,
:c:func:`SUNMemoryHelper_Dealloc`, and :c:func:`SUNMemoryHelper_Copy` have been
updated to accept an opaque handle as the last input. At a minimum, user-defined
:c:type:`SUNMemoryHelper` implementations will need to update these functions to
accept the additional argument. Typically, this handle is the execution stream
(e.g., a CUDA/HIP stream or SYCL queue) for the operation. The :ref:`CUDA
<SUNMemory.CUDA>`, :ref:`HIP <SUNMemory.HIP>`, and :ref:`SYCL <SUNMemory.SYCL>`
implementations have been updated accordingly. Additionally, the constructor
:c:func:`SUNMemoryHelper_Sycl` has been updated to remove the SYCL queue as an
input.

**NVector**

Two new optional vector operations, :c:func:`N_VDotProdMultiLocal` and
:c:func:`N_VDotProdMultiAllReduce`, have been added to support
low-synchronization methods for Anderson acceleration.

The CUDA, HIP, and SYCL execution policies have been moved from the ``sundials``
namespace to the ``sundials::cuda``, ``sundials::hip``, and ``sundials::sycl``
namespaces respectively. Accordingly, the prefixes "Cuda", "Hip", and "Sycl"
have been removed from the execution policy classes and methods.

The ``Sundials`` namespace used by the Trilinos Tpetra NVector has been replaced
with the ``sundials::trilinos::nvector_tpetra`` namespace.

The serial, PThreads, PETSc, *hypre*, Parallel, OpenMP_DEV, and OpenMP vector
functions ``N_VCloneVectorArray_*`` and ``N_VDestroyVectorArray_*`` have been
deprecated. The generic :c:func:`N_VCloneVectorArray` and
:c:func:`N_VDestroyVectorArray` functions should be used instead.

The previously deprecated constructor ``N_VMakeWithManagedAllocator_Cuda`` and
the function ``N_VSetCudaStream_Cuda`` have been removed and replaced with
:c:func:`N_VNewWithMemHelp_Cuda` and :c:func:`N_VSetKerrnelExecPolicy_Cuda`
respectively.

The previously deprecated macros ``PVEC_REAL_MPI_TYPE`` and
``PVEC_INTEGER_MPI_TYPE`` have been removed and replaced with
``MPI_SUNREALTYPE`` and ``MPI_SUNINDEXTYPE`` respectively.

**SUNLinearSolver**

The following previously deprecated functions have been removed:

+-----------------------------+------------------------------------------+
| Removed                     | Replacement                              |
+=============================+==========================================+
| ``SUNBandLinearSolver``     | :c:func:`SUNLinSol_Band`                 |
+-----------------------------+------------------------------------------+
| ``SUNDenseLinearSolver``    | :c:func:`SUNLinSol_Dense`                |
+-----------------------------+------------------------------------------+
| ``SUNKLU``                  | :c:func:`SUNLinSol_KLU`                  |
+-----------------------------+------------------------------------------+
| ``SUNKLUReInit``            | :c:func:`SUNLinSol_KLUReInit`            |
+-----------------------------+------------------------------------------+
| ``SUNKLUSetOrdering``       | :c:func:`SUNLinSol_KLUSetOrdering`       |
+-----------------------------+------------------------------------------+
| ``SUNLapackBand``           | :c:func:`SUNLinSol_LapackBand`           |
+-----------------------------+------------------------------------------+
| ``SUNLapackDense``          | :c:func:`SUNLinSol_LapackDense`          |
+-----------------------------+------------------------------------------+
| ``SUNPCG``                  | :c:func:`SUNLinSol_PCG`                  |
+-----------------------------+------------------------------------------+
| ``SUNPCGSetPrecType``       | :c:func:`SUNLinSol_PCGSetPrecType`       |
+-----------------------------+------------------------------------------+
| ``SUNPCGSetMaxl``           | :c:func:`SUNLinSol_PCGSetMaxl`           |
+-----------------------------+------------------------------------------+
| ``SUNSPBCGS``               | :c:func:`SUNLinSol_SPBCGS`               |
+-----------------------------+------------------------------------------+
| ``SUNSPBCGSSetPrecType``    | :c:func:`SUNLinSol_SPBCGSSetPrecType`    |
+-----------------------------+------------------------------------------+
| ``SUNSPBCGSSetMaxl``        | :c:func:`SUNLinSol_SPBCGSSetMaxl`        |
+-----------------------------+------------------------------------------+
| ``SUNSPFGMR``               | :c:func:`SUNLinSol_SPFGMR`               |
+-----------------------------+------------------------------------------+
| ``SUNSPFGMRSetPrecType``    | :c:func:`SUNLinSol_SPFGMRSetPrecType`    |
+-----------------------------+------------------------------------------+
| ``SUNSPFGMRSetGSType``      | :c:func:`SUNLinSol_SPFGMRSetGSType`      |
+-----------------------------+------------------------------------------+
| ``SUNSPFGMRSetMaxRestarts`` | :c:func:`SUNLinSol_SPFGMRSetMaxRestarts` |
+-----------------------------+------------------------------------------+
| ``SUNSPGMR``                | :c:func:`SUNLinSol_SPGMR`                |
+-----------------------------+------------------------------------------+
| ``SUNSPGMRSetPrecType``     | :c:func:`SUNLinSol_SPGMRSetPrecType`     |
+-----------------------------+------------------------------------------+
| ``SUNSPGMRSetGSType``       | :c:func:`SUNLinSol_SPGMRSetGSType`       |
+-----------------------------+------------------------------------------+
| ``SUNSPGMRSetMaxRestarts``  | :c:func:`SUNLinSol_SPGMRSetMaxRestarts`  |
+-----------------------------+------------------------------------------+
| ``SUNSPTFQMR``              | :c:func:`SUNLinSol_SPTFQMR`              |
+-----------------------------+------------------------------------------+
| ``SUNSPTFQMRSetPrecType``   | :c:func:`SUNLinSol_SPTFQMRSetPrecType`   |
+-----------------------------+------------------------------------------+
| ``SUNSPTFQMRSetMaxl``       | :c:func:`SUNLinSol_SPTFQMRSetMaxl`       |
+-----------------------------+------------------------------------------+
| ``SUNSuperLUMT``            | :c:func:`SUNLinSol_SuperLUMT`            |
+-----------------------------+------------------------------------------+
| ``SUNSuperLUMTSetOrdering`` | :c:func:`SUNLinSol_SuperLUMTSetOrdering` |
+-----------------------------+------------------------------------------+

**CVODES**

Added a new function :c:func:`CVodeGetLinSolveStats` to get the CVODES linear
solver statistics as a group.

Added a new function, :c:func:`CVodeSetMonitorFn`, that takes a user-function
to be called by CVODES after every `nst` successfully completed time-steps.
This is intended to provide a way of monitoring the CVODES statistics
throughout the simulation.

The previously deprecated function ``CVodeSetMaxStepsBetweenJac`` has been
removed and replaced with :c:func:`CVodeSetJacEvalFrequency`.

**Deprecations**

In addition to the deprecations noted elsewhere, many constants, types, and
functions have been renamed so that they are properly namespaced. The old names
have been deprecated and will be removed in SUNDIALS v7.0.0.

The following constants, macros, and typedefs are now deprecated:

+------------------------------+-------------------------------------+
| Deprecated Name              | New Name                            |
+==============================+=====================================+
| ``realtype``                 | ``sunrealtype``                     |
+------------------------------+-------------------------------------+
| ``booleantype``              | ``sunbooleantype``                  |
+------------------------------+-------------------------------------+
| ``RCONST``                   | ``SUN_RCONST``                      |
+------------------------------+-------------------------------------+
| ``BIG_REAL``                 | ``SUN_BIG_REAL``                    |
+------------------------------+-------------------------------------+
| ``SMALL_REAL``               | ``SUN_SMALL_REAL``                  |
+------------------------------+-------------------------------------+
| ``UNIT_ROUNDOFF``            | ``SUN_UNIT_ROUNDOFF``               |
+------------------------------+-------------------------------------+
| ``PREC_NONE``                | ``SUN_PREC_NONE``                   |
+------------------------------+-------------------------------------+
| ``PREC_LEFT``                | ``SUN_PREC_LEFT``                   |
+------------------------------+-------------------------------------+
| ``PREC_RIGHT``               | ``SUN_PREC_RIGHT``                  |
+------------------------------+-------------------------------------+
| ``PREC_BOTH``                | ``SUN_PREC_BOTH``                   |
+------------------------------+-------------------------------------+
| ``MODIFIED_GS``              | ``SUN_MODIFIED_GS``                 |
+------------------------------+-------------------------------------+
| ``CLASSICAL_GS``             | ``SUN_CLASSICAL_GS``                |
+------------------------------+-------------------------------------+
| ``ATimesFn``                 | ``SUNATimesFn``                     |
+------------------------------+-------------------------------------+
| ``PSetupFn``                 | ``SUNPSetupFn``                     |
+------------------------------+-------------------------------------+
| ``PSolveFn``                 | ``SUNPSolveFn``                     |
+------------------------------+-------------------------------------+
| ``DlsMat``                   | ``SUNDlsMat``                       |
+------------------------------+-------------------------------------+
| ``DENSE_COL``                | ``SUNDLS_DENSE_COL``                |
+------------------------------+-------------------------------------+
| ``DENSE_ELEM``               | ``SUNDLS_DENSE_ELEM``               |
+------------------------------+-------------------------------------+
| ``BAND_COL``                 | ``SUNDLS_BAND_COL``                 |
+------------------------------+-------------------------------------+
| ``BAND_COL_ELEM``            | ``SUNDLS_BAND_COL_ELEM``            |
+------------------------------+-------------------------------------+
| ``BAND_ELEM``                | ``SUNDLS_BAND_ELEM``                |
+------------------------------+-------------------------------------+

In addition, the following functions are now deprecated (compile-time warnings
will be thrown if supported by the compiler):

+---------------------------------+--------------------------------+
| Deprecated Name                 | New Name                       |
+=================================+================================+
| ``CVSpilsSetLinearSolver``      | ``CVodeSetLinearSolver``       |
+---------------------------------+--------------------------------+
| ``CVSpilsSetEpsLin``            | ``CVodeSetEpsLin``             |
+---------------------------------+--------------------------------+
| ``CVSpilsSetPreconditioner``    | ``CVodeSetPreconditioner``     |
+---------------------------------+--------------------------------+
| ``CVSpilsSetJacTimes``          | ``CVodeSetJacTimes``           |
+---------------------------------+--------------------------------+
| ``CVSpilsGetWorkSpace``         | ``CVodeGetLinWorkSpace``       |
+---------------------------------+--------------------------------+
| ``CVSpilsGetNumPrecEvals``      | ``CVodeGetNumPrecEvals``       |
+---------------------------------+--------------------------------+
| ``CVSpilsGetNumPrecSolves``     | ``CVodeGetNumPrecSolves``      |
+---------------------------------+--------------------------------+
| ``CVSpilsGetNumLinIters``       | ``CVodeGetNumLinIters``        |
+---------------------------------+--------------------------------+
| ``CVSpilsGetNumConvFails``      | ``CVodeGetNumConvFails``       |
+---------------------------------+--------------------------------+
| ``CVSpilsGetNumJTSetupEvals``   | ``CVodeGetNumJTSetupEvals``    |
+---------------------------------+--------------------------------+
| ``CVSpilsGetNumJtimesEvals``    | ``CVodeGetNumJtimesEvals``     |
+---------------------------------+--------------------------------+
| ``CVSpilsGetNumRhsEvals``       | ``CVodeGetNumLinRhsEvals``     |
+---------------------------------+--------------------------------+
| ``CVSpilsGetLastFlag``          | ``CVodeGetLastLinFlag``        |
+---------------------------------+--------------------------------+
| ``CVSpilsGetReturnFlagName``    | ``CVodeGetLinReturnFlagName``  |
+---------------------------------+--------------------------------+
| ``CVSpilsSetLinearSolverB``     | ``CVodeSetLinearSolverB``      |
+---------------------------------+--------------------------------+
| ``CVSpilsSetEpsLinB``           | ``CVodeSetEpsLinB``            |
+---------------------------------+--------------------------------+
| ``CVSpilsSetPreconditionerB``   | ``CVodeSetPreconditionerB``    |
+---------------------------------+--------------------------------+
| ``CVSpilsSetPreconditionerBS``  | ``CVodeSetPreconditionerBS``   |
+---------------------------------+--------------------------------+
| ``CVSpilsSetJacTimesB``         | ``CVodeSetJacTimesB``          |
+---------------------------------+--------------------------------+
| ``CVSpilsSetJacTimesBS``        | ``CVodeSetJacTimesBS``         |
+---------------------------------+--------------------------------+
| ``CVDlsSetLinearSolver``        | ``CVodeSetLinearSolver``       |
+---------------------------------+--------------------------------+
| ``CVDlsSetJacFn``               | ``CVodeSetJacFn``              |
+---------------------------------+--------------------------------+
| ``CVDlsGetWorkSpace``           | ``CVodeGetLinWorkSpace``       |
+---------------------------------+--------------------------------+
| ``CVDlsGetNumJacEvals``         | ``CVodeGetNumJacEvals``        |
+---------------------------------+--------------------------------+
| ``CVDlsGetNumRhsEvals``         | ``CVodeGetNumLinRhsEvals``     |
+---------------------------------+--------------------------------+
| ``CVDlsGetLastFlag``            | ``CVodeGetLastLinFlag``        |
+---------------------------------+--------------------------------+
| ``CVDlsGetReturnFlagName``      | ``CVodeGetLinReturnFlagName``  |
+---------------------------------+--------------------------------+
| ``CVDlsSetLinearSolverB``       | ``CVodeSetLinearSolverB``      |
+---------------------------------+--------------------------------+
| ``CVDlsSetJacFnB``              | ``CVodeSetJacFnB``             |
+---------------------------------+--------------------------------+
| ``CVDlsSetJacFnBS``             | ``CVodeSetJacFnBS``            |
+---------------------------------+--------------------------------+
| ``DenseGETRF``                  | ``SUNDlsMat_DenseGETRF``       |
+---------------------------------+--------------------------------+
| ``DenseGETRS``                  | ``SUNDlsMat_DenseGETRS``       |
+---------------------------------+--------------------------------+
| ``denseGETRF``                  | ``SUNDlsMat_denseGETRF``       |
+---------------------------------+--------------------------------+
| ``denseGETRS``                  | ``SUNDlsMat_denseGETRS``       |
+---------------------------------+--------------------------------+
| ``DensePOTRF``                  | ``SUNDlsMat_DensePOTRF``       |
+---------------------------------+--------------------------------+
| ``DensePOTRS``                  | ``SUNDlsMat_DensePOTRS``       |
+---------------------------------+--------------------------------+
| ``densePOTRF``                  | ``SUNDlsMat_densePOTRF``       |
+---------------------------------+--------------------------------+
| ``densePOTRS``                  | ``SUNDlsMat_densePOTRS``       |
+---------------------------------+--------------------------------+
| ``DenseGEQRF``                  | ``SUNDlsMat_DenseGEQRF``       |
+---------------------------------+--------------------------------+
| ``DenseORMQR``                  | ``SUNDlsMat_DenseORMQR``       |
+---------------------------------+--------------------------------+
| ``denseGEQRF``                  | ``SUNDlsMat_denseGEQRF``       |
+---------------------------------+--------------------------------+
| ``denseORMQR``                  | ``SUNDlsMat_denseORMQR``       |
+---------------------------------+--------------------------------+
| ``DenseCopy``                   | ``SUNDlsMat_DenseCopy``        |
+---------------------------------+--------------------------------+
| ``denseCopy``                   | ``SUNDlsMat_denseCopy``        |
+---------------------------------+--------------------------------+
| ``DenseScale``                  | ``SUNDlsMat_DenseScale``       |
+---------------------------------+--------------------------------+
| ``denseScale``                  | ``SUNDlsMat_denseScale``       |
+---------------------------------+--------------------------------+
| ``denseAddIdentity``            | ``SUNDlsMat_denseAddIdentity`` |
+---------------------------------+--------------------------------+
| ``DenseMatvec``                 | ``SUNDlsMat_DenseMatvec``      |
+---------------------------------+--------------------------------+
| ``denseMatvec``                 | ``SUNDlsMat_denseMatvec``      |
+---------------------------------+--------------------------------+
| ``BandGBTRF``                   | ``SUNDlsMat_BandGBTRF``        |
+---------------------------------+--------------------------------+
| ``bandGBTRF``                   | ``SUNDlsMat_bandGBTRF``        |
+---------------------------------+--------------------------------+
| ``BandGBTRS``                   | ``SUNDlsMat_BandGBTRS``        |
+---------------------------------+--------------------------------+
| ``bandGBTRS``                   | ``SUNDlsMat_bandGBTRS``        |
+---------------------------------+--------------------------------+
| ``BandCopy``                    | ``SUNDlsMat_BandCopy``         |
+---------------------------------+--------------------------------+
| ``bandCopy``                    | ``SUNDlsMat_bandCopy``         |
+---------------------------------+--------------------------------+
| ``BandScale``                   | ``SUNDlsMat_BandScale``        |
+---------------------------------+--------------------------------+
| ``bandScale``                   | ``SUNDlsMat_bandScale``        |
+---------------------------------+--------------------------------+
| ``bandAddIdentity``             | ``SUNDlsMat_bandAddIdentity``  |
+---------------------------------+--------------------------------+
| ``BandMatvec``                  | ``SUNDlsMat_BandMatvec``       |
+---------------------------------+--------------------------------+
| ``bandMatvec``                  | ``SUNDlsMat_bandMatvec``       |
+---------------------------------+--------------------------------+
| ``ModifiedGS``                  | ``SUNModifiedGS``              |
+---------------------------------+--------------------------------+
| ``ClassicalGS``                 | ``SUNClassicalGS``             |
+---------------------------------+--------------------------------+
| ``QRfact``                      | ``SUNQRFact``                  |
+---------------------------------+--------------------------------+
| ``QRsol``                       | ``SUNQRsol``                   |
+---------------------------------+--------------------------------+
| ``DlsMat_NewDenseMat``          | ``SUNDlsMat_NewDenseMat``      |
+---------------------------------+--------------------------------+
| ``DlsMat_NewBandMat``           | ``SUNDlsMat_NewBandMat``       |
+---------------------------------+--------------------------------+
| ``DestroyMat``                  | ``SUNDlsMat_DestroyMat``       |
+---------------------------------+--------------------------------+
| ``NewIntArray``                 | ``SUNDlsMat_NewIntArray``      |
+---------------------------------+--------------------------------+
| ``NewIndexArray``               | ``SUNDlsMat_NewIndexArray``    |
+---------------------------------+--------------------------------+
| ``NewRealArray``                | ``SUNDlsMat_NewRealArray``     |
+---------------------------------+--------------------------------+
| ``DestroyArray``                | ``SUNDlsMat_DestroyArray``     |
+---------------------------------+--------------------------------+
| ``AddIdentity``                 | ``SUNDlsMat_AddIdentity``      |
+---------------------------------+--------------------------------+
| ``SetToZero``                   | ``SUNDlsMat_SetToZero``        |
+---------------------------------+--------------------------------+
| ``PrintMat``                    | ``SUNDlsMat_PrintMat``         |
+---------------------------------+--------------------------------+
| ``newDenseMat``                 | ``SUNDlsMat_newDenseMat``      |
+---------------------------------+--------------------------------+
| ``newBandMat``                  | ``SUNDlsMat_newBandMat``       |
+---------------------------------+--------------------------------+
| ``destroyMat``                  | ``SUNDlsMat_destroyMat``       |
+---------------------------------+--------------------------------+
| ``newIntArray``                 | ``SUNDlsMat_newIntArray``      |
+---------------------------------+--------------------------------+
| ``newIndexArray``               | ``SUNDlsMat_newIndexArray``    |
+---------------------------------+--------------------------------+
| ``newRealArray``                | ``SUNDlsMat_newRealArray``     |
+---------------------------------+--------------------------------+
| ``destroyArray``                | ``SUNDlsMat_destroyArray``     |
+---------------------------------+--------------------------------+

In addition, the entire ``sundials_lapack.h`` header file is now deprecated for
removal in SUNDIALS v7.0.0. Note, this header file is not needed to use the
SUNDIALS LAPACK linear solvers.

### IDA Changes in v6.0.0

**SUNContext**

SUNDIALS v6.0.0 introduces a new :c:type:`SUNContext` object on which all other
SUNDIALS objects depend. As such, the constructors for all SUNDIALS packages,
vectors, matrices, linear solvers, nonlinear solvers, and memory helpers have
been updated to accept a context as the last input. Users upgrading to SUNDIALS
v6.0.0 will need to call :c:func:`SUNContext_Create` to create a context object
with before calling any other SUNDIALS library function, and then provide this
object to other SUNDIALS constructors. The context object has been introduced to
allow SUNDIALS to provide new features, such as the profiling/instrumentation
also introduced in this release, while maintaining thread-safety. See the
documentation section on the :c:type:`SUNContext` for more details.

A script ``upgrade-to-sundials-6-from-5.sh`` has been provided with the release
(obtainable from the GitHub release page) to help ease the transition to
SUNDIALS v6.0.0. The script will add a ``SUNCTX_PLACEHOLDER`` argument to all of
the calls to SUNDIALS constructors that now require a ``SUNContext`` object. It
can also update deprecated SUNDIALS constants/types to the new names. It can be
run like this:

.. code-block::

   > ./upgrade-to-sundials-6-from-5.sh <files to update>

**SUNProfiler**

A capability to profile/instrument SUNDIALS library code has been added. This
can be enabled with the CMake option :cmakeop:`SUNDIALS_BUILD_WITH_PROFILING`. A
built-in profiler will be used by default, but the `Caliper
<https://github.com/LLNL/Caliper>`_ library can also be used instead with the
CMake option :cmakeop:`ENABLE_CALIPER`. See the documentation section on
profiling for more details.  **WARNING**: Profiling will impact performance, and
should be enabled judiciously.

**SUNMemoryHelper**

The :c:type:`SUNMemoryHelper` functions :c:func:`SUNMemoryHelper_Alloc`,
:c:func:`SUNMemoryHelper_Dealloc`, and :c:func:`SUNMemoryHelper_Copy` have been
updated to accept an opaque handle as the last input. At a minimum, user-defined
:c:type:`SUNMemoryHelper` implementations will need to update these functions to
accept the additional argument. Typically, this handle is the execution stream
(e.g., a CUDA/HIP stream or SYCL queue) for the operation. The :ref:`CUDA
<SUNMemory.CUDA>`, :ref:`HIP <SUNMemory.HIP>`, and :ref:`SYCL <SUNMemory.SYCL>`
implementations have been updated accordingly. Additionally, the constructor
:c:func:`SUNMemoryHelper_Sycl` has been updated to remove the SYCL queue as an
input.

**NVector**

Two new optional vector operations, :c:func:`N_VDotProdMultiLocal` and
:c:func:`N_VDotProdMultiAllReduce`, have been added to support
low-synchronization methods for Anderson acceleration.

The CUDA, HIP, and SYCL execution policies have been moved from the ``sundials``
namespace to the ``sundials::cuda``, ``sundials::hip``, and ``sundials::sycl``
namespaces respectively. Accordingly, the prefixes "Cuda", "Hip", and "Sycl"
have been removed from the execution policy classes and methods.

The ``Sundials`` namespace used by the Trilinos Tpetra NVector has been replaced
with the ``sundials::trilinos::nvector_tpetra`` namespace.

The serial, PThreads, PETSc, *hypre*, Parallel, OpenMP_DEV, and OpenMP vector
functions ``N_VCloneVectorArray_*`` and ``N_VDestroyVectorArray_*`` have been
deprecated. The generic :c:func:`N_VCloneVectorArray` and
:c:func:`N_VDestroyVectorArray` functions should be used instead.

The previously deprecated constructor ``N_VMakeWithManagedAllocator_Cuda`` and
the function ``N_VSetCudaStream_Cuda`` have been removed and replaced with
:c:func:`N_VNewWithMemHelp_Cuda` and :c:func:`N_VSetKerrnelExecPolicy_Cuda`
respectively.

The previously deprecated macros ``PVEC_REAL_MPI_TYPE`` and
``PVEC_INTEGER_MPI_TYPE`` have been removed and replaced with
``MPI_SUNREALTYPE`` and ``MPI_SUNINDEXTYPE`` respectively.

**SUNLinearSolver**

The following previously deprecated functions have been removed:

+-----------------------------+------------------------------------------+
| Removed                     | Replacement                              |
+=============================+==========================================+
| ``SUNBandLinearSolver``     | :c:func:`SUNLinSol_Band`                 |
+-----------------------------+------------------------------------------+
| ``SUNDenseLinearSolver``    | :c:func:`SUNLinSol_Dense`                |
+-----------------------------+------------------------------------------+
| ``SUNKLU``                  | :c:func:`SUNLinSol_KLU`                  |
+-----------------------------+------------------------------------------+
| ``SUNKLUReInit``            | :c:func:`SUNLinSol_KLUReInit`            |
+-----------------------------+------------------------------------------+
| ``SUNKLUSetOrdering``       | :c:func:`SUNLinSol_KLUSetOrdering`       |
+-----------------------------+------------------------------------------+
| ``SUNLapackBand``           | :c:func:`SUNLinSol_LapackBand`           |
+-----------------------------+------------------------------------------+
| ``SUNLapackDense``          | :c:func:`SUNLinSol_LapackDense`          |
+-----------------------------+------------------------------------------+
| ``SUNPCG``                  | :c:func:`SUNLinSol_PCG`                  |
+-----------------------------+------------------------------------------+
| ``SUNPCGSetPrecType``       | :c:func:`SUNLinSol_PCGSetPrecType`       |
+-----------------------------+------------------------------------------+
| ``SUNPCGSetMaxl``           | :c:func:`SUNLinSol_PCGSetMaxl`           |
+-----------------------------+------------------------------------------+
| ``SUNSPBCGS``               | :c:func:`SUNLinSol_SPBCGS`               |
+-----------------------------+------------------------------------------+
| ``SUNSPBCGSSetPrecType``    | :c:func:`SUNLinSol_SPBCGSSetPrecType`    |
+-----------------------------+------------------------------------------+
| ``SUNSPBCGSSetMaxl``        | :c:func:`SUNLinSol_SPBCGSSetMaxl`        |
+-----------------------------+------------------------------------------+
| ``SUNSPFGMR``               | :c:func:`SUNLinSol_SPFGMR`               |
+-----------------------------+------------------------------------------+
| ``SUNSPFGMRSetPrecType``    | :c:func:`SUNLinSol_SPFGMRSetPrecType`    |
+-----------------------------+------------------------------------------+
| ``SUNSPFGMRSetGSType``      | :c:func:`SUNLinSol_SPFGMRSetGSType`      |
+-----------------------------+------------------------------------------+
| ``SUNSPFGMRSetMaxRestarts`` | :c:func:`SUNLinSol_SPFGMRSetMaxRestarts` |
+-----------------------------+------------------------------------------+
| ``SUNSPGMR``                | :c:func:`SUNLinSol_SPGMR`                |
+-----------------------------+------------------------------------------+
| ``SUNSPGMRSetPrecType``     | :c:func:`SUNLinSol_SPGMRSetPrecType`     |
+-----------------------------+------------------------------------------+
| ``SUNSPGMRSetGSType``       | :c:func:`SUNLinSol_SPGMRSetGSType`       |
+-----------------------------+------------------------------------------+
| ``SUNSPGMRSetMaxRestarts``  | :c:func:`SUNLinSol_SPGMRSetMaxRestarts`  |
+-----------------------------+------------------------------------------+
| ``SUNSPTFQMR``              | :c:func:`SUNLinSol_SPTFQMR`              |
+-----------------------------+------------------------------------------+
| ``SUNSPTFQMRSetPrecType``   | :c:func:`SUNLinSol_SPTFQMRSetPrecType`   |
+-----------------------------+------------------------------------------+
| ``SUNSPTFQMRSetMaxl``       | :c:func:`SUNLinSol_SPTFQMRSetMaxl`       |
+-----------------------------+------------------------------------------+
| ``SUNSuperLUMT``            | :c:func:`SUNLinSol_SuperLUMT`            |
+-----------------------------+------------------------------------------+
| ``SUNSuperLUMTSetOrdering`` | :c:func:`SUNLinSol_SuperLUMTSetOrdering` |
+-----------------------------+------------------------------------------+

**IDA**

The IDA Fortran 77 interface has been removed. See :numref:`SUNDIALS.Fortran`
and the F2003 example programs for more details using the SUNDIALS Fortran 2003
module interfaces.

**Deprecations**

In addition to the deprecations noted elsewhere, many constants, types, and
functions have been renamed so that they are properly namespaced. The old names
have been deprecated and will be removed in SUNDIALS v7.0.0.

The following constants, macros, and typedefs are now deprecated:

+------------------------------+-------------------------------------+
| Deprecated Name              | New Name                            |
+==============================+=====================================+
| ``realtype``                 | ``sunrealtype``                     |
+------------------------------+-------------------------------------+
| ``booleantype``              | ``sunbooleantype``                  |
+------------------------------+-------------------------------------+
| ``RCONST``                   | ``SUN_RCONST``                      |
+------------------------------+-------------------------------------+
| ``BIG_REAL``                 | ``SUN_BIG_REAL``                    |
+------------------------------+-------------------------------------+
| ``SMALL_REAL``               | ``SUN_SMALL_REAL``                  |
+------------------------------+-------------------------------------+
| ``UNIT_ROUNDOFF``            | ``SUN_UNIT_ROUNDOFF``               |
+------------------------------+-------------------------------------+
| ``PREC_NONE``                | ``SUN_PREC_NONE``                   |
+------------------------------+-------------------------------------+
| ``PREC_LEFT``                | ``SUN_PREC_LEFT``                   |
+------------------------------+-------------------------------------+
| ``PREC_RIGHT``               | ``SUN_PREC_RIGHT``                  |
+------------------------------+-------------------------------------+
| ``PREC_BOTH``                | ``SUN_PREC_BOTH``                   |
+------------------------------+-------------------------------------+
| ``MODIFIED_GS``              | ``SUN_MODIFIED_GS``                 |
+------------------------------+-------------------------------------+
| ``CLASSICAL_GS``             | ``SUN_CLASSICAL_GS``                |
+------------------------------+-------------------------------------+
| ``ATimesFn``                 | ``SUNATimesFn``                     |
+------------------------------+-------------------------------------+
| ``PSetupFn``                 | ``SUNPSetupFn``                     |
+------------------------------+-------------------------------------+
| ``PSolveFn``                 | ``SUNPSolveFn``                     |
+------------------------------+-------------------------------------+
| ``DlsMat``                   | ``SUNDlsMat``                       |
+------------------------------+-------------------------------------+
| ``DENSE_COL``                | ``SUNDLS_DENSE_COL``                |
+------------------------------+-------------------------------------+
| ``DENSE_ELEM``               | ``SUNDLS_DENSE_ELEM``               |
+------------------------------+-------------------------------------+
| ``BAND_COL``                 | ``SUNDLS_BAND_COL``                 |
+------------------------------+-------------------------------------+
| ``BAND_COL_ELEM``            | ``SUNDLS_BAND_COL_ELEM``            |
+------------------------------+-------------------------------------+
| ``BAND_ELEM``                | ``SUNDLS_BAND_ELEM``                |
+------------------------------+-------------------------------------+

In addition, the following functions are now deprecated (compile-time warnings
will be thrown if supported by the compiler):

+---------------------------------+--------------------------------+
| Deprecated Name                 | New Name                       |
+=================================+================================+
| ``IDASpilsSetLinearSolver``     | ``IDASetLinearSolver``         |
+---------------------------------+--------------------------------+
| ``IDASpilsSetPreconditioner``   | ``IDASetPreconditioner``       |
+---------------------------------+--------------------------------+
| ``IDASpilsSetJacTimes``         | ``IDASetJacTimes``             |
+---------------------------------+--------------------------------+
| ``IDASpilsSetEpsLin``           | ``IDASetEpsLin``               |
+---------------------------------+--------------------------------+
| ``IDASpilsSetIncrementFactor``  | ``IDASetIncrementFactor``      |
+---------------------------------+--------------------------------+
| ``IDASpilsGetWorkSpace``        | ``IDAGetLinWorkSpace``         |
+---------------------------------+--------------------------------+
| ``IDASpilsGetNumPrecEvals``     | ``IDAGetNumPrecEvals``         |
+---------------------------------+--------------------------------+
| ``IDASpilsGetNumPrecSolves``    | ``IDAGetNumPrecSolves``        |
+---------------------------------+--------------------------------+
| ``IDASpilsGetNumLinIters``      | ``IDAGetNumLinIters``          |
+---------------------------------+--------------------------------+
| ``IDASpilsGetNumConvFails``     | ``IDAGetNumLinConvFails``      |
+---------------------------------+--------------------------------+
| ``IDASpilsGetNumJTSetupEvals``  | ``IDAGetNumJTSetupEvals``      |
+---------------------------------+--------------------------------+
| ``IDASpilsGetNumJtimesEvals``   | ``IDAGetNumJtimesEvals``       |
+---------------------------------+--------------------------------+
| ``IDASpilsGetNumResEvals``      | ``IDAGetNumLinResEvals``       |
+---------------------------------+--------------------------------+
| ``IDASpilsGetLastFlag``         | ``IDAGetLastLinFlag``          |
+---------------------------------+--------------------------------+
| ``IDASpilsGetReturnFlagName``   | ``IDAGetLinReturnFlagName``    |
+---------------------------------+--------------------------------+
| ``IDADlsSetLinearSolver``       | ``IDASetLinearSolver``         |
+---------------------------------+--------------------------------+
| ``IDADlsSetJacFn``              | ``IDASetJacFn``                |
+---------------------------------+--------------------------------+
| ``IDADlsGetWorkSpace``          | ``IDAGetLinWorkSpace``         |
+---------------------------------+--------------------------------+
| ``IDADlsGetNumJacEvals``        | ``IDAGetNumJacEvals``          |
+---------------------------------+--------------------------------+
| ``IDADlsGetNumResEvals``        | ``IDAGetNumLinResEvals``       |
+---------------------------------+--------------------------------+
| ``IDADlsGetLastFlag``           | ``IDAGetLastLinFlag``          |
+---------------------------------+--------------------------------+
| ``IDADlsGetReturnFlagName``     | ``IDAGetLinReturnFlagName``    |
+---------------------------------+--------------------------------+
| ``DenseGETRF``                  | ``SUNDlsMat_DenseGETRF``       |
+---------------------------------+--------------------------------+
| ``DenseGETRS``                  | ``SUNDlsMat_DenseGETRS``       |
+---------------------------------+--------------------------------+
| ``denseGETRF``                  | ``SUNDlsMat_denseGETRF``       |
+---------------------------------+--------------------------------+
| ``denseGETRS``                  | ``SUNDlsMat_denseGETRS``       |
+---------------------------------+--------------------------------+
| ``DensePOTRF``                  | ``SUNDlsMat_DensePOTRF``       |
+---------------------------------+--------------------------------+
| ``DensePOTRS``                  | ``SUNDlsMat_DensePOTRS``       |
+---------------------------------+--------------------------------+
| ``densePOTRF``                  | ``SUNDlsMat_densePOTRF``       |
+---------------------------------+--------------------------------+
| ``densePOTRS``                  | ``SUNDlsMat_densePOTRS``       |
+---------------------------------+--------------------------------+
| ``DenseGEQRF``                  | ``SUNDlsMat_DenseGEQRF``       |
+---------------------------------+--------------------------------+
| ``DenseORMQR``                  | ``SUNDlsMat_DenseORMQR``       |
+---------------------------------+--------------------------------+
| ``denseGEQRF``                  | ``SUNDlsMat_denseGEQRF``       |
+---------------------------------+--------------------------------+
| ``denseORMQR``                  | ``SUNDlsMat_denseORMQR``       |
+---------------------------------+--------------------------------+
| ``DenseCopy``                   | ``SUNDlsMat_DenseCopy``        |
+---------------------------------+--------------------------------+
| ``denseCopy``                   | ``SUNDlsMat_denseCopy``        |
+---------------------------------+--------------------------------+
| ``DenseScale``                  | ``SUNDlsMat_DenseScale``       |
+---------------------------------+--------------------------------+
| ``denseScale``                  | ``SUNDlsMat_denseScale``       |
+---------------------------------+--------------------------------+
| ``denseAddIdentity``            | ``SUNDlsMat_denseAddIdentity`` |
+---------------------------------+--------------------------------+
| ``DenseMatvec``                 | ``SUNDlsMat_DenseMatvec``      |
+---------------------------------+--------------------------------+
| ``denseMatvec``                 | ``SUNDlsMat_denseMatvec``      |
+---------------------------------+--------------------------------+
| ``BandGBTRF``                   | ``SUNDlsMat_BandGBTRF``        |
+---------------------------------+--------------------------------+
| ``bandGBTRF``                   | ``SUNDlsMat_bandGBTRF``        |
+---------------------------------+--------------------------------+
| ``BandGBTRS``                   | ``SUNDlsMat_BandGBTRS``        |
+---------------------------------+--------------------------------+
| ``bandGBTRS``                   | ``SUNDlsMat_bandGBTRS``        |
+---------------------------------+--------------------------------+
| ``BandCopy``                    | ``SUNDlsMat_BandCopy``         |
+---------------------------------+--------------------------------+
| ``bandCopy``                    | ``SUNDlsMat_bandCopy``         |
+---------------------------------+--------------------------------+
| ``BandScale``                   | ``SUNDlsMat_BandScale``        |
+---------------------------------+--------------------------------+
| ``bandScale``                   | ``SUNDlsMat_bandScale``        |
+---------------------------------+--------------------------------+
| ``bandAddIdentity``             | ``SUNDlsMat_bandAddIdentity``  |
+---------------------------------+--------------------------------+
| ``BandMatvec``                  | ``SUNDlsMat_BandMatvec``       |
+---------------------------------+--------------------------------+
| ``bandMatvec``                  | ``SUNDlsMat_bandMatvec``       |
+---------------------------------+--------------------------------+
| ``ModifiedGS``                  | ``SUNModifiedGS``              |
+---------------------------------+--------------------------------+
| ``ClassicalGS``                 | ``SUNClassicalGS``             |
+---------------------------------+--------------------------------+
| ``QRfact``                      | ``SUNQRFact``                  |
+---------------------------------+--------------------------------+
| ``QRsol``                       | ``SUNQRsol``                   |
+---------------------------------+--------------------------------+
| ``DlsMat_NewDenseMat``          | ``SUNDlsMat_NewDenseMat``      |
+---------------------------------+--------------------------------+
| ``DlsMat_NewBandMat``           | ``SUNDlsMat_NewBandMat``       |
+---------------------------------+--------------------------------+
| ``DestroyMat``                  | ``SUNDlsMat_DestroyMat``       |
+---------------------------------+--------------------------------+
| ``NewIntArray``                 | ``SUNDlsMat_NewIntArray``      |
+---------------------------------+--------------------------------+
| ``NewIndexArray``               | ``SUNDlsMat_NewIndexArray``    |
+---------------------------------+--------------------------------+
| ``NewRealArray``                | ``SUNDlsMat_NewRealArray``     |
+---------------------------------+--------------------------------+
| ``DestroyArray``                | ``SUNDlsMat_DestroyArray``     |
+---------------------------------+--------------------------------+
| ``AddIdentity``                 | ``SUNDlsMat_AddIdentity``      |
+---------------------------------+--------------------------------+
| ``SetToZero``                   | ``SUNDlsMat_SetToZero``        |
+---------------------------------+--------------------------------+
| ``PrintMat``                    | ``SUNDlsMat_PrintMat``         |
+---------------------------------+--------------------------------+
| ``newDenseMat``                 | ``SUNDlsMat_newDenseMat``      |
+---------------------------------+--------------------------------+
| ``newBandMat``                  | ``SUNDlsMat_newBandMat``       |
+---------------------------------+--------------------------------+
| ``destroyMat``                  | ``SUNDlsMat_destroyMat``       |
+---------------------------------+--------------------------------+
| ``newIntArray``                 | ``SUNDlsMat_newIntArray``      |
+---------------------------------+--------------------------------+
| ``newIndexArray``               | ``SUNDlsMat_newIndexArray``    |
+---------------------------------+--------------------------------+
| ``newRealArray``                | ``SUNDlsMat_newRealArray``     |
+---------------------------------+--------------------------------+
| ``destroyArray``                | ``SUNDlsMat_destroyArray``     |
+---------------------------------+--------------------------------+

In addition, the entire ``sundials_lapack.h`` header file is now deprecated for
removal in SUNDIALS v7.0.0. Note, this header file is not needed to use the
SUNDIALS LAPACK linear solvers.

### IDAS Changes in v5.0.0

**SUNContext**

SUNDIALS v6.0.0 introduces a new :c:type:`SUNContext` object on which all other
SUNDIALS objects depend. As such, the constructors for all SUNDIALS packages,
vectors, matrices, linear solvers, nonlinear solvers, and memory helpers have
been updated to accept a context as the last input. Users upgrading to SUNDIALS
v6.0.0 will need to call :c:func:`SUNContext_Create` to create a context object
with before calling any other SUNDIALS library function, and then provide this
object to other SUNDIALS constructors. The context object has been introduced to
allow SUNDIALS to provide new features, such as the profiling/instrumentation
also introduced in this release, while maintaining thread-safety. See the
documentation section on the :c:type:`SUNContext` for more details.

A script ``upgrade-to-sundials-6-from-5.sh`` has been provided with the release
(obtainable from the GitHub release page) to help ease the transition to
SUNDIALS v6.0.0. The script will add a ``SUNCTX_PLACEHOLDER`` argument to all of
the calls to SUNDIALS constructors that now require a ``SUNContext`` object. It
can also update deprecated SUNDIALS constants/types to the new names. It can be
run like this:

.. code-block::

   > ./upgrade-to-sundials-6-from-5.sh <files to update>

**SUNProfiler**

A capability to profile/instrument SUNDIALS library code has been added. This
can be enabled with the CMake option :cmakeop:`SUNDIALS_BUILD_WITH_PROFILING`. A
built-in profiler will be used by default, but the `Caliper
<https://github.com/LLNL/Caliper>`_ library can also be used instead with the
CMake option :cmakeop:`ENABLE_CALIPER`. See the documentation section on
profiling for more details.  **WARNING**: Profiling will impact performance, and
should be enabled judiciously.

**SUNMemoryHelper**

The :c:type:`SUNMemoryHelper` functions :c:func:`SUNMemoryHelper_Alloc`,
:c:func:`SUNMemoryHelper_Dealloc`, and :c:func:`SUNMemoryHelper_Copy` have been
updated to accept an opaque handle as the last input. At a minimum, user-defined
:c:type:`SUNMemoryHelper` implementations will need to update these functions to
accept the additional argument. Typically, this handle is the execution stream
(e.g., a CUDA/HIP stream or SYCL queue) for the operation. The :ref:`CUDA
<SUNMemory.CUDA>`, :ref:`HIP <SUNMemory.HIP>`, and :ref:`SYCL <SUNMemory.SYCL>`
implementations have been updated accordingly. Additionally, the constructor
:c:func:`SUNMemoryHelper_Sycl` has been updated to remove the SYCL queue as an
input.

**NVector**

Two new optional vector operations, :c:func:`N_VDotProdMultiLocal` and
:c:func:`N_VDotProdMultiAllReduce`, have been added to support
low-synchronization methods for Anderson acceleration.

The CUDA, HIP, and SYCL execution policies have been moved from the ``sundials``
namespace to the ``sundials::cuda``, ``sundials::hip``, and ``sundials::sycl``
namespaces respectively. Accordingly, the prefixes "Cuda", "Hip", and "Sycl"
have been removed from the execution policy classes and methods.

The ``Sundials`` namespace used by the Trilinos Tpetra NVector has been replaced
with the ``sundials::trilinos::nvector_tpetra`` namespace.

The serial, PThreads, PETSc, *hypre*, Parallel, OpenMP_DEV, and OpenMP vector
functions ``N_VCloneVectorArray_*`` and ``N_VDestroyVectorArray_*`` have been
deprecated. The generic :c:func:`N_VCloneVectorArray` and
:c:func:`N_VDestroyVectorArray` functions should be used instead.

The previously deprecated constructor ``N_VMakeWithManagedAllocator_Cuda`` and
the function ``N_VSetCudaStream_Cuda`` have been removed and replaced with
:c:func:`N_VNewWithMemHelp_Cuda` and :c:func:`N_VSetKerrnelExecPolicy_Cuda`
respectively.

The previously deprecated macros ``PVEC_REAL_MPI_TYPE`` and
``PVEC_INTEGER_MPI_TYPE`` have been removed and replaced with
``MPI_SUNREALTYPE`` and ``MPI_SUNINDEXTYPE`` respectively.

**SUNLinearSolver**

The following previously deprecated functions have been removed:

+-----------------------------+------------------------------------------+
| Removed                     | Replacement                              |
+=============================+==========================================+
| ``SUNBandLinearSolver``     | :c:func:`SUNLinSol_Band`                 |
+-----------------------------+------------------------------------------+
| ``SUNDenseLinearSolver``    | :c:func:`SUNLinSol_Dense`                |
+-----------------------------+------------------------------------------+
| ``SUNKLU``                  | :c:func:`SUNLinSol_KLU`                  |
+-----------------------------+------------------------------------------+
| ``SUNKLUReInit``            | :c:func:`SUNLinSol_KLUReInit`            |
+-----------------------------+------------------------------------------+
| ``SUNKLUSetOrdering``       | :c:func:`SUNLinSol_KLUSetOrdering`       |
+-----------------------------+------------------------------------------+
| ``SUNLapackBand``           | :c:func:`SUNLinSol_LapackBand`           |
+-----------------------------+------------------------------------------+
| ``SUNLapackDense``          | :c:func:`SUNLinSol_LapackDense`          |
+-----------------------------+------------------------------------------+
| ``SUNPCG``                  | :c:func:`SUNLinSol_PCG`                  |
+-----------------------------+------------------------------------------+
| ``SUNPCGSetPrecType``       | :c:func:`SUNLinSol_PCGSetPrecType`       |
+-----------------------------+------------------------------------------+
| ``SUNPCGSetMaxl``           | :c:func:`SUNLinSol_PCGSetMaxl`           |
+-----------------------------+------------------------------------------+
| ``SUNSPBCGS``               | :c:func:`SUNLinSol_SPBCGS`               |
+-----------------------------+------------------------------------------+
| ``SUNSPBCGSSetPrecType``    | :c:func:`SUNLinSol_SPBCGSSetPrecType`    |
+-----------------------------+------------------------------------------+
| ``SUNSPBCGSSetMaxl``        | :c:func:`SUNLinSol_SPBCGSSetMaxl`        |
+-----------------------------+------------------------------------------+
| ``SUNSPFGMR``               | :c:func:`SUNLinSol_SPFGMR`               |
+-----------------------------+------------------------------------------+
| ``SUNSPFGMRSetPrecType``    | :c:func:`SUNLinSol_SPFGMRSetPrecType`    |
+-----------------------------+------------------------------------------+
| ``SUNSPFGMRSetGSType``      | :c:func:`SUNLinSol_SPFGMRSetGSType`      |
+-----------------------------+------------------------------------------+
| ``SUNSPFGMRSetMaxRestarts`` | :c:func:`SUNLinSol_SPFGMRSetMaxRestarts` |
+-----------------------------+------------------------------------------+
| ``SUNSPGMR``                | :c:func:`SUNLinSol_SPGMR`                |
+-----------------------------+------------------------------------------+
| ``SUNSPGMRSetPrecType``     | :c:func:`SUNLinSol_SPGMRSetPrecType`     |
+-----------------------------+------------------------------------------+
| ``SUNSPGMRSetGSType``       | :c:func:`SUNLinSol_SPGMRSetGSType`       |
+-----------------------------+------------------------------------------+
| ``SUNSPGMRSetMaxRestarts``  | :c:func:`SUNLinSol_SPGMRSetMaxRestarts`  |
+-----------------------------+------------------------------------------+
| ``SUNSPTFQMR``              | :c:func:`SUNLinSol_SPTFQMR`              |
+-----------------------------+------------------------------------------+
| ``SUNSPTFQMRSetPrecType``   | :c:func:`SUNLinSol_SPTFQMRSetPrecType`   |
+-----------------------------+------------------------------------------+
| ``SUNSPTFQMRSetMaxl``       | :c:func:`SUNLinSol_SPTFQMRSetMaxl`       |
+-----------------------------+------------------------------------------+
| ``SUNSuperLUMT``            | :c:func:`SUNLinSol_SuperLUMT`            |
+-----------------------------+------------------------------------------+
| ``SUNSuperLUMTSetOrdering`` | :c:func:`SUNLinSol_SuperLUMTSetOrdering` |
+-----------------------------+------------------------------------------+

**Deprecations**

In addition to the deprecations noted elsewhere, many constants, types, and
functions have been renamed so that they are properly namespaced. The old names
have been deprecated and will be removed in SUNDIALS v7.0.0.

The following constants, macros, and typedefs are now deprecated:

+------------------------------+-------------------------------------+
| Deprecated Name              | New Name                            |
+==============================+=====================================+
| ``realtype``                 | ``sunrealtype``                     |
+------------------------------+-------------------------------------+
| ``booleantype``              | ``sunbooleantype``                  |
+------------------------------+-------------------------------------+
| ``RCONST``                   | ``SUN_RCONST``                      |
+------------------------------+-------------------------------------+
| ``BIG_REAL``                 | ``SUN_BIG_REAL``                    |
+------------------------------+-------------------------------------+
| ``SMALL_REAL``               | ``SUN_SMALL_REAL``                  |
+------------------------------+-------------------------------------+
| ``UNIT_ROUNDOFF``            | ``SUN_UNIT_ROUNDOFF``               |
+------------------------------+-------------------------------------+
| ``PREC_NONE``                | ``SUN_PREC_NONE``                   |
+------------------------------+-------------------------------------+
| ``PREC_LEFT``                | ``SUN_PREC_LEFT``                   |
+------------------------------+-------------------------------------+
| ``PREC_RIGHT``               | ``SUN_PREC_RIGHT``                  |
+------------------------------+-------------------------------------+
| ``PREC_BOTH``                | ``SUN_PREC_BOTH``                   |
+------------------------------+-------------------------------------+
| ``MODIFIED_GS``              | ``SUN_MODIFIED_GS``                 |
+------------------------------+-------------------------------------+
| ``CLASSICAL_GS``             | ``SUN_CLASSICAL_GS``                |
+------------------------------+-------------------------------------+
| ``ATimesFn``                 | ``SUNATimesFn``                     |
+------------------------------+-------------------------------------+
| ``PSetupFn``                 | ``SUNPSetupFn``                     |
+------------------------------+-------------------------------------+
| ``PSolveFn``                 | ``SUNPSolveFn``                     |
+------------------------------+-------------------------------------+
| ``DlsMat``                   | ``SUNDlsMat``                       |
+------------------------------+-------------------------------------+
| ``DENSE_COL``                | ``SUNDLS_DENSE_COL``                |
+------------------------------+-------------------------------------+
| ``DENSE_ELEM``               | ``SUNDLS_DENSE_ELEM``               |
+------------------------------+-------------------------------------+
| ``BAND_COL``                 | ``SUNDLS_BAND_COL``                 |
+------------------------------+-------------------------------------+
| ``BAND_COL_ELEM``            | ``SUNDLS_BAND_COL_ELEM``            |
+------------------------------+-------------------------------------+
| ``BAND_ELEM``                | ``SUNDLS_BAND_ELEM``                |
+------------------------------+-------------------------------------+

In addition, the following functions are now deprecated (compile-time warnings
will be thrown if supported by the compiler):

+---------------------------------+--------------------------------+
| Deprecated Name                 | New Name                       |
+=================================+================================+
| ``IDASpilsSetLinearSolver``     | ``IDASetLinearSolver``         |
+---------------------------------+--------------------------------+
| ``IDASpilsSetPreconditioner``   | ``IDASetPreconditioner``       |
+---------------------------------+--------------------------------+
| ``IDASpilsSetJacTimes``         | ``IDASetJacTimes``             |
+---------------------------------+--------------------------------+
| ``IDASpilsSetEpsLin``           | ``IDASetEpsLin``               |
+---------------------------------+--------------------------------+
| ``IDASpilsSetIncrementFactor``  | ``IDASetIncrementFactor``      |
+---------------------------------+--------------------------------+
| ``IDASpilsGetWorkSpace``        | ``IDAGetLinWorkSpace``         |
+---------------------------------+--------------------------------+
| ``IDASpilsGetNumPrecEvals``     | ``IDAGetNumPrecEvals``         |
+---------------------------------+--------------------------------+
| ``IDASpilsGetNumPrecSolves``    | ``IDAGetNumPrecSolves``        |
+---------------------------------+--------------------------------+
| ``IDASpilsGetNumLinIters``      | ``IDAGetNumLinIters``          |
+---------------------------------+--------------------------------+
| ``IDASpilsGetNumConvFails``     | ``IDAGetNumLinConvFails``      |
+---------------------------------+--------------------------------+
| ``IDASpilsGetNumJTSetupEvals``  | ``IDAGetNumJTSetupEvals``      |
+---------------------------------+--------------------------------+
| ``IDASpilsGetNumJtimesEvals``   | ``IDAGetNumJtimesEvals``       |
+---------------------------------+--------------------------------+
| ``IDASpilsGetNumResEvals``      | ``IDAGetNumLinResEvals``       |
+---------------------------------+--------------------------------+
| ``IDASpilsGetLastFlag``         | ``IDAGetLastLinFlag``          |
+---------------------------------+--------------------------------+
| ``IDASpilsGetReturnFlagName``   | ``IDAGetLinReturnFlagName``    |
+---------------------------------+--------------------------------+
| ``IDASpilsSetLinearSolverB``    | ``IDASetLinearSolverB``        |
+---------------------------------+--------------------------------+
| ``IDASpilsSetEpsLinB``          | ``IDASetEpsLinB``              |
+---------------------------------+--------------------------------+
| ``IDASpilsSetIncrementFactorB`` | ``IDASetIncrementFactorB``     |
+---------------------------------+--------------------------------+
| ``IDASpilsSetPreconditionerB``  | ``IDASetPreconditionerB``      |
+---------------------------------+--------------------------------+
| ``IDASpilsSetPreconditionerBS`` | ``IDASetPreconditionerBS``     |
+---------------------------------+--------------------------------+
| ``IDASpilsSetJacTimesB``        | ``IDASetJacTimesB``            |
+---------------------------------+--------------------------------+
| ``IDASpilsSetJacTimesBS``       | ``IDASetJacTimesBS``           |
+---------------------------------+--------------------------------+
| ``IDADlsSetLinearSolver``       | ``IDASetLinearSolver``         |
+---------------------------------+--------------------------------+
| ``IDADlsSetJacFn``              | ``IDASetJacFn``                |
+---------------------------------+--------------------------------+
| ``IDADlsGetWorkSpace``          | ``IDAGetLinWorkSpace``         |
+---------------------------------+--------------------------------+
| ``IDADlsGetNumJacEvals``        | ``IDAGetNumJacEvals``          |
+---------------------------------+--------------------------------+
| ``IDADlsGetNumResEvals``        | ``IDAGetNumLinResEvals``       |
+---------------------------------+--------------------------------+
| ``IDADlsGetLastFlag``           | ``IDAGetLastLinFlag``          |
+---------------------------------+--------------------------------+
| ``IDADlsGetReturnFlagName``     | ``IDAGetLinReturnFlagName``    |
+---------------------------------+--------------------------------+
| ``IDADlsSetLinearSolverB``      | ``IDASetLinearSolverB``        |
+---------------------------------+--------------------------------+
| ``IDADlsSetJacFnB``             | ``IDASetJacFnB``               |
+---------------------------------+--------------------------------+
| ``IDADlsSetJacFnBS``            | ``IDASetJacFnBS``              |
+---------------------------------+--------------------------------+
| ``DenseGETRF``                  | ``SUNDlsMat_DenseGETRF``       |
+---------------------------------+--------------------------------+
| ``DenseGETRS``                  | ``SUNDlsMat_DenseGETRS``       |
+---------------------------------+--------------------------------+
| ``denseGETRF``                  | ``SUNDlsMat_denseGETRF``       |
+---------------------------------+--------------------------------+
| ``denseGETRS``                  | ``SUNDlsMat_denseGETRS``       |
+---------------------------------+--------------------------------+
| ``DensePOTRF``                  | ``SUNDlsMat_DensePOTRF``       |
+---------------------------------+--------------------------------+
| ``DensePOTRS``                  | ``SUNDlsMat_DensePOTRS``       |
+---------------------------------+--------------------------------+
| ``densePOTRF``                  | ``SUNDlsMat_densePOTRF``       |
+---------------------------------+--------------------------------+
| ``densePOTRS``                  | ``SUNDlsMat_densePOTRS``       |
+---------------------------------+--------------------------------+
| ``DenseGEQRF``                  | ``SUNDlsMat_DenseGEQRF``       |
+---------------------------------+--------------------------------+
| ``DenseORMQR``                  | ``SUNDlsMat_DenseORMQR``       |
+---------------------------------+--------------------------------+
| ``denseGEQRF``                  | ``SUNDlsMat_denseGEQRF``       |
+---------------------------------+--------------------------------+
| ``denseORMQR``                  | ``SUNDlsMat_denseORMQR``       |
+---------------------------------+--------------------------------+
| ``DenseCopy``                   | ``SUNDlsMat_DenseCopy``        |
+---------------------------------+--------------------------------+
| ``denseCopy``                   | ``SUNDlsMat_denseCopy``        |
+---------------------------------+--------------------------------+
| ``DenseScale``                  | ``SUNDlsMat_DenseScale``       |
+---------------------------------+--------------------------------+
| ``denseScale``                  | ``SUNDlsMat_denseScale``       |
+---------------------------------+--------------------------------+
| ``denseAddIdentity``            | ``SUNDlsMat_denseAddIdentity`` |
+---------------------------------+--------------------------------+
| ``DenseMatvec``                 | ``SUNDlsMat_DenseMatvec``      |
+---------------------------------+--------------------------------+
| ``denseMatvec``                 | ``SUNDlsMat_denseMatvec``      |
+---------------------------------+--------------------------------+
| ``BandGBTRF``                   | ``SUNDlsMat_BandGBTRF``        |
+---------------------------------+--------------------------------+
| ``bandGBTRF``                   | ``SUNDlsMat_bandGBTRF``        |
+---------------------------------+--------------------------------+
| ``BandGBTRS``                   | ``SUNDlsMat_BandGBTRS``        |
+---------------------------------+--------------------------------+
| ``bandGBTRS``                   | ``SUNDlsMat_bandGBTRS``        |
+---------------------------------+--------------------------------+
| ``BandCopy``                    | ``SUNDlsMat_BandCopy``         |
+---------------------------------+--------------------------------+
| ``bandCopy``                    | ``SUNDlsMat_bandCopy``         |
+---------------------------------+--------------------------------+
| ``BandScale``                   | ``SUNDlsMat_BandScale``        |
+---------------------------------+--------------------------------+
| ``bandScale``                   | ``SUNDlsMat_bandScale``        |
+---------------------------------+--------------------------------+
| ``bandAddIdentity``             | ``SUNDlsMat_bandAddIdentity``  |
+---------------------------------+--------------------------------+
| ``BandMatvec``                  | ``SUNDlsMat_BandMatvec``       |
+---------------------------------+--------------------------------+
| ``bandMatvec``                  | ``SUNDlsMat_bandMatvec``       |
+---------------------------------+--------------------------------+
| ``ModifiedGS``                  | ``SUNModifiedGS``              |
+---------------------------------+--------------------------------+
| ``ClassicalGS``                 | ``SUNClassicalGS``             |
+---------------------------------+--------------------------------+
| ``QRfact``                      | ``SUNQRFact``                  |
+---------------------------------+--------------------------------+
| ``QRsol``                       | ``SUNQRsol``                   |
+---------------------------------+--------------------------------+
| ``DlsMat_NewDenseMat``          | ``SUNDlsMat_NewDenseMat``      |
+---------------------------------+--------------------------------+
| ``DlsMat_NewBandMat``           | ``SUNDlsMat_NewBandMat``       |
+---------------------------------+--------------------------------+
| ``DestroyMat``                  | ``SUNDlsMat_DestroyMat``       |
+---------------------------------+--------------------------------+
| ``NewIntArray``                 | ``SUNDlsMat_NewIntArray``      |
+---------------------------------+--------------------------------+
| ``NewIndexArray``               | ``SUNDlsMat_NewIndexArray``    |
+---------------------------------+--------------------------------+
| ``NewRealArray``                | ``SUNDlsMat_NewRealArray``     |
+---------------------------------+--------------------------------+
| ``DestroyArray``                | ``SUNDlsMat_DestroyArray``     |
+---------------------------------+--------------------------------+
| ``AddIdentity``                 | ``SUNDlsMat_AddIdentity``      |
+---------------------------------+--------------------------------+
| ``SetToZero``                   | ``SUNDlsMat_SetToZero``        |
+---------------------------------+--------------------------------+
| ``PrintMat``                    | ``SUNDlsMat_PrintMat``         |
+---------------------------------+--------------------------------+
| ``newDenseMat``                 | ``SUNDlsMat_newDenseMat``      |
+---------------------------------+--------------------------------+
| ``newBandMat``                  | ``SUNDlsMat_newBandMat``       |
+---------------------------------+--------------------------------+
| ``destroyMat``                  | ``SUNDlsMat_destroyMat``       |
+---------------------------------+--------------------------------+
| ``newIntArray``                 | ``SUNDlsMat_newIntArray``      |
+---------------------------------+--------------------------------+
| ``newIndexArray``               | ``SUNDlsMat_newIndexArray``    |
+---------------------------------+--------------------------------+
| ``newRealArray``                | ``SUNDlsMat_newRealArray``     |
+---------------------------------+--------------------------------+
| ``destroyArray``                | ``SUNDlsMat_destroyArray``     |
+---------------------------------+--------------------------------+

In addition, the entire ``sundials_lapack.h`` header file is now deprecated for
removal in SUNDIALS v7.0.0. Note, this header file is not needed to use the
SUNDIALS LAPACK linear solvers.

### KINSOL Changes in v6.0.0

**SUNContext**

SUNDIALS v6.0.0 introduces a new :c:type:`SUNContext` object on which all other
SUNDIALS objects depend. As such, the constructors for all SUNDIALS packages,
vectors, matrices, linear solvers, nonlinear solvers, and memory helpers have
been updated to accept a context as the last input. Users upgrading to SUNDIALS
v6.0.0 will need to call :c:func:`SUNContext_Create` to create a context object
with before calling any other SUNDIALS library function, and then provide this
object to other SUNDIALS constructors. The context object has been introduced to
allow SUNDIALS to provide new features, such as the profiling/instrumentation
also introduced in this release, while maintaining thread-safety. See the
documentation section on the :c:type:`SUNContext` for more details.

A script ``upgrade-to-sundials-6-from-5.sh`` has been provided with the release
(obtainable from the GitHub release page) to help ease the transition to
SUNDIALS v6.0.0. The script will add a ``SUNCTX_PLACEHOLDER`` argument to all of
the calls to SUNDIALS constructors that now require a ``SUNContext`` object. It
can also update deprecated SUNDIALS constants/types to the new names. It can be
run like this:

.. code-block::

   > ./upgrade-to-sundials-6-from-5.sh <files to update>

**SUNProfiler**

A capability to profile/instrument SUNDIALS library code has been added. This
can be enabled with the CMake option :cmakeop:`SUNDIALS_BUILD_WITH_PROFILING`. A
built-in profiler will be used by default, but the `Caliper
<https://github.com/LLNL/Caliper>`_ library can also be used instead with the
CMake option :cmakeop:`ENABLE_CALIPER`. See the documentation section on
profiling for more details.  **WARNING**: Profiling will impact performance, and
should be enabled judiciously.

**SUNMemoryHelper**

The :c:type:`SUNMemoryHelper` functions :c:func:`SUNMemoryHelper_Alloc`,
:c:func:`SUNMemoryHelper_Dealloc`, and :c:func:`SUNMemoryHelper_Copy` have been
updated to accept an opaque handle as the last input. At a minimum, user-defined
:c:type:`SUNMemoryHelper` implementations will need to update these functions to
accept the additional argument. Typically, this handle is the execution stream
(e.g., a CUDA/HIP stream or SYCL queue) for the operation. The :ref:`CUDA
<SUNMemory.CUDA>`, :ref:`HIP <SUNMemory.HIP>`, and :ref:`SYCL <SUNMemory.SYCL>`
implementations have been updated accordingly. Additionally, the constructor
:c:func:`SUNMemoryHelper_Sycl` has been updated to remove the SYCL queue as an
input.

**NVector**

Two new optional vector operations, :c:func:`N_VDotProdMultiLocal` and
:c:func:`N_VDotProdMultiAllReduce`, have been added to support
low-synchronization methods for Anderson acceleration.

The CUDA, HIP, and SYCL execution policies have been moved from the ``sundials``
namespace to the ``sundials::cuda``, ``sundials::hip``, and ``sundials::sycl``
namespaces respectively. Accordingly, the prefixes "Cuda", "Hip", and "Sycl"
have been removed from the execution policy classes and methods.

The ``Sundials`` namespace used by the Trilinos Tpetra NVector has been replaced
with the ``sundials::trilinos::nvector_tpetra`` namespace.

The serial, PThreads, PETSc, *hypre*, Parallel, OpenMP_DEV, and OpenMP vector
functions ``N_VCloneVectorArray_*`` and ``N_VDestroyVectorArray_*`` have been
deprecated. The generic :c:func:`N_VCloneVectorArray` and
:c:func:`N_VDestroyVectorArray` functions should be used instead.

The previously deprecated constructor ``N_VMakeWithManagedAllocator_Cuda`` and
the function ``N_VSetCudaStream_Cuda`` have been removed and replaced with
:c:func:`N_VNewWithMemHelp_Cuda` and :c:func:`N_VSetKerrnelExecPolicy_Cuda`
respectively.

The previously deprecated macros ``PVEC_REAL_MPI_TYPE`` and
``PVEC_INTEGER_MPI_TYPE`` have been removed and replaced with
``MPI_SUNREALTYPE`` and ``MPI_SUNINDEXTYPE`` respectively.

**SUNLinearSolver**

The following previously deprecated functions have been removed:

+-----------------------------+------------------------------------------+
| Removed                     | Replacement                              |
+=============================+==========================================+
| ``SUNBandLinearSolver``     | :c:func:`SUNLinSol_Band`                 |
+-----------------------------+------------------------------------------+
| ``SUNDenseLinearSolver``    | :c:func:`SUNLinSol_Dense`                |
+-----------------------------+------------------------------------------+
| ``SUNKLU``                  | :c:func:`SUNLinSol_KLU`                  |
+-----------------------------+------------------------------------------+
| ``SUNKLUReInit``            | :c:func:`SUNLinSol_KLUReInit`            |
+-----------------------------+------------------------------------------+
| ``SUNKLUSetOrdering``       | :c:func:`SUNLinSol_KLUSetOrdering`       |
+-----------------------------+------------------------------------------+
| ``SUNLapackBand``           | :c:func:`SUNLinSol_LapackBand`           |
+-----------------------------+------------------------------------------+
| ``SUNLapackDense``          | :c:func:`SUNLinSol_LapackDense`          |
+-----------------------------+------------------------------------------+
| ``SUNPCG``                  | :c:func:`SUNLinSol_PCG`                  |
+-----------------------------+------------------------------------------+
| ``SUNPCGSetPrecType``       | :c:func:`SUNLinSol_PCGSetPrecType`       |
+-----------------------------+------------------------------------------+
| ``SUNPCGSetMaxl``           | :c:func:`SUNLinSol_PCGSetMaxl`           |
+-----------------------------+------------------------------------------+
| ``SUNSPBCGS``               | :c:func:`SUNLinSol_SPBCGS`               |
+-----------------------------+------------------------------------------+
| ``SUNSPBCGSSetPrecType``    | :c:func:`SUNLinSol_SPBCGSSetPrecType`    |
+-----------------------------+------------------------------------------+
| ``SUNSPBCGSSetMaxl``        | :c:func:`SUNLinSol_SPBCGSSetMaxl`        |
+-----------------------------+------------------------------------------+
| ``SUNSPFGMR``               | :c:func:`SUNLinSol_SPFGMR`               |
+-----------------------------+------------------------------------------+
| ``SUNSPFGMRSetPrecType``    | :c:func:`SUNLinSol_SPFGMRSetPrecType`    |
+-----------------------------+------------------------------------------+
| ``SUNSPFGMRSetGSType``      | :c:func:`SUNLinSol_SPFGMRSetGSType`      |
+-----------------------------+------------------------------------------+
| ``SUNSPFGMRSetMaxRestarts`` | :c:func:`SUNLinSol_SPFGMRSetMaxRestarts` |
+-----------------------------+------------------------------------------+
| ``SUNSPGMR``                | :c:func:`SUNLinSol_SPGMR`                |
+-----------------------------+------------------------------------------+
| ``SUNSPGMRSetPrecType``     | :c:func:`SUNLinSol_SPGMRSetPrecType`     |
+-----------------------------+------------------------------------------+
| ``SUNSPGMRSetGSType``       | :c:func:`SUNLinSol_SPGMRSetGSType`       |
+-----------------------------+------------------------------------------+
| ``SUNSPGMRSetMaxRestarts``  | :c:func:`SUNLinSol_SPGMRSetMaxRestarts`  |
+-----------------------------+------------------------------------------+
| ``SUNSPTFQMR``              | :c:func:`SUNLinSol_SPTFQMR`              |
+-----------------------------+------------------------------------------+
| ``SUNSPTFQMRSetPrecType``   | :c:func:`SUNLinSol_SPTFQMRSetPrecType`   |
+-----------------------------+------------------------------------------+
| ``SUNSPTFQMRSetMaxl``       | :c:func:`SUNLinSol_SPTFQMRSetMaxl`       |
+-----------------------------+------------------------------------------+
| ``SUNSuperLUMT``            | :c:func:`SUNLinSol_SuperLUMT`            |
+-----------------------------+------------------------------------------+
| ``SUNSuperLUMTSetOrdering`` | :c:func:`SUNLinSol_SuperLUMTSetOrdering` |
+-----------------------------+------------------------------------------+

**KINSOL**

New orthogonalization methods were added for use within the KINSOL Anderson
acceleration routine. See :numref:`Anderson_QR` and :c:func:`KINSetOrthAA`
for more details.

The KINSOL Fortran 77 interface has been removed. See :numref:`SUNDIALS.Fortran`
and the F2003 example programs for more details using the SUNDIALS Fortran 2003
module interfaces.

**Deprecations**

In addition to the deprecations noted elsewhere, many constants, types, and
functions have been renamed so that they are properly namespaced. The old names
have been deprecated and will be removed in SUNDIALS v7.0.0.

The following constants, macros, and typedefs are now deprecated:

+------------------------------+-------------------------------------+
| Deprecated Name              | New Name                            |
+==============================+=====================================+
| ``realtype``                 | ``sunrealtype``                     |
+------------------------------+-------------------------------------+
| ``booleantype``              | ``sunbooleantype``                  |
+------------------------------+-------------------------------------+
| ``RCONST``                   | ``SUN_RCONST``                      |
+------------------------------+-------------------------------------+
| ``BIG_REAL``                 | ``SUN_BIG_REAL``                    |
+------------------------------+-------------------------------------+
| ``SMALL_REAL``               | ``SUN_SMALL_REAL``                  |
+------------------------------+-------------------------------------+
| ``UNIT_ROUNDOFF``            | ``SUN_UNIT_ROUNDOFF``               |
+------------------------------+-------------------------------------+
| ``PREC_NONE``                | ``SUN_PREC_NONE``                   |
+------------------------------+-------------------------------------+
| ``PREC_LEFT``                | ``SUN_PREC_LEFT``                   |
+------------------------------+-------------------------------------+
| ``PREC_RIGHT``               | ``SUN_PREC_RIGHT``                  |
+------------------------------+-------------------------------------+
| ``PREC_BOTH``                | ``SUN_PREC_BOTH``                   |
+------------------------------+-------------------------------------+
| ``MODIFIED_GS``              | ``SUN_MODIFIED_GS``                 |
+------------------------------+-------------------------------------+
| ``CLASSICAL_GS``             | ``SUN_CLASSICAL_GS``                |
+------------------------------+-------------------------------------+
| ``ATimesFn``                 | ``SUNATimesFn``                     |
+------------------------------+-------------------------------------+
| ``PSetupFn``                 | ``SUNPSetupFn``                     |
+------------------------------+-------------------------------------+
| ``PSolveFn``                 | ``SUNPSolveFn``                     |
+------------------------------+-------------------------------------+
| ``DlsMat``                   | ``SUNDlsMat``                       |
+------------------------------+-------------------------------------+
| ``DENSE_COL``                | ``SUNDLS_DENSE_COL``                |
+------------------------------+-------------------------------------+
| ``DENSE_ELEM``               | ``SUNDLS_DENSE_ELEM``               |
+------------------------------+-------------------------------------+
| ``BAND_COL``                 | ``SUNDLS_BAND_COL``                 |
+------------------------------+-------------------------------------+
| ``BAND_COL_ELEM``            | ``SUNDLS_BAND_COL_ELEM``            |
+------------------------------+-------------------------------------+
| ``BAND_ELEM``                | ``SUNDLS_BAND_ELEM``                |
+------------------------------+-------------------------------------+

In addition, the following functions are now deprecated (compile-time warnings
will be thrown if supported by the compiler):

+---------------------------------+--------------------------------+
| Deprecated Name                 | New Name                       |
+=================================+================================+
| ``KINDlsSetLinearSolver``       | ``KINSetLinearSolver``         |
+---------------------------------+--------------------------------+
| ``KINDlsSetJacFn``              | ``KINSetJacFn``                |
+---------------------------------+--------------------------------+
| ``KINDlsGetWorkSpace``          | ``KINGetLinWorkSpace``         |
+---------------------------------+--------------------------------+
| ``KINDlsGetNumJacEvals``        | ``KINGetNumJacEvals``          |
+---------------------------------+--------------------------------+
| ``KINDlsGetNumFuncEvals``       | ``KINGetNumLinFuncEvals``      |
+---------------------------------+--------------------------------+
| ``KINDlsGetLastFlag``           | ``KINGetLastLinFlag``          |
+---------------------------------+--------------------------------+
| ``KINDlsGetReturnFlagName``     | ``KINGetLinReturnFlagName``    |
+---------------------------------+--------------------------------+
| ``KINSpilsSetLinearSolver``     | ``KINSetLinearSolver``         |
+---------------------------------+--------------------------------+
| ``KINSpilsSetPreconditioner``   | ``KINSetPreconditioner``       |
+---------------------------------+--------------------------------+
| ``KINSpilsSetJacTimesVecFn``    | ``KINSetJacTimesVecFn``        |
+---------------------------------+--------------------------------+
| ``KINSpilsGetWorkSpace``        | ``KINGetLinWorkSpace``         |
+---------------------------------+--------------------------------+
| ``KINSpilsGetNumPrecEvals``     | ``KINGetNumPrecEvals``         |
+---------------------------------+--------------------------------+
| ``KINSpilsGetNumPrecSolves``    | ``KINGetNumPrecSolves``        |
+---------------------------------+--------------------------------+
| ``KINSpilsGetNumLinIters``      | ``KINGetNumLinIters``          |
+---------------------------------+--------------------------------+
| ``KINSpilsGetNumConvFails``     | ``KINGetNumLinConvFails``      |
+---------------------------------+--------------------------------+
| ``KINSpilsGetNumJtimesEvals``   | ``KINGetNumJtimesEvals``       |
+---------------------------------+--------------------------------+
| ``KINSpilsGetNumFuncEvals``     | ``KINGetNumLinFuncEvals``      |
+---------------------------------+--------------------------------+
| ``KINSpilsGetLastFlag``         | ``KINGetLastLinFlag``          |
+---------------------------------+--------------------------------+
| ``KINSpilsGetReturnFlagName``   | ``KINGetLinReturnFlagName``    |
+---------------------------------+--------------------------------+
| ``DenseGETRF``                  | ``SUNDlsMat_DenseGETRF``       |
+---------------------------------+--------------------------------+
| ``DenseGETRS``                  | ``SUNDlsMat_DenseGETRS``       |
+---------------------------------+--------------------------------+
| ``denseGETRF``                  | ``SUNDlsMat_denseGETRF``       |
+---------------------------------+--------------------------------+
| ``denseGETRS``                  | ``SUNDlsMat_denseGETRS``       |
+---------------------------------+--------------------------------+
| ``DensePOTRF``                  | ``SUNDlsMat_DensePOTRF``       |
+---------------------------------+--------------------------------+
| ``DensePOTRS``                  | ``SUNDlsMat_DensePOTRS``       |
+---------------------------------+--------------------------------+
| ``densePOTRF``                  | ``SUNDlsMat_densePOTRF``       |
+---------------------------------+--------------------------------+
| ``densePOTRS``                  | ``SUNDlsMat_densePOTRS``       |
+---------------------------------+--------------------------------+
| ``DenseGEQRF``                  | ``SUNDlsMat_DenseGEQRF``       |
+---------------------------------+--------------------------------+
| ``DenseORMQR``                  | ``SUNDlsMat_DenseORMQR``       |
+---------------------------------+--------------------------------+
| ``denseGEQRF``                  | ``SUNDlsMat_denseGEQRF``       |
+---------------------------------+--------------------------------+
| ``denseORMQR``                  | ``SUNDlsMat_denseORMQR``       |
+---------------------------------+--------------------------------+
| ``DenseCopy``                   | ``SUNDlsMat_DenseCopy``        |
+---------------------------------+--------------------------------+
| ``denseCopy``                   | ``SUNDlsMat_denseCopy``        |
+---------------------------------+--------------------------------+
| ``DenseScale``                  | ``SUNDlsMat_DenseScale``       |
+---------------------------------+--------------------------------+
| ``denseScale``                  | ``SUNDlsMat_denseScale``       |
+---------------------------------+--------------------------------+
| ``denseAddIdentity``            | ``SUNDlsMat_denseAddIdentity`` |
+---------------------------------+--------------------------------+
| ``DenseMatvec``                 | ``SUNDlsMat_DenseMatvec``      |
+---------------------------------+--------------------------------+
| ``denseMatvec``                 | ``SUNDlsMat_denseMatvec``      |
+---------------------------------+--------------------------------+
| ``BandGBTRF``                   | ``SUNDlsMat_BandGBTRF``        |
+---------------------------------+--------------------------------+
| ``bandGBTRF``                   | ``SUNDlsMat_bandGBTRF``        |
+---------------------------------+--------------------------------+
| ``BandGBTRS``                   | ``SUNDlsMat_BandGBTRS``        |
+---------------------------------+--------------------------------+
| ``bandGBTRS``                   | ``SUNDlsMat_bandGBTRS``        |
+---------------------------------+--------------------------------+
| ``BandCopy``                    | ``SUNDlsMat_BandCopy``         |
+---------------------------------+--------------------------------+
| ``bandCopy``                    | ``SUNDlsMat_bandCopy``         |
+---------------------------------+--------------------------------+
| ``BandScale``                   | ``SUNDlsMat_BandScale``        |
+---------------------------------+--------------------------------+
| ``bandScale``                   | ``SUNDlsMat_bandScale``        |
+---------------------------------+--------------------------------+
| ``bandAddIdentity``             | ``SUNDlsMat_bandAddIdentity``  |
+---------------------------------+--------------------------------+
| ``BandMatvec``                  | ``SUNDlsMat_BandMatvec``       |
+---------------------------------+--------------------------------+
| ``bandMatvec``                  | ``SUNDlsMat_bandMatvec``       |
+---------------------------------+--------------------------------+
| ``ModifiedGS``                  | ``SUNModifiedGS``              |
+---------------------------------+--------------------------------+
| ``ClassicalGS``                 | ``SUNClassicalGS``             |
+---------------------------------+--------------------------------+
| ``QRfact``                      | ``SUNQRFact``                  |
+---------------------------------+--------------------------------+
| ``QRsol``                       | ``SUNQRsol``                   |
+---------------------------------+--------------------------------+
| ``DlsMat_NewDenseMat``          | ``SUNDlsMat_NewDenseMat``      |
+---------------------------------+--------------------------------+
| ``DlsMat_NewBandMat``           | ``SUNDlsMat_NewBandMat``       |
+---------------------------------+--------------------------------+
| ``DestroyMat``                  | ``SUNDlsMat_DestroyMat``       |
+---------------------------------+--------------------------------+
| ``NewIntArray``                 | ``SUNDlsMat_NewIntArray``      |
+---------------------------------+--------------------------------+
| ``NewIndexArray``               | ``SUNDlsMat_NewIndexArray``    |
+---------------------------------+--------------------------------+
| ``NewRealArray``                | ``SUNDlsMat_NewRealArray``     |
+---------------------------------+--------------------------------+
| ``DestroyArray``                | ``SUNDlsMat_DestroyArray``     |
+---------------------------------+--------------------------------+
| ``AddIdentity``                 | ``SUNDlsMat_AddIdentity``      |
+---------------------------------+--------------------------------+
| ``SetToZero``                   | ``SUNDlsMat_SetToZero``        |
+---------------------------------+--------------------------------+
| ``PrintMat``                    | ``SUNDlsMat_PrintMat``         |
+---------------------------------+--------------------------------+
| ``newDenseMat``                 | ``SUNDlsMat_newDenseMat``      |
+---------------------------------+--------------------------------+
| ``newBandMat``                  | ``SUNDlsMat_newBandMat``       |
+---------------------------------+--------------------------------+
| ``destroyMat``                  | ``SUNDlsMat_destroyMat``       |
+---------------------------------+--------------------------------+
| ``newIntArray``                 | ``SUNDlsMat_newIntArray``      |
+---------------------------------+--------------------------------+
| ``newIndexArray``               | ``SUNDlsMat_newIndexArray``    |
+---------------------------------+--------------------------------+
| ``newRealArray``                | ``SUNDlsMat_newRealArray``     |
+---------------------------------+--------------------------------+
| ``destroyArray``                | ``SUNDlsMat_destroyArray``     |
+---------------------------------+--------------------------------+

In addition, the entire ``sundials_lapack.h`` header file is now deprecated for
removal in SUNDIALS v7.0.0. Note, this header file is not needed to use the
SUNDIALS LAPACK linear solvers.




## Changes to SUNDIALS in release 5.8.0

The RAJA NVECTOR implementation has been updated to support the SYCL backend
in addition to the CUDA and HIP backend. Users can choose the backend when
configuring SUNDIALS by using the `SUNDIALS_RAJA_BACKENDS` CMake variable. This
module remains experimental and is subject to change from version to version.

A new SUNMatrix and SUNLinearSolver implementation were added to interface
with the Intel oneAPI Math Kernel Library (oneMKL). Both the matrix and the
linear solver support general dense linear systems as well as block diagonal
linear systems. This module is experimental and is subject to change from
version to version.

Added a new *optional* function to the SUNLinearSolver API,
`SUNLinSolSetZeroGuess`, to indicate that the next call to `SUNlinSolSolve` will
be made with a zero initial guess. SUNLinearSolver implementations that do not
use the `SUNLinSolNewEmpty` constructor will, at a minimum, need set the
`setzeroguess` function pointer in the linear solver `ops` structure to
`NULL`. The SUNDIALS iterative linear solver implementations have been updated
to leverage this new set function to remove one dot product per solve.

The time integrator packages (ARKODE, CVODE(S), and IDA(S)) all now support a
new "matrix-embedded" SUNLinearSolver type.  This type supports user-supplied
SUNLinearSolver implementations that set up and solve the specified linear
system at each linear solve call.  Any matrix-related data structures are held
internally to the linear solver itself, and are not provided by the SUNDIALS
package.

Added functions to ARKODE and CVODE(S) for supplying an alternative right-hand
side function and to IDA(S) for supplying an alternative residual for use within
nonlinear system function evaluations.

Support for user-defined inner (fast) integrators has been to the MRIStep module
in ARKODE. See the "MRIStep Custom Inner Steppers" section in the user guide for
more information on providing a user-defined integration method.

Added specialized fused HIP kernels to CVODE which may offer better
performance on smaller problems when using CVODE with the `NVECTOR_HIP`
module. See the optional input function `CVodeSetUseIntegratorFusedKernels`
for more information. As with other SUNDIALS HIP features, this is
feature is experimental and may change from version to version.

New KINSOL options have been added to apply a constant damping factor in the
fixed point and Picard iterations (see `KINSetDamping`), to delay the start of
Anderson acceleration with the fixed point and Picard iterations (see
`KINSetDelayAA`), and to return the newest solution with the fixed point
iteration (see `KINSetReturnNewest`).

The installed SUNDIALSConfig.cmake file now supports the `COMPONENTS` option
to `find_package`. The exported targets no longer have IMPORTED_GLOBAL set.

A bug was fixed in `SUNMatCopyOps` where the matrix-vector product setup
function pointer was not copied.

A bug was fixed in the SPBCGS and SPTFQMR solvers for the case where a non-zero
initial guess and a solution scaling vector are provided. This fix only impacts
codes using SPBCGS or SPTFQMR as standalone solvers as all SUNDIALS packages
utilize a zero initial guess.

A bug was fixed in the ARKODE stepper modules where the stop time may be passed
after resetting the integrator.

A bug was fixed in `IDASetJacTimesResFn` in IDAS where the supplied function was
used in the dense finite difference Jacobian computation rather than the finite
difference Jacobian-vector product approximation.

A bug was fixed in the KINSOL Picard iteration where the value of
`KINSetMaxSetupCalls` would be ignored.

### ARKODE Changes in v4.8.0

The RAJA NVECTOR implementation has been updated to support the SYCL backend
in addition to the CUDA and HIP backend. Users can choose the backend when
configuring SUNDIALS by using the ``SUNDIALS_RAJA_BACKENDS`` CMake variable.
This module remains experimental and is subject to change from version to
version.

A new SUNMatrix and SUNLinearSolver implementation were added to interface with
the Intel oneAPI Math Kernel Library (oneMKL). Both the matrix and the linear
solver support general dense linear systems as well as block diagonal linear
systems. See :numref:`SUNLinSol.OneMklDense` for more details. This module is
experimental and is subject to change from version to version.

Added a new *optional* function to the SUNLinearSolver API,
:c:func:`SUNLinSolSetZeroGuess`, to indicate that the next call to
:c:func:`SUNLinSolSolve` will be made with a zero initial guess. SUNLinearSolver
implementations that do not use the :c:func:`SUNLinSolNewEmpty` constructor
will, at a minimum, need set the ``setzeroguess`` function pointer in the linear
solver ``ops`` structure to ``NULL``. The SUNDIALS iterative linear solver
implementations have been updated to leverage this new set function to remove
one dot product per solve.

ARKODE now supports a new "matrix-embedded" SUNLinearSolver type.  This type
supports user-supplied SUNLinearSolver implementations that set up and solve
the specified linear system at each linear solve call.  Any matrix-related data
structures are held internally to the linear solver itself, and are not
provided by the SUNDIALS package.

Support for user-defined inner (fast) integrators has been to the MRIStep
module. See :numref:`ARKODE.Usage.MRIStep.CustomInnerStepper` for more information on providing
a user-defined integration method.

Added the functions :c:func:`ARKStepSetNlsRhsFn()` and
:c:func:`MRIStepSetNlsRhsFn()` to supply an alternative implicit right-hand side
function for use within nonlinear system function evaluations.

The installed SUNDIALSConfig.cmake file now supports the ``COMPONENTS`` option
to ``find_package``. The exported targets no longer have IMPORTED_GLOBAL set.

A bug was fixed in :c:func:`SUNMatCopyOps` where the matrix-vector product setup
function pointer was not copied.

A bug was fixed in the SPBCGS and SPTFQMR solvers for the case where a non-zero
initial guess and a solution scaling vector are provided. This fix only impacts
codes using SPBCGS or SPTFQMR as standalone solvers as all SUNDIALS packages
utilize a zero initial guess.

A bug was fixed in the ARKODE stepper modules where the stop time may be passed
after resetting the integrator.


### CVODE Changes in v5.8.0

The :ref:`RAJA N_Vector <NVectors.RAJA>` implementation has been updated to
support the SYCL backend in addition to the CUDA and HIP backend. Users can
choose the backend when configuring SUNDIALS by using the
``SUNDIALS_RAJA_BACKENDS`` CMake variable.  This module remains experimental
and is subject to change from version to version.

New :c:type:`SUNMatrix` and :c:type:`SUNLinearSolver` implementations were added to
interface with the Intel oneAPI Math Kernel Library (oneMKL). Both the matrix
and the linear solver support general dense linear systems as well as block
diagonal linear systems. See :numref:`SUNLinSol.OneMklDense` for more details.
This module is experimental and is subject to change from version to version.

Added a new *optional* function to the SUNLinearSolver API,
:c:func:`SUNLinSolSetZeroGuess()`, to indicate that the next call to
:c:func:`SUNlinSolSolve()` will be made with a zero initial guess. SUNLinearSolver
implementations that do not use the :c:func:`SUNLinSolNewEmpty` constructor will,
at a minimum, need set the ``setzeroguess`` function pointer in the linear solver
``ops`` structure to ``NULL``. The SUNDIALS iterative linear solver
implementations have been updated to leverage this new set function to remove
one dot product per solve.

CVODE now supports a new "matrix-embedded" :c:type:`SUNLinearSolver` type.  This
type supports user-supplied :c:type:`SUNLinearSolver` implementations that set up
and solve the specified linear system at each linear solve call.  Any
matrix-related data structures are held internally to the linear solver itself,
and are not provided by the SUNDIALS package.

Added specialized fused HIP kernels to CVODE which may offer better
performance on smaller problems when using CVODE with the
:ref:`N_Vector HIP <NVectors.HIP>` module. See the optional input function
:c:func:`CVodeSetUseIntegratorFusedKernels()` for more information. As with
other SUNDIALS HIP features, this capability is considered experimental and may
change from version to version.

Added the function :c:func:`CVodeSetNlsRhsFn()` to supply an alternative right-hand
side function for use within nonlinear system function evaluations.

The installed ``SUNDIALSConfig.cmake`` file now supports the ``COMPONENTS`` option
to ``find_package``. The exported targets no longer have ``IMPORTED_GLOBAL``
set.

A bug was fixed in :c:func:`SUNMatCopyOps()` where the matrix-vector product setup
function pointer was not copied.

A bug was fixed in the :c:ref:`SPBCGS <SUNLinSol.SPBCGS>` and
:c:ref:`SPTFQMR <SUNLinSol.SPTFQMR>` solvers for the case where a non-zero
initial guess and a solution scaling vector are provided. This fix only impacts
codes using :c:ref:`SPBCGS <SUNLinSol.SPBCGS>` or :c:ref:`SPTFQMR <SUNLinSol.SPTFQMR>`
as standalone solvers as all SUNDIALS packages utilize a zero initial guess.

### CVODES Changes in v5.8.0

The RAJA ``N_Vector`` implementation has been updated to support the SYCL
backend in addition to the CUDA and HIP backend. Users can choose the backend
when configuring SUNDIALS by using the ``SUNDIALS_RAJA_BACKENDS`` CMake
variable. This module remains experimental and is subject to change from version
to version.

A new ``SUNMatrix`` and ``SUNLinearSolver`` implementation were added to
interface with the Intel oneAPI Math Kernel Library (oneMKL). Both the matrix
and the linear solver support general dense linear systems as well as block
diagonal linear systems. See Chapter :numref:`SUNLinSol.OneMklDense` for more
details. This module is experimental and is subject to change from version to
version.

Added a new *optional* function to the SUNLinearSolver API,
``SUNLinSolSetZeroGuess``, to indicate that the next call to ``SUNlinSolSolve``
will be made with a zero initial guess. SUNLinearSolver implementations that do
not use the ``SUNLinSolNewEmpty`` constructor will, at a minimum, need set the
``setzeroguess`` function pointer in the linear solver ``ops`` structure to
``NULL``. The SUNDIALS iterative linear solver implementations have been updated
to leverage this new set function to remove one dot product per solve.

CVODES now supports a new "matrix-embedded" ``SUNLinearSolver`` type. This type
supports user-supplied ``SUNLinearSolver`` implementations that set up and solve
the specified linear system at each linear solve call. Any matrix-related data
structures are held internally to the linear solver itself, and are not provided
by the SUNDIALS package.

Added the function ``CVodeSetNlsRhsFn`` to supply an alternative right-hand side
function for use within nonlinear system function evaluations.

The installed SUNDIALSConfig.cmake file now supports the ``COMPONENTS`` option
to ``find_package``. The exported targets no longer have ``IMPORTED_GLOBAL``
set.

A bug was fixed in ``SUNMatCopyOps`` where the matrix-vector product setup
function pointer was not copied.

A bug was fixed in the SPBCGS and SPTFQMR solvers for the case where a non-zero
initial guess and a solution scaling vector are provided. This fix only impacts
codes using SPBCGS or SPTFQMR as standalone solvers as all SUNDIALS packages
utilize a zero initial guess.

### IDA Changes in v5.8.0

The :ref:`RAJA N_Vector <NVectors.RAJA>` implementation has been updated to
support the SYCL backend in addition to the CUDA and HIP backends. Users can
choose the backend when configuring SUNDIALS by using the
:cmakeop:`SUNDIALS_RAJA_BACKENDS` CMake variable. This module remains
experimental and is subject to change from version to version.

A new ``SUNMatrix`` and ``SUNLinearSolver`` implementation were added to
interface with the Intel oneAPI Math Kernel Library (oneMKL). Both the matrix
and the linear solver support general dense linear systems as well as block
diagonal linear systems. See :numref:`SUNLinSol.OneMklDense` for more
details. This module is experimental and is subject to change from version to
version.

Added a new *optional* function to the ``SUNLinearSolver`` API,
:c:func:`SUNLinSolSetZeroGuess`, to indicate that the next call to
:c:func:`SUNLinSolSolve` will be made with a zero initial guess.
``SUNLinearSolver`` implementations that do not use the
:c:func:`SUNLinSolNewEmpty` constructor will, at a minimum, need set the
``setzeroguess`` function pointer in the linear solver ``ops`` structure to
``NULL``. The SUNDIALS iterative linear solver implementations have been updated
to leverage this new set function to remove one dot product per solve.

IDA now supports a new "matrix-embedded" ``SUNLinearSolver`` type. This type
supports user-supplied ``SUNLinearSolver`` implementations that set up and solve
the specified linear system at each linear solve call. Any matrix-related data
structures are held internally to the linear solver itself, and are not provided
by the SUNDIALS package.

Added the function :c:func:`IDASetNlsResFn` to supply an alternative residual
side function for use within nonlinear system function evaluations.

The installed ``SUNDIALSConfig.cmake`` file now supports the ``COMPONENTS``
option to ``find_package``.

A bug was fixed in :c:func:`SUNMatCopyOps` where the matrix-vector product setup
function pointer was not copied.

A bug was fixed in the SPBCGS and SPTFQMR solvers for the case where a non-zero
initial guess and a solution scaling vector are provided. This fix only impacts
codes using SPBCGS or SPTFQMR as standalone solvers as all SUNDIALS packages
utilize a zero initial guess.

### IDAS Changes in v4.8.0

The :ref:`RAJA N_Vector <NVectors.RAJA>` implementation has been updated to
support the SYCL backend in addition to the CUDA and HIP backends. Users can
choose the backend when configuring SUNDIALS by using the
:cmakeop:`SUNDIALS_RAJA_BACKENDS` CMake variable. This module remains
experimental and is subject to change from version to version.

A new ``SUNMatrix`` and ``SUNLinearSolver`` implementation were added to
interface with the Intel oneAPI Math Kernel Library (oneMKL). Both the matrix
and the linear solver support general dense linear systems as well as block
diagonal linear systems. See :numref:`SUNLinSol.OneMklDense` for more
details. This module is experimental and is subject to change from version to
version.

Added a new *optional* function to the ``SUNLinearSolver`` API,
:c:func:`SUNLinSolSetZeroGuess`, to indicate that the next call to
:c:func:`SUNLinSolSolve` will be made with a zero initial guess.
``SUNLinearSolver`` implementations that do not use the
:c:func:`SUNLinSolNewEmpty` constructor will, at a minimum, need set the
``setzeroguess`` function pointer in the linear solver ``ops`` structure to
``NULL``. The SUNDIALS iterative linear solver implementations have been updated
to leverage this new set function to remove one dot product per solve.

IDAS now supports a new "matrix-embedded" ``SUNLinearSolver`` type. This type
supports user-supplied ``SUNLinearSolver`` implementations that set up and solve
the specified linear system at each linear solve call. Any matrix-related data
structures are held internally to the linear solver itself, and are not provided
by the SUNDIALS package.

Added the function :c:func:`IDASetNlsResFn` to supply an alternative residual
side function for use within nonlinear system function evaluations.

The installed ``SUNDIALSConfig.cmake`` file now supports the ``COMPONENTS``
option to ``find_package``.

A bug was fixed in :c:func:`SUNMatCopyOps` where the matrix-vector product setup
function pointer was not copied.

A bug was fixed in the SPBCGS and SPTFQMR solvers for the case where a non-zero
initial guess and a solution scaling vector are provided. This fix only impacts
codes using SPBCGS or SPTFQMR as standalone solvers as all SUNDIALS packages
utilize a zero initial guess.

### KINSOL Changes in v5.8.0

The RAJA ``N_Vector`` implementation has been updated to support the SYCL backend in addition to the CUDA and HIP
backend. Users can choose the backend when configuring SUNDIALS by using the ``SUNDIALS_RAJA_BACKENDS`` CMake variable.
This module remains experimental and is subject to change from version to version.

A new ``SUNMatrix`` and ``SUNLinearSolver`` implementation were added to interface with the Intel oneAPI Math Kernel
Library (oneMKL). Both the matrix and the linear solver support general dense linear systems as well as block diagonal
linear systems. See :numref:`SUNLinSol.OneMklDense`  for more details. This module is
experimental and is subject to change from version to version.

Added a new *optional* function to the SUNLinearSolver API, ``SUNLinSolSetZeroGuess``, to indicate that the next call to
``SUNlinSolSolve`` will be made with a zero initial guess. SUNLinearSolver implementations that do not use the
``SUNLinSolNewEmpty`` constructor will, at a minimum, need set the ``setzeroguess`` function pointer in the linear
solver ``ops`` structure to ``NULL``. The SUNDIALS iterative linear solver implementations have been updated to leverage
this new set function to remove one dot product per solve.

New KINSOL options have been added to apply a constant damping in the fixed point and Picard iterations (see
``KINSetDamping``), to delay the start of Anderson acceleration with the fixed point and Picard iterations (see
``KINSetDelayAA``), and to return the newest solution with the fixed point iteration (see ``KINSetReturnNewest``).

The installed SUNDIALSConfig.cmake file now supports the ``COMPONENTS`` option to ``find_package``. The exported targets
no longer have ``IMPORTED_GLOBAL`` set.

A bug was fixed in ``SUNMatCopyOps`` where the matrix-vector product setup function pointer was not copied.

A bug was fixed in the SPBCGS and SPTFQMR solvers for the case where a non-zero initial guess and a solution scaling
vector are provided. This fix only impacts codes using SPBCGS or SPTFQMR as standalone solvers as all SUNDIALS
packages utilize a zero initial guess.

A bug was fixed in the Picard iteration where the value of ``KINSetMaxSetupCalls`` would be ignored.


## Changes to SUNDIALS in release 5.7.0

A new NVECTOR implementation based on the SYCL abstraction layer has been added
targeting Intel GPUs. At present the only SYCL compiler supported is the DPC++
(Intel oneAPI) compiler. See the SYCL NVECTOR section in the user guide for more
details. This module is considered experimental and is subject to major changes
even in minor releases.

A new SUNMatrix and SUNLinearSolver implementation were added to interface
with the MAGMA linear algebra library. Both the matrix and the linear solver
support general dense linear systems as well as block diagonal linear systems,
and both are targeted at GPUs (AMD or NVIDIA).

### ARKODE Changes in v4.7.0

A new NVECTOR implementation based on the SYCL abstraction layer has been added
targeting Intel GPUs. At present the only SYCL compiler supported is the DPC++
(Intel oneAPI) compiler. See :numref:`NVectors.SYCL` for more details. This module
is considered experimental and is subject to major changes even in minor
releases.

A new SUNMatrix and SUNLinearSolver implementation were added to interface
with the MAGMA linear algebra library. Both the matrix and the linear solver
support general dense linear systems as well as block diagonal linear systems,
and both are targeted at GPUs (AMD or NVIDIA). See :numref:`SUNLinSol.MagmaDense`
for more details.

### CVODE Changes in v5.7.0

A new :c:type:`N_Vector` implementation based on the SYCL abstraction layer
has been added targeting Intel GPUs. At present the only SYCL
compiler supported is the DPC++ (Intel oneAPI) compiler. See
:numref:`NVectors.sycl` for more details. This module is
considered experimental and is subject to major changes even in minor
releases.

New ``SUNMatrix`` and ``SUNLinearSolver`` implementations were added to
interface with the MAGMA linear algebra library. Both the matrix and the
linear solver support general dense linear systems as well as block
diagonal linear systems, and both are targeted at GPUs (AMD or NVIDIA).
See :numref:`SUNLinSol.MagmaDense` for more details.

### CVODES Changes in v5.7.0

A new ``N_Vector`` implementation based on the SYCL abstraction layer has been added targeting Intel GPUs. At
present the only SYCL compiler supported is the DPC++ (Intel oneAPI) compiler. See Section
:numref:`NVectors.sycl` for more details. This module is considered experimental and is subject to major
changes even in minor releases.

A new ``SUNMatrix`` and ``SUNLinearSolver`` implementation were added to interface with the MAGMA linear algebra library.
Both the matrix and the linear solver support general dense linear systems as well as block diagonal linear systems, and
both are targeted at GPUs (AMD or NVIDIA). See Section :numref:`SUNLinSol.magmadense` for more
details.

### IDA Changes in v5.7.0

A new ``N_Vector`` implementation based on the SYCL abstraction layer has been
added targeting Intel GPUs. At present the only SYCL compiler supported is the
DPC++ (Intel oneAPI) compiler. See :numref:`NVectors.SYCL` for more details.
This module is considered experimental and is subject to major changes even in
minor releases.

A new ``SUNMatrix`` and ``SUNLinearSolver`` implementation were added to
interface with the MAGMA linear algebra library. Both the matrix and the linear
solver support general dense linear systems as well as block diagonal linear
systems, and both are targeted at GPUs (AMD or NVIDIA). See
:numref:`SUNLinSol.MagmaDense` for more details.

### IDAS Changes in v4.7.0

A new ``N_Vector`` implementation based on the SYCL abstraction layer has been
added targeting Intel GPUs. At present the only SYCL compiler supported is the
DPC++ (Intel oneAPI) compiler. See :numref:`NVectors.SYCL` for more details.
This module is considered experimental and is subject to major changes even in
minor releases.

A new ``SUNMatrix`` and ``SUNLinearSolver`` implementation were added to
interface with the MAGMA linear algebra library. Both the matrix and the linear
solver support general dense linear systems as well as block diagonal linear
systems, and both are targeted at GPUs (AMD or NVIDIA). See
:numref:`SUNLinSol.MagmaDense` for more details.

### KINSOL Changes in v5.7.0

A new ``N_Vector`` implementation based on the SYCL abstraction layer has been added targeting Intel GPUs. At
present the only SYCL compiler supported is the DPC++ (Intel oneAPI) compiler. See :numref:`NVectors.SYCL` for more details. This module is considered experimental and is subject to major
changes even in minor releases.

A new ``SUNMatrix`` and ``SUNLinearSolver`` implementation were added to interface with the MAGMA linear algebra library.
Both the matrix and the linear solver support general dense linear systems as well as block diagonal linear systems, and
both are targeted at GPUs (AMD or NVIDIA). See :numref:`SUNLinSol.MagmaDense` for more
details.



## Changes to SUNDIALS in release 5.6.1

Fixed a bug in the SUNDIALS CMake which caused an error if the
`CMAKE_CXX_STANDARD` and `SUNDIALS_RAJA_BACKENDS` options were not provided.

Fixed some compiler warnings when using the IBM XL compilers.

### ARKODE Changes in v4.6.1

Fixed a bug in the SUNDIALS CMake which caused an error
if the CMAKE_CXX_STANDARD and SUNDIALS_RAJA_BACKENDS options
were not provided.

Fixed some compiler warnings when using the IBM XL compilers.

### CVODE Changes in v5.6.1

Fixed a bug in the SUNDIALS CMake which caused an error if the
``CMAKE_CXX_STANDARD`` and ``SUNDIALS_RAJA_BACKENDS`` options were not provided.

Fixed some compiler warnings when using the IBM XL compilers.

### CVODES Changes in v5.6.1

Fixed a bug in the SUNDIALS CMake which caused an error if the CMAKE_CXX_STANDARD and SUNDIALS_RAJA_BACKENDS
options were not provided.

Fixed some compiler warnings when using the IBM XL compilers.

### IDA Changes in v5.6.1

Fixed a bug in the SUNDIALS CMake which caused an error if the
:cmakeop:`CMAKE_CXX_STANDARD` and :cmakeop:`SUNDIALS_RAJA_BACKENDS` options were
not provided.

Fixed some compiler warnings when using the IBM XL compilers.

### IDAS Changes in v4.6.1

Fixed a bug in the SUNDIALS CMake which caused an error if the
:cmakeop:`CMAKE_CXX_STANDARD` and :cmakeop:`SUNDIALS_RAJA_BACKENDS` options were
not provided.

Fixed some compiler warnings when using the IBM XL compilers.

### KINSOL Changes in v5.6.1

Fixed a bug in the SUNDIALS CMake which caused an error if the CMAKE_CXX_STANDARD and SUNDIALS_RAJA_BACKENDS
options were not provided.

Fixed some compiler warnings when using the IBM XL compilers.



## Changes to SUNDIALS in release 5.6.0

A new NVECTOR implementation based on the AMD ROCm HIP platform has been added.
This vector can target NVIDIA or AMD GPUs. See HIP NVECTOR section in the user
guide for more details. This module is considered experimental and is subject to
change from version to version.

The RAJA NVECTOR implementation has been updated to support the HIP backend
in addition to the CUDA backend. Users can choose the backend when configuring
SUNDIALS by using the `SUNDIALS_RAJA_BACKENDS` CMake variable. This module
remains experimental and is subject to change from version to version.

A new optional operation, `N_VGetDeviceArrayPointer`, was added to the N_Vector
API. This operation is useful for N_Vectors that utilize dual memory spaces,
e.g. the native SUNDIALS CUDA N_Vector.

The SUNMATRIX_CUSPARSE and SUNLINEARSOLVER_CUSOLVERSP_BATCHQR implementations
no longer require the SUNDIALS CUDA N_Vector. Instead, they require that the
vector utilized provides the `N_VGetDeviceArrayPointer` operation, and that the
pointer returned by `N_VGetDeviceArrayPointer` is a valid CUDA device pointer.

### ARKODE Changes in v4.6.0

A new NVECTOR implementation based on the AMD ROCm HIP platform has been added.
This vector can target NVIDIA or AMD GPUs. See :numref:`NVectors.HIP` for more
details. This module is considered experimental and is subject to change from
version to version.

The RAJA NVECTOR implementation has been updated to support the HIP backend
in addition to the CUDA backend. Users can choose the backend when configuring
SUNDIALS by using the ``SUNDIALS_RAJA_BACKENDS`` CMake variable. This module
remains experimental and is subject to change from version to version.

A new optional operation, :c:func:`N_VGetDeviceArrayPointer`, was added to the
N_Vector API. This operation is useful for N_Vectors that utilize dual memory
spaces, e.g. the native SUNDIALS CUDA N_Vector.

The SUNMATRIX_CUSPARSE and SUNLINEARSOLVER_CUSOLVERSP_BATCHQR implementations
no longer require the SUNDIALS CUDA N_Vector. Instead, they require that the
vector utilized provides the :c:func:`N_VGetDeviceArrayPointer` operation, and
that the pointer returned by :c:func:`N_VGetDeviceArrayPointer` is a valid CUDA
device pointer.


### CVODE Changes in v5.6.0

A new :c:type:`N_Vector` implementation based on the AMD ROCm HIP platform has
been added. This vector can target NVIDIA or AMD GPUs. See
:numref:`NVectors.hip` for more details. This module is
considered experimental and is subject to change from version to
version.

The :ref:`RAJA N_Vector <NVectors.RAJA>` implementation has been updated to support the HIP
backend in addition to the CUDA backend. Users can choose the backend
when configuring SUNDIALS by using the ``SUNDIALS_RAJA_BACKENDS`` CMake variable. This module
remains experimental and is subject to change from version to version.

A new optional operation, :c:func:`N_VGetDeviceArrayPointer`, was added to the N_Vector API. This
operation is useful for N_Vectors that utilize dual memory spaces, e.g.
the native SUNDIALS CUDA N_Vector.

The :ref:`SUNMATRIX_CUSPARSE <SUNMatrix.cuSparse>` and
:ref:`SUNLINEARSOLVER_CUSOLVERSP_BATCHQR <SUNLinSol.cuSolverSp>`
implementations no longer require the SUNDIALS :ref:`CUDA N_Vector <NVectors.CUDA>`. Instead,
they require that the vector utilized provides the :c:func:`N_VGetDeviceArrayPointer` operation, and that
the pointer returned by :c:func:`N_VGetDeviceArrayPointer` is a valid CUDA device pointer.

### CVODES Changes in v5.6.0

A new ``N_Vector`` implementation based on the AMD ROCm HIP platform has been added. This vector can target NVIDIA or
AMD GPUs. See :numref:`NVectors.hip` for more details. This module is considered experimental and is subject
to change from version to version.

The RAJA ``N_Vector`` implementation has been updated to support the HIP backend in addition to the CUDA backend. Users
can choose the backend when configuring SUNDIALS by using the ``SUNDIALS_RAJA_BACKENDS`` CMake variable. This module
remains experimental and is subject to change from version to version.

A new optional operation, ``N_VGetDeviceArrayPointer``, was added to the N_Vector API. This operation is useful for
N_Vectors that utilize dual memory spaces, e.g. the native SUNDIALS CUDA N_Vector.

The SUNMATRIX_CUSPARSE and SUNLINEARSOLVER_CUSOLVERSP_BATCHQR implementations no longer require the SUNDIALS CUDA
N_Vector. Instead, they require that the vector utilized provides the ``N_VGetDeviceArrayPointer`` operation, and that
the pointer returned by ``N_VGetDeviceArrayPointer`` is a valid CUDA device pointer.

### IDA Changes in v5.6.0

A new ``N_Vector`` implementation based on the AMD ROCm HIP platform has been
added. This vector can target NVIDIA or AMD GPUs. See :numref:`NVectors.Hip` for
more details. This module is considered experimental and is subject to change
from version to version.

The :ref:`NVECTOR_RAJA <NVectors.RAJA>` implementation has been updated to
support the HIP backend in addition to the CUDA backend. Users can choose the
backend when configuring SUNDIALS by using the :cmakeop:`SUNDIALS_RAJA_BACKENDS`
CMake variable. This module remains experimental and is subject to change from
version to version.

A new optional operation, :c:func:`N_VGetDeviceArrayPointer`, was added to the
``N_Vector`` API. This operation is useful for :c:type:`N_Vectors` that utilize
dual memory spaces, e.g. the native SUNDIALS CUDA ``N_Vector``.

The :ref:`SUNMATRIX_CUSPARSE <SUNMatrix.cuSparse>` and
:ref:`SUNLINEARSOLVER_CUSOLVERSP_BATCHQR <SUNLinSol.cuSolverSp>` implementations
no longer require the SUNDIALS CUDA ``N_Vector``. Instead, they require that the
vector utilized provides the :c:func:`N_VGetDeviceArrayPointer` operation, and
that the pointer returned by :c:func:`N_VGetDeviceArrayPointer` is a valid CUDA
device pointer.

### IDAS Changes in v4.6.0

A new ``N_Vector`` implementation based on the AMD ROCm HIP platform has been
added. This vector can target NVIDIA or AMD GPUs. See :numref:`NVectors.Hip` for
more details. This module is considered experimental and is subject to change
from version to version.

The :ref:`NVECTOR_RAJA <NVectors.RAJA>` implementation has been updated to
support the HIP backend in addition to the CUDA backend. Users can choose the
backend when configuring SUNDIALS by using the :cmakeop:`SUNDIALS_RAJA_BACKENDS`
CMake variable. This module remains experimental and is subject to change from
version to version.

A new optional operation, :c:func:`N_VGetDeviceArrayPointer`, was added to the
``N_Vector`` API. This operation is useful for :c:type:`N_Vectors` that utilize
dual memory spaces, e.g. the native SUNDIALS CUDA ``N_Vector``.

The :ref:`SUNMATRIX_CUSPARSE <SUNMatrix.cuSparse>` and
:ref:`SUNLINEARSOLVER_CUSOLVERSP_BATCHQR <SUNLinSol.cuSolverSp>` implementations
no longer require the SUNDIALS CUDA ``N_Vector``. Instead, they require that the
vector utilized provides the :c:func:`N_VGetDeviceArrayPointer` operation, and
that the pointer returned by :c:func:`N_VGetDeviceArrayPointer` is a valid CUDA
device pointer.

### KINSOL Changes in v5.6.0

A new ``N_Vector`` implementation based on the AMD ROCm HIP platform has been added. This vector can target NVIDIA or
AMD GPUs. See :numref:`NVectors.HIP` for more details. This module is considered experimental and is subject
to change from version to version.

The RAJA ``N_Vector`` implementation has been updated to support the HIP backend in addition to the CUDA backend. Users
can choose the backend when configuring SUNDIALS by using the ``SUNDIALS_RAJA_BACKENDS`` CMake variable. This module
remains experimental and is subject to change from version to version.

A new optional operation, ``N_VGetDeviceArrayPointer``, was added to the N_Vector API. This operation is useful for
N_Vectors that utilize dual memory spaces, e.g. the native SUNDIALS CUDA N_Vector.

The SUNMATRIX_CUSPARSE and SUNLINEARSOLVER_CUSOLVERSP_BATCHQR implementations no longer require the SUNDIALS CUDA
N_Vector. Instead, they require that the vector utilized provides the ``N_VGetDeviceArrayPointer`` operation, and that
the pointer returned by ``N_VGetDeviceArrayPointer`` is a valid CUDA device pointer.



## Changes to SUNDIALS in release 5.5.0

Refactored the SUNDIALS build system. CMake 3.12.0 or newer is now required.
Users will likely see deprecation warnings, but otherwise the changes
should be fully backwards compatible for almost all users. SUNDIALS
now exports CMake targets and installs a `SUNDIALSConfig.cmake` file.

Added support for SuperLU DIST 6.3.0 or newer.

### ARKODE Changes in v4.5.0

Refactored the SUNDIALS build system. CMake 3.12.0 or newer is now required.
Users will likely see deprecation warnings, but otherwise the changes
should be fully backwards compatible for almost all users. SUNDIALS
now exports CMake targets and installs a SUNDIALSConfig.cmake file.

Added support for SuperLU DIST 6.3.0 or newer.

### CVODE Changes in v5.5.0

Refactored the SUNDIALS build system. CMake 3.12.0 or newer is now
required. Users will likely see deprecation warnings, but otherwise the
changes should be fully backwards compatible for almost all users.
SUNDIALS now exports CMake targets and installs a
SUNDIALSConfig.cmake file.

Added support for SuperLU DIST 6.3.0 or newer.

### CVODES Changes in v5.5.0

Refactored the SUNDIALS build system. CMake 3.12.0 or newer is now required. Users will likely see deprecation
warnings, but otherwise the changes should be fully backwards compatible for almost all users. SUNDIALS now
exports CMake targets and installs a SUNDIALSConfig.cmake file.

Added support for SuperLU DIST 6.3.0 or newer.

### IDA Changes in v5.5.0

Refactored the SUNDIALS build system. CMake 3.12.0 or newer is now required.
Users will likely see deprecation warnings, but otherwise the changes should be
fully backwards compatible for almost all users. SUNDIALS now exports CMake
targets and installs a ``SUNDIALSConfig.cmake`` file.

Added support for SuperLU_DIST 6.3.0 or newer.

### IDAS Changes in v4.5.0

Refactored the SUNDIALS build system. CMake 3.12.0 or newer is now required.
Users will likely see deprecation warnings, but otherwise the changes should be
fully backwards compatible for almost all users. SUNDIALS now exports CMake
targets and installs a ``SUNDIALSConfig.cmake`` file.

Added support for SuperLU_DIST 6.3.0 or newer.

### KINSOL Changes in v5.5.0

Refactored the SUNDIALS build system. CMake 3.12.0 or newer is now required. Users will likely see deprecation
warnings, but otherwise the changes should be fully backwards compatible for almost all users. SUNDIALS now
exports CMake targets and installs a SUNDIALSConfig.cmake file.

Added support for SuperLU DIST 6.3.0 or newer.



## Changes to SUNDIALS in release 5.4.0

Added full support for time-dependent mass matrices in ARKStep, and expanded
existing non-identity mass matrix infrastructure to support use of the
fixed point nonlinear solver.  Fixed bug for ERK method integration with
static mass matrices.

An interface between ARKStep and the XBraid multigrid reduction in time (MGRIT)
library has been added to enable parallel-in-time integration. See the ARKStep
documentation and examples for more details. This interface required the
addition of three new N_Vector operations to exchange vector data between
computational nodes, see `N_VBufSize`, `N_VBufPack`, and `N_VBufUnpack`. These
N_Vector operations are only used within the XBraid interface and need not be
implemented for any other context.

Updated the MRIStep time-stepping module in ARKODE to support higher-order
MRI-GARK methods [Sandu, SIAM J. Numer. Anal., 57, 2019], including methods that
involve solve-decoupled, diagonally-implicit treatment of the slow time scale.

A new API, `SUNMemoryHelper`, was added to support **GPU users** who have
complex memory management needs such as using memory pools. This is paired with
new constructors for the `NVECTOR_CUDA` and `NVECTOR_RAJA` modules that accept a
`SUNMemoryHelper` object. Refer to "The SUNMemoryHelper API", "NVECTOR CUDA" and
"NVECTOR RAJA" sections in the documentation for more information.

The `NVECTOR_RAJA` module has been updated to mirror the `NVECTOR_CUDA` module.
Notably, the update adds managed memory support to the `NVECTOR_RAJA` module.
Users of the module will need to update any calls to the `N_VMake_Raja` function
because that signature was changed. This module remains experimental and is
subject to change from version to version.

Added new `SetLSNormFactor` functions to CVODE(S), ARKODE, and IDA(S) to
to specify the factor for converting between integrator tolerances (WRMS norm)
and linear solver tolerances (L2 norm) i.e., `tol_L2 = nrmfac * tol_WRMS`.

Added new reset functions `ARKStepReset`, `ERKStepReset`, and
`MRIStepReset` to reset the stepper time and state vector to user-provided
values for continuing the integration from that point while retaining the
integration history. These function complement the reinitialization functions
`ARKStepReInit`, `ERKStepReInit`, and `MRIStepReInit` which reinitialize
the stepper so that the problem integration should resume as if started from
scratch.

Added new functions for advanced users providing a custom `SUNNonlinSolSysFn`.

The expected behavior of `SUNNonlinSolGetNumIters` and
`SUNNonlinSolGetNumConvFails` in the SUNNonlinearSolver API have been updated to
specify that they should return the number of nonlinear solver iterations and
convergence failures in the most recent solve respectively rather than the
cumulative number of iterations and failures across all solves respectively. The
API documentation and SUNDIALS provided SUNNonlinearSolver implementations and
have been updated accordingly. As before, the cumulative number of nonlinear
iterations and failures may be retrieved by calling the integrator provided get
functions.

**This change may cause a runtime error in existing user code**.
In IDAS and CVODES, the functions for forward integration with checkpointing
(`IDASolveF`, `CVodeF`) are now subject to a restriction on the number of time
steps allowed to reach the output time. This is the same restriction applied to
the `IDASolve` and `CVode` functions. The default maximum number of steps is
500, but this may be changed using the `<IDA|CVode>SetMaxNumSteps` function.
This change fixes a bug that could cause an infinite loop in the `IDASolveF`
and `CVodeF` and functions.

A minor inconsistency in CVODE(S) and a bug ARKODE when checking the Jacobian
evaluation frequency has been fixed. As a result codes using using a
non-default Jacobian update frequency through a call to
`CVodeSetMaxStepsBetweenJac` or `ARKStepSetMaxStepsBetweenJac` will need to
increase the provided value by 1 to achieve the same behavior as before. For
greater clarity the functions `CVodeSetMaxStepsBetweenJac`,
`ARKStepSetMaxStepsBetweenJac`, and `ARKStepSetMaxStepsBetweenLSet` have been
deprecated and replaced with `CVodeSetJacEvalFrequency`,
`ARKStepSetJacEvalFrequency`, and `ARKStepSetLSetupFrequency` respectively.
Additionally, the function `CVodeSetLSetupFrequency` has been added to CVODE(S)
to set the frequency of calls to the linear solver setup function.

The `NVECTOR_TRILINOS` module has been updated to work with Trilinos 12.18+.
This update changes the local ordinal type to always be an `int`.

Added support for CUDA v11.

### ARKODE Changes in v4.4.0

Added full support for time-dependent mass matrices in ARKStep, and expanded
existing non-identity mass matrix infrastructure to support use of the
fixed point nonlinear solver. Fixed bug for ERK method integration with
static mass matrices.

An interface between ARKStep and the XBraid multigrid reduction in time (MGRIT)
library :cite:p:`xbraid` has been added to enable parallel-in-time integration. See the
:numref:`ARKODE.Usage.ARKStep.XBraid` section for more information and the example
codes in ``examples/arkode/CXX_xbraid``. This interface required the addition of
three new N_Vector operations to exchange vector data between computational
nodes, see :c:func:`N_VBufSize()`, :c:func:`N_VBufPack()`, and
:c:func:`N_VBufUnpack()`.  These N_Vector operations are only used within the
XBraid interface and need not be implemented for any other context.

Updated the MRIStep time-stepping module in ARKODE to support
higher-order MRI-GARK methods :cite:p:`Sandu:19`, including methods that
involve solve-decoupled, diagonally-implicit treatment of the
slow time scale.

Added the functions :c:func:`ARKStepSetLSNormFactor()`,
:c:func:`ARKStepSetMassLSNormFactor()`, and :c:func:`MRIStepSetLSNormFactor()`
to specify the factor for converting between integrator tolerances (WRMS norm)
and linear solver tolerances (L2 norm) i.e.,
``tol_L2 = nrmfac * tol_WRMS``.

Added new reset functions :c:func:`ARKStepReset()`, :c:func:`ERKStepReset()`,
and :c:func:`MRIStepReset()` to reset the stepper time and state vector to
user-provided values for continuing the integration from that point while
retaining the integration history. These function complement the
reinitialization functions :c:func:`ARKStepReInit()`, :c:func:`ERKStepReInit()`,
and :c:func:`MRIStepReInit()` which reinitialize the stepper so that the problem
integration should resume as if started from scratch.

Added new functions :c:func:`ARKStepComputeState`,
:c:func:`ARKStepGetNonlinearSystemData`, :c:func:`MRIStepComputeState`, and
:c:func:`MRIStepGetNonlinearSystemData` which advanced users might find useful
if providing a custom :c:func:`SUNNonlinSolSysFn`.

The expected behavior of :c:func:`SUNNonlinSolGetNumIters()` and
:c:func:`SUNNonlinSolGetNumConvFails()` in the SUNNonlinearSolver API have been
updated to specify that they should return the number of nonlinear solver
iterations and convergence failures in the most recent solve respectively rather
than the cumulative number of iterations and failures across all solves
respectively. The API documentation and SUNDIALS provided SUNNonlinearSolver
implementations have been updated accordingly. As before, the cumulative number
of nonlinear iterations may be retrieved by calling
:c:func:`ARKStepGetNumNonlinSolvIters()`, the cumulative number of failures with
:c:func:`ARKStepGetNumNonlinSolvConvFails()`, or both with
:c:func:`ARKStepGetNonlinSolvStats()`.

A minor bug in checking the Jacobian evaluation frequency has been fixed. As a
result codes using using a non-default Jacobian update frequency through a call
to :c:func:`ARKStepSetMaxStepsBetweenJac()` will need to increase the provided
value by 1 to achieve the same behavior as before. Additionally, for greater
clarity the functions :c:func:`ARKStepSetMaxStepsBetweenLSet()` and
:c:func:`ARKStepSetMaxStepsBetweenJac()` have been deprecated and replaced with
:c:func:`ARKStepSetLSetupFrequency()` and :c:func:`ARKStepSetJacEvalFrequency()`
respectively.

The ``NVECTOR_RAJA`` module has been updated to mirror the ``NVECTOR_CUDA`` module.
Notably, the update adds managed memory support to the ``NVECTOR_RAJA`` module.
Users of the module will need to update any calls to the ``N_VMake_Raja`` function
because that signature was changed. This module remains experimental and is
subject to change from version to version.

The ``NVECTOR_TRILINOS`` module has been updated to work with Trilinos 12.18+.
This update changes the local ordinal type to always be an ``int``.

Added support for CUDA v11.

### CVODE Changes in v5.4.0

Added new functions :c:func:`CVodeComputeState`, and
:c:func:`CVodeGetNonlinearSystemData` which advanced users might find useful if
providing a custom :c:type:`SUNNonlinSolSysFn`.

Added the function :c:func:`CVodeSetLSNormFactor` to specify the factor for
converting between integrator tolerances (WRMS norm) and linear solver
tolerances (L2 norm) i.e., ``tol_L2 = nrmfac * tol_WRMS``.

The expected behavior of :c:func:`SUNNonlinSolGetNumIters` and
:c:func:`SUNNonlinSolGetNumConvFails` in the :c:type:`SUNNonlinearSolver` API have
been updated to specify that they should return the number of nonlinear solver
iterations and convergence failures in the most recent solve respectively rather
than the cumulative number of iterations and failures across all solves
respectively. The API documentation and SUNDIALS provided :c:type:`SUNNonlinearSolver`
implementations have been updated accordingly. As before, the cumulative number
of nonlinear iterations may be retreived by calling
:c:func:`CVodeGetNumNonlinSolvIters`, the cumulative number of failures with
:c:func:`CVodeGetNumNonlinSolvConvFails`, or both with
:c:func:`CVodeGetNonlinSolvStats`.

A minor inconsistency in checking the Jacobian evaluation frequency has
been fixed. As a result codes using using a non-default Jacobian update
frequency through a call to :c:func:`CVodeSetMaxStepsBetweenJac` will need to increase the provided value by
1 to achieve the same behavior as before. For greater clarity the
function has been deprecated and replaced with :c:func:`CVodeSetJacEvalFrequency`. Additionally, the
function :c:func:`CVodeSetLSetupFrequency` has been added to set the frequency of calls to the linear
solver setup function.

A new class, :ref:`SUNMemoryHelper <SUNMemory>`, was added to support **GPU
users** who have complex memory management needs such as using memory pools.
This is paired with new constructors for the ``NVECTOR_CUDA`` and
``NVECTOR_RAJA`` modules that accept a ``SUNMemoryHelper`` object. Refer to
:numref:`SUNDIALS.GPU`, :numref:`SUNMemory`, :numref:`NVectors.cuda` and
:numref:`NVectors.raja` for more information.

The ``NVECTOR_RAJA`` vector implementation has been updated to mirror the
``NVECTOR_CUDA`` implementation. Notably, the update adds managed memory
support. Users of the vector will need to update any calls to the function
because that signature was changed. This vector remains experimental and is
subject to change from version to version.

The ``NVECTOR_TRILINOS`` vector implementation has been updated to work with
Trilinos 12.18+. This update changes the local ordinal type to always be an
``int``.

### CVODES Changes in v5.4.0

Added the function ``CVodeSetLSNormFactor`` to specify the factor for converting between integrator tolerances (WRMS
norm) and linear solver tolerances (L2 norm) i.e., ``tol_L2 = nrmfac * tol_WRMS``.

Added new functions ``CVodeComputeState``, and ``CVodeGetNonlinearSystemData`` which advanced users might find useful if
providing a custom ``SUNNonlinSolSysFn``.

**This change may cause an error in existing user code**. The ``CVodeF`` function for forward integration with
checkpointing is now subject to a restriction on the number of time steps allowed to reach the output time. This is the
same restriction applied to the ``CVode`` function. The default maximum number of steps is 500, but this may be changed
using the ``CVodeSetMaxNumSteps`` function. This change fixes a bug that could cause an infinite loop in the ``CVodeF``
function.

The expected behavior of ``SUNNonlinSolGetNumIters`` and ``SUNNonlinSolGetNumConvFails`` in the ``SUNNonlinearSolver`` API
have been updated to specify that they should return the number of nonlinear solver iterations and convergence failures
in the most recent solve respectively rather than the cumulative number of iterations and failures across all solves
respectively. The API documentation and SUNDIALS provided ``SUNNonlinearSolver`` implementations have been updated
accordingly. As before, the cumulative number of nonlinear iterations may be retreived by calling
``CVodeGetNumNonlinSolvIters``, ``CVodeGetSensNumNonlinSolvIters``, ``CVodeGetStgrSensNumNonlinSolvIters``, the
cumulative number of failures with ``CVodeGetNumNonlinSolvConvFails``, ``CVodeGetSensNumNonlinSolvConvFails``,
``CVodeGetStgrSensNumNonlinSolvConvFails``, or both with ``CVodeGetNonlinSolvStats``, ``CVodeGetSensNonlinSolvStats``,
``CVodeGetStgrSensNonlinSolvStats``.

A minor inconsistency in checking the Jacobian evaluation frequency has been fixed. As a result codes using using a
non-default Jacobian update frequency through a call to ``CVodeSetMaxStepsBetweenJac`` will need to increase the
provided value by 1 to achieve the same behavior as before. For greater clarity the function
``CVodeSetMaxStepsBetweenJac`` has been deprecated and replaced with ``CVodeSetJacEvalFrequency``. Additionally, the
function ``CVodeSetLSetupFrequency`` has been added to set the frequency of calls to the linear solver setup function.

A new API, ``SUNMemoryHelper``, was added to support **GPU users** who have complex memory management needs such as
using memory pools. This is paired with new constructors for the ``NVECTOR_CUDA`` and ``NVECTOR_RAJA`` modules that accept a
``SUNMemoryHelper`` object. Refer to :numref:`SUNDIALS.GPU.Model`, :numref:`SUNMemory`,
:numref:`NVectors.cuda` and :numref:`NVectors.raja` for more information.

The ``NVECTOR_RAJA`` module has been updated to mirror the ``NVECTOR_CUDA`` module. Notably, the update adds managed
memory support to the ``NVECTOR_RAJA`` module. Users of the module will need to update any calls to the ``N_VMake_Raja``
function because that signature was changed. This module remains experimental and is subject to change from version to
version.

The ``NVECTOR_TRILINOS`` module has been updated to work with Trilinos 12.18+. This update changes the local ordinal
type to always be an ``int``.

Added support for CUDA v11.

### IDA Changes in v5.4.0

Added the function :c:func:`IDASetLSNormFactor` to specify the factor for
converting between integrator tolerances (WRMS norm) and linear solver
tolerances (L2 norm) i.e., ``tol_L2 = nrmfac * tol_WRMS``.

The expected behavior of :c:func:`SUNNonlinSolGetNumIters` and
:c:func:`SUNNonlinSolGetNumConvFails` in the ``SUNNonlinearSolver`` API have
been updated to specify that they should return the number of nonlinear solver
iterations and convergence failures in the most recent solve respectively rather
than the cumulative number of iterations and failures across all solves
respectively. The API documentation and SUNDIALS provided ``SUNNonlinearSolver``
implementations have been updated accordingly. As before, the cumulative number
of nonlinear iterations may be retreived by calling
:c:func:`IDAGetNumNonlinSolvIters`, the cumulative number of failures with
:c:func:`IDAGetNumNonlinSolvConvFails`, or both with
:c:func:`IDAGetNonlinSolvStats`.

A new API, ``SUNMemoryHelper``, was added to support **GPU users** who have
complex memory management needs such as using memory pools. This is paired with
new constructors for the :ref:`NVECTOR_CUDA <NVectors.CUDA>` and
:ref:`NVECTOR_RAJA <NVectors.RAJA>` modules that accept a ``SUNMemoryHelper``
object. Refer to :numref:`SUNDIALS.GPU` and :numref:`SUNMemory` for more
information.

The :ref:`NVECTOR_RAJA <NVectors.RAJA>` module has been updated to mirror the
:ref:`NVECTOR_CUDA <NVectors.CUDA>` module.  Notably, the update adds managed
memory support to the :ref:`NVECTOR_RAJA <NVectors.RAJA>` module.  Users of the
module will need to update any calls to the :c:func:`N_VMake_Raja` function
because that signature was changed. This module remains experimental and is
subject to change from version to version.

The :ref:`NVECTOR_TRILINOS <NVectors.NVTrilinos>` module has been updated to
work with Trilinos 12.18+. This update changes the local ordinal type to always
be an ``int``.

Added support for CUDA v11.

### IDAS Changes in v4.4.0

Added the function :c:func:`IDASetLSNormFactor` to specify the factor for
converting between integrator tolerances (WRMS norm) and linear solver
tolerances (L2 norm) i.e., ``tol_L2 = nrmfac * tol_WRMS``.

Added a new function :c:func:`IDAGetNonlinearSystemData` which advanced users might
find useful if providing a custom :c:type:`SUNNonlinSolSysFn`.

**This change may cause an error in existing user code**. The :c:func:`IDASolveF` function for forward integration with
checkpointing is now subject to a restriction on the number of time steps
allowed to reach the output time. This is the same restriction applied to the
:c:func:`IDASolve` function. The default maximum number of steps is 500, but
this may be changed using the :c:func:`IDASetMaxNumSteps` function. This change
fixes a bug that could cause an infinite loop in the :c:func:`IDASolveF`
function.


The expected behavior of :c:func:`SUNNonlinSolGetNumIters` and
:c:func:`SUNNonlinSolGetNumConvFails` in the ``SUNNonlinearSolver`` API have
been updated to specify that they should return the number of nonlinear solver
iterations and convergence failures in the most recent solve respectively rather
than the cumulative number of iterations and failures across all solves
respectively. The API documentation and SUNDIALS provided ``SUNNonlinearSolver``
implementations have been updated accordingly. As before, the cumulative number
of nonlinear iterations may be retreived by calling
:c:func:`IDAGetNumNonlinSolvIters`, the cumulative number of failures with
:c:func:`IDAGetNumNonlinSolvConvFails`, or both with
:c:func:`IDAGetNonlinSolvStats`.

A new API, ``SUNMemoryHelper``, was added to support **GPU users** who have
complex memory management needs such as using memory pools. This is paired with
new constructors for the :ref:`NVECTOR_CUDA <NVectors.CUDA>` and
:ref:`NVECTOR_RAJA <NVectors.RAJA>` modules that accept a ``SUNMemoryHelper``
object. Refer to :numref:`SUNDIALS.GPU` and :numref:`SUNMemory` for more
information.

The :ref:`NVECTOR_RAJA <NVectors.RAJA>` module has been updated to mirror the
:ref:`NVECTOR_CUDA <NVectors.CUDA>` module.  Notably, the update adds managed
memory support to the :ref:`NVECTOR_RAJA <NVectors.RAJA>` module.  Users of the
module will need to update any calls to the :c:func:`N_VMake_Raja` function
because that signature was changed. This module remains experimental and is
subject to change from version to version.

The :ref:`NVECTOR_TRILINOS <NVectors.NVTrilinos>` module has been updated to
work with Trilinos 12.18+. This update changes the local ordinal type to always
be an ``int``.

Added support for CUDA v11.

### KINSOL Changes in v5.4.0

A new API, ``SUNMemoryHelper``, was added to support **GPU users** who have complex memory management needs such as
using memory pools. This is paired with new constructors for the ``NVECTOR_CUDA`` and ``NVECTOR_RAJA`` modules that accept a
``SUNMemoryHelper`` object. Refer to :numref:`SUNDIALS.GPU.Model`, :numref:`NVectors.CUDA`, :numref:`NVectors.RAJA`, and :numref:`SUNMemory` for more information.

The ``NVECTOR_RAJA`` module has been updated to mirror the ``NVECTOR_CUDA`` module. Notably, the update adds managed
memory support to the ``NVECTOR_RAJA`` module. Users of the module will need to update any calls to the ``N_VMake_Raja``
function because that signature was changed. This module remains experimental and is subject to change from version to
version.

The ``NVECTOR_TRILINOS`` module has been updated to work with Trilinos 12.18+. This update changes the local ordinal
type to always be an ``int``.

Added support for CUDA v11.



## Changes to SUNDIALS in release 5.3.0

Fixed a bug in ARKODE where the prototypes for `ERKStepSetMinReduction` and
`ARKStepSetMinReduction` were not included in `arkode_erkstep.h` and
`arkode_arkstep.h` respectively.

Fixed a bug in ARKODE where inequality constraint checking would need to be
disabled and then re-enabled to update the inequality constraint values after
resizing a problem. Resizing a problem will now disable constraints and a call
to `ARKStepSetConstraints` or `ERKStepSetConstraints` is required to re-enable
constraint checking for the new problem size.

Fixed a bug in the iterative linear solver modules where an error is not
returned if the Atimes function is `NULL` or, if preconditioning is enabled, the
PSolve function is `NULL`.

Added specialized fused CUDA kernels to CVODE which may offer better
performance on smaller problems when using CVODE with the `NVECTOR_CUDA`
module. See the optional input function `CVodeSetUseIntegratorFusedKernels`
for more information. As with other SUNDIALS CUDA features, this is
feature is experimental and may change from version to version.

Added the ability to control the CUDA kernel launch parameters for the
`NVECTOR_CUDA` and `SUNMATRIX_CUSPARSE` modules. These modules remain
experimental and are subject to change from version to version.
In addition, the `NVECTOR_CUDA` kernels were rewritten to be more flexible.
Most users should see equivalent performance or some improvement, but a select
few may observe minor performance degradation with the default settings. Users
are encouraged to contact the SUNDIALS team about any performance changes
that they notice.

Added new capabilities for monitoring the solve phase in the
`SUNNONLINSOL_NEWTON` and `SUNNONLINSOL_FIXEDPOINT` modules, and the SUNDIALS
iterative linear solver modules. SUNDIALS must be built with the CMake option
`SUNDIALS_BUILD_WITH_MONITORING` to use these capabilities.

Added a new function, `CVodeSetMonitorFn`, that takes a user-function
to be called by CVODE after every `nst` successfully completed time-steps.
This is intended to provide a way of monitoring the CVODE statistics
throughout the simulation.

Added a new function `CVodeGetLinSolveStats` to get the CVODE linear solver
statistics as a group.

Added optional set functions to provide an alternative ODE right-hand side
function (ARKODE and CVODE(S)), DAE residual function (IDA(S)), or nonlinear
system function (KINSOL) for use when computing Jacobian-vector products with
the internal difference quotient approximation.

Added support to CVODE for integrating IVPs with constraints using BDF methods
and projecting the solution onto the constraint manifold with a user defined
projection function. This implementation is accompanied by additions to the
CVODE user documentation and examples.

### ARKODE Changes in v4.3.0

Fixed a bug in ARKODE where the prototypes for :c:func:`ERKStepSetMinReduction()`
and :c:func:`ARKStepSetMinReduction()` were not included in ``arkode_erkstep.h``
and ``arkode_arkstep.h`` respectively.

Fixed a bug where inequality constraint checking would need to be disabled and
then re-enabled to update the inequality constraint values after resizing a
problem. Resizing a problem will now disable constraints and a call to
:c:func:`ARKStepSetConstraints()` or :c:func:`ERKStepSetConstraints()` is
required to re-enable constraint checking for the new problem size.

Fixed a bug in the iterative linear solver modules where an error is not
returned if the Atimes function is ``NULL`` or, if preconditioning is enabled,
the PSolve function is ``NULL``.

Added the ability to control the CUDA kernel launch parameters for the
``NVECTOR_CUDA`` and ``SUNMATRIX_CUSPARSE`` modules. These modules remain
experimental and are subject to change from version to version.
In addition, the ``NVECTOR_CUDA`` kernels were rewritten to be more flexible.
Most users should see equivalent performance or some improvement, but a select
few may observe minor performance degradation with the default settings. Users
are encouraged to contact the SUNDIALS team about any perfomance changes
that they notice.

Added the optional function :c:func:`ARKStepSetJacTimesRhsFn()` to specify an
alternative implicit right-hand side function for computing Jacobian-vector
products with the internal difference quotient approximation.

Added new capabilities for monitoring the solve phase in the ``SUNNONLINSOL_NEWTON``
and ``SUNNONLINSOL_FIXEDPOINT`` modules, and the SUNDIALS iterative linear solver
modules. SUNDIALS must be built with the CMake option
``SUNDIALS_BUILD_WITH_MONITORING`` to use these capabilties.


### CVODE Changes in v5.3.0

Fixed a bug in the iterative linear solver modules where an error is not
returned if the Atimes function is ``NULL`` or, if preconditioning is enabled,
the PSolve function is ``NULL``.

Added specialized fused CUDA kernels to CVODE which may offer
better performance on smaller problems when using CVODE with the
``NVECTOR_CUDA`` module. See the optional input function for more
information. As with other SUNDIALS CUDA features, this
capability is considered experimental and may change from version to
version.

Added the ability to control the CUDA kernel launch parameters for the
``NVECTOR_CUDA`` and ``SUNMATRIX_CUSPARSE`` modules. These modules remain
experimental and are subject to change from version to version. In addition, the
kernels were rewritten to be more flexible. Most users should see equivalent
performance or some improvement, but a select few may observe minor performance
degradation with the default settings. Users are encouraged to contact the
SUNDIALS team about any perfomance changes that they notice.

Added new capabilities for monitoring the solve phase in the
``SUNNONLINSOL_NEWTON`` and ``SUNNONLINSOL_FIXEDPOINT`` modules, and the
SUNDIALS iterative linear solver modules. SUNDIALS must be built
with the ``SUNDIALS_BUILD_WITH_MONITORING`` CMake option set to ``TRUE`` to use these capabilties.

Added a new function, :c:func:`CVodeSetMonitorFn`, that takes a user-function to be called by
CVODE after every :math:`nst` succesfully completed time-steps. This
is intended to provide a way of monitoring the CVODE statistics
throughout the simulation.

Added a new function :c:func:`CVodeGetLinSolveStats` to get the CVODE linear solver statistics as a
group.

Added the optional function :c:func:`CVodeSetJacTimsRhsFn` to specify an alternative right-hand side
function for computing Jacobian-vector products with the internal
difference quotient approximation.

Added support for integrating IVPs with constraints using BDF methods
and projecting the solution onto the constraint manifold with a user
defined projection function. This implementation is accompanied by
additions to user documentation and CVODE examples. See
:c:func:`CVodeSetProjFn` for more information.

Added support for CUDA v11.

### CVODES Changes in v5.3.0

Fixed a bug in the iterative linear solver modules where an error is not returned if the Atimes function is ``NULL`` or,
if preconditioning is enabled, the PSolve function is ``NULL``.

Added the ability to control the CUDA kernel launch parameters for the ``NVECTOR_CUDA`` and ``SUNMATRIX_CUSPARSE``
modules. These modules remain experimental and are subject to change from version to version. In addition, the
``NVECTOR_CUDA`` kernels were rewritten to be more flexible. Most users should see equivalent performance or some
improvement, but a select few may observe minor performance degradation with the default settings. Users are encouraged
to contact the SUNDIALS team about any perfomance changes that they notice.

Added new capabilities for monitoring the solve phase in the ``SUNNONLINSOL_NEWTON`` and ``SUNNONLINSOL_FIXEDPOINT``
modules, and the SUNDIALS iterative linear solver modules. SUNDIALS must be built with the CMake option
``SUNDIALS_BUILD_WITH_MONITORING`` to use these capabilties.

Added the optional functions ``CVodeSetJacTimesRhsFn`` and ``CVodeSetJacTimesRhsFnB`` to specify an alternative
right-hand side function for computing Jacobian-vector products with the internal difference quotient approximation.

### IDA Changes in v5.3.0

Fixed a bug in the iterative linear solver modules where an error is not
returned if the ATimes function is ``NULL`` or, if preconditioning is enabled,
the PSolve function is ``NULL``.

Added a new function :c:func:`IDAGetNonlinearSystemData` which advanced users
might find useful if providing a custom :c:type:`SUNNonlinSolSysFn`.

Added the ability to control the CUDA kernel launch parameters for the
:ref:`NVECTOR_CUDA <NVectors.CUDA>` and
:ref:`SUNMATRIX_CUSPARSE <SUNMatrix.cuSparse>` modules. These modules remain
experimental and are subject to change from version to version.  In addition,
the :ref:`NVECTOR_CUDA <NVectors.CUDA>` kernels were rewritten to be more
flexible. Most users should see equivalent performance or some improvement, but
a select few may observe minor performance degradation with the default
settings. Users are encouraged to contact the SUNDIALS team about any
performance changes that they notice.

Added new capabilities for monitoring the solve phase in the
:ref:`SUNNONLINSOL_NEWTON <SUNNonlinSol.Newton>`
and :ref:`SUNNONLINSOL_FIXEDPOINT <SUNNonlinSol.FixedPoint>` modules, and the
SUNDIALS iterative linear solver modules. SUNDIALS must be built with the CMake
option :cmakeop:`SUNDIALS_BUILD_WITH_MONITORING` to use these capabilities.

Added the optional function :c:func:`IDASetJacTimesResFn` to specify an
alternative residual function for computing Jacobian-vector products with the
internal difference quotient approximation.

### IDAS Changes in v4.3.0

Fixed a bug in the iterative linear solver modules where an error is not
returned if the ATimes function is ``NULL`` or, if preconditioning is enabled,
the PSolve function is ``NULL``.

Added a new function :c:func:`IDAGetNonlinearSystemData` which advanced users
might find useful if providing a custom :c:type:`SUNNonlinSolSysFn`.

Added the ability to control the CUDA kernel launch parameters for the
:ref:`NVECTOR_CUDA <NVectors.CUDA>` and
:ref:`SUNMATRIX_CUSPARSE <SUNMatrix.cuSparse>` modules. These modules remain
experimental and are subject to change from version to version.  In addition,
the :ref:`NVECTOR_CUDA <NVectors.CUDA>` kernels were rewritten to be more
flexible. Most users should see equivalent performance or some improvement, but
a select few may observe minor performance degradation with the default
settings. Users are encouraged to contact the SUNDIALS team about any
performance changes that they notice.

Added new capabilities for monitoring the solve phase in the
:ref:`SUNNONLINSOL_NEWTON <SUNNonlinSol.Newton>`
and :ref:`SUNNONLINSOL_FIXEDPOINT <SUNNonlinSol.FixedPoint>` modules, and the
SUNDIALS iterative linear solver modules. SUNDIALS must be built with the CMake
option :cmakeop:`SUNDIALS_BUILD_WITH_MONITORING` to use these capabilities.

Added the optional functions :c:func:`IDASetJacTimesResFn` and
:c:func:`IDASetJacTimesResFnB` to specify an alternative residual function for
computing Jacobian-vector products with the internal difference quotient
approximation.

### KINSOL Changes in v5.3.0

Fixed a bug in the iterative linear solver modules where an error is not returned if the Atimes function is ``NULL`` or,
if preconditioning is enabled, the PSolve function is ``NULL``.

Added the ability to control the CUDA kernel launch parameters for the ``NVECTOR_CUDA`` and ``SUNMATRIX_CUSPARSE``
modules. These modules remain experimental and are subject to change from version to version. In addition, the
``NVECTOR_CUDA`` kernels were rewritten to be more flexible. Most users should see equivalent performance or some
improvement, but a select few may observe minor performance degradation with the default settings. Users are encouraged
to contact the SUNDIALS team about any perfomance changes that they notice.

Added new capabilities for monitoring the solve phase in the ``SUNNONLINSOL_NEWTON`` and ``SUNNONLINSOL_FIXEDPOINT``
modules, and the SUNDIALS iterative linear solver modules. SUNDIALS must be built with the CMake option
``SUNDIALS_BUILD_WITH_MONITORING`` to use these capabilties.

Added the optional function ``KINSetJacTimesVecSysFn`` to specify an alternative system function for computing
Jacobian-vector products with the internal difference quotient approximation.



## Changes to SUNDIALS in release 5.2.0

Fixed a build system bug related to the Fortran 2003 interfaces when using the
IBM XL compiler. When building the Fortran 2003 interfaces with an XL compiler
it is recommended to set `CMAKE_Fortran_COMPILER` to `f2003`, `xlf2003`, or
`xlf2003_r`.

Fixed a bug in how ARKODE interfaces with a user-supplied, iterative, unscaled
linear solver. In this case, ARKODE adjusts the linear solver tolerance in an
attempt to account for the lack of support for left/right scaling matrices.
Previously, ARKODE computed this scaling factor using the error weight vector,
`ewt`; this fix changes that to the residual weight vector, `rwt`, that can
differ from `ewt` when solving problems with non-identity mass matrix.

Fixed a linkage bug affecting Windows users that stemmed from
dllimport/dllexport attribute missing on some SUNDIALS API functions.

Fixed a memory leak in CVODES and IDAS from not deallocating the `atolSmin0` and
`atolQSmin0` arrays.

Fixed a bug where a non-default value for the maximum allowed growth factor
after the first step would be ignored.

Functions were added to each of the time integration packages to enable or
disable the scaling applied to linear system solutions with matrix-based linear
solvers to account for lagged matrix information.

Added two new functions, `ARKStepSetMinReduction` and `ERKStepSetMinReduction`
to change the minimum allowed step size reduction factor after an error test
failure.

Added a new `SUNMatrix` implementation, `SUNMATRIX_CUSPARSE`, that interfaces to
the sparse matrix implementation from the NVIDIA cuSPARSE library. In addition,
the `SUNLINSOL_CUSOLVER_BATCHQR` linear solver has been updated to use this
matrix, therefore, users of this module will need to update their code. These
modules are still considered to be experimental, thus they are subject to
breaking changes even in minor releases.

Added a new "stiff" interpolation module to ARKODE, based on Lagrange polynomial
interpolation, that is accessible to each of the ARKStep, ERKStep and MRIStep
time-stepping modules. This module is designed to provide increased
interpolation accuracy when integrating stiff problems, as opposed to the ARKODE
standard Hermite interpolation module that can suffer when the IVP right-hand
side has large Lipschitz constant. While the Hermite module remains the default,
the new Lagrange module may be enabled using one of the routines
`ARKStepSetInterpolantType`, `ERKStepSetInterpolantType`, or
`MRIStepSetInterpolantType`. The serial example problem `ark_brusselator.c` has
been converted to use this Lagrange interpolation module. Created accompanying
routines `ARKStepSetInterpolantDegree`, `ARKStepSetInterpolantDegree` and
`ARKStepSetInterpolantDegree` to provide user control over these interpolating
polynomials. While the routines `ARKStepSetDenseOrder`, `ARKStepSetDenseOrder`
and `ARKStepSetDenseOrder` still exist, these have been deprecated and will be
removed in a future release.

### ARKODE Changes in v4.2.0

Fixed a build system bug related to the Fortran 2003 interfaces when using the
IBM XL compiler. When building the Fortran 2003 interfaces with an XL compiler
it is recommended to set ``CMAKE_Fortran_COMPILER`` to ``f2003``, ``xlf2003``,
or ``xlf2003_r``.

Fixed a bug in how ARKODE interfaces with a user-supplied, iterative, unscaled linear solver.
In this case, ARKODE adjusts the linear solver tolerance in an attempt to account for the
lack of support for left/right scaling matrices.  Previously, ARKODE computed this scaling
factor using the error weight vector, ``ewt``; this fix changes that to the residual weight vector,
``rwt``, that can differ from ``ewt`` when solving problems with non-identity mass matrix.

Fixed a similar bug in how ARKODE interfaces with scaled linear solvers when solving problems
with non-identity mass matrices.  Here, the left scaling matrix should correspond with ``rwt``
and the right scaling matrix with ``ewt``; these were reversed but are now correct.

Fixed a bug where a non-default value for the maximum allowed growth factor
after the first step would be ignored.

The function :c:func:`ARKStepSetLinearSolutionScaling()` was added to
enable or disable the scaling applied to linear system solutions with
matrix-based linear solvers to account for a lagged value of :math:`\gamma` in
the linear system matrix e.g., :math:`M - \gamma J` or :math:`I - \gamma J`.
Scaling is enabled by default when using a matrix-based linear solver.

Added two new functions, :c:func:`ARKStepSetMinReduction()` and
:c:func:`ERKStepSetMinReduction()`, to change the minimum allowed step size
reduction factor after an error test failure.

Added a new ``SUNMatrix`` implementation, :numref:`SUNMatrix.cuSparse`, that interfaces
to the sparse matrix implementation from the NVIDIA cuSPARSE library. In addition,
the :numref:`SUNLinSol.cuSolverSp` ``SUNLinearSolver`` has been updated to
use this matrix, as such, users of this module will need to update their code.
These modules are still considered to be experimental, thus they are subject to
breaking changes even in minor releases.

Added a new "stiff" interpolation module, based on Lagrange polynomial interpolation,
that is accessible to each of the ARKStep, ERKStep and MRIStep time-stepping modules.
This module is designed to provide increased interpolation accuracy when integrating
stiff problems, as opposed to the ARKODE-standard Hermite interpolation module that
can suffer when the IVP right-hand side has large Lipschitz constant.  While the
Hermite module remains the default, the new Lagrange module may be enabled using one
of the routines :c:func:`ARKStepSetInterpolantType()`, :c:func:`ERKStepSetInterpolantType()`,
or :c:func:`MRIStepSetInterpolantType()`.  The serial example problem ``ark_brusselator.c``
has been converted to use this Lagrange interpolation module.  Created accompanying routines
:c:func:`ARKStepSetInterpolantDegree()`, :c:func:`ARKStepSetInterpolantDegree()` and
:c:func:`ARKStepSetInterpolantDegree()` to provide user control over these
interpolating polynomials.  While the routines :c:func:`ARKStepSetDenseOrder()`,
:c:func:`ARKStepSetDenseOrder()` and :c:func:`ARKStepSetDenseOrder()` still exist,
these have been deprecated and will be removed in a future release.


### CVODE Changes in v5.2.0

Fixed a build system bug related to the Fortran 2003 interfaces when using the
IBM XL compiler. When building the Fortran 2003 interfaces with an XL compiler
it is recommended to set ``CMAKE_Fortran_COMPILER`` to ``f2003``, ``xlf2003``,
or ``xlf2003_r``.

Fixed a linkage bug affecting Windows users that stemmed from
dllimport/dllexport attributes missing on some SUNDIALS API functions.

Added a new ``SUNMatrix`` implementation, ``SUNMATRIX_CUSPARSE``, that
interfaces to the sparse matrix implementation from the NVIDIA cuSPARSE library.
In addition, the linear solver has been updated to use this matrix, therefore,
users of this module will need to update their code. These modules are still
considered to be experimental, thus they are subject to breaking changes even in
minor releases.

The function :c:func:`CVodeSetLinearSolutionScaling` was added to enable or
disable the scaling applied to linear system solutions with matrix-based linear
solvers to account for a lagged value of :math:`\gamma` in the linear system
matrix :math:`I - \gamma J`. Scaling is enabled by default when using a
matrix-based linear solver with BDF methods.

### CVODES Changes in v5.2.0

Fixed a build system bug related to the Fortran 2003 interfaces when using the IBM XL compiler. When building the
Fortran 2003 interfaces with an XL compiler it is recommended to set ``CMAKE_Fortran_COMPILER`` to ``f2003``,
``xlf2003``, or ``xlf2003_r``.

Fixed a linkage bug affecting Windows users that stemmed from dllimport/dllexport attributes missing on some
SUNDIALS API functions.

Fixed a memory leak from not deallocating the ``atolSmin0`` and ``atolQSmin0`` arrays.

Added a new ``SUNMatrix`` implementation, ``SUNMATRIX_CUSPARSE``, that interfaces to the sparse matrix implementation
from the NVIDIA cuSPARSE library. In addition, the ``SUNLINSOL_CUSOLVER_BATCHQR`` linear solver has been updated to use
this matrix, therefore, users of this module will need to update their code. These modules are still considered to be
experimental, thus they are subject to breaking changes even in minor releases.

The functions ``CVodeSetLinearSolutionScaling`` and ``CVodeSetLinearSolutionScalingB`` were added to enable or disable
the scaling applied to linear system solutions with matrix-based linear solvers to account for a lagged value of
:math:`\gamma` in the linear system matrix :math:`I - \gamma J`. Scaling is enabled by default when using a matrix-based
linear solver with BDF methods.

### IDA Changes in v5.2.0

Fixed a build system bug related to the Fortran 2003 interfaces when using the
IBM XL compiler. When building the Fortran 2003 interfaces with an XL compiler
it is recommended to set :cmakeop:`CMAKE_Fortran_COMPILER` to ``f2003``,
``xlf2003``, or ``xlf2003_r``.

Fixed a linkage bug affecting Windows users that stemmed from
dllimport/dllexport attributes missing on some SUNDIALS API functions.

Added a new ``SUNMatrix`` implementation, :ref:`SUNMATRIX_CUSPARSE
<SUNMatrix.cuSparse>`, that interfaces to the sparse matrix implementation from
the NVIDIA cuSPARSE library. In addition, the :ref:`SUNLINSOL_CUSOLVER_BATCHQR
<SUNLinSol.cuSolverSp>` linear solver has been updated to use this matrix,
therefore, users of this module will need to update their code.  These modules
are still considered to be experimental, thus they are subject to breaking
changes even in minor releases.

The function :c:func:`IDASetLinearSolutionScaling` was added to enable or
disable the scaling applied to linear system solutions with matrix-based linear
solvers to account for a lagged value of :math:`\alpha` in the linear system
matrix :math:`J = \frac{\partial F}{\partial y} + \alpha\frac{\partial
F}{\partial \dot{y}}`.  Scaling is enabled by default when using a matrix-based
linear solver.

### IDAS Changes in v4.2.0

Fixed a build system bug related to the Fortran 2003 interfaces when using the
IBM XL compiler. When building the Fortran 2003 interfaces with an XL compiler
it is recommended to set :cmakeop:`CMAKE_Fortran_COMPILER` to ``f2003``,
``xlf2003``, or ``xlf2003_r``.

Fixed a linkage bug affecting Windows users that stemmed from
dllimport/dllexport attributes missing on some SUNDIALS API functions.

Added a new ``SUNMatrix`` implementation, :ref:`SUNMATRIX_CUSPARSE
<SUNMatrix.cuSparse>`, that interfaces to the sparse matrix implementation from
the NVIDIA cuSPARSE library. In addition, the :ref:`SUNLINSOL_CUSOLVER_BATCHQR
<SUNLinSol.cuSolverSp>` linear solver has been updated to use this matrix,
therefore, users of this module will need to update their code.  These modules
are still considered to be experimental, thus they are subject to breaking
changes even in minor releases.

The function :c:func:`IDASetLinearSolutionScaling` and
``IDASetLinearSolutionScalingB`` was added to enable or disable the scaling
applied to linear system solutions with matrix-based linear solvers to account
for a lagged value of :math:`\alpha` in the linear system matrix
:math:`J = \frac{\partial F}{\partial y} + \alpha\frac{\partial F}{\partial \dot{y}}`.
Scaling is enabled by default when using a matrix-based linear solver.

### KINSOL Changes in v5.2.0

Fixed a build system bug related to the Fortran 2003 interfaces when using the IBM XL compiler. When building the
Fortran 2003 interfaces with an XL compiler it is recommended to set ``CMAKE_Fortran_COMPILER`` to ``f2003``,
``xlf2003``, or ``xlf2003_r``.

Fixed a linkage bug affecting Windows users that stemmed from dllimport/dllexport attributes missing on some
SUNDIALS API functions.

Added a new ``SUNMatrix`` implementation, ``SUNMATRIX_CUSPARSE``, that interfaces to the sparse matrix implementation
from the NVIDIA cuSPARSE library. In addition, the ``SUNLINSOL_CUSOLVER_BATCHQR`` linear solver has been updated to use
this matrix, therefore, users of this module will need to update their code. These modules are still considered to be
experimental, thus they are subject to breaking changes even in minor releases.



## Changes to SUNDIALS in release 5.1.0

Added support for a user-supplied function to update the prediction for each
implicit stage solution in ARKStep.  If supplied, this routine will be called
*after* any existing ARKStep predictor algorithm completes, so that the
predictor may be modified by the user as desired.  The new user-supplied routine
has type `ARKStepStagePredictFn`, and may be set by calling
`ARKStepSetStagePredictFn`.

The MRIStep module has been updated to support attaching different user data
pointers to the inner and outer integrators. If applicable, user codes will
need to add a call to `ARKStepSetUserData` to attach their user data
pointer to the inner integrator memory as `MRIStepSetUserData` will
not set the pointer for both the inner and outer integrators. The MRIStep
examples have been updated to reflect this change.

Added support for damping when using Anderson acceleration in KINSOL. See the
mathematical considerations section of the user guide and the description of the
`KINSetDampingAA` function for more details.

Added support for damping to the `SUNNonlinearSolver_FixedPoint` module when
using Anderson acceleration. See the `SUNNonlinearSolver_FixedPoint` section in
the user guides and the description of the `SUNNonlinSolSetDamping_FixedPoint`
function for more details.

Fixed a build system bug related to finding LAPACK/BLAS.

Fixed a build system bug related to checking if the KLU library works.

Fixed a build system bug related to finding PETSc when using the CMake
variables `PETSC_INCLUDES` and `PETSC_LIBRARIES` instead of `PETSC_DIR`.

Added a new build system option, `CUDA_ARCH`, to specify the CUDA architecture
to target.

Fixed a bug in the Fortran 2003 interfaces to the ARKODE Butcher table routines
and structure. This includes changing the `ARKodeButcherTable` type to be a
`type(c_ptr)` in Fortran.

Added two utility functions, `SUNDIALSFileOpen` and `SUNDIALSFileClose` for
creating/destroying file pointers. These are useful when using the Fortran 2003
interfaces.

### ARKODE Changes in v4.1.0

Fixed a build system bug related to finding LAPACK/BLAS.

Fixed a build system bug related to checking if the KLU library works.

Fixed a build system bug related to finding PETSc when using the CMake
variables ``PETSC_INCLUDES`` and ``PETSC_LIBRARIES`` instead of
``PETSC_DIR``.

Added a new build system option, ``CUDA_ARCH``, that can be used to specify
the CUDA architecture to compile for.

Fixed a bug in the Fortran 2003 interfaces to the ARKODE Butcher table routines and structure.
This includes changing the ``ARKodeButcherTable`` type to be a ``type(c_ptr)`` in Fortran.

Added two utility functions, ``SUNDIALSFileOpen`` and ``SUNDIALSFileClose``
for creating/destroying file pointers that are useful when using the Fortran
2003 interfaces.

Added support for a user-supplied function to update the prediction for each
implicit stage solution in ARKStep.  If supplied, this routine will be called
*after* any existing ARKStep predictor algorithm completes, so that the
predictor may be modified by the user as desired.  The new user-supplied routine
has type :c:type:`ARKStepStagePredictFn`, and may be set by calling
:c:func:`ARKStepSetStagePredictFn()`.

The MRIStep module has been updated to support attaching different user data
pointers to the inner and outer integrators. If applicable, user codes will
need to add a call to :c:func:`ARKStepSetUserData()` to attach their user data
pointer to the inner integrator memory as :c:func:`MRIStepSetUserData()` will
not set the pointer for both the inner and outer integrators. The MRIStep
examples have been updated to reflect this change.

Added support for constant damping to the ``SUNNonlinearSolver_FixedPoint``
module when using Anderson acceleration. See :numref:`SUNNonlinSol.FixedPoint.Math`
and the :c:func:`SUNNonlinSolSetDamping_FixedPoint()` for more details.


### CVODE Changes in v5.1.0

Fixed a build system bug related to finding LAPACK/BLAS.

Fixed a build system bug related to checking if the KLU library works.

Fixed a build system bug related to finding PETSc when using the CMake
variables ``PETSC_INCLUDES`` and ``PETSC_LIBRARIES`` instead of ``PETSC_DIR``.

Added a new build system option, ``CUDA_ARCH``, that can be used to specify the CUDA
architecture to compile for.

Added two utility functions, :c:func:`SUNDIALSFileOpen` and
:c:func:`SUNDIALSFileClose` for creating/destroying file pointers that are
useful when using the Fortran 2003 interfaces.

Added support for constant damping to the :ref:`SUNNonlinearSolver_FixedPoint
<SUNNonlinSol.FixedPoint>` module when using Anderson acceleration.

### CVODES Changes in v5.1.0

Fixed a build system bug related to finding LAPACK/BLAS.

Fixed a build system bug related to checking if the KLU library works.

Fixed a build system bug related to finding PETSc when using the CMake
variables ``PETSC_INCLUDES`` and ``PETSC_LIBRARIES`` instead of ``PETSC_DIR``.

Added a new build system option, ``CUDA_ARCH``, that can be used to specify the CUDA
architecture to compile for.

Added two utility functions, :c:func:`SUNDIALSFileOpen` and
:c:func:`SUNDIALSFileClose` for creating/destroying file pointers that are
useful when using the Fortran 2003 interfaces.

Added support for constant damping to the :ref:`SUNNonlinearSolver_FixedPoint
<SUNNonlinSol.FixedPoint>` module when using Anderson acceleration.

### IDA Changes in v5.1.0

Fixed a build system bug related to finding LAPACK/BLAS.

Fixed a build system bug related to checking if the KLU library works.

Fixed a build system bug related to finding PETSc when using the CMake variables
:cmakeop:`PETSC_INCLUDES` and :cmakeop:`PETSC_LIBRARIES` instead of
:cmakeop:`PETSC_DIR`.

Added a new build system option, :cmakeop:`CUDA_ARCH`, that can be used to
specify the CUDA architecture to compile for.

Added two utility functions, :f:func:`FSUNDIALSFileOpen` and
:f:subr:`FSUNDIALSFileClose` for creating/destroying file pointers that are
useful when using the Fortran 2003 interfaces.

### IDAS Changes in v4.1.0

Fixed a build system bug related to finding LAPACK/BLAS.

Fixed a build system bug related to checking if the KLU library works.

Fixed a build system bug related to finding PETSc when using the CMake variables
:cmakeop:`PETSC_INCLUDES` and :cmakeop:`PETSC_LIBRARIES` instead of
:cmakeop:`PETSC_DIR`.

Added a new build system option, :cmakeop:`CUDA_ARCH`, that can be used to
specify the CUDA architecture to compile for.

Added two utility functions, :f:func:`FSUNDIALSFileOpen` and
:f:subr:`FSUNDIALSFileClose` for creating/destroying file pointers that are
useful when using the Fortran 2003 interfaces.

### KINSOL Changes in v5.1.0

Fixed a build system bug related to finding LAPACK/BLAS.

Fixed a build system bug related to checking if the KLU library works.

Fixed a build system bug related to finding PETSc when using the CMake variables ``PETSC_INCLUDES`` and
``PETSC_LIBRARIES`` instead of ``PETSC_DIR``.

Added a new build system option, ``CUDA_ARCH``, that can be used to specify the CUDA architecture to compile for.

Added two utility functions, ``SUNDIALSFileOpen`` and ``SUNDIALSFileClose`` for creating/destroying file pointers that
are useful when using the Fortran 2003 interfaces.

Added support for constant damping when using Anderson acceleration. See :numref:`KINSOL.Mathematics` and the
description of the ``KINSetDampingAA`` function for more details.



## Changes to SUNDIALS in release 5.0.0

### Build System

Increased the minimum required CMake version to 3.5 for most SUNDIALS
configurations, and 3.10 when CUDA or OpenMP with device offloading are enabled.

The CMake option `BLAS_ENABLE` and the variable `BLAS_LIBRARIES` have been
removed to simplify builds as SUNDIALS packages do not use BLAS directly. For
third party libraries that require linking to BLAS, the path to the BLAS
library should be included in the `_LIBRARIES` variable for the third party
library e.g., `SUPERLUDIST_LIBRARIES` when enabling SuperLU_DIST.

Fixed a bug in the build system that prevented the PThreads NVECTOR module from
being built.

### NVector

Two new functions were added to aid in creating custom NVECTOR objects. The
constructor `N_VNewEmpty` allocates an "empty" generic NVECTOR with the object's
content pointer and the function pointers in the operations structure
initialized to NULL. When used in the constructor for custom objects this
function will ease the introduction of any new optional operations to the
NVECTOR API by ensuring only required operations need to be set. Additionally,
the function `N_VCopyOps(w, v)` has been added to copy the operation function
pointers between vector objects. When used in clone routines for custom vector
objects these functions also will ease the introduction of any new optional
operations to the NVECTOR API by ensuring all operations are copied when cloning
objects.

Two new NVECTOR implementations, NVECTOR_MANYVECTOR and NVECTOR_MPIMANYVECTOR,
have been created to support flexible partitioning of solution data among
different processing elements (e.g., CPU + GPU) or for multi-physics problems
that couple distinct MPI-based simulations together (see the NVECTOR_MANYVECTOR
and NVECTOR_MPIMANYVECTOR sections in the user guides for more details). This
implementation is accompanied by additions to user documentation and SUNDIALS
examples.

An additional NVECTOR implementation, NVECTOR_MPIPLUSX, has been created to
support the MPI+X paradigm where X is a type of on-node parallelism (e.g.,
OpenMP, CUDA). The implementation is accompanied by additions to user
documentation and SUNDIALS examples.

One new required vector operation and ten new optional vector operations have
been added to the NVECTOR API. The new required operation, `N_VGetLength`,
returns the global length of an N_Vector. The optional operations have been
added to support the new NVECTOR_MPIMANYVECTOR implementation. The operation
`N_VGetCommunicator` must be implemented by subvectors that are combined to
create an NVECTOR_MPIMANYVECTOR, but is not used outside of this context. The
remaining nine operations are optional local reduction operations intended to
eliminate unnecessary latency when performing vector reduction operations
(norms, etc.) on distributed memory systems. The optional local reduction vector
operations are `N_VDotProdLocal`, `N_VMaxNormLocal`, `N_VMinLocal`,
`N_VL1NormLocal`, `N_VWSqrSumLocal`, `N_VWSqrSumMaskLocal`, `N_VInvTestLocal`,
`N_VConstrMaskLocal`, and `N_VMinQuotientLocal`. If an NVECTOR implementation
defines any of the local operations as NULL, then the NVECTOR_MPIMANYVECTOR will
call standard NVECTOR operations to complete the computation.

The `*_MPICuda` and `*_MPIRaja` functions have been removed from the
NVECTOR_CUDA and NVECTOR_RAJA implementations respectively. Accordingly, the
`nvector_mpicuda.h`, `nvector_mpiraja.h`, `libsundials_nvecmpicuda.lib`, and
`libsundials_nvecmpicudaraja.lib` files have been removed. Users should use the
NVECTOR_MPIPLUSX module in conjunction with the NVECTOR_CUDA or NVECTOR_RAJA
modules to replace the functionality. The necessary changes are minimal and
should require few code modifications.

Fixed a memory leak in the NVECTOR_PETSC clone function.

Made performance improvements to the CUDA NVECTOR. Users who utilize a
non-default stream should no longer see default stream synchronizations after
memory transfers.

Added a new constructor to the CUDA NVECTOR that allows a user to provide
custom allocate and free functions for the vector data array and internal
reduction buffer.

Added new Fortran 2003 interfaces for most NVECTOR modules. See NEVTOR section
in the user guides for more details on how to use the interfaces.

Added three new NVECTOR utility functions, `FN_VGetVecAtIndexVectorArray`,
`FN_VSetVecAtIndexVectorArray`, and `FN_VNewVectorArray`, for working with
`N_Vector` arrays when using the Fortran 2003 interfaces.

### SUNMatrix

Two new functions were added to aid in creating custom SUNMATRIX objects. The
constructor `SUNMatNewEmpty` allocates an "empty" generic SUNMATRIX with the
object's content pointer and the function pointers in the operations structure
initialized to NULL. When used in the constructor for custom objects this
function will ease the introduction of any new optional operations to the
SUNMATRIX API by ensuring only required operations need to be set. Additionally,
the function `SUNMatCopyOps(A, B)` has been added to copy the operation function
pointers between matrix objects. When used in clone routines for custom matrix
objects these functions also will ease the introduction of any new optional
operations to the SUNMATRIX API by ensuring all operations are copied when
cloning objects.

A new operation, `SUNMatMatvecSetup`, was added to the SUNMatrix API. Users
who have implemented custom SUNMatrix modules will need to at least update
their code to set the corresponding ops structure member, matvecsetup, to NULL.

The generic SUNMatrix API now defines error codes to be returned by SUNMatrix
operations. Operations which return an integer flag indiciating success/failure
may return different values than previously.

A new SUNMatrix (and SUNLinearSolver) implementation was added to facilitate
the use of the SuperLU_DIST library with SUNDIALS.

Added new Fortran 2003 interfaces for most SUNMATRIX modules. See SUNMATRIX
section in the user guides for more details on how to use the interfaces.

### SUNLinearSolver

A new function was added to aid in creating custom SUNLINEARSOLVER objects. The
constructor `SUNLinSolNewEmpty` allocates an "empty" generic SUNLINEARSOLVER
with the object's content pointer and the function pointers in the operations
structure initialized to NULL. When used in the constructor for custom objects
this function will ease the introduction of any new optional operations to the
SUNLINEARSOLVER API by ensuring only required operations need to be set.

The return type of the SUNLinearSolver API function `SUNLinSolLastFlag` has
changed from `long int` to `sunindextype` to be consistent with the type
used to store row indices in dense and banded linear solver modules.

Added a new optional operation to the SUNLINEARSOLVER API, `SUNLinSolGetID`,
that returns a `SUNLinearSolver_ID` for identifying the linear solver module.

The SUNLinearSolver API has been updated to make the initialize and setup
functions optional.

A new SUNLinearSolver (and SUNMatrix) implementation was added to facilitate
the use of the SuperLU_DIST library with SUNDIALS.

Added a new SUNLinearSolver implementation,
`SUNLinearSolver_cuSolverSp_batchQR`, which leverages the NVIDIA cuSOLVER sparse
batched QR method for efficiently solving block diagonal linear systems on
NVIDIA GPUs.

Added three new accessor functions to the SUNLinSol_KLU module,
`SUNLinSol_KLUGetSymbolic`, `SUNLinSol_KLUGetNumeric`, and
`SUNLinSol_KLUGetCommon`, to provide user access to the underlying
KLU solver structures.

Added new Fortran 2003 interfaces for most SUNLINEARSOLVER modules. See
SUNLINEARSOLVER section in the user guides for more details on how to use
the interfaces.

### SUNNonlinearSolver

A new function was added to aid in creating custom SUNNONLINEARSOLVER objects.
The constructor `SUNNonlinSolNewEmpty` allocates an "empty" generic
SUNNONLINEARSOLVER with the object's content pointer and the function pointers
in the operations structure initialized to NULL. When used in the constructor
for custom objects this function will ease the introduction of any new optional
operations to the SUNNONLINEARSOLVER API by ensuring only required operations
need to be set.

To facilitate the use of user supplied nonlinear solver convergence test
functions the `SUNNonlinSolSetConvTestFn` function in the SUNNonlinearSolver API
has been updated to take a `void*` data pointer as input. The supplied data
pointer will be passed to the nonlinear solver convergence test function on each
call.

The inputs values passed to the first two inputs of the `SUNNonlinSolSolve`
function in the SUNNONLINEARSOLVER have been changed to be the predicted
state and the initial guess for the correction to that state. Additionally,
the definitions of `SUNNonlinSolLSetupFn` and `SUNNonlinSolLSolveFn` in the
SUNNonlinearSolver API have been updated to remove unused input parameters.
For more information on the nonlinear system formulation and the API functions
see the SUNNONLINEARSOLVER chapter in the user guides.

Added a new `SUNNonlinearSolver` implementation, `SUNNonlinsol_PetscSNES`,
which interfaces to the PETSc SNES nonlinear solver API.

Added new Fortran 2003 interfaces for most SUNNONLINEARSOLVER modules. See
SUNNONLINEARSOLVER section in the user guides for more details on how to use
the interfaces.

### CVODE and CVODES

Fixed a bug in the CVODE and CVODES constraint handling where the step size
could be set below the minimum step size.

Fixed a bug in the CVODE and CVODES nonlinear solver interfaces where the norm
of the accumulated correction was not updated when using a non-default
convergence test function.

Fixed a bug in the CVODES `cvRescale` function where the loops to compute the
array of scalars for the fused vector scale operation stopped one iteration
early.

Fixed a bug in CVODES where CVodeF would return the wrong flag under certain
cirumstances.

Fixed a bug in CVODES where CVodeF would not return a root in NORMAL_STEP mode
if the root occurred after the desired output time.

Fixed a memeory leak in FCVODE when not using the default nonlinear solver.

Removed extraneous calls to `N_VMin` for simulations where the scalar valued
absolute tolerance, or all entries of the vector-valued absolute tolerance
array, are strictly positive. In this scenario CVODE and CVODES will remove
at least one global reduction per time step.

The CVLS interface has been updated to only zero the Jacobian matrix before
calling a user-supplied Jacobian evaluation function when the attached linear
solver has type `SUNLINEARSOLVER_DIRECT`.

A new linear solver interface function, `CVLsLinSysFn`, was added as an
alternative method for evaluating the linear systems I - gamma J.

Added functions to get the current state and gamma value to CVODE and CVODES.
These functions may be useful to users who chose to provide their own nonlinear
solver implementation.

Added New Fortran 2003 interfaces to CVODE and CVODES were added. These new
interfaces were generated with SWIG-Fortran and provide a user an idiomatic
Fortran 2003 interface to most of the SUNDIALS C API. The existing CVODE F2003
interface, and all module implementations with existing Fortran 2003 interfaces
were updated accordingly. See the section "Using CVODE for Fortran
Applications" and "Using CVODES for Fortran Applications" in the appropriate
user guide for more details on how to use the interfaces.

### ARKODE

The MRIStep module has been updated to support explicit, implicit, or IMEX
methods as the fast integrator using the ARKStep module. As a result some
function signatures have been changed including MRIStepCreate which now
takes an ARKStep memory structure for the fast integration as an input.

Fixed a bug in the ARKStep time-stepping module in ARKODE that would result in
an infinite loop if the nonlinear solver failed to converge more than the
maximum allowed times during a single step.

Fixed a bug in ARKODE that would result in a "too much accuracy requested" error
when using fixed time step sizes with explicit methods in some cases.

Fixed a bug in ARKStep where the mass matrix linear solver setup function was
not called in the Matrix-free case.

Fixed a minor bug in ARKStep where an incorrect flag is reported when an
error occurs in the mass matrix setup or Jacobian-vector product setup
functions.

Fixed a memeory leak in FARKODE when not using the default nonlinear solver.

The reinitialization functions `ERKStepReInit`, `ARKStepReInit`, and
`MRIStepReInit` have been updated to retain the minimum and maxiumum step
size values from before reinitialization rather than resetting them to the
default values.

Removed extraneous calls to `N_VMin` for simulations where the scalar valued
absolute tolerance, or all entries of the vector-valued absolute tolerance
array, are strictly positive. In this scenario ARKODE steppers will remove
at least one global reduction per time step.

The ARKLS interface has been updated to only zero the Jacobian matrix before
calling a user-supplied Jacobian evaluation function when the attached linear
solver has type `SUNLINEARSOLVER_DIRECT`.

A new linear solver interface function, `ARKLsLinSysFn`, was added as an
alternative method for evaluating the linear systems M - gamma J and
I - gamma J.

Added two new embedded ARK methods of orders 4 and 5 to ARKODE (from
Kennedy & Carpenter, Appl. Numer. Math., 136:183--205, 2019).

Support for optional inequality constraints on individual components of the
solution vector has been added the ARKODE ERKStep and ARKStep modules. See
the descriptions of `ERKStepSetConstraints` and `ARKStepSetConstraints` for
more details. Note that enabling constraint handling requires the NVECTOR
operations `N_VMinQuotient`, `N_VConstrMask`, and `N_VCompare` that were not
previously required by ARKODE.

Added functions to get the current state and gamma value to the ARKStep module.
These functions may be useful to users who chose to provide their own nonlinear
solver implementation.

Add two new 'Set' functions to MRIStep, `MRIStepSetPreInnerFn` and
`MRIStepSetPostInnerFn` for performing communication or memory
transfers needed before or after the inner integration.

Added new Fortran 2003 interfaces to all ARKODE stepper modules. These new
interfaces were generated with SWIG-Fortran and provide a user an idiomatic
Fortran 2003 interface to most of the SUNDIALS C API. See the section "Using
ARKODE for Fortran Applications" in the user guide for more details on how
to use the interfaces.

### IDA and IDAS

A bug was fixed in the IDA and IDAS linear solver interfaces where an incorrect
Jacobian-vector product increment was used with iterative solvers other than
SPGMR and SPFGMR.

Fixed a bug in IDAS where IDASolveF would return the wrong flag under certain
cirumstances.

Fixed a bug in IDAS where IDASolveF would not return a root in NORMAL_STEP mode
if the root occurred after the desired output time.

Fixed a bug the IDAS IDAQuadReInitB function where an incorrect memory structure
was passed to IDAQuadReInit.

Fixed a memeory leak in FIDA when not using the default nonlinear solver.

Removed extraneous calls to `N_VMin` for simulations where the scalar valued
absolute tolerance, or all entries of the vector-valued absolute tolerance
array, are strictly positive. In this scenario IDA and IDAS will remove
at least one global reduction per time step.

The IDALS interface has been updated to only zero the Jacobian matrix before
calling a user-supplied Jacobian evaluation function when the attached linear
solver has type SUNLINEARSOLVER_DIRECT.

Added new Fortran 2003 interfaces to IDA and IDAS. These new interfaces were
generated with SWIG-Fortran and provide a user an idiomatic Fortran 2003
interface to most of the SUNDIALS C API.  See the section "Using IDA for Fortran
Applications" and "Using IDAS for Fortran Applications" in the appropriate
user guide for more details on how to use the interfaces.

### KINSOL

Fixed a bug in the KINSOL linear solver interface where the auxiliary scalar
`sJpnorm` was not computed when necessary with the Picard iteration and the
auxiliary scalar `sFdotJp` was unnecessarily computed in some cases.

The KINLS interface has been updated to only zero the Jacobian matrix before
calling a user-supplied Jacobian evaluation function when the attached linear
solver has type SUNLINEARSOLVER_DIRECT.

Added new Fortran 2003 interfaces to KINSOL. These new interfaces were
generated with SWIG-Fortran and provide a user an idiomatic Fortran 2003
interface to most of the SUNDIALS C API.  See the section "Using KINSOL for
Fortran Applications" for more details on how to use the interfaces.

### ARKODE Changes in v4.0.0

**Build system changes**

Increased the minimum required CMake version to 3.5 for most SUNDIALS
configurations, and 3.10 when CUDA or OpenMP with device offloading are enabled.

The CMake option ``BLAS_ENABLE`` and the variable ``BLAS_LIBRARIES`` have been
removed to simplify builds as SUNDIALS packages do not use BLAS directly. For
third party libraries that require linking to BLAS, the path to the BLAS
library should be included in the ``_LIBRARIES`` variable for the third party
library e.g., ``SUPERLUDIST_LIBRARIES`` when enabling SuperLU_DIST.

Fixed a bug in the build system that prevented the PThreads NVECTOR module from
being built.

**NVECTOR module changes**

Two new functions were added to aid in creating custom NVECTOR objects. The
constructor :c:func:`N_VNewEmpty` allocates an "empty" generic NVECTOR with
the object's content pointer and the function pointers in the operations
structure initialized to ``NULL``. When used in the constructor for custom
objects this function will ease the introduction of any new optional operations
to the NVECTOR API by ensuring only required operations need to be set.
Additionally, the function :c:func:`N_VCopyOps()` has been added to copy the
operation function pointers between vector objects. When used in clone routines
for custom vector objects these functions also will ease the introduction of
any new optional operations to the NVECTOR API by ensuring all operations
are copied when cloning objects.

Two new NVECTOR implementations, NVECTOR_MANYVECTOR and
NVECTOR_MPIMANYVECTOR, have been created to support flexible partitioning
of solution data among different processing elements (e.g., CPU + GPU) or for
multi-physics problems that couple distinct MPI-based simulations together. This
implementation is accompanied by additions to user documentation and SUNDIALS
examples.

One new required vector operation and ten new optional vector operations have
been added to the NVECTOR API. The new required operation, :c:func:`N_VGetLength()`,
returns the global length of an ``N_Vector``. The optional operations have
been added to support the new NVECTOR_MPIMANYVECTOR implementation. The
operation :c:func:`N_VGetCommunicator()` must be implemented by subvectors that are
combined to create an NVECTOR_MPIMANYVECTOR, but is not used outside of
this context. The remaining nine operations are optional local reduction
operations intended to eliminate unnecessary latency when performing vector
reduction operations (norms, etc.) on distributed memory systems. The optional
local reduction vector operations are
:c:func:`N_VDotProdLocal`,
:c:func:`N_VMaxNormLocal`,
:c:func:`N_VMinLocal`,
:c:func:`N_VL1NormLocal`,
:c:func:`N_VWSqrSumLocal`,
:c:func:`N_VWSqrSumMaskLocal`,
:c:func:`N_VInvTestLocal`,
:c:func:`N_VConstrMaskLocal`, and
:c:func:`N_VMinQuotientLocal`.
If an NVECTOR implementation defines any of the local operations as
``NULL``, then the NVECTOR_MPIMANYVECTOR will call standard NVECTOR
operations to complete the computation.

An additional NVECTOR implementation, NVECTOR_MPIPLUSX, has been created to
support the MPI+X paradigm where X is a type of on-node parallelism
(*e.g.*, OpenMP, CUDA). The implementation is accompanied by additions to
user documentation and SUNDIALS examples.

The ``*_MPICuda`` and ``*_MPIRaja`` functions have been removed from the
NVECTOR_CUDA and NVECTOR_RAJA implementations respectively. Accordingly, the
``nvector_mpicuda.h``, ``nvector_mpiraja.h``, ``libsundials_nvecmpicuda.lib``,
and ``libsundials_nvecmpicudaraja.lib`` files have been removed. Users should
use the NVECTOR_MPIPLUSX module coupled in conjunction with the NVECTOR_CUDA
or NVECTOR_RAJA modules to replace the functionality. The necessary changes are
minimal and should require few code modifications. See the programs in
``examples/ida/mpicuda`` and ``examples/ida/mpiraja`` for examples of how to
use the NVECTOR_MPIPLUSX module with the NVECTOR_CUDA and NVECTOR_RAJA modules
respectively.

Fixed a memory leak in the NVECTOR_PETSC module clone function.

Made performance improvements to the NVECTOR_CUDA module. Users who utilize a
non-default stream should no longer see default stream synchronizations
after memory transfers.

Added a new constructor to the NVECTOR_CUDA module that allows a user to provide
custom allocate and free functions for the vector data array and internal
reduction buffer.

Added new Fortran 2003 interfaces for most NVECTOR modules. See the
:numref:`SUNDIALS.Fortran` section for more details.

Added three new NVECTOR utility functions,
:c:func:`N_VGetVecAtIndexVectorArray()`
:c:func:`N_VSetVecAtIndexVectorArray()`, and
:c:func:`N_VNewVectorArray`,
for working with ``N_Vector`` arrays when using the Fortran 2003 interfaces.

**SUNMatrix module changes**

Two new functions were added to aid in creating custom SUNMATRIX objects. The
constructor :c:func:`SUNMatNewEmpty` allocates an "empty" generic SUNMATRIX with
the object's content pointer and the function pointers in the operations
structure initialized to ``NULL``. When used in the constructor for custom
objects this function will ease the introduction of any new optional operations
to the SUNMATRIX API by ensuring only required operations need to be set.
Additionally, the function :c:func:`SUNMatCopyOps()` has been added to copy the
operation function pointers between matrix objects. When used in clone routines
for custom matrix objects these functions also will ease the introduction of any
new optional operations to the SUNMATRIX API by ensuring all operations are
copied when cloning objects.

A new operation, :c:func:`SUNMatMatvecSetup()`, was added to the SUNMATRIX API.
Users who have implemented custom SUNMATRIX modules will need to at least
update their code to set the corresponding ``ops`` structure member,
``matvecsetup``, to ``NULL``.

A new operation, :c:func:`SUNMatMatvecSetup()`, was added to the SUNMATRIX API
to perform any setup necessary for computing a matrix-vector product. This
operation is useful for SUNMATRIX implementations which need to prepare the
matrix itself, or communication structures before performing the matrix-vector
product. Users who have implemented custom SUNMATRIX modules will need to at
least update their code to set the corresponding ``ops`` structure member,
``matvecsetup``, to ``NULL``.

The generic SUNMATRIX API now defines error codes to be returned by
SUNMATRIX operations. Operations which return an integer flag indiciating
success/failure may return different values than previously.

A new SUNMATRIX (and SUNLINEARSOLVER) implementation was added to
facilitate the use of the SuperLU_DIST library with SUNDIALS.

Added new Fortran 2003 interfaces for most SUNMATRIX modules. See the
:numref:`SUNDIALS.Fortran` section for more details.

**SUNLinearSolver module changes**

A new function was added to aid in creating custom SUNLINEARSOLVER objects.
The constructor :c:func:`SUNLinSolNewEmpty` allocates an "empty" generic
SUNLINEARSOLVER with the object's content pointer and the function pointers
in the operations structure initialized to ``NULL``. When used in the
constructor for custom objects this function will ease the introduction of any
new optional operations to the SUNLINEARSOLVER API by ensuring only required
operations need to be set.

The return type of the SUNLINEARSOLVER API function :c:func:`SUNLinSolLastFlag()`
has changed from ``long int`` to ``sunindextype`` to be consistent with the
type used to store row indices in dense and banded linear solver modules.

Added a new optional operation to the SUNLINEARSOLVER API,
:c:func:`SUNLinSolGetID`, that returns a ``SUNLinearSolver_ID`` for identifying
the linear solver module.

The SUNLINEARSOLVER API has been updated to make the initialize and setup
functions optional.

A new SUNLINEARSOLVER (and SUNMATRIX) implementation was added to
facilitate the use of the SuperLU_DIST library with SUNDIALS.

Added a new SUNLinearSolver implementation, ``SUNLinearSolver_cuSolverSp_batchQR``,
which leverages the NVIDIA cuSOLVER sparse batched QR method for efficiently
solving block diagonal linear systems on NVIDIA GPUs.


Added three new accessor functions to the SUNLinSol_KLU module,
:c:func:`SUNLinSol_KLUGetSymbolic()`, :c:func:`SUNLinSol_KLUGetNumeric()`, and
:c:func:`SUNLinSol_KLUGetCommon()`, to provide user access to the underlying
KLU solver structures.

Added new Fortran 2003 interfaces for most SUNLINEARSOLVER modules. See the
:numref:`SUNDIALS.Fortran` section for more details.

**SUNNonlinearSolver module changes**

A new function was added to aid in creating custom SUNNONLINEARSOLVER
objects. The constructor :c:func:`SUNNonlinSolNewEmpty` allocates an "empty"
generic SUNNONLINEARSOLVER with the object's content pointer and the function
pointers in the operations structure initialized to ``NULL``. When used in the
constructor for custom objects this function will ease the introduction of any
new optional operations to the SUNNONLINEARSOLVER API by ensuring only
required operations need to be set.

To facilitate the use of user supplied nonlinear solver convergence test
functions the :c:func:`SUNNonlinSolSetConvTestFn()` function in the
SUNNONLINEARSOLVER API has been updated to take a ``void*`` data pointer as
input. The supplied data pointer will be passed to the nonlinear solver
convergence test function on each call.

The inputs values passed to the first two inputs of the :c:func:`SUNNonlinSolSolve()`
function in the SUNNONLINEARSOLVER have been changed to be the predicted
state and the initial guess for the correction to that state. Additionally,
the definitions of :c:type:`SUNNonlinSolLSetupFn` and :c:type:`SUNNonlinSolLSolveFn`
in the SUNNONLINEARSOLVER API have been updated to remove unused input
parameters.

Added a new ``SUNNonlinearSolver`` implementation, ``SUNNonlinsol_PetscSNES``,
which interfaces to the PETSc SNES nonlinear solver API.

Added new Fortran 2003 interfaces for most SUNNONLINEARSOLVER modules. See the
:numref:`SUNDIALS.Fortran` section for more details.

**ARKODE changes**

The MRIStep module has been updated to support explicit, implicit, or ImEx
methods as the fast integrator using the ARKStep module. As a result some
function signatures have been changed including :c:func:`MRIStepCreate` which
now takes an ARKStep memory structure for the fast integration as an input.

Fixed a bug in the ARKStep time-stepping module that would result in an infinite
loop if the nonlinear solver failed to converge more than the maximum allowed times
during a single step.

Fixed a bug that would result in a "too much accuracy requested" error when
using fixed time step sizes with explicit methods in some cases.

Fixed a bug in ARKStep where the mass matrix linear solver setup function was
not called in the Matrix-free case.

Fixed a minor bug in ARKStep where an incorrect flag is reported when an
error occurs in the mass matrix setup or Jacobian-vector product setup
functions.

Fixed a memeory leak in FARKODE when not using the default nonlinear solver.

The reinitialization functions :c:func:`ERKStepReInit()`,
:c:func:`ARKStepReInit()`, and :c:func:`MRIStepReInit()` have been updated to
retain the minimum and maxiumum step size values from before reinitialization
rather than resetting them to the default values.

Removed extraneous calls to :c:func:`N_VMin()` for simulations where
the scalar valued absolute tolerance, or all entries of the
vector-valued absolute tolerance array, are strictly positive.  In
this scenario, ARKODE will remove at least one global reduction per
time step.

The ARKLS interface has been updated to only zero the Jacobian matrix before
calling a user-supplied Jacobian evaluation function when the attached linear
solver has type ``SUNLINEARSOLVER_DIRECT``.

A new linear solver interface function :c:func:`ARKLsLinSysFn` was added as an
alternative method for evaluating the linear system :math:`A = M - \gamma J`.

Added two new embedded ARK methods of orders 4 and 5 to ARKODE (from :cite:p:`KenCarp:19`).

Support for optional inequality constraints on individual components of the
solution vector has been added the ARKODE ERKStep and ARKStep modules. See
the descriptions of :c:func:`ERKStepSetConstraints()` and
:c:func:`ARKStepSetConstraints()` for more details. Note that enabling
constraint handling requires the NVECTOR operations :c:func:`N_VMinQuotient()`,
:c:func:`N_VConstrMask()`, and :c:func:`N_VCompare()` that were not previously
required by ARKODE.

Added two new 'Get' functions to ARKStep, :c:func:`ARKStepGetCurrentGamma()`,
and :c:func:`ARKStepGetCurrentState`, that may be useful to users who choose
to provide their own nonlinear solver implementation.

Add two new 'Set' functions to MRIStep, :c:func:`MRIStepSetPreInnerFn()` and
:c:func:`MRIStepSetPostInnerFn()` for performing communication or memory
transfers needed before or after the inner integration.

A new Fortran 2003 interface to ARKODE was added. This includes Fortran 2003 interfaces
to the ARKStep, ERKStep, and MRIStep time-stepping modules. See the
:numref:`SUNDIALS.Fortran` section for more details.


### CVODE Changes in v5.0.0

**Build system changes**

-  Increased the minimum required CMake version to 3.5 for most
   SUNDIALS configurations, and 3.10 when CUDA or OpenMP with device
   offloading are enabled.

-  The CMake option ``BLAS_ENABLE`` and the variable ``BLAS_LIBRARIES`` have been removed to simplify
   builds as SUNDIALS packages do not use BLAS directly. For third
   party libraries that require linking to BLAS, the path to the BLAS
   library should be included in the variable for the third party
   library *e.g.*, ``SUPERLUDIST_LIBRARIES`` when enabling SuperLU_DIST.

-  Fixed a bug in the build system that prevented the ``NVECTOR_PTHREADS``
   module from being built.

**NVECTOR module changes**

-  Two new functions were added to aid in creating custom ``N_Vector``
   objects. The constructor :c:func:`N_VNewEmpty` allocates an "empty" generic ``N_Vector``
   with the object's content pointer and the function pointers in the
   operations structure initialized to  ``NULL``. When used in the constructor
   for custom objects this function will ease the introduction of any
   new optional operations to the ``N_Vector`` API by ensuring only
   required operations need to be set. Additionally, the function :c:func:`N_VCopyOps` has
   been added to copy the operation function pointers between vector
   objects. When used in clone routines for custom vector objects these
   functions also will ease the introduction of any new optional
   operations to the ``N_Vector`` API by ensuring all operations are
   copied when cloning objects. See :numref:`NVectors.Description.custom_implementation` for more details.

-  Two new ``N_Vector`` implementations, ``NVECTOR_MANYVECTOR`` and
   ``NVECTOR_MPIMANYVECTOR``, have been created to support flexible
   partitioning of solution data among different processing elements
   (e.g., CPU + GPU) or for multi-physics problems that couple distinct
   MPI-based simulations together. This implementation is accompanied by
   additions to user documentation and SUNDIALS examples. See
   :numref:`NVectors.manyvector` and :numref:`NVectors.mpimanyvector` for more
   details.

-  One new required vector operation and ten new optional vector
   operations have been added to the ``N_Vector`` API. The new required
   operation, , returns the global length of an . The optional operations have
   been added to support the new ``NVECTOR_MPIMANYVECTOR`` implementation. The
   operation must be implemented by subvectors that are combined to create an
   ``NVECTOR_MPIMANYVECTOR``, but is not used outside of this context. The
   remaining nine operations are optional local reduction operations intended to
   eliminate unnecessary latency when performing vector reduction operations
   (norms, etc.) on distributed memory systems. The optional local reduction
   vector operations are :c:func:`N_VDotProdLocal`, :c:func:`N_VMaxNormLocal`,
   :c:func:`N_VL1NormLocal`, :c:func:`N_VWSqrSumLocal`,
   :c:func:`N_VWSqrSumMaskLocal`, :c:func:`N_VInvTestLocal`,
   :c:func:`N_VConstrMaskLocal`, :c:func:`N_VMinLocal`, and
   :c:func:`N_VMinQuotientLocal`. If an ``N_Vector`` implementation defines any
   of the local operations as , then the ``NVECTOR_MPIMANYVECTOR`` will call
   standard ``N_Vector`` operations to complete the computation.

-  An additional ``N_Vector`` implementation, ``NVECTOR_MPIPLUSX``, has been
   created to support the MPI+X paradigm where X is a type of on-node
   parallelism (*e.g.*, OpenMP, CUDA). The implementation is accompanied
   by additions to user documentation and SUNDIALS examples. See
   :numref:`NVectors.mpiplusx` for more details.

-  The and functions have been removed from the ``NVECTOR_CUDA`` and
   ``NVECTOR_RAJA`` implementations respectively. Accordingly, the
   ``nvector_mpicuda.h``, ``libsundials_nvecmpicuda.lib``,
   ``libsundials_nvecmpicudaraja.lib``, and files have been removed. Users
   should use the ``NVECTOR_MPIPLUSX`` module coupled in conjunction with the
   ``NVECTOR_CUDA`` or ``NVECTOR_RAJA`` modules to replace the functionality.
   The necessary changes are minimal and should require few code modifications.
   See the programs in and for examples of how to use the ``NVECTOR_MPIPLUSX``
   module with the ``NVECTOR_CUDA`` and ``NVECTOR_RAJA`` modules respectively.

-  Fixed a memory leak in the ``NVECTOR_PETSC`` module clone function.

-  Made performance improvements to the ``NVECTOR_CUDA`` module. Users who
   utilize a non-default stream should no longer see default stream
   synchronizations after memory transfers.

-  Added a new constructor to the ``NVECTOR_CUDA`` module that allows a user
   to provide custom allocate and free functions for the vector data
   array and internal reduction buffer. See :numref:`NVectors.Cuda` for more details.

-  Added new Fortran 2003 interfaces for most ``N_Vector`` modules. See
   :numref:`NVectors` for more details on how to use
   the interfaces.

-  Added three new ``N_Vector`` utility functions :c:func:`N_VGetVecAtIndexVectorArray`,
   :c:func:`N_VSetVecAtIndexVectorArray`, and :c:func:`N_VNewVectorArray` for working
   with arrays when using the Fortran 2003 interfaces.

**SUNMatrix module changes**

-  Two new functions were added to aid in creating custom ``SUNMatrix``
   objects. The constructor :c:func:`SUNMatNewEmpty` allocates an "empty" generic ``SUNMatrix``
   with the object's content pointer and the function pointers in the
   operations structure initialized to . When used in the constructor
   for custom objects this function will ease the introduction of any
   new optional operations to the ``SUNMatrix`` API by ensuring only
   required operations need to be set. Additionally, the function :c:func:`SUNMatCopyOps` has
   been added to copy the operation function pointers between matrix
   objects. When used in clone routines for custom matrix objects these
   functions also will ease the introduction of any new optional
   operations to the ``SUNMatrix`` API by ensuring all operations are
   copied when cloning objects. See :numref:`SUNMatrix` for more
   details.
-  A new operation, :c:func:`SUNMatMatvecSetup`, was added to the ``SUNMatrix`` API to perform any
   setup necessary for computing a matrix-vector product. This operation
   is useful for ``SUNMatrix`` implementations which need to prepare the
   matrix itself, or communication structures before performing the
   matrix-vector product. Users who have implemented custom
   ``SUNMatrix`` modules will need to at least update their code to set
   the corresponding structure member to ``NULL``. See :numref:`SUNMatrix.Ops`
   for more details.
-  The generic ``SUNMatrix`` API now defines error codes to be returned
   by ``SUNMatrix`` operations. Operations which return an integer flag
   indiciating success/failure may return different values than
   previously.
-  A new ``SUNMatrix`` (and ``SUNLinearSolver``) implementation was added to
   facilitate the use of the SuperLU_DIST library with SUNDIALS. See
   :numref:`SUNMatrix.SLUNRloc` for more details.
-  Added new Fortran 2003 interfaces for most ``SUNMatrix`` modules. See
   :numref:`SUNMatrix` for more details on how to
   use the interfaces.

**SUNLinearSolver module changes**

-  A new function was added to aid in creating custom ``SUNLinearSolver``
   objects. The constructor allocates an "empty" generic ``SUNLinearSolver``
   with the object's content pointer and the function pointers in the operations
   structure initialized to . When used in the constructor for custom objects
   this function will ease the introduction of any new optional operations to
   the ``SUNLinearSolver`` API by ensuring only required operations need to be
   set. See :numref:`SUNLinSol.API.Custom` for more details.
-  The return type of the ``SUNLinearSolver`` API function has changed from to
   to be consistent with the type used to store row indices in dense and banded
   linear solver modules.
-  Added a new optional operation to the ``SUNLinearSolver`` API,
   :c:func:`SUNLinSolLastFlag`, that returns a for identifying the linear solver module.
-  The ``SUNLinearSolver`` API has been updated to make the initialize and
   setup functions optional.
-  A new ``SUNLinearSolver`` (and ``SUNMatrix``) implementation was added to
   facilitate the use of the SuperLU_DIST library with SUNDIALS. See
   :numref:`SUNLinSol.SuperLUDIST` for more details.
-  Added a new ``SUNLinearSolver`` implementation, :ref:`SUNLINEARSOLVER_CUSOLVERSP <SUNLinSol.cuSolverSp>`,
   which leverages the NVIDIA cuSOLVER sparse batched QR method for efficiently solving block
   diagonal linear systems on NVIDIA GPUs.
-  Added three new accessor functions to the ``SUNLINSOL_KLU`` module, :c:func:`SUNLinSol_KLUGetSymbolic`,
   , :c:func:`SUNLinSol_KLUGetNumeric` and :c:func:`SUNLinSol_KLUGetCommon`, to
   provide user access to the underlying KLU solver structures. See
   :numref:`SUNLinSol.KLU` for more details.
-  Added new Fortran 2003 interfaces for most ``SUNLinearSolver`` modules. See
   :numref:`SUNLinSol` for more details on how to use the interfaces.

**SUNNonlinearSolver module changes**

-  A new function was added to aid in creating custom ``SUNNonlinearSolver``
   objects. The constructor :c:func:`SUNNonlinSolSetConvTestFN` allocates an
   "empty" generic ``SUNNonlinearSolver`` with the object's content pointer and
   the function pointers in the operations structure initialized to . When used
   in the constructor for custom objects this function will ease the
   introduction of any new optional operations to the ``SUNNonlinearSolver`` API
   by ensuring only required operations need to be set. See
   :numref:`SUNNonlinSol.API.Custom` for more details.
-  To facilitate the use of user supplied nonlinear solver convergence
   test functions the function in the ``SUNNonlinearSolver`` API has been
   updated to take a data pointer as input. The supplied data pointer will be
   passed to the nonlinear solver convergence test function on each call.
-  The inputs values passed to the first two inputs of the function
   :c:func:`SUNNonlinSolSolve` in the ``SUNNonlinearSolver`` have been changed to
   be the predicted state and the initial guess for the correction to that state. Additionally, the
   definitions of :c:func:`SUNNonlinSolLSetupFn` and
   :c:func:`SUNNonlinSolLSolveFn` in the ``SUNNonlinearSolver`` API have been
   updated to remove unused input parameters. For more information on the
   nonlinear system formulation see :numref:`SUNNonlinSol.CVODE` and for more
   details on the API functions see :numref:`SUNNonlinSol`.
-  Added a new ``SUNNonlinearSolver`` implementation, ``SUNNONLINSOL_PETSC``,
   which interfaces to the PETSc SNES nonlinear solver API. See
   :numref:`SUNNonlinSol.PetscSNES` for more details.
-  Added new Fortran 2003 interfaces for most ``SUNNonlinearSolver`` modules.
   See :numref:`SUNDIALS.Fortran` for more details on how to use the
   interfaces.

**CVODE changes**

-  Fixed a bug in the CVODE constraint handling where the step size
   could be set below the minimum step size.
-  Fixed a bug in the CVODE nonlinear solver interface where the
   norm of the accumulated correction was not updated when using a non-default
   convergence test function.
-  Fixed a memeory leak in FCVODE when not using the default
   nonlinear solver.
-  Removed extraneous calls to for simulations where the scalar valued
   absolute tolerance, or all entries of the vector-valued absolute tolerance
   array, are strictly positive. In this scenario, CVODE will remove at least
   one global reduction per time step.
-  The CVLS interface has been updated to only zero the Jacobian matrix
   before calling a user-supplied Jacobian evaluation function when the attached
   linear solver has type ``SUNLINEARSOLVER_DIRECT``.
-  A new linear solver interface function :c:func:`CVLsLinSysFn` was added as an
   alternative method for evaluating the linear system :math:`M = I - \gamma J`.
-  Added two new functions, :c:func:`CVodeGetCurrentGamma` and :c:func:`CVodeGetCurrentState`, which may be useful to users who
   choose to provide their own nonlinear solver implementations.
-  The CVODE Fortran 2003 interface was completely redone to be more
   sustainable and to allow users to write more idiomatic Fortran. See
   :numref:`SUNDIALS.Fortran` for more details.

### CVODES Changes in v5.0.0

**Build system changes**

-  Increased the minimum required CMake version to 3.5 for most
   SUNDIALS configurations, and 3.10 when CUDA or OpenMP with device
   offloading are enabled.

-  The CMake option ``BLAS_ENABLE`` and the variable ``BLAS_LIBRARIES`` have been removed to simplify
   builds as SUNDIALS packages do not use BLAS directly. For third
   party libraries that require linking to BLAS, the path to the BLAS
   library should be included in the variable for the third party
   library *e.g.*, ``SUPERLUDIST_LIBRARIES`` when enabling SuperLU_DIST.

-  Fixed a bug in the build system that prevented the ``NVECTOR_PTHREADS``
   module from being built.

**NVECTOR module changes**

-  Two new functions were added to aid in creating custom ``N_Vector``
   objects. The constructor :c:func:`N_VNewEmpty` allocates an "empty" generic ``N_Vector``
   with the object's content pointer and the function pointers in the
   operations structure initialized to  ``NULL``. When used in the constructor
   for custom objects this function will ease the introduction of any
   new optional operations to the ``N_Vector`` API by ensuring only
   required operations need to be set. Additionally, the function :c:func:`N_VCopyOps` has
   been added to copy the operation function pointers between vector
   objects. When used in clone routines for custom vector objects these
   functions also will ease the introduction of any new optional
   operations to the ``N_Vector`` API by ensuring all operations are
   copied when cloning objects. See :numref:`NVectors.Description.custom_implementation` for more details.

-  Two new ``N_Vector`` implementations, ``NVECTOR_MANYVECTOR`` and
   ``NVECTOR_MPIMANYVECTOR``, have been created to support flexible
   partitioning of solution data among different processing elements
   (e.g., CPU + GPU) or for multi-physics problems that couple distinct
   MPI-based simulations together. This implementation is accompanied by
   additions to user documentation and SUNDIALS examples. See
   :numref:`NVectors.manyvector` and :numref:`NVectors.mpimanyvector` for more
   details.

-  One new required vector operation and ten new optional vector
   operations have been added to the ``N_Vector`` API. The new required
   operation, , returns the global length of an . The optional operations have
   been added to support the new ``NVECTOR_MPIMANYVECTOR`` implementation. The
   operation must be implemented by subvectors that are combined to create an
   ``NVECTOR_MPIMANYVECTOR``, but is not used outside of this context. The
   remaining nine operations are optional local reduction operations intended to
   eliminate unnecessary latency when performing vector reduction operations
   (norms, etc.) on distributed memory systems. The optional local reduction
   vector operations are :c:func:`N_VDotProdLocal`, :c:func:`N_VMaxNormLocal`,
   :c:func:`N_VL1NormLocal`, :c:func:`N_VWSqrSumLocal`,
   :c:func:`N_VWSqrSumMaskLocal`, :c:func:`N_VInvTestLocal`,
   :c:func:`N_VConstrMaskLocal`, :c:func:`N_VMinLocal`, and
   :c:func:`N_VMinQuotientLocal`. If an ``N_Vector`` implementation defines any
   of the local operations as , then the ``NVECTOR_MPIMANYVECTOR`` will call
   standard ``N_Vector`` operations to complete the computation.

-  An additional ``N_Vector`` implementation, ``NVECTOR_MPIPLUSX``, has been
   created to support the MPI+X paradigm where X is a type of on-node
   parallelism (*e.g.*, OpenMP, CUDA). The implementation is accompanied
   by additions to user documentation and SUNDIALS examples. See
   :numref:`NVectors.mpiplusx` for more details.

-  The and functions have been removed from the ``NVECTOR_CUDA`` and
   ``NVECTOR_RAJA`` implementations respectively. Accordingly, the
   ``nvector_mpicuda.h``, ``libsundials_nvecmpicuda.lib``,
   ``libsundials_nvecmpicudaraja.lib``, and files have been removed. Users
   should use the ``NVECTOR_MPIPLUSX`` module coupled in conjunction with the
   ``NVECTOR_CUDA`` or ``NVECTOR_RAJA`` modules to replace the functionality.
   The necessary changes are minimal and should require few code modifications.
   See the programs in and for examples of how to use the ``NVECTOR_MPIPLUSX``
   module with the ``NVECTOR_CUDA`` and ``NVECTOR_RAJA`` modules respectively.

-  Fixed a memory leak in the ``NVECTOR_PETSC`` module clone function.

-  Made performance improvements to the ``NVECTOR_CUDA`` module. Users who
   utilize a non-default stream should no longer see default stream
   synchronizations after memory transfers.

-  Added a new constructor to the ``NVECTOR_CUDA`` module that allows a user
   to provide custom allocate and free functions for the vector data
   array and internal reduction buffer. See :numref:`NVectors.Cuda` for more details.

-  Added new Fortran 2003 interfaces for most ``N_Vector`` modules. See
   :numref:`NVectors` for more details on how to use
   the interfaces.

-  Added three new ``N_Vector`` utility functions :c:func:`N_VGetVecAtIndexVectorArray`,
   :c:func:`N_VSetVecAtIndexVectorArray`, and :c:func:`N_VNewVectorArray` for working
   with arrays when using the Fortran 2003 interfaces.

**SUNMatrix module changes**

-  Two new functions were added to aid in creating custom ``SUNMatrix``
   objects. The constructor :c:func:`SUNMatNewEmpty` allocates an "empty" generic ``SUNMatrix``
   with the object's content pointer and the function pointers in the
   operations structure initialized to . When used in the constructor
   for custom objects this function will ease the introduction of any
   new optional operations to the ``SUNMatrix`` API by ensuring only
   required operations need to be set. Additionally, the function :c:func:`SUNMatCopyOps` has
   been added to copy the operation function pointers between matrix
   objects. When used in clone routines for custom matrix objects these
   functions also will ease the introduction of any new optional
   operations to the ``SUNMatrix`` API by ensuring all operations are
   copied when cloning objects. See :numref:`SUNMatrix` for more
   details.
-  A new operation, :c:func:`SUNMatMatvecSetup`, was added to the ``SUNMatrix`` API to perform any
   setup necessary for computing a matrix-vector product. This operation
   is useful for ``SUNMatrix`` implementations which need to prepare the
   matrix itself, or communication structures before performing the
   matrix-vector product. Users who have implemented custom
   ``SUNMatrix`` modules will need to at least update their code to set
   the corresponding structure member to ``NULL``. See :numref:`SUNMatrix.Ops`
   for more details.
-  The generic ``SUNMatrix`` API now defines error codes to be returned
   by ``SUNMatrix`` operations. Operations which return an integer flag
   indiciating success/failure may return different values than
   previously.
-  A new ``SUNMatrix`` (and ``SUNLinearSolver``) implementation was added to
   facilitate the use of the SuperLU_DIST library with SUNDIALS. See
   :numref:`SUNMatrix.SLUNRloc` for more details.
-  Added new Fortran 2003 interfaces for most ``SUNMatrix`` modules. See
   :numref:`SUNMatrix` for more details on how to
   use the interfaces.

**SUNLinearSolver module changes**

-  A new function was added to aid in creating custom ``SUNLinearSolver``
   objects. The constructor allocates an "empty" generic ``SUNLinearSolver``
   with the object's content pointer and the function pointers in the operations
   structure initialized to . When used in the constructor for custom objects
   this function will ease the introduction of any new optional operations to
   the ``SUNLinearSolver`` API by ensuring only required operations need to be
   set. See :numref:`SUNLinSol.API.Custom` for more details.
-  The return type of the ``SUNLinearSolver`` API function has changed from to
   to be consistent with the type used to store row indices in dense and banded
   linear solver modules.
-  Added a new optional operation to the ``SUNLinearSolver`` API,
   :c:func:`SUNLinSolLastFlag`, that returns a for identifying the linear solver module.
-  The ``SUNLinearSolver`` API has been updated to make the initialize and
   setup functions optional.
-  A new ``SUNLinearSolver`` (and ``SUNMatrix``) implementation was added to
   facilitate the use of the SuperLU_DIST library with SUNDIALS. See
   :numref:`SUNLinSol.SuperLUDIST` for more details.
-  Added a new ``SUNLinearSolver`` implementation, :ref:`SUNLINEARSOLVER_CUSOLVERSP <SUNLinSol.cuSolverSp>`,
   which leverages the NVIDIA cuSOLVER sparse batched QR method for efficiently solving block
   diagonal linear systems on NVIDIA GPUs.
-  Added three new accessor functions to the ``SUNLINSOL_KLU`` module, :c:func:`SUNLinSol_KLUGetSymbolic`,
   , :c:func:`SUNLinSol_KLUGetNumeric` and :c:func:`SUNLinSol_KLUGetCommon`, to
   provide user access to the underlying KLU solver structures. See
   :numref:`SUNLinSol.KLU` for more details.
-  Added new Fortran 2003 interfaces for most ``SUNLinearSolver`` modules. See
   :numref:`SUNLinSol` for more details on how to use the interfaces.

**SUNNonlinearSolver module changes**

-  A new function was added to aid in creating custom ``SUNNonlinearSolver``
   objects. The constructor :c:func:`SUNNonlinSolSetConvTestFN` allocates an
   "empty" generic ``SUNNonlinearSolver`` with the object's content pointer and
   the function pointers in the operations structure initialized to . When used
   in the constructor for custom objects this function will ease the
   introduction of any new optional operations to the ``SUNNonlinearSolver`` API
   by ensuring only required operations need to be set. See
   :numref:`SUNNonlinSol.API.Custom` for more details.
-  To facilitate the use of user supplied nonlinear solver convergence
   test functions the function in the ``SUNNonlinearSolver`` API has been
   updated to take a data pointer as input. The supplied data pointer will be
   passed to the nonlinear solver convergence test function on each call.
-  The inputs values passed to the first two inputs of the function
   :c:func:`SUNNonlinSolSolve` in the ``SUNNonlinearSolver`` have been changed to
   be the predicted state and the initial guess for the correction to that state. Additionally, the
   definitions of :c:func:`SUNNonlinSolLSetupFn` and
   :c:func:`SUNNonlinSolLSolveFn` in the ``SUNNonlinearSolver`` API have been
   updated to remove unused input parameters. For more information on the
   nonlinear system formulation see :numref:`SUNNonlinSol.CVODES` and for more
   details on the API functions see :numref:`SUNNonlinSol`.
-  Added a new ``SUNNonlinearSolver`` implementation, ``SUNNONLINSOL_PETSC``,
   which interfaces to the PETSc SNES nonlinear solver API. See
   :numref:`SUNNonlinSol.PetscSNES` for more details.
-  Added new Fortran 2003 interfaces for most ``SUNNonlinearSolver`` modules.
   See :numref:`SUNDIALS.Fortran` for more details on how to use the
   interfaces.


CVODES changes
^^^^^^^^^^^^^^

-  Fixed a bug in the CVODES constraint handling where the step size could be
   set below the minimum step size.

-  Fixed a bug in the CVODES nonlinear solver interface where the norm of the
   accumulated correction was not updated when using a non-default convergence
   test function.

-  Fixed a bug in the CVODES ``cvRescale`` function where the loops to compute
   the array of scalars for the fused vector scale operation stopped one
   iteration early.

-  Fixed a bug where the ``CVodeF`` function would return the wrong flag under
   certrain cirumstances.

-  Fixed a bug where the ``CVodeF`` function would not return a root in
   ``CV_NORMAL_STEP`` mode if the root occurred after the desired output time.

-  Removed extraneous calls to ``N_VMin`` for simulations where the scalar
   valued absolute tolerance, or all entries of the vector-valued absolute
   tolerance array, are strictly positive. In this scenario, CVODES will remove
   at least one global reduction per time step.

-  The CVLS interface has been updated to only zero the Jacobian matrix before
   calling a user-supplied Jacobian evaluation function when the attached linear
   solver has type ``SUNLINEARSOLVER_DIRECT``.

-  A new linear solver interface function ``CVLsLinSysFn`` was added as an
   alternative method for evaluating the linear system :math:`M = I - \gamma J`.

-  Added new functions, ``CVodeGetCurrentGamma``, ``CVodeGetCurrentState``,
   ``CVodeGetCurrentStateSens``, and ``CVodeGetCurrentSensSolveIndex`` which may
   be useful to users who choose to provide their own nonlinear solver
   implementations.

-  Added a Fortran 2003 interface to CVODES. See
   Chapter :numref:`SUNDIALS.Fortran` for more details.

### IDA Changes in v5.0.0

Build system changes
^^^^^^^^^^^^^^^^^^^^

* Increased the minimum required CMake version to 3.5 for most SUNDIALS
  configurations, and 3.10 when CUDA or OpenMP with device offloading are
  enabled.

* The CMake option ``BLAS_ENABLE`` and the variable ``BLAS_LIBRARIES`` have been
  removed to simplify builds as SUNDIALS packages do not use BLAS directly. For
  third party libraries that require linking to BLAS, the path to the BLAS
  library should be included in the ``*_LIBRARIES`` variable for the third party
  library *e.g.*, :cmakeop:`SUPERLUDIST_LIBRARIES` when enabling SuperLU_DIST.

* Fixed a bug in the build system that prevented the
  :ref:`NVECTOR_PTHREADS <NVectors.Pthreads>` module from being built.

NVECTOR module changes
^^^^^^^^^^^^^^^^^^^^^^

* Two new functions were added to aid in creating custom ``N_Vector``
  objects. The constructor :c:func:`N_VNewEmpty` allocates an "empty" generic
  ``N_Vector`` with the object's content pointer and the function pointers in
  the operations structure initialized to ``NULL``. When used in the constructor
  for custom objects this function will ease the introduction of any new
  optional operations to the ``N_Vector`` API by ensuring only required
  operations need to be set.  Additionally, the function :c:func:`N_VCopyOps`
  has been added to copy the operation function pointers between vector
  objects. When used in clone routines for custom vector objects these functions
  also will ease the introduction of any new optional operations to the
  ``N_Vector`` API by ensuring all operations are copied when cloning
  objects. See :numref:`NVectors.Description.utilities` for more details.

* Two new ``N_Vector`` implementations,
  :ref:`NVECTOR_MANYVECTOR <NVectors.ManyVector>` and
  :ref:`NVECTOR_MPIMANYVECTOR <NVectors.MPIManyVector>`, have been created to
  support flexible partitioning of solution data among different processing
  elements (e.g., CPU + GPU) or for multi-physics problems that couple distinct
  MPI-based simulations together. This implementation is accompanied by
  additions to user documentation and SUNDIALS examples. See
  :numref:`NVectors.ManyVector` and :numref:`NVectors.MPIManyVector` for more
  details.

* One new required vector operation and ten new optional vector operations have
  been added to the ``N_Vector`` API. The new required operation,
  :c:func:`N_VGetLength`, returns the global length of an ``N_Vector``.  The
  optional operations have been added to support the new
  :ref:`NVECTOR_MPIMANYVECTOR <NVectors.MPIManyVector>` implementation. The
  operation :c:func:`N_VGetCommunicator` must be implemented by subvectors that
  are combined to create an
  :ref:`NVECTOR_MPIMANYVECTOR <NVectors.MPIManyVector>`, but is not used outside
  of this context. The remaining nine operations are optional local reduction
  operations intended to eliminate unnecessary latency when performing vector
  reduction operations (norms, etc.) on distributed memory systems. The optional
  local reduction vector operations are :c:func:`N_VDotProdLocal`,
  :c:func:`N_VMaxNormLocal`, :c:func:`N_VMinLocal`, :c:func:`N_VL1NormLocal`,
  :c:func:`N_VWSqrSumLocal`, :c:func:`N_VWSqrSumMaskLocal`,
  :c:func:`N_VInvTestLocal`, :c:func:`N_VConstrMaskLocal`, and
  :c:func:`N_VMinQuotientLocal`. If an ``N_Vector`` implementation defines any
  of the local operations as ``NULL``, then the
  :ref:`NVECTOR_MPIMANYVECTOR <NVectors.MPIManyVector>` will call standard
  ``N_Vector`` operations to complete the computation. See
  :numref:`NVectors.Ops.Local` for more details.

* An additional ``N_Vector`` implementation, :ref:`NVECTOR_MPIPLUSX
  <NVectors.MPIPlusX>`, has been created to support the MPI+X paradigm where X
  is a type of on-node parallelism (*e.g.*, OpenMP, CUDA). The implementation is
  accompanied by additions to user documentation and SUNDIALS examples. See
  :numref:`NVectors.MPIPlusX` for more details.

* The ``*_MPICuda`` and ``*_MPIRaja`` functions have been removed from the
  :ref:`NVECTOR_CUDA <NVectors.CUDA>` and :ref:`NVECTOR_RAJA <NVectors.RAJA>`
  implementations respectively. Accordingly, the ``nvector_mpicuda.h``,
  ``nvector_mpiraja.h``, ``libsundials_nvecmpicuda.lib``, and
  ``libsundials_nvecmpicudaraja.lib`` files have been removed. Users should use
  the :ref:`NVECTOR_MPIPLUSX <NVectors.MPIPlusX>` module coupled in conjunction
  with the :ref:`NVECTOR_CUDA <NVectors.CUDA>` or :ref:`NVECTOR_RAJA
  <NVectors.RAJA>` modules to replace the functionality. The necessary changes
  are minimal and should require few code modifications. See the programs in
  ``examples/ida/mpicuda`` and ``examples/ida/mpiraja`` for examples of how to
  use the :ref:`NVECTOR_MPIPLUSX <NVectors.MPIPlusX>` module with the
  :ref:`NVECTOR_CUDA <NVectors.CUDA>` and :ref:`NVECTOR_RAJA <NVectors.RAJA>`
  modules respectively.

* Fixed a memory leak in the :ref:`NVECTOR_PETSC <NVectors.NVPETSc>` module
  clone function.

* Made performance improvements to the :ref:`NVECTOR_CUDA <NVectors.CUDA>`
  module. Users who utilize a non-default stream should no longer see default
  stream synchronizations after memory transfers.

* Added a new constructor to the :ref:`NVECTOR_CUDA <NVectors.CUDA>` module that
  allows a user to provide custom allocate and free functions for the vector
  data array and internal reduction buffer. See :numref:`NVectors.CUDA` for more
  details.

* Added new Fortran 2003 interfaces for most ``N_Vector`` modules. See
  :numref:`NVectors` for more details on how to use the interfaces.

* Added three new ``N_Vector`` utility functions,
  :c:func:`FN_VGetVecAtIndexVectorArray`,
  :c:func:`FN_VSetVecAtIndexVectorArray`, and :c:func:`FN_VNewVectorArray`, for
  working with ``N_Vector`` arrays when using the Fortran 2003 interfaces.  See
  :numref:`NVectors.Description.utilities` for more details.

SUNMatrix module changes
^^^^^^^^^^^^^^^^^^^^^^^^

* Two new functions were added to aid in creating custom ``SUNMatrix``
  objects. The constructor :c:func:`SUNMatNewEmpty` allocates an "empty" generic
  ``SUNMatrix`` with the object's content pointer and the function pointers in
  the operations structure initialized to ``NULL``. When used in the constructor
  for custom objects this function will ease the introduction of any new
  optional operations to the ``SUNMatrix`` API by ensuring only required
  operations need to be set.  Additionally, the function :c:func:`SUNMatCopyOps`
  has been added to copy the operation function pointers between matrix
  objects. When used in clone routines for custom matrix objects these functions
  also will ease the introduction of any new optional operations to the
  ``SUNMatrix`` API by ensuring all operations are copied when cloning
  objects. See :numref:`SUNMatrix.Description` for more details.

* A new operation, :c:func:`SUNMatMatvecSetup`, was added to the ``SUNMatrix``
  API to perform any setup necessary for computing a matrix-vector product. This
  operation is useful for ``SUNMatrix`` implementations which need to prepare
  the matrix itself, or communication structures before performing the
  matrix-vector product. Users who have implemented custom ``SUNMatrix`` modules
  will need to at least update their code to set the corresponding ``ops``
  structure member, ``matvecsetup``, to ``NULL``. See
  :numref:`SUNMatrix.Description` for more details.

* The generic ``SUNMatrix`` API now defines error codes to be returned by
  ``SUNMatrix`` operations. Operations which return an integer flag indicating
  success/failure may return different values than previously.

* A new ``SUNMatrix`` (and ``SUNLinearSolver``) implementation was added to
  facilitate the use of the SuperLU_DIST library with SUNDIALS. See
  :numref:`SUNMatrix.SLUNRloc` for more details.

* Added new Fortran 2003 interfaces for most ``SUNMatrix`` modules. See
  :numref:`SUNMatrix` for more details on how to use the interfaces.

SUNLinearSolver module changes
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* A new function was added to aid in creating custom ``SUNLinearSolver``
  objects.  The constructor :c:func:`SUNLinSolNewEmpty` allocates an "empty"
  generic ``SUNLinearSolver`` with the object's content pointer and the function
  pointers in the operations structure initialized to ``NULL``. When used in the
  constructor for custom objects this function will ease the introduction of any
  new optional operations to the ``SUNLinearSolver`` API by ensuring only
  required operations need to be set. See :numref:`SUNLinSol.API.Custom` for
  more details.

* The return type of the ``SUNLinearSolver`` API function
  :c:func:`SUNLinSolLastFlag` has changed from ``long int`` to ``sunindextype``
  to be consistent with the type used to store row indices in dense and banded
  linear solver modules.

* Added a new optional operation to the ``SUNLinearSolver`` API,
  :c:func:`SUNLinSolGetID`, that returns a ``SUNLinearSolver_ID`` for
  identifying the linear solver module.

* The ``SUNLinearSolver`` API has been updated to make the initialize and setup
  functions optional.

* A new ``SUNLinearSolver`` (and ``SUNMatrix``) implementation was added to
  facilitate the use of the SuperLU_DIST library with SUNDIALS. See
  :numref:`SUNLinSol.SuperLUDIST` for more details.

* Added a new ``SUNLinearSolver`` implementation,
  SUNLinearSolver_cuSolverSp_batchQR, which leverages the NVIDIA cuSOLVER sparse
  batched QR method for efficiently solving block diagonal linear systems on
  NVIDIA GPUs. See :numref:`SUNLinSol.cuSolverSp` for more details.

* Added three new accessor functions to the SUNLINSOL_KLU module,
  :c:func:`SUNLinSol_KLUGetSymbolic`, :c:func:`SUNLinSol_KLUGetNumeric`, and
  :c:func:`SUNLinSol_KLUGetCommon`, to provide user access to the underlying KLU
  solver structures. See :numref:`SUNLinSol.KLU` for more details.

* Added new Fortran 2003 interfaces for most ``SUNLinearSolver`` modules.  See
  :numref:`SUNLinSol` for more details on how to use the interfaces.

SUNNonlinearSolver module changes
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* A new function was added to aid in creating custom ``SUNNonlinearSolver``
  objects. The constructor :c:func:`SUNNonlinSolNewEmpty` allocates an "empty"
  generic ``SUNNonlinearSolver`` with the object's content pointer and the
  function pointers in the operations structure initialized to ``NULL``. When
  used in the constructor for custom objects this function will ease the
  introduction of any new optional operations to the ``SUNNonlinearSolver`` API
  by ensuring only required operations need to be set. See
  :numref:`SUNNonlinSol.API.Custom` for more details.

* To facilitate the use of user supplied nonlinear solver convergence test
  functions the :c:type:`SUNNonlinSolSetConvTestFn` function in the
  ``SUNNonlinearSolver`` API has been updated to take a ``void*`` data pointer
  as input. The supplied data pointer will be passed to the nonlinear solver
  convergence test function on each call.

* The inputs values passed to the first two inputs of the
  :c:func:`SUNNonlinSolSolve` function in the ``SUNNonlinearSolver`` have been
  changed to be the predicted state and the initial guess for the correction to
  that state. Additionally, the definitions of :c:type:`SUNNonlinSolLSetupFn`
  and :c:type:`SUNNonlinSolLSolveFn` in the ``SUNNonlinearSolver`` API have been
  updated to remove unused input parameters. For more information see
  :numref:`SUNNonlinSol`.

* Added a new ``SUNNonlinearSolver`` implementation,
  :ref:`SUNNONLINSOL_PETSC <SUNNonlinSol.PetscSNES>`, which interfaces to the
  PETSc SNES nonlinear solver API. See :numref:`SUNNonlinSol.PetscSNES` for more
  details.

* Added new Fortran 2003 interfaces for most ``SUNNonlinearSolver`` modules. See
  :numref:`SUNNonlinSol` for more details on how to use the interfaces.

IDA changes
^^^^^^^^^^^

* A bug was fixed in the IDA linear solver interface where an incorrect
  Jacobian-vector product increment was used with iterative solvers other than
  :ref:`SUNLINSOL_SPGMR <SUNLinSol.SPGMR>` and
  :ref:`SUNLINSOL_SPFGMR <SUNLinSol.SPFGMR>`.

* Fixed a memeory leak in FIDA when not using the default nonlinear solver.

* Removed extraneous calls to :c:func:`N_VMin` for simulations where the scalar
  valued absolute tolerance, or all entries of the vector-valued absolute
  tolerance array, are strictly positive. In this scenario, IDA will remove at
  least one global reduction per time step.

* The IDALS interface has been updated to only zero the Jacobian matrix before
  calling a user-supplied Jacobian evaluation function when the attached linear
  solver has type ``SUNLINEARSOLVER_DIRECT``.

* Added the new functions, :c:func:`IDAGetCurrentCj`, :c:func:`IDAGetCurrentY`,
  :c:func:`IDAGetCurrentYp`, :c:func:`IDAComputeY`, and :c:func:`IDAComputeYp`
  which may be useful to users who choose to provide their own nonlinear solver
  implementations.

* Added a Fortran 2003 interface to IDA. See :numref:`SUNDIALS.Fortran` for more
  details.

### IDAS Changes in v4.0.0

Build system changes
^^^^^^^^^^^^^^^^^^^^

* Increased the minimum required CMake version to 3.5 for most SUNDIALS
  configurations, and 3.10 when CUDA or OpenMP with device offloading are
  enabled.

* The CMake option ``BLAS_ENABLE`` and the variable ``BLAS_LIBRARIES`` have been
  removed to simplify builds as SUNDIALS packages do not use BLAS directly. For
  third party libraries that require linking to BLAS, the path to the BLAS
  library should be included in the ``*_LIBRARIES`` variable for the third party
  library *e.g.*, :cmakeop:`SUPERLUDIST_LIBRARIES` when enabling SuperLU_DIST.

* Fixed a bug in the build system that prevented the
  :ref:`NVECTOR_PTHREADS <NVectors.Pthreads>` module from being built.

NVECTOR module changes
^^^^^^^^^^^^^^^^^^^^^^

* Two new functions were added to aid in creating custom ``N_Vector``
  objects. The constructor :c:func:`N_VNewEmpty` allocates an "empty" generic
  ``N_Vector`` with the object's content pointer and the function pointers in
  the operations structure initialized to ``NULL``. When used in the constructor
  for custom objects this function will ease the introduction of any new
  optional operations to the ``N_Vector`` API by ensuring only required
  operations need to be set.  Additionally, the function :c:func:`N_VCopyOps`
  has been added to copy the operation function pointers between vector
  objects. When used in clone routines for custom vector objects these functions
  also will ease the introduction of any new optional operations to the
  ``N_Vector`` API by ensuring all operations are copied when cloning
  objects. See :numref:`NVectors.Description.utilities` for more details.

* Two new ``N_Vector`` implementations,
  :ref:`NVECTOR_MANYVECTOR <NVectors.ManyVector>` and
  :ref:`NVECTOR_MPIMANYVECTOR <NVectors.MPIManyVector>`, have been created to
  support flexible partitioning of solution data among different processing
  elements (e.g., CPU + GPU) or for multi-physics problems that couple distinct
  MPI-based simulations together. This implementation is accompanied by
  additions to user documentation and SUNDIALS examples. See
  :numref:`NVectors.ManyVector` and :numref:`NVectors.MPIManyVector` for more
  details.

* One new required vector operation and ten new optional vector operations have
  been added to the ``N_Vector`` API. The new required operation,
  :c:func:`N_VGetLength`, returns the global length of an ``N_Vector``.  The
  optional operations have been added to support the new
  :ref:`NVECTOR_MPIMANYVECTOR <NVectors.MPIManyVector>` implementation. The
  operation :c:func:`N_VGetCommunicator` must be implemented by subvectors that
  are combined to create an
  :ref:`NVECTOR_MPIMANYVECTOR <NVectors.MPIManyVector>`, but is not used outside
  of this context. The remaining nine operations are optional local reduction
  operations intended to eliminate unnecessary latency when performing vector
  reduction operations (norms, etc.) on distributed memory systems. The optional
  local reduction vector operations are :c:func:`N_VDotProdLocal`,
  :c:func:`N_VMaxNormLocal`, :c:func:`N_VMinLocal`, :c:func:`N_VL1NormLocal`,
  :c:func:`N_VWSqrSumLocal`, :c:func:`N_VWSqrSumMaskLocal`,
  :c:func:`N_VInvTestLocal`, :c:func:`N_VConstrMaskLocal`, and
  :c:func:`N_VMinQuotientLocal`. If an ``N_Vector`` implementation defines any
  of the local operations as ``NULL``, then the
  :ref:`NVECTOR_MPIMANYVECTOR <NVectors.MPIManyVector>` will call standard
  ``N_Vector`` operations to complete the computation. See
  :numref:`NVectors.Ops.Local` for more details.

* An additional ``N_Vector`` implementation, :ref:`NVECTOR_MPIPLUSX
  <NVectors.MPIPlusX>`, has been created to support the MPI+X paradigm where X
  is a type of on-node parallelism (*e.g.*, OpenMP, CUDA). The implementation is
  accompanied by additions to user documentation and SUNDIALS examples. See
  :numref:`NVectors.MPIPlusX` for more details.

* The ``*_MPICuda`` and ``*_MPIRaja`` functions have been removed from the
  :ref:`NVECTOR_CUDA <NVectors.CUDA>` and :ref:`NVECTOR_RAJA <NVectors.RAJA>`
  implementations respectively. Accordingly, the ``nvector_mpicuda.h``,
  ``nvector_mpiraja.h``, ``libsundials_nvecmpicuda.lib``, and
  ``libsundials_nvecmpicudaraja.lib`` files have been removed. Users should use
  the :ref:`NVECTOR_MPIPLUSX <NVectors.MPIPlusX>` module coupled in conjunction
  with the :ref:`NVECTOR_CUDA <NVectors.CUDA>` or :ref:`NVECTOR_RAJA
  <NVectors.RAJA>` modules to replace the functionality. The necessary changes
  are minimal and should require few code modifications. See the programs in
  ``examples/ida/mpicuda`` and ``examples/ida/mpiraja`` for examples of how to
  use the :ref:`NVECTOR_MPIPLUSX <NVectors.MPIPlusX>` module with the
  :ref:`NVECTOR_CUDA <NVectors.CUDA>` and :ref:`NVECTOR_RAJA <NVectors.RAJA>`
  modules respectively.

* Fixed a memory leak in the :ref:`NVECTOR_PETSC <NVectors.NVPETSc>` module
  clone function.

* Made performance improvements to the :ref:`NVECTOR_CUDA <NVectors.CUDA>`
  module. Users who utilize a non-default stream should no longer see default
  stream synchronizations after memory transfers.

* Added a new constructor to the :ref:`NVECTOR_CUDA <NVectors.CUDA>` module that
  allows a user to provide custom allocate and free functions for the vector
  data array and internal reduction buffer. See :numref:`NVectors.CUDA` for more
  details.

* Added new Fortran 2003 interfaces for most ``N_Vector`` modules. See
  :numref:`NVectors` for more details on how to use the interfaces.

* Added three new ``N_Vector`` utility functions,
  :c:func:`FN_VGetVecAtIndexVectorArray`,
  :c:func:`FN_VSetVecAtIndexVectorArray`, and :c:func:`FN_VNewVectorArray`, for
  working with ``N_Vector`` arrays when using the Fortran 2003 interfaces.  See
  :numref:`NVectors.Description.utilities` for more details.

SUNMatrix module changes
^^^^^^^^^^^^^^^^^^^^^^^^

* Two new functions were added to aid in creating custom ``SUNMatrix``
  objects. The constructor :c:func:`SUNMatNewEmpty` allocates an "empty" generic
  ``SUNMatrix`` with the object's content pointer and the function pointers in
  the operations structure initialized to ``NULL``. When used in the constructor
  for custom objects this function will ease the introduction of any new
  optional operations to the ``SUNMatrix`` API by ensuring only required
  operations need to be set.  Additionally, the function :c:func:`SUNMatCopyOps`
  has been added to copy the operation function pointers between matrix
  objects. When used in clone routines for custom matrix objects these functions
  also will ease the introduction of any new optional operations to the
  ``SUNMatrix`` API by ensuring all operations are copied when cloning
  objects. See :numref:`SUNMatrix.Description` for more details.

* A new operation, :c:func:`SUNMatMatvecSetup`, was added to the ``SUNMatrix``
  API to perform any setup necessary for computing a matrix-vector product. This
  operation is useful for ``SUNMatrix`` implementations which need to prepare
  the matrix itself, or communication structures before performing the
  matrix-vector product. Users who have implemented custom ``SUNMatrix`` modules
  will need to at least update their code to set the corresponding ``ops``
  structure member, ``matvecsetup``, to ``NULL``. See
  :numref:`SUNMatrix.Description` for more details.

* The generic ``SUNMatrix`` API now defines error codes to be returned by
  ``SUNMatrix`` operations. Operations which return an integer flag indicating
  success/failure may return different values than previously.

* A new ``SUNMatrix`` (and ``SUNLinearSolver``) implementation was added to
  facilitate the use of the SuperLU_DIST library with SUNDIALS. See
  :numref:`SUNMatrix.SLUNRloc` for more details.

* Added new Fortran 2003 interfaces for most ``SUNMatrix`` modules. See
  :numref:`SUNMatrix` for more details on how to use the interfaces.

SUNLinearSolver module changes
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* A new function was added to aid in creating custom ``SUNLinearSolver``
  objects.  The constructor :c:func:`SUNLinSolNewEmpty` allocates an "empty"
  generic ``SUNLinearSolver`` with the object's content pointer and the function
  pointers in the operations structure initialized to ``NULL``. When used in the
  constructor for custom objects this function will ease the introduction of any
  new optional operations to the ``SUNLinearSolver`` API by ensuring only
  required operations need to be set. See :numref:`SUNLinSol.API.Custom` for
  more details.

* The return type of the ``SUNLinearSolver`` API function
  :c:func:`SUNLinSolLastFlag` has changed from ``long int`` to ``sunindextype``
  to be consistent with the type used to store row indices in dense and banded
  linear solver modules.

* Added a new optional operation to the ``SUNLinearSolver`` API,
  :c:func:`SUNLinSolGetID`, that returns a ``SUNLinearSolver_ID`` for
  identifying the linear solver module.

* The ``SUNLinearSolver`` API has been updated to make the initialize and setup
  functions optional.

* A new ``SUNLinearSolver`` (and ``SUNMatrix``) implementation was added to
  facilitate the use of the SuperLU_DIST library with SUNDIALS. See
  :numref:`SUNLinSol.SuperLUDIST` for more details.

* Added a new ``SUNLinearSolver`` implementation,
  SUNLinearSolver_cuSolverSp_batchQR, which leverages the NVIDIA cuSOLVER sparse
  batched QR method for efficiently solving block diagonal linear systems on
  NVIDIA GPUs. See :numref:`SUNLinSol.cuSolverSp` for more details.

* Added three new accessor functions to the SUNLINSOL_KLU module,
  :c:func:`SUNLinSol_KLUGetSymbolic`, :c:func:`SUNLinSol_KLUGetNumeric`, and
  :c:func:`SUNLinSol_KLUGetCommon`, to provide user access to the underlying KLU
  solver structures. See :numref:`SUNLinSol.KLU` for more details.

* Added new Fortran 2003 interfaces for most ``SUNLinearSolver`` modules.  See
  :numref:`SUNLinSol` for more details on how to use the interfaces.

SUNNonlinearSolver module changes
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* A new function was added to aid in creating custom ``SUNNonlinearSolver``
  objects. The constructor :c:func:`SUNNonlinSolNewEmpty` allocates an "empty"
  generic ``SUNNonlinearSolver`` with the object's content pointer and the
  function pointers in the operations structure initialized to ``NULL``. When
  used in the constructor for custom objects this function will ease the
  introduction of any new optional operations to the ``SUNNonlinearSolver`` API
  by ensuring only required operations need to be set. See
  :numref:`SUNNonlinSol.API.Custom` for more details.

* To facilitate the use of user supplied nonlinear solver convergence test
  functions the :c:type:`SUNNonlinSolSetConvTestFn` function in the
  ``SUNNonlinearSolver`` API has been updated to take a ``void*`` data pointer
  as input. The supplied data pointer will be passed to the nonlinear solver
  convergence test function on each call.

* The inputs values passed to the first two inputs of the
  :c:func:`SUNNonlinSolSolve` function in the ``SUNNonlinearSolver`` have been
  changed to be the predicted state and the initial guess for the correction to
  that state. Additionally, the definitions of :c:type:`SUNNonlinSolLSetupFn`
  and :c:type:`SUNNonlinSolLSolveFn` in the ``SUNNonlinearSolver`` API have been
  updated to remove unused input parameters. For more information see
  :numref:`SUNNonlinSol`.

* Added a new ``SUNNonlinearSolver`` implementation,
  :ref:`SUNNONLINSOL_PETSC <SUNNonlinSol.PetscSNES>`, which interfaces to the
  PETSc SNES nonlinear solver API. See :numref:`SUNNonlinSol.PetscSNES` for more
  details.

* Added new Fortran 2003 interfaces for most ``SUNNonlinearSolver`` modules. See
  :numref:`SUNNonlinSol` for more details on how to use the interfaces.

IDAS changes
^^^^^^^^^^^^

* A bug was fixed in the IDAS linear solver interface where an incorrect
  Jacobian-vector product increment was used with iterative solvers other than
  :ref:`SUNLINSOL_SPGMR <SUNLinSol.SPGMR>` and
  :ref:`SUNLINSOL_SPFGMR <SUNLinSol.SPFGMR>`.

* Fixed a memeory leak in FIDA when not using the default nonlinear solver.

* Fixed a bug where the :c:func:`IDASolveF` function would not return a root in
  ``IDA_NORMAL_STEP`` mode if the root occurred
  after the desired output time.

* Fixed a bug where the :c:func:`IDASolveF` function would return the wrong flag
  under certrain cirumstances.

* Fixed a bug in :c:func:`IDAQuadReInitB` where an incorrect memory structure was
  passed to :c:func:`IDAQuadReInit`.

* Removed extraneous calls to :c:func:`N_VMin` for simulations where the scalar
  valued absolute tolerance, or all entries of the vector-valued absolute
  tolerance array, are strictly positive. In this scenario, IDAS will remove at
  least one global reduction per time step.

* The IDALS interface has been updated to only zero the Jacobian matrix before
  calling a user-supplied Jacobian evaluation function when the attached linear
  solver has type ``SUNLINEARSOLVER_DIRECT``.

* Added the new functions, :c:func:`IDAGetCurrentCj`, :c:func:`IDAGetCurrentY`,
  :c:func:`IDAGetCurrentYp`, :c:func:`IDAComputeY`, and :c:func:`IDAComputeYp`
  which may be useful to users who choose to provide their own nonlinear solver
  implementations.

* Added a Fortran 2003 interface to IDAS. See :numref:`SUNDIALS.Fortran` for more
  details.

### KINSOL Changes in v5.0.0

Build system changes
^^^^^^^^^^^^^^^^^^^^

-  Increased the minimum required CMake version to 3.5 for most SUNDIALS configurations, and 3.10 when CUDA or
   OpenMP with device offloading are enabled.

-  The CMake option ``BLAS_ENABLE`` and the variable ``BLAS_LIBRARIES`` have been removed to simplify builds as
   SUNDIALS packages do not use BLAS directly. For third party libraries that require linking to BLAS, the path to
   the BLAS library should be included in the ``_LIBRARIES`` variable for the third party library *e.g.*,
   ``SUPERLUDIST_LIBRARIES`` when enabling SuperLU_DIST.

-  Fixed a bug in the build system that prevented the ``NVECTOR_PTHREADS`` module from being built.

NVECTOR module changes
^^^^^^^^^^^^^^^^^^^^^^

-  Two new functions were added to aid in creating custom ``N_Vector`` objects. The constructor ``N_VNewEmpty``
   allocates an "empty" generic ``N_Vector`` with the object's content pointer and the function pointers in the
   operations structure initialized to ``NULL``. When used in the constructor for custom objects this function will ease
   the introduction of any new optional operations to the ``N_Vector`` API by ensuring only required operations need to
   be set. Additionally, the function ``N_VCopyOps(w, v)`` has been added to copy the operation function pointers
   between vector objects. When used in clone routines for custom vector objects these functions also will ease the
   introduction of any new optional operations to the ``N_Vector`` API by ensuring all operations are copied when
   cloning objects. See :numref:`NVectors.Description.utilities` for more details.

-  Two new ``N_Vector`` implementations, ``NVECTOR_MANYVECTOR`` and ``NVECTOR_MPIMANYVECTOR``, have been created to support
   flexible partitioning of solution data among different processing elements (e.g., CPU + GPU) or for multi-physics
   problems that couple distinct MPI-based simulations together. This implementation is accompanied by additions to user
   documentation and SUNDIALS examples. See :numref:`NVectors.ManyVector` and :numref:`NVectors.MPIManyVector` for more details.

-  One new required vector operation and ten new optional vector operations have been added to the ``N_Vector`` API.
   The new required operation, ``N_VGetLength``, returns the global length of an ``N_Vector``. The optional operations
   have been added to support the new ``NVECTOR_MPIMANYVECTOR`` implementation. The operation ``N_VGetCommunicator`` must
   be implemented by subvectors that are combined to create an ``NVECTOR_MPIMANYVECTOR``, but is not used outside of this
   context. The remaining nine operations are optional local reduction operations intended to eliminate unnecessary
   latency when performing vector reduction operations (norms, etc.) on distributed memory systems. The optional local
   reduction vector operations are ``N_VDotProdLocal``, ``N_VMaxNormLocal``, ``N_VMinLocal``, ``N_VL1NormLocal``,
   ``N_VWSqrSumLocal``, ``N_VWSqrSumMaskLocal``, ``N_VInvTestLocal``, ``N_VConstrMaskLocal``, and
   ``N_VMinQuotientLocal``. If an ``N_Vector`` implementation defines any of the local operations as ``NULL``, then the
   ``NVECTOR_MPIMANYVECTOR`` will call standard ``N_Vector`` operations to complete the computation. See
   :numref:`NVectors.Ops.Local` for more details.

-  An additional ``N_Vector`` implementation, ``NVECTOR_MPIPLUSX``, has been created to support the MPI+X paradigm where
   X is a type of on-node parallelism (*e.g.*, OpenMP, CUDA). The implementation is accompanied by additions to user
   documentation and SUNDIALS examples. See :numref:`NVectors.MPIPlusX` for more details.

-  The ``*_MPICuda`` and ``*_MPIRaja`` functions have been removed from the ``NVECTOR_CUDA`` and ``NVECTOR_RAJA``
   implementations respectively. Accordingly, the ``nvector_mpicuda.h``, ``nvector_mpiraja.h``,
   ``libsundials_nvecmpicuda.lib``, and ``libsundials_nvecmpicudaraja.lib`` files have been removed. Users should use
   the ``NVECTOR_MPIPLUSX`` module coupled in conjunction with the ``NVECTOR_CUDA`` or ``NVECTOR_RAJA`` modules to replace the
   functionality. The necessary changes are minimal and should require few code modifications. See the programs in
   ``examples/ida/mpicuda`` and ``examples/ida/mpiraja`` for examples of how to use the ``NVECTOR_MPIPLUSX`` module with
   the ``NVECTOR_CUDA`` and ``NVECTOR_RAJA`` modules respectively.

-  Fixed a memory leak in the ``NVECTOR_PETSC`` module clone function.

-  Made performance improvements to the ``NVECTOR_CUDA`` module. Users who utilize a non-default stream should no longer
   see default stream synchronizations after memory transfers.

-  Added a new constructor to the ``NVECTOR_CUDA`` module that allows a user to provide custom allocate and free functions
   for the vector data array and internal reduction buffer. See
   :numref:`NVectors.CUDA.Functions` for more details.

-  Added new Fortran 2003 interfaces for most ``N_Vector`` modules. See Chapter :numref:`NVectors` for more
   details on how to use the interfaces.

-  Added three new ``N_Vector`` utility functions, ``FN_VGetVecAtIndexVectorArray``, ``FN_VSetVecAtIndexVectorArray``,
   and ``FN_VNewVectorArray``, for working with ``N_Vector`` arrays when using the Fortran
   2003 interfaces. See :numref:`NVectors.Description.utilities` for more details.

SUNMatrix module changes
^^^^^^^^^^^^^^^^^^^^^^^^

-  Two new functions were added to aid in creating custom ``SUNMatrix`` objects. The constructor ``SUNMatNewEmpty``
   allocates an "empty" generic ``SUNMatrix`` with the object's content pointer and the function pointers in the
   operations structure initialized to ``NULL``. When used in the constructor for custom objects this function will ease
   the introduction of any new optional operations to the ``SUNMatrix`` API by ensuring only required operations need
   to be set. Additionally, the function ``SUNMatCopyOps(A, B)`` has been added to copy the operation function pointers
   between matrix objects. When used in clone routines for custom matrix objects these functions also will ease the
   introduction of any new optional operations to the ``SUNMatrix`` API by ensuring all operations are copied when
   cloning objects. See :numref:`SUNMatrix.Description` for more details.

-  A new operation, ``SUNMatMatvecSetup``, was added to the ``SUNMatrix`` API to perform any setup necessary for
   computing a matrix-vector product. This operation is useful for ``SUNMatrix`` implementations which need to prepare
   the matrix itself, or communication structures before performing the matrix-vector product. Users who have
   implemented custom ``SUNMatrix`` modules will need to at least update their code to set the corresponding ``ops``
   structure member, ``matvecsetup``, to ``NULL``. See :numref:`SUNMatrix.Ops` for
   more details.

-  The generic ``SUNMatrix`` API now defines error codes to be returned by ``SUNMatrix`` operations. Operations
   which return an integer flag indiciating success/failure may return different values than previously. See
   "SUNMatrix Error Codes" for more details.

-  A new ``SUNMatrix`` (and ``SUNLinearSolver``) implementation was added to facilitate the use of the SuperLU_DIST
   library with SUNDIALS. See :numref:`SUNMatrix.SLUNRloc` for more details.

-  Added new Fortran 2003 interfaces for most ``SUNMatrix`` modules. See Chapter :numref:`SUNMatrix` for
   more details on how to use the interfaces.

SUNLinearSolver module changes
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

-  A new function was added to aid in creating custom ``SUNLinearSolver`` objects. The constructor ``SUNLinSolNewEmpty``
   allocates an "empty" generic ``SUNLinearSolver`` with the object's content pointer and the function pointers in the
   operations structure initialized to ``NULL``. When used in the constructor for custom objects this function will ease
   the introduction of any new optional operations to the ``SUNLinearSolver`` API by ensuring only required operations need
   to be set. See :numref:`SUNLinSol.API.Custom` for more details.

-  The return type of the ``SUNLinearSolver`` API function ``SUNLinSolLastFlag`` has changed from ``long int`` to
   ``sunindextype`` to be consistent with the type used to store row indices in dense and banded linear solver modules.

-  Added a new optional operation to the ``SUNLinearSolver`` API, ``SUNLinSolGetID``, that returns a ``SUNLinearSolver_ID``
   for identifying the linear solver module.

-  The ``SUNLinearSolver`` API has been updated to make the initialize and setup functions optional.

-  A new ``SUNLinearSolver`` (and ``SUNMatrix``) implementation was added to facilitate the use of the SuperLU_DIST
   library with SUNDIALS. See :numref:`SUNLinSol.SuperLUDIST` for more details.

-  Added a new ``SUNLinearSolver`` implementation, ``SUNLinearSolver_cuSolverSp_batchQR``, which leverages the NVIDIA
   cuSOLVER sparse batched QR method for efficiently solving block diagonal linear systems on NVIDIA GPUs. See :numref:`SUNLinSol.cuSolverSp` for more details.

-  Added three new accessor functions to the ``SUNLINSOL_KLU`` module, ``SUNLinSol_KLUGetSymbolic``,
   ``SUNLinSol_KLUGetNumeric``, and ``SUNLinSol_KLUGetCommon``, to provide user access to the underlying KLU solver
   structures. See :numref:`SUNLinSol.KLU.Usage` for more details.

-  Added new Fortran 2003 interfaces for most ``SUNLinearSolver`` modules. See Chapter :numref:`SUNLinSol` for
   more details on how to use the interfaces.

KINSOL changes
^^^^^^^^^^^^^^

-  Fixed a bug in the KINSOL linear solver interface where the auxiliary scalar ``sJpnorm`` was not computed when
   necessary with the Picard iteration and the auxiliary scalar ``sFdotJp`` was unnecessarily computed in some cases.

-  The KINLS interface has been updated to only zero the Jacobian matrix before calling a user-supplied Jacobian
   evaluation function when the attached linear solver has type ``SUNLINEARSOLVER_DIRECT``.

-  Added a Fortran 2003 interface to KINSOL. See :numref:`SUNDIALS.Fortran` for more details.


## Changes to SUNDIALS in release 4.1.0

An additional N_Vector implementation was added for Tpetra vector from
Trilinos library to facilitate interoperability between SUNDIALS and Trilinos.
This implementation is accompanied by additions to user documentation and
SUNDIALS examples.

A bug was fixed where a nonlinear solver object could be freed twice in some use
cases.

The EXAMPLES_ENABLE_RAJA CMake option has been removed. The option
`EXAMPLES_ENABLE_CUDA` enables all examples that use CUDA including the RAJA
examples with a CUDA back end (if the RAJA NVECTOR is enabled).

The implementation header files (e.g. `arkode_impl.h`) are no longer installed.
This means users who are directly manipulating package memory structures will
need to update their code to use the package's public API.

Python is no longer required to run `make test` and `make test_install`.

Fixed a bug in `ARKodeButcherTable_Write` when printing a Butcher table
without an embedding.

### ARKODE Changes in v3.1.0

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
to use ARKODE's public API.

Python is no longer required to run ``make test`` and ``make test_install``.

Fixed a bug in ``ARKodeButcherTable_Write`` when printing a Butcher table
without an embedding.


### CVODE Changes in v4.1.0

An additional ``N_Vector`` implementation was added for the Tpetra
vector from the Trilinos library to facilitate interoperability
between SUNDIALS and Trilinos. This implementation is accompanied
by additions to user documentation and SUNDIALS examples.

A bug was fixed where a nonlinear solver object could be freed twice in
some use cases.

The CMake option ``EXAMPLES_ENABLE_RAJA`` has been removed. The option enables all examples that
use CUDA including the RAJA examples with a CUDA back end (if the RAJA
``N_Vector`` is enabled).

The implementation header file is no longer installed. This means users
who are directly manipulating the structure will need to update their
code to use CVODE's public API.

Python is no longer required to run ``make test`` and ``make test_install``.

### CVODES Changes in v4.1.0

An additional ``N_Vector`` implementation was added for the TPETRA vector from
the Trilinos library to facilitate interoperability between SUNDIALS and
Trilinos. This implementation is accompanied by additions to user documentation
and SUNDIALS examples.

A bug was fixed where a nonlinear solver object could be freed twice in some use
cases.

The ``EXAMPLES_ENABLE_RAJA`` CMake option has been removed. The option
``EXAMPLES_ENABLE_CUDA`` enables all examples that use CUDA including the RAJA
examples with a CUDA back end (if the RAJA ``N_Vector`` is enabled).

The implementation header file ``cvodes_impl.h`` is no longer installed. This
means users who are directly manipulating the ``CVodeMem`` structure will need
to update their code to use CVODES's public API.

Python is no longer required to run ``make test`` and ``make test_install``.

### IDA Changes in v4.1.0

An additional ``N_Vector`` implementation was added for the TPETRA vector from
the TRILINOS library to facilitate interoperability between SUNDIALS and
TRILINOS. This implementation is accompanied by additions to user documentation
and SUNDIALS examples.

A bug was fixed where a nonlinear solver object could be freed twice in some use
cases.

The ``EXAMPLES_ENABLE_RAJA`` CMake option has been removed. The option
:cmakeop:`EXAMPLES_ENABLE_CUDA` enables all examples that use CUDA including the
RAJA examples with a CUDA back end (if the RAJA ``N_Vector`` is enabled).

The implementation header file ``ida_impl.h`` is no longer installed. This means
users who are directly manipulating the ``IDAMem`` structure will need to update
their code to use IDA's public API.

Python is no longer required to run ``make test`` and ``make test_install``.

### IDAS Changes in v3.1.0

An additional ``N_Vector`` implementation was added for the TPETRA vector from
the TRILINOS library to facilitate interoperability between SUNDIALS and
TRILINOS. This implementation is accompanied by additions to user documentation
and SUNDIALS examples.

A bug was fixed where a nonlinear solver object could be freed twice in some use
cases.

The ``EXAMPLES_ENABLE_RAJA`` CMake option has been removed. The option
:cmakeop:`EXAMPLES_ENABLE_CUDA` enables all examples that use CUDA including the
RAJA examples with a CUDA back end (if the RAJA ``N_Vector`` is enabled).

The implementation header file ``idas_impl.h`` is no longer installed. This means
users who are directly manipulating the ``IDAMem`` structure will need to update
their code to use IDAS's public API.

Python is no longer required to run ``make test`` and ``make test_install``.

### KINSOL Changes in v4.1.0

An additional ``N_Vector`` implementation was added for the TPetra vector from the Trilinos library to
facilitate interoperability between SUNDIALS and Trilinos. This implementation is accompanied by additions
to user documentation and SUNDIALS examples.

The ``EXAMPLES_ENABLE_RAJA`` CMake option has been removed. The option ``EXAMPLES_ENABLE_CUDA`` enables all examples
that use CUDA including the RAJA examples with a CUDA back end (if the RAJA ``N_Vector`` is enabled).

The implementation header file ``kin_impl.h`` is no longer installed. This means users who are directly manipulating the
``KINMem`` structure will need to update their code to use KINSOL's public API.

Python is no longer required to run ``make test`` and ``make test_install``.


## Changes to SUNDIALS in release 4.0.2

Added information on how to contribute to SUNDIALS and a contributing agreement.

Moved definitions of DLS and SPILS backwards compatibility functions to a source
file. The symbols are now included in the appropriate package library, e.g.
`libsundials_cvode.lib`.

### ARKODE Changes in v3.0.2

Added information on how to contribute to SUNDIALS and a contributing agreement.

### CVODE Changes in v4.0.2

Added information on how to contribute to SUNDIALS and a
contributing agreement.

Moved definitions of DLS and SPILS backwards compatibility functions to
a source file. The symbols are now included in the CVODE library, ``libsundials_cvode``.

### CVODES Changes in v4.0.2

Added information on how to contribute to SUNDIALS and a contributing agreement.

Moved definitions of DLS and SPILS backwards compatibility functions to a source
file. The symbols are now included in the CVODES library,
``libsundials_cvodes``.

### IDA Changes in v4.0.2

Added information on how to contribute to SUNDIALS and a contributing agreement.

Moved definitions of DLS and SPILS backwards compatibility functions to a source
file. The symbols are now included in the IDA library, ``libsundials_ida``.

### IDAS Changes in v3.0.2

Added information on how to contribute to SUNDIALS and a contributing agreement.

Moved definitions of DLS and SPILS backwards compatibility functions to a source
file. The symbols are now included in the IDAS library, ``libsundials_idas``.

### KINSOL Changes in v4.0.2

Added information on how to contribute to SUNDIALS and a contributing agreement.

Moved definitions of DLS and SPILS backwards compatibility functions to a source file. The symbols are now included in
the KINSOL library, ``libsundials_kinsol``.


## Changes to SUNDIALS in release 4.0.1

A bug in ARKODE where single precision builds would fail to compile has been
fixed.

### ARKODE Changes in v3.0.1

A bug in ARKODE where single precision builds would fail to compile has been fixed.

### CVODE Changes in v4.0.1

No changes were made in this release.

### CVODES Changes in v4.0.1

No changes were made in this release.

### IDA Changes in v4.0.1

No changes were made in this release.

### IDAS Changes in v3.0.1

No changes were made in this release.

### KINSOL Changes in v4.0.1

No changes were made in this release.


## Changes to SUNDIALS in release 4.0.0

The direct and iterative linear solver interfaces in all SUNDIALS packages have
been merged into a single unified linear solver interface to support any valid
SUNLINSOL module. This includes the DIRECT and ITERATIVE types as well as the
new MATRIX_ITERATIVE type. Details regarding how SUNDIALS packages utilize
linear solvers of each type as well as discussion regarding intended use cases
for user-supplied SUNLINSOL implementations are included in the SUNLINSOL
chapter of the user guides. All example programs have been updated to use the
new unified interfaces.

The unified interface is very similar to the previous DLS and SPILS interfaces.
To minimize challenges in user migration to the unified linear solver interface,
the previous DLS and SPILS routines for all packages may still be used; these
will be deprecated in future releases, so we recommend that users migrate to the
new names soon. Additionally, we note that Fortran users will need to enlarge
their iout array of optional integer outputs, and update the indices that they
query for certain linear-solver-related statistics.

The names of all constructor routines for SUNDIALS-provided SUNLinSol
implementations have been updated to follow the naming convention SUNLinSol_*
where * is the name of the linear solver e.g., Dense, KLU, SPGMR, PCG, etc.
Solver-specific "set" routine names have been similarly standardized. To
minimize challenges in user migration to the new names, the previous routine
names may still be used; these will be deprecated in future releases, so we
recommend that users migrate to the new names soon. All example programs have
been updated to used the new naming convention.

The SUNBandMatrix constructor has been simplified to remove the storage upper
bandwidth argument.

SUNDIALS integrators (ARKODE, CVODE, CVODES, IDA, and IDAS) have been updated to
utilize generic nonlinear solver modules through the SUNNONLINSOL API. This API
will ease the addition of new nonlinear solver options and allow for external or
user-supplied nonlinear solvers. The SUNNONLINSOL API and provided SUNNONLINSOL
modules are described in a new user guide chapter and follow the same object
oriented design and implementation used by the NVECTOR, SUNMATRIX, and
SUNLINSOL modules. All integrator example programs have also been updated to
used the new nonlinear solver API.

Three fused vector operations and seven vector array operations have been added
to the NVECTOR API. These optional operations are disabled by default and may be
activated by calling vector specific routines after creating an NVECTOR. See the
NVECTOR chapter in the user guides for more information on the new operations.

Added a new NVECTOR (NVECTOR_OPENMPDEV) which leverages OpenMP 4.5+ device
offloading.

Multiple updates to the CUDA NVECTOR were made:

* Changed the `N_VMake_Cuda` function to take a host data pointer and a device
  data pointer instead of an `N_VectorContent_Cuda` object.

* Changed `N_VGetLength_Cuda` to return the global vector length instead of
  the local vector length.

* Added `N_VGetLocalLength_Cuda` to return the local vector length.

* Added `N_VGetMPIComm_Cuda` to return the MPI communicator used.

* Removed the accessor functions in the namespace suncudavec.

* Added the ability to set the `cudaStream_t` used for execution of the CUDA
  NVECTOR kernels. See the function `N_VSetCudaStreams_Cuda`.

* Added `N_VNewManaged_Cuda`, `N_VMakeManaged_Cuda`, and
  `N_VIsManagedMemory_Cuda` functions to accommodate using managed memory with
  the CUDA NVECTOR.

Multiple updates to the RAJA NVECTOR were made:

* Changed `N_VGetLength_Raja` to return the global vector length instead of
  the local vector length.

* Added `N_VGetLocalLength_Raja` to return the local vector length.

* Added `N_VGetMPIComm_Raja` to return the MPI communicator used.

* Removed the accessor functions in the namespace sunrajavec.

Two changes were made in the CVODE/CVODES/ARKODE initial step size algorithm:

  * Fixed an efficiency bug where an extra call to the RHS function was made.

  * Changed the behavior of the algorithm if the max-iterations case is hit.
    Before the algorithm would exit with the step size calculated on the
    penultimate iteration. Now it will exit with the step size calculated
    on the final iteration.

Fortran 2003 interfaces to CVODE, the fixed-point and Newton nonlinear solvers,
the dense, band, KLU, PCG, SPBCGS, SPFGMR, SPGMR, and SPTFQMR linear solvers,
and the serial, PThreads, and OpenMP NVECTORs have been added.

The ARKODE library has been entirely rewritten to support a modular approach to
one-step methods, which should allow for rapid research and development of novel
integration methods without affecting existing solver functionality.

A new ARKODE stepper, MRIStep, has been added for two rate explicit-explicit
multirate infinitesimal step methods.

ARKODE's dense output infrastructure has been improved to support higher-degree
Hermite polynomial interpolants (up to degree 5) over the last successful time
step.

### ARKODE Changes in v3.0.0

The ARKODE library has been entirely rewritten to support a modular
approach to one-step methods, which should allow rapid research and
development of novel integration methods without affecting existing
solver functionality.  To support this, the existing ARK-based methods
have been encapsulated inside the new ``ARKStep`` time-stepping
module. Two new time-stepping modules have been added:

* The ``ERKStep`` module provides an optimized implementation for explicit
  Runge--Kutta methods with reduced storage and number of calls to the ODE
  right-hand side function.

* The ``MRIStep`` module implements two-rate explicit-explicit multirate
  infinitesimal step methods utilizing different step sizes for slow
  and fast processes in an additive splitting.

This restructure has resulted in numerous small changes to the user
interface, particularly the suite of "Set" routines for user-provided
solver parameters and "Get" routines to access solver statistics,
that are now prefixed with the name of time-stepping module (e.g., ``ARKStep``
or ``ERKStep``) instead of ``ARKODE``.  Aside from affecting the names of these
routines, user-level changes have been kept to a minimum.  However, we recommend
that users consult both this documentation and the ARKODE example programs for
further details on the updated infrastructure.

As part of the ARKODE restructuring an :c:type:`ARKodeButcherTable` structure
has been added for storing Butcher tables. Functions for creating new Butcher
tables and checking their analytic order are provided along with other utility
routines. For more details see :numref:`ARKodeButcherTable`.

Two changes were made in the initial step size algorithm:

* Fixed an efficiency bug where an extra call to the right hand side function was made.

* Changed the behavior of the algorithm if the max-iterations case is hit.
  Before the algorithm would exit with the step size calculated on the
  penultimate iteration. Now it will exit with the step size calculated
  on the final iteration.

ARKODE's dense output infrastructure has been improved to support
higher-degree Hermite polynomial interpolants (up to degree 5) over
the last successful time step.

ARKODE's previous direct and iterative linear solver interfaces, ARKDLS and
ARKSPILS, have been merged into a single unified linear solver interface, ARKLS,
to support any valid SUNLINSOL module. This includes ``DIRECT`` and
``ITERATIVE`` types as well as the new ``MATRIX_ITERATIVE`` type. Details
regarding how ARKLS utilizes linear solvers of each type as well as discussion
regarding intended use cases for user-supplied SUNLinSol implementations are
included in the chapter :numref:`SUNLinSol`. All ARKODE examples programs and the
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
new names soon. All ARKODE example programs and the standalone linear solver
examples have been updated to use the new naming convention.

The ``SUNBandMatrix`` constructor has been simplified to remove the
storage upper bandwidth argument.

SUNDIALS integrators have been updated to utilize generic nonlinear solver
modules defined through the SUNNONLINSOL API. This API will ease the addition of
new nonlinear solver options and allow for external or user-supplied nonlinear
solvers. The SUNNONLINSOL API and SUNDIALS provided modules are described in
:numref:`SUNNonlinSol` and follow the same object oriented design and
implementation used by the NVector, SUNMatrix, and SUNLinSol modules. Currently
two SUNNONLINSOL implementations are provided, SUNNonlinSol_Newton and
SUNNonlinSol_FixedPoint. These replicate the previous integrator specific
implementations of a Newton iteration and an accelerated fixed-point iteration,
respectively. Example programs using each of these nonlinear solver modules in a
standalone manner have been added and all ARKODE example programs have been
updated to use generic SUNNonlinSol modules.

As with previous versions, ARKODE will use the Newton solver (now
provided by SUNNonlinSol_Newton) by default.  Use of the
:c:func:`ARKStepSetLinear()` routine (previously named
``ARKodeSetLinear``) will indicate that the problem is
linearly-implicit, using only a single Newton iteration per implicit
stage.  Users wishing to switch to the accelerated fixed-point solver
are now required to create a SUNNonlinSol_FixedPoint object and attach
that to ARKODE, instead of calling the previous
``ARKodeSetFixedPoint`` routine.  See the documentation sections
:numref:`ARKODE.Usage.ARKStep.Skeleton`,
:numref:`ARKODE.Usage.ARKStep.NonlinearSolvers`, and
:numref:`SUNNonlinSol.FixedPoint` for further details, or the serial C
example program ``ark_brusselator_fp.c`` for an example.

Three fused vector operations and seven vector array operations have been added
to the NVECTOR API. These *optional* operations are disabled by default and may
be activated by calling vector specific routines after creating an NVector (see
:numref:`NVectors.Description` for more details). The new operations are intended
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
been added, NVECTOR_OpenMPDEV. See :numref:`NVectors.OpenMPDEV` for more details.



### CVODE Changes in v4.0.0

CVODE's previous direct and iterative linear solver interfaces,
CVDLS and CVSPILS, have been merged into a single unified linear
solver interface, CVLS, to support any valid ``SUNLinearSolver`` module.
This includes the "DIRECT" and "ITERATIVE" types as well as the new
"MATRIX_ITERATIVE" type. Details regarding how CVLS utilizes linear
solvers of each type as well as discussion regarding intended use cases
for user-supplied ``SUNLinearSolver`` implementations are included in
:numref:`SUNLinSol`. All CVODE example programs
and the standalone linear solver examples have been updated to use the
unified linear solver interface.

The unified interface for the new CVLS module is very similar to the
previous CVDLS and CVSPILS interfaces. To minimize challenges in
user migration to the new names, the previous C and Fortran routine
names may still be used; these will be deprecated in future releases, so
we recommend that users migrate to the new names soon. Additionally, we
note that Fortran users, however, may need to enlarge their array of
optional integer outputs, and update the indices that they query for
certain linear-solver-related statistics.

The names of all constructor routines for SUNDIALS-provided ``SUNLinearSolver``
implementations have been updated to follow the naming convention
``SUNLinSol_*`` where is the name of the linear solver. Solver-specific "set"
routine names have been similarly standardized. To minimize challenges in user
migration to the new names, the previous routine names may still be used; these
will be deprecated in future releases, so we recommend that users migrate to the
new names soon. All CVODE example programs and the standalone linear solver
examples have been updated to use the new naming convention.

The :ref:`SUNMATRIX_BAND <SUNMatrix.Band>` constructor has been simplified to remove the storage upper
bandwidth argument.

SUNDIALS integrators have been updated to utilize generic nonlinear
solver modules defined through the ``SUNNonlinearSolver`` API. This API will
ease the addition of new nonlinear solver options and allow for external
or user-supplied nonlinear solvers. The ``SUNNonlinearSolver`` API and
SUNDIALS provided modules are described in
:numref:`SUNNonlinSol` and follow the same
object oriented design and implementation used by the ``N_Vector``,
``SUNMatrix``, and ``SUNLinearSolver`` modules. Currently two ``SUNNonlinearSolver``
implementations are provided, ``SUNNONLINSOL_NEWTON`` and
``SUNNONLINSOL_FIXEDPOINT``. These replicate the previous integrator
specific implementations of a Newton iteration and a fixed-point
iteration (previously referred to as a functional iteration),
respectively. Note the ``SUNNONLINSOL_FIXEDPOINT`` module can optionally
utilize Anderson's method to accelerate convergence. Example programs
using each of these nonlinear solver modules in a standalone manner have
been added and all CVODE example programs have been updated to use
generic ``SUNNonlinearSolver`` modules.

With the introduction of ``SUNNonlinearSolver`` modules, the ``iter`` input
parameter to :c:func:`CVodeCreate` has been removed along with the function
:c:func:`CVodeSetIterType` and the constants ``CV_NEWTON`` and
``CV_FUNCTIONAL``. Similarly, the parameter has been removed from the Fortran
interface function ``FCVMALLOC``. Instead of specifying the nonlinear iteration
type when creating the CVODE memory structure, CVODE uses the
``SUNNONLINSOL_NEWTON`` module implementation of a Newton iteration by default.
For details on using a non-default or user-supplied nonlinear solver see
:numref:CVODE.Usage.CC. CVODE functions for setting the nonlinear solver options
(e.g., :c:func:`CVodeSetMaxNonlinIters`) or getting nonlinear solver statistics
(e.g., :c:func:`CVodeGetNumNonlinSolvIters`) remain unchanged and internally
call generic ``SUNNonlinearSolver`` functions as needed.

Three fused vector operations and seven vector array operations have been added
to the ``N_Vector`` API. These *optional* operations are disabled by default and
may be activated by calling vector specific routines after creating an
``N_Vector`` (see :numref:`NVectors` for more details). The new operations
are intended to increase data reuse in vector operations, reduce parallel
communication on distributed memory systems, and lower the number of kernel
launches on systems with accelerators. The fused operations are
:c:func:`N_VLinearCombination`,  :c:func:`N_VScaleAddMulti`, and
:c:func:`N_VDotProdMulti` and the vector array operations are
:c:func:`N_VLinearCombinationVectorArray`, :c:func:`N_VScaleVectorArray`,
:c:func:`N_VConstVectorArray`, :c:func:`N_VWrmsNormVectorArray`,
:c:func:`N_VWrmsNormMaskVectorArray`, and :c:func:`N_VScaleAddMultiVectorArray`.
If an ``N_Vector`` implementation defines any of these operations as, then
standard ``N_Vector`` operations will automatically be called as necessary to
complete the computation.

Multiple updates to ``NVECTOR_CUDA`` were made:

* Changed to return the global vector length instead of the local
  vector length.
* Added to return the local vector length.
* Added to return the MPI communicator used.
* Removed the accessor functions in the namespace suncudavec.
* Changed the function to take a host data pointer and a device data
  pointer instead of an object.
* Added the ability to set the used for execution of the ``NVECTOR_CUDA``
  kernels. See the function :c:func:`N_VSetCudaStream_Cuda()`.
* Added :c:func:`N_VNewManaged_Cuda`, :c:func:`N_VMakeManaged_Cuda`, and
  :c:func:`N_VIsManagedMemory_Cuda()` functions to accommodate using managed
  memory with ``NVECTOR_CUDA``.

Multiple changes to ``NVECTOR_RAJA`` were made:

   - Changed to return the global vector length instead of the local vector length.
   - Added to return the local vector length.
   - Added to return the MPI communicator used.
   - Removed the accessor functions in the namespace suncudavec.
   - A new ``N_Vector`` implementation for leveraging OpenMP 4.5+ device
     offloading has been added, ``NVECTOR_OPENMPDEV``.
   - Two changes were made in the CVODE/CVODES/ARKODE initial step size algorithm:

     - Fixed an efficiency bug where an extra call to the right hand side function was made.
     - Changed the behavior of the algorithm if the max-iterations case is hit. Before the algorithm would exit with the step size calculated on the penultimate iteration. Now it will exit with the step size calculated on the final iteration.

A Fortran 2003 interface to CVODE has been added along with Fortran 2003
interfaces to the following shared SUNDIALS modules:

   -  ``SUNNONLINSOL_FIXEDPOINT`` and ``SUNNONLINSOL_NEWTON`` nonlinear solver modules
   -  ``SUNLINSOL_BAND``, ``SUNLINSOL_DENSE``, ``SUNLINSOL_KLU``, ``SUNLINSOL_PCG``, ``SUNLINSOL_SPBCGS``, ``SUNLINSOL_SPFGMR``, ``SUNLINSOL_SPGMR``, and ``SUNLINSOL_SPTFQMR`` linear solver modules
   -  ``NVECTOR_SERIAL``, ``NVECTOR_PTHREADS``, and ``NVECTOR_OPENMP`` vector modules

### CVODES Changes in v4.0.0

CVODES' previous direct and iterative linear solver interfaces, CVDLS and
CVSPILS, have been merged into a single unified linear solver interface,
CVLS, to support any valid ``SUNLinearSolver`` module. This includes the
"DIRECT" and "ITERATIVE" types as well as the new "MATRIX_ITERATIVE" type.
Details regarding how CVLS utilizes linear solvers of each type as well as
discussion regarding intended use cases for user-supplied ``SUNLinearSolver``
implementations are included in Chapter :numref:`SUNLinSol`. All CVODES example
programs and the standalone linear solver examples have been updated to use the
unified linear solver interface.

The unified interface for the new CVLS module is very similar to the previous
CVDLS and CVSPILS interfaces. To minimize challenges in user
migration to the new names, the previous C routine names may still be used;
these will be deprecated in future releases, so we recommend that users migrate
to the new names soon.

The names of all constructor routines for SUNDIALS-provided ``SUNLinearSolver``
implementations have been updated to follow the naming convention
``SUNLinSol_*`` where ``*`` is the name of the linear solver. The new names are
``SUNLinSol_Band``, ``SUNLinSol_Dense``, ``SUNLinSol_KLU``,
``SUNLinSol_LapackBand``, ``SUNLinSol_LapackDense``, ``SUNLinSol_PCG``,
``SUNLinSol_SPBCGS``, ``SUNLinSol_SPFGMR``, ``SUNLinSol_SPGMR``,
``SUNLinSol_SPTFQMR``, and ``SUNLinSol_SuperLUMT``. Solver-specific "set"
routine names have been similarly standardized. To minimize challenges in user
migration to the new names, the previous routine names may still be used; these
will be deprecated in future releases, so we recommend that users migrate to the
new names soon. All CVODES example programs and the standalone linear solver
examples have been updated to use the new naming convention.

The ``SUNBandMatrix`` constructor has been simplified to remove the storage
upper bandwidth argument.

SUNDIALS integrators have been updated to utilize generic nonlinear solver
modules defined through the ``SUNNonlinearSolver`` API. This API will ease the
addition of new nonlinear solver options and allow for external or user-supplied
nonlinear solvers. The ``SUNNonlinearSolver`` API and SUNDIALS provided modules
are described in Chapter :numref:`SUNNonlinSol` and follow the same object
oriented design and implementation used by the ``N_Vector``, ``SUNMatrix``, and
``SUNLinearSolver`` modules. Currently two ``SUNNonlinearSolver``
implementations are provided, ``SUNNONLINSOL_NEWTON`` and
``SUNNONLINSOL_FIXEDPOINT``. These replicate the previous integrator specific
implementations of a Newton iteration and a fixed-point iteration (previously
referred to as a functional iteration), respectively. Note the
``SUNNONLINSOL_FIXEDPOINT`` module can optionally utilize Anderson's method to
accelerate convergence. Example programs using each of these nonlinear solver
modules in a standalone manner have been added and all CVODES example programs
have been updated to use generic ``SUNNonlinearSolver`` modules.

With the introduction of ``SUNNonlinearSolver`` modules, the input parameter
``iter`` to ``CVodeCreate`` has been removed along with the function
``CVodeSetIterType`` and the constants ``CV_NEWTON`` and ``CV_FUNCTIONAL``.
Instead of specifying the nonlinear iteration type when creating the CVODES
memory structure, CVODES uses the ``SUNNONLINSOL_NEWTON`` module implementation
of a Newton iteration by default. For details on using a non-default or
user-supplied nonlinear solver see Chapters :numref:`CVODES.Usage.SIM`,
:numref:`CVODES.Usage.FSA`, and :numref:`CVODES.Usage.ADJ`. CVODES functions for
setting the nonlinear solver options (e.g., ``CVodeSetMaxNonlinIters``) or
getting nonlinear solver statistics (e.g., ``CVodeGetNumNonlinSolvIters``)
remain unchanged and internally call generic ``SUNNonlinearSolver`` functions as
needed.

Three fused vector operations and seven vector array operations have been added
to the ``N_Vector`` API. These *optional* operations are disabled by default and
may be activated by calling vector specific routines after creating an
``N_Vector`` (see Chapter :numref:`NVectors` for more details). The
new operations are intended to increase data reuse in vector operations, reduce
parallel communication on distributed memory systems, and lower the number of
kernel launches on systems with accelerators. The fused operations are
``N_VLinearCombination``, ``N_VScaleAddMulti``, and ``N_VDotProdMulti`` and the
vector array operations are ``N_VLinearCombinationVectorArray``,
``N_VScaleVectorArray``, ``N_VConstVectorArray``, ``N_VWrmsNormVectorArray``,
``N_VWrmsNormMaskVectorArray``, ``N_VScaleAddMultiVectorArray``, and
``N_VLinearCombinationVectorArray``. If an ``N_Vector`` implementation defines
any of these operations as ``NULL``, then standard ``N_Vector`` operations will
automatically be called as necessary to complete the computation. Multiple
updates to ``NVECTOR_CUDA`` were made:

-  Changed ``N_VGetLength_Cuda`` to return the global vector length instead of the local vector length.

-  Added ``N_VGetLocalLength_Cuda`` to return the local vector length.

-  Added ``N_VGetMPIComm_Cuda`` to return the MPI communicator used.

-  Removed the accessor functions in the namespace suncudavec.

-  Changed the ``N_VMake_Cuda`` function to take a host data pointer and a device data pointer instead of an
   ``N_VectorContent_Cuda`` object.

-  Added the ability to set the ``cudaStream_t`` used for execution of the ``NVECTOR_CUDA`` kernels. See the function
   ``N_VSetCudaStreams_Cuda``.

-  Added ``N_VNewManaged_Cuda``, ``N_VMakeManaged_Cuda``, and ``N_VIsManagedMemory_Cuda`` functions to accommodate using
   managed memory with the ``NVECTOR_CUDA``.

Multiple changes to ``NVECTOR_RAJA`` were made:

-  Changed ``N_VGetLength_Raja`` to return the global vector length instead of the local vector length.

-  Added ``N_VGetLocalLength_Raja`` to return the local vector length.

-  Added ``N_VGetMPIComm_Raja`` to return the MPI communicator used.

-  Removed the accessor functions in the namespace suncudavec.

A new ``N_Vector`` implementation for leveraging OpenMP 4.5+ device offloading
has been added, ``NVECTOR_OPENMPDEV``. See :numref:`NVectors.openmpdev` for more
details. Two changes were made in the CVODE/CVODES/ARKODE initial step size
algorithm:

#. Fixed an efficiency bug where an extra call to the right hand side function was made.

#. Changed the behavior of the algorithm if the max-iterations case is hit. Before the algorithm would exit with the
   step size calculated on the penultimate iteration. Now it will exit with the step size calculated on the final
   iteration.

### IDA Changes in v4.0.0

IDA's previous direct and iterative linear solver interfaces, IDADLS and
IDASPILS, have been merged into a single unified linear solver interface, IDALS,
to support any valid ``SUNLinearSolver`` module.  This includes the "DIRECT" and
"ITERATIVE" types as well as the new "MATRIX_ITERATIVE" type. Details regarding
how IDALS utilizes linear solvers of each type as well as discussion regarding
intended use cases for user-supplied ``SUNLinearSolver`` implementations are
included in :numref:`SUNLinSol`. All IDA example programs and the standalone
linear solver examples have been updated to use the unified linear solver
interface.

The unified interface for the new IDALS module is very similar to the previous
IDADLS and IDASPILS interfaces. To minimize challenges in user migration to the
new names, the previous C and Fortran routine names may still be used; these
will be deprecated in future releases, so we recommend that users migrate to the
new names soon. Additionally, we note that Fortran users, however, may need to
enlarge their ``iout`` array of optional integer outputs, and update the indices
that they query for certain linear-solver-related statistics.

The names of all constructor routines for SUNDIALS-provided ``SUNLinearSolver``
implementations have been updated to follow the naming convention ``SUNLinSol_``
where ``*`` is the name of the linear solver. The new names are
:c:func:`SUNLinSol_Band`, :c:func:`SUNLinSol_Dense`, :c:func:`SUNLinSol_KLU`,
:c:func:`SUNLinSol_LapackBand`, :c:func:`SUNLinSol_LapackDense`,
:c:func:`SUNLinSol_PCG`, :c:func:`SUNLinSol_SPBCGS`, :c:func:`SUNLinSol_SPFGMR`,
:c:func:`SUNLinSol_SPGMR`, :c:func:`SUNLinSol_SPTFQMR`, and
:c:func:`SUNLinSol_SuperLUMT`. Solver-specific "set" routine names have been
similarly standardized. To minimize challenges in user migration to the new
names, the previous routine names may still be used; these will be deprecated in
future releases, so we recommend that users migrate to the new names soon. All
IDA example programs and the standalone linear solver examples have been updated
to use the new naming convention.

The ``SUNBandMatrix`` constructor has been simplified to remove the storage
upper bandwidth argument.

SUNDIALS integrators have been updated to utilize generic nonlinear solver
modules defined through the ``SUNNonlinearSolver`` API. This API will ease the
addition of new nonlinear solver options and allow for external or user-supplied
nonlinear solvers. The ``SUNNonlinearSolver`` API and SUNDIALS provided modules
are described in :numref:`SUNNonlinSol` and follow the same object oriented
design and implementation used by the ``N_Vector``, ``SUNMatrix``, and
``SUNLinearSolver`` modules. Currently two ``SUNNonlinearSolver``
implementations are provided, :ref:`SUNNONLINSOL_NEWTON <SUNNonlinSol.Newton>`
and
:ref:`SUNNONLINSOL_FIXEDPOINT <SUNNonlinSol.FixedPoint>`. These replicate the
previous integrator specific implementations of a Newton iteration and a
fixed-point iteration (previously referred to as a functional iteration),
respectively. Note the :ref:`SUNNONLINSOL_FIXEDPOINT <SUNNonlinSol.FixedPoint>`
module can optionally utilize Anderson's method to accelerate
convergence. Example programs using each of these nonlinear solver modules in a
standalone manner have been added and all IDA example programs have been updated
to use generic ``SUNNonlinearSolver`` modules.

By default IDA uses the :ref:`SUNNONLINSOL_NEWTON <SUNNonlinSol.Newton>`
module. Since IDA previously only used an internal implementation of a Newton
iteration no changes are required to user programs and functions for setting the
nonlinear solver options (e.g., :c:func:`IDASetMaxNonlinIters`) or getting
nonlinear solver statistics (e.g., :c:func:`IDAGetNumNonlinSolvIters`) remain
unchanged and internally call generic ``SUNNonlinearSolver`` functions as
needed. While SUNDIALS includes a fixed-point nonlinear solver module, it is not
currently supported in IDA. For details on attaching a user-supplied nonlinear
solver to IDA see :numref:IDA.Usage.CC. Additionally, the example program
``idaRoberts_dns.c`` explicitly creates an attaches a :ref:`SUNNONLINSOL_NEWTON
<SUNNonlinSol.Newton>` object to demonstrate the process of creating and
attaching a nonlinear solver module (note this is not necessary in general as
IDA uses the :ref:`SUNNONLINSOL_NEWTON <SUNNonlinSol.Newton>` module by
default).

Three fused vector operations and seven vector array operations have been added
to the ``N_Vector`` API. These *optional* operations are disabled by default and
may be activated by calling vector specific routines after creating an
``N_Vector`` (see :numref:`NVectors` for more details). The new operations are
intended to increase data reuse in vector operations, reduce parallel
communication on distributed memory systems, and lower the number of kernel
launches on systems with accelerators. The fused operations are
:c:func:`N_VLinearCombination`, :c:func:`N_VScaleAddMulti`, and
:c:func:`N_VDotProdMulti` and the vector array operations are
:c:func:`N_VLinearCombinationVectorArray`, :c:func:`N_VScaleVectorArray`,
:c:func:`N_VConstVectorArray`, :c:func:`N_VWrmsNormVectorArray`,
:c:func:`N_VWrmsNormMaskVectorArray`, :c:func:`N_VScaleAddMultiVectorArray`, and
:c:func:`N_VLinearCombinationVectorArray`.

If an ``N_Vector`` implementation defines any of these operations as ``NULL``,
then standard ``N_Vector`` operations will automatically be called as necessary
to complete the computation.

Multiple updates to :ref:`NVECTOR_CUDA <NVectors.CUDA>` were made:

* Changed :c:func:`N_VGetLength_Cuda` to return the global vector length instead
  of the local vector length.

* Added :c:func:`N_VGetLocalLength_Cuda` to return the local vector length.

* Added :c:func:`N_VGetMPIComm_Cuda` to return the MPI communicator used.

* Removed the accessor functions in the namespace ``suncudavec``.

* Changed the :c:func:`N_VMake_Cuda` function to take a host data pointer and a
  device data pointer instead of an ``N_VectorContent_Cuda`` object.

* Added the ability to set the ``cudaStream_t`` used for execution of the
  :ref:`NVECTOR_CUDA <NVectors.CUDA>` kernels. See the function
  :c:func:`N_VSetCudaStreams_Cuda`.

* Added :c:func:`N_VNewManaged_Cuda`, :c:func:`N_VMakeManaged_Cuda`, and
  :c:func:`N_VIsManagedMemory_Cuda` functions to accommodate using managed
  memory with the :ref:`NVECTOR_CUDA <NVectors.CUDA>`.

Multiple changes to :ref:`NVECTOR_RAJA <NVectors.RAJA>` were made:

* Changed :c:func:`N_VGetLength_Raja` to return the global vector length instead
  of the local vector length.

* Added :c:func:`N_VGetLocalLength_Raja` to return the local vector length.

* Added :c:func:`N_VGetMPIComm_Raja` to return the MPI communicator used.

* Removed the accessor functions in the namespace ``suncudavec``.

A new ``N_Vector`` implementation for leveraging OpenMP 4.5+ device offloading
has been added, :ref:`NVECTOR_OPENMPDEV <NVectors.OpenMPDEV>`. See
:numref:`NVectors.OpenMPDEV` for more details.

### IDAS Changes in v3.0.0

IDA's previous direct and iterative linear solver interfaces, IDADLS and
IDASPILS, have been merged into a single unified linear solver interface, IDALS,
to support any valid ``SUNLinearSolver`` module.  This includes the "DIRECT" and
"ITERATIVE" types as well as the new "MATRIX_ITERATIVE" type. Details regarding
how IDALS utilizes linear solvers of each type as well as discussion regarding
intended use cases for user-supplied ``SUNLinearSolver`` implementations are
included in :numref:`SUNLinSol`. All IDAS example programs and the standalone
linear solver examples have been updated to use the unified linear solver
interface.

The unified interface for the new IDALS module is very similar to the previous
IDADLS and IDASPILS interfaces. To minimize challenges in user migration to the
new names, the previous C and Fortran routine names may still be used; these
will be deprecated in future releases, so we recommend that users migrate to the
new names soon. Additionally, we note that Fortran users, however, may need to
enlarge their ``iout`` array of optional integer outputs, and update the indices
that they query for certain linear-solver-related statistics.

The names of all constructor routines for SUNDIALS-provided ``SUNLinearSolver``
implementations have been updated to follow the naming convention ``SUNLinSol_``
where ``*`` is the name of the linear solver. The new names are
:c:func:`SUNLinSol_Band`, :c:func:`SUNLinSol_Dense`, :c:func:`SUNLinSol_KLU`,
:c:func:`SUNLinSol_LapackBand`, :c:func:`SUNLinSol_LapackDense`,
:c:func:`SUNLinSol_PCG`, :c:func:`SUNLinSol_SPBCGS`, :c:func:`SUNLinSol_SPFGMR`,
:c:func:`SUNLinSol_SPGMR`, :c:func:`SUNLinSol_SPTFQMR`, and
:c:func:`SUNLinSol_SuperLUMT`. Solver-specific "set" routine names have been
similarly standardized. To minimize challenges in user migration to the new
names, the previous routine names may still be used; these will be deprecated in
future releases, so we recommend that users migrate to the new names soon. All
IDAS example programs and the standalone linear solver examples have been updated
to use the new naming convention.

The ``SUNBandMatrix`` constructor has been simplified to remove the storage
upper bandwidth argument.

SUNDIALS integrators have been updated to utilize generic nonlinear solver
modules defined through the ``SUNNonlinearSolver`` API. This API will ease the
addition of new nonlinear solver options and allow for external or user-supplied
nonlinear solvers. The ``SUNNonlinearSolver`` API and SUNDIALS provided modules
are described in :numref:`SUNNonlinSol` and follow the same object oriented
design and implementation used by the ``N_Vector``, ``SUNMatrix``, and
``SUNLinearSolver`` modules. Currently two ``SUNNonlinearSolver``
implementations are provided, :ref:`SUNNONLINSOL_NEWTON <SUNNonlinSol.Newton>`
and
:ref:`SUNNONLINSOL_FIXEDPOINT <SUNNonlinSol.FixedPoint>`. These replicate the
previous integrator specific implementations of a Newton iteration and a
fixed-point iteration (previously referred to as a functional iteration),
respectively. Note the :ref:`SUNNONLINSOL_FIXEDPOINT <SUNNonlinSol.FixedPoint>`
module can optionally utilize Anderson's method to accelerate
convergence. Example programs using each of these nonlinear solver modules in a
standalone manner have been added and all IDAS example programs have been updated
to use generic ``SUNNonlinearSolver`` modules.

By default IDAS uses the :ref:`SUNNONLINSOL_NEWTON <SUNNonlinSol.Newton>`
module. Since IDAS previously only used an internal implementation of a Newton
iteration no changes are required to user programs and functions for setting the
nonlinear solver options (e.g., :c:func:`IDASetMaxNonlinIters`) or getting
nonlinear solver statistics (e.g., :c:func:`IDAGetNumNonlinSolvIters`) remain
unchanged and internally call generic ``SUNNonlinearSolver`` functions as
needed. While SUNDIALS includes a fixed-point nonlinear solver module, it is not
currently supported in IDAS. For details on attaching a user-supplied nonlinear
solver to IDAS see :numref:`IDAS.Usage`. Additionally, the example program
``idaRoberts_dns.c`` explicitly creates an attaches a :ref:`SUNNONLINSOL_NEWTON
<SUNNonlinSol.Newton>` object to demonstrate the process of creating and
attaching a nonlinear solver module (note this is not necessary in general as
IDAS uses the :ref:`SUNNONLINSOL_NEWTON <SUNNonlinSol.Newton>` module by
default).

Three fused vector operations and seven vector array operations have been added
to the ``N_Vector`` API. These *optional* operations are disabled by default and
may be activated by calling vector specific routines after creating an
``N_Vector`` (see :numref:`NVectors` for more details). The new operations are
intended to increase data reuse in vector operations, reduce parallel
communication on distributed memory systems, and lower the number of kernel
launches on systems with accelerators. The fused operations are
:c:func:`N_VLinearCombination`, :c:func:`N_VScaleAddMulti`, and
:c:func:`N_VDotProdMulti` and the vector array operations are
:c:func:`N_VLinearCombinationVectorArray`, :c:func:`N_VScaleVectorArray`,
:c:func:`N_VConstVectorArray`, :c:func:`N_VWrmsNormVectorArray`,
:c:func:`N_VWrmsNormMaskVectorArray`, :c:func:`N_VScaleAddMultiVectorArray`, and
:c:func:`N_VLinearCombinationVectorArray`.

If an ``N_Vector`` implementation defines any of these operations as ``NULL``,
then standard ``N_Vector`` operations will automatically be called as necessary
to complete the computation.

Multiple updates to :ref:`NVECTOR_CUDA <NVectors.CUDA>` were made:

* Changed :c:func:`N_VGetLength_Cuda` to return the global vector length instead
  of the local vector length.

* Added :c:func:`N_VGetLocalLength_Cuda` to return the local vector length.

* Added :c:func:`N_VGetMPIComm_Cuda` to return the MPI communicator used.

* Removed the accessor functions in the namespace ``suncudavec``.

* Changed the :c:func:`N_VMake_Cuda` function to take a host data pointer and a
  device data pointer instead of an ``N_VectorContent_Cuda`` object.

* Added the ability to set the ``cudaStream_t`` used for execution of the
  :ref:`NVECTOR_CUDA <NVectors.CUDA>` kernels. See the function
  :c:func:`N_VSetCudaStreams_Cuda`.

* Added :c:func:`N_VNewManaged_Cuda`, :c:func:`N_VMakeManaged_Cuda`, and
  :c:func:`N_VIsManagedMemory_Cuda` functions to accommodate using managed
  memory with the :ref:`NVECTOR_CUDA <NVectors.CUDA>`.

Multiple changes to :ref:`NVECTOR_RAJA <NVectors.RAJA>` were made:

* Changed :c:func:`N_VGetLength_Raja` to return the global vector length instead
  of the local vector length.

* Added :c:func:`N_VGetLocalLength_Raja` to return the local vector length.

* Added :c:func:`N_VGetMPIComm_Raja` to return the MPI communicator used.

* Removed the accessor functions in the namespace ``suncudavec``.

A new ``N_Vector`` implementation for leveraging OpenMP 4.5+ device offloading
has been added, :ref:`NVECTOR_OPENMPDEV <NVectors.OpenMPDEV>`. See
:numref:`NVectors.OpenMPDEV` for more details.

### KINSOL Changes in v4.0.0

KINSOL's previous direct and iterative linear solver interfaces, KINDls and KINSpils, have been merged
into a single unified linear solver interface, KINLs, to support any valid ``SUNLinearSolver`` module. This includes
the "DIRECT" and "ITERATIVE" types as well as the new "MATRIX_ITERATIVE" type. Details regarding how KINLs
utilizes linear solvers of each type as well as discussion regarding intended use cases for user-supplied
``SUNLinearSolver`` implementations are included in Chapter :numref:`SUNLinSol`. All KINSOL example
programs and the standalone linear solver examples have been updated to use the unified linear solver interface.

The unified interface for the new KINLs module is very similar to the previous KINDls and KINSpils
interfaces. To minimize challenges in user migration to the new names, the previous C and Fortran routine names
may still be used; these will be deprecated in future releases, so we recommend that users migrate to the new names
soon. Additionally, we note that Fortran users, however, may need to enlarge their ``iout`` array of optional integer
outputs, and update the indices that they query for certain linear-solver-related statistics.

The names of all constructor routines for SUNDIALS-provided ``SUNLinearSolver`` implementations have been updated to
follow the naming convention ``SUNLinSol_*`` where ``*`` is the name of the linear solver. The new names are
``SUNLinSol_Band``, ``SUNLinSol_Dense``, ``SUNLinSol_KLU``, ``SUNLinSol_LapackBand``, ``SUNLinSol_LapackDense``,
``SUNLinSol_PCG``, ``SUNLinSol_SPBCGS``, ``SUNLinSol_SPFGMR``, ``SUNLinSol_SPGMR``, ``SUNLinSol_SPTFQMR``, and
``SUNLinSol_SuperLUMT``. Solver-specific "set" routine names have been similarly standardized. To minimize challenges in
user migration to the new names, the previous routine names may still be used; these will be deprecated in future
releases, so we recommend that users migrate to the new names soon. All KINSOL example programs and the standalone
linear solver examples have been updated to use the new naming convention.

The ``SUNBandMatrix`` constructor has been simplified to remove the storage upper bandwidth argument.

Three fused vector operations and seven vector array operations have been added to the ``N_Vector`` API.
These *optional* operations are disabled by default and may be activated by calling vector specific
routines after creating an ``N_Vector`` (see Chapter :numref:`NVectors` for more details). The new
operations are intended to increase data reuse in vector operations, reduce parallel communication on
distributed memory systems, and lower the number of kernel launches on systems with accelerators. The
fused operations are ``N_VLinearCombination``, ``N_VScaleAddMulti``, and ``N_VDotProdMulti`` and the
vector array operations are ``N_VLinearCombinationVectorArray``, ``N_VScaleVectorArray``,
``N_VConstVectorArray``, ``N_VWrmsNormVectorArray``, ``N_VWrmsNormMaskVectorArray``,
``N_VScaleAddMultiVectorArray``, and ``N_VLinearCombinationVectorArray``. If an ``N_Vector``
implementation defines any of these operations as ``NULL``, then standard ``N_Vector`` operations will
automatically be called as necessary to complete the computation. Multiple updates to ``NVECTOR_CUDA``
were made:

-  Changed ``N_VGetLength_Cuda`` to return the global vector length instead of the local vector length.

-  Added ``N_VGetLocalLength_Cuda`` to return the local vector length.

-  Added ``N_VGetMPIComm_Cuda`` to return the MPI communicator used.

-  Removed the accessor functions in the namespace suncudavec.

-  Changed the ``N_VMake_Cuda`` function to take a host data pointer and a device data pointer instead of an
   ``N_VectorContent_Cuda`` object.

-  Added the ability to set the ``cudaStream_t`` used for execution of the ``NVECTOR_CUDA`` kernels. See the function
   ``N_VSetCudaStreams_Cuda``.

-  Added ``N_VNewManaged_Cuda``, ``N_VMakeManaged_Cuda``, and ``N_VIsManagedMemory_Cuda`` functions to accommodate using
   managed memory with the ``NVECTOR_CUDA``.

Multiple changes to ``NVECTOR_RAJA`` were made:

-  Changed ``N_VGetLength_Raja`` to return the global vector length instead of the local vector length.

-  Added ``N_VGetLocalLength_Raja`` to return the local vector length.

-  Added ``N_VGetMPIComm_Raja`` to return the MPI communicator used.

-  Removed the accessor functions in the namespace suncudavec.

A new ``N_Vector`` implementation for leveraging OpenMP 4.5+ device offloading has been added, ``NVECTOR_OPENMPDEV``. See
:numref:`NVectors.OpenMPDEV` for more details.


## Changes to SUNDIALS in release 3.2.1

Fixed a bug in the CUDA NVECTOR where the `N_VInvTest` operation could write
beyond the allocated vector data.

Fixed library installation path for multiarch systems. This fix changes the
default library installation path to `CMAKE_INSTALL_PREFIX/CMAKE_INSTALL_LIBDIR`
from `CMAKE_INSTALL_PREFIX/lib`. `CMAKE_INSTALL_LIBDIR` is automatically set,
but is available as a CMAKE option that can modified.

### ARKODE Changes in v2.2.1

Fixed a bug in the CUDA NVECTOR where the ``N_VInvTest`` operation could
write beyond the allocated vector data.

Fixed library installation path for multiarch systems. This fix changes the default
library installation path to ``CMAKE_INSTALL_PREFIX/CMAKE_INSTALL_LIBDIR``
from ``CMAKE_INSTALL_PREFIX/lib``. ``CMAKE_INSTALL_LIBDIR`` is automatically
set, but is available as a CMAKE option that can modified.


### CVODE Changes in v3.2.1

The changes in this minor release include the following:

-  Fixed a bug in the CUDA ``N_Vector`` where the operation could
   write beyond the allocated vector data.

-  Fixed library installation path for multiarch systems. This fix
   changes the default library installation path to
   ``CMAKE_INSTALL_PREFIX/CMAKE_INSTALL_LIBDIR`` from
   ``CMAKE_INSTALL_PREFIX/lib``. ``CMAKE_INSTALL_LIBDIR`` is automatically set,
   but is available as a CMake option that can modified.

### CVODES Changes in v3.2.1

The changes in this minor release include the following:

-  Fixed a bug in the CUDA ``N_Vector`` where the ``N_VInvTest`` operation could
   write beyond the allocated vector data.

-  Fixed library installation path for multiarch systems. This fix changes the default library installation path to
   ``CMAKE_INSTALL_PREFIX/CMAKE_INSTALL_LIBDIR`` from
   ``CMAKE_INSTALL_PREFIX/lib``. ``CMAKE_INSTALL_LIBDIR`` is automatically set,
   but is available as a CMake option that can modified.

### IDA Changes in v3.2.1

The changes in this minor release include the following:

* Fixed a bug in the :ref:`CUDA N_Vector <NVectors.CUDA>` where the
  :c:func:`N_VInvTest` operation could write beyond the allocated vector data.

* Fixed library installation path for multiarch systems. This fix changes the
  default library installation path to
  ``CMAKE_INSTALL_PREFIX/CMAKE_INSTALL_LIBDIR`` from
  ``CMAKE_INSTALL_PREFIX/lib``. Note :cmakeop:`CMAKE_INSTALL_LIBDIR` is
  automatically set, but is available as a CMake option that can be modified.

### IDAS Changes in v2.2.1

The changes in this minor release include the following:

* Fixed a bug in the :ref:`CUDA N_Vector <NVectors.CUDA>` where the
  :c:func:`N_VInvTest` operation could write beyond the allocated vector data.

* Fixed library installation path for multiarch systems. This fix changes the
  default library installation path to
  ``CMAKE_INSTALL_PREFIX/CMAKE_INSTALL_LIBDIR`` from
  ``CMAKE_INSTALL_PREFIX/lib``. Note :cmakeop:`CMAKE_INSTALL_LIBDIR` is
  automatically set, but is available as a CMake option that can be modified.

### KINSOL Changes in v3.2.1

The changes in this minor release include the following:

-  Fixed a bug in the CUDA ``N_Vector`` where the ``N_VInvTest`` operation could write beyond the allocated
   vector data.

-  Fixed library installation path for multiarch systems. This fix changes the default library
   installation path to ``CMAKE_INSTALL_PREFIX/CMAKE_INSTALL_LIBDIR`` from ``CMAKE_INSTALL_PREFIX/lib``.
   ``CMAKE_INSTALL_LIBDIR`` is automatically set, but is available as a CMake option that can modified.


## Changes to SUNDIALS in release 3.2.0

Fixed problem with index types which would occur with some compilers (e.g.
armclang) that did not define `__STDC_VERSION__`. The fix includes a
depcrecation of the current behavior of the `SUNDIALS_INDEX_TYPE` CMake option.

Fixed a thread-safety issue in CVODES and IDAS when using adjoint sensitivity
analysis.

Added hybrid MPI/CUDA and MPI/RAJA vectors to allow use of more than one MPI
rank when using a GPU system. The vectors assume one GPU device per MPI rank.

Changed the name of the RAJA nvector library to `libsundials_nveccudaraja.lib`
from `libsundials_nvecraja.lib` to better reflect that we only support CUDA as a
backend for RAJA currently.

Increased CMake minimum version to 3.1.3

Add constraint handling feature to CVODE and CVODES.

Fixed a bug in IDAS where the saved residual value used in the nonlinear solve
for consistent initial conditions was passed as temporary workspace and could be
overwritten.

Several changes were made to the build system. If MPI is enabled and MPI
compiler wrappers are not set, the build system will check if
`CMAKE_<language>_COMPILER` can compile MPI programs before trying to locate and
use an MPI installation. The native CMake FindMPI module is now used to locate
an MPI installation. The options for setting MPI compiler wrappers and the
executable for running MPI programs have been updated to align with those in
native CMake FindMPI module. This included changing `MPI_MPICC` to
`MPI_C_COMPILER`, `MPI_MPICXX` to `MPI_CXX_COMPILER` combining `MPI_MPIF77` and
`MPI_MPIF90` to `MPI_Fortran_COMPILER`, and changing `MPI_RUN_COMMAND` to
`MPIEXEC_EXECUTABLE`. When a Fortran name-mangling scheme is needed (e.g.,
`LAPACK_ENABLE` is `ON`) the build system will infer the scheme from the Fortran
compiler. If a Fortran compiler is not available or the inferred or default
scheme needs to be overridden, the advanced options `SUNDIALS_F77_FUNC_CASE` and
`SUNDIALS_F77_FUNC_UNDERSCORES` can be used to manually set the name-mangling
scheme and bypass trying to infer the scheme. Additionally, parts of the main
`CMakeLists.txt` file were moved to new files in the src and example directories
to make the CMake configuration file structure more modular.

### ARKODE Changes in v2.2.0

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

* When a Fortran name-mangling scheme is needed (e.g., ``ENABLE_LAPACK``
  is ``ON``) the build system will infer the scheme from the Fortran
  compiler. If a Fortran compiler is not available or the inferred or default
  scheme needs to be overridden, the advanced options
  ``SUNDIALS_F77_FUNC_CASE`` and ``SUNDIALS_F77_FUNC_UNDERSCORES`` can
  be used to manually set the name-mangling scheme and bypass trying to infer
  the scheme.

* Parts of the main CMakeLists.txt file were moved to new files in the
  ``src`` and ``example`` directories to make the CMake configuration file
  structure more modular.


### CVODE Changes in v3.2.0

Support for optional inequality constraints on individual components of the
solution vector has been added to CVODE and CVODES. See
:numref:`CVODE.Mathematics` and the description of in :numref:`CVODE.Usage.CC.optional_input` for
more details. Use of :c:func:`CVodeSetConstraints` requires the ``N_Vector``
operations :c:func:`N_VMinQuotient`, :c:func:`N_VConstMask`, and
:c:func:`N_VCompare` that were not previously required by CVODE and CVODES.

Fixed a problem with setting which would occur with some compilers (e.g.
armclang) that did not define ``__STDC_VERSION__``.

Added hybrid MPI/CUDA and MPI/RAJA vectors to allow use of more than one MPI
rank when using a GPU system. The vectors assume one GPU device per MPI rank.

Changed the name of the RAJA ``N_Vector`` library to from to better reflect that
we only support CUDA as a backend for RAJA currently.

Several changes were made to the build system:

   - CMake 3.1.3 is now the minimum required CMake version.
   - Deprecate the behavior of the CMake option and added the CMake option to select the integer size.
   - The native CMake FindMPI module is now used to locate an MPI installation.
   - If MPI is enabled and MPI compiler wrappers are not set, the build system will
     check if can compile MPI programs before trying to locate and use an MPI
     installation.
   - The previous options for setting MPI compiler wrappers and the executable for
     running MPI programs have been have been depreated. The new options that align
     with those used in native CMake FindMPI module are ``MPI_C_COMPILER``,
     ``MPO_CXX_COMPILER``, ``MPI_Fortran_COMPILER``, and ``MPIEXEC_EXECUTABLE``.
   - When a Fortran name-mangling scheme is needed (e.g., is ) the build system will
     infer the scheme from the Fortran compiler. If a Fortran compiler is not
     available or the inferred or default scheme needs to be overridden, the advanced
     options and can be used to manually set the name-mangling scheme and bypass
     trying to infer the scheme.
   - Parts of the main CMakeLists.txt file were moved to new files in the and
     directories to make the CMake configuration file structure more modular.

### CVODES Changes in v3.2.0

Support for optional inequality constraints on individual components of the
solution vector has been added to CVODE and CVODES. See Chapter
:numref:`CVODES.Mathematics` and the description of
:c:func:`CVodeSetConstraints` for more details. Use of ``CVodeSetConstraints``
requires the ``N_Vector`` operations ``N_MinQuotient``, ``N_VConstrMask``, and
``N_VCompare`` that were not previously required by CVODE and CVODES.

Fixed a thread-safety issue when using ajdoint sensitivity analysis.

Fixed a problem with setting ``sunindextype`` which would occur with some
compilers (e.g. armclang) that did not define ``__STDC_VERSION__``.

Added hybrid MPI/CUDA and MPI/RAJA vectors to allow use of more than one MPI
rank when using a GPU system. The vectors assume one GPU device per MPI rank.

Changed the name of the RAJA ``N_Vector`` library to
``libsundials_nveccudaraja.lib`` from ``libsundials_nvecraja.lib`` to better
reflect that we only support CUDA as a backend for RAJA currently.


Several changes were made to the build system:

-  CMake 3.1.3 is now the minimum required CMake version.

-  Deprecate the behavior of the ``SUNDIALS_INDEX_TYPE`` CMake option and added the ``SUNDIALS_INDEX_SIZE`` CMake option
   to select the ``sunindextype`` integer size.

-  The native CMake FindMPI module is now used to locate an MPI installation.

-  If MPI is enabled and MPI compiler wrappers are not set, the build system will check if ``CMAKE_<language>_COMPILER``
   can compile MPI programs before trying to locate and use an MPI installation.

-  The previous options for setting MPI compiler wrappers and the executable for running MPI programs have been have
   been depreated. The new options that align with those used in native CMake FindMPI module are ``MPI_C_COMPILER``,
   ``MPI_CXX_COMPILER``, ``MPI_Fortran_COMPILER``, and ``MPIEXEC_EXECUTABLE``.

-  When a Fortran name-mangling scheme is needed (e.g., ``ENABLE_LAPACK`` is ``ON``) the build system will infer the
   scheme from the Fortran compiler. If a Fortran compiler is not available or the inferred or default scheme needs to
   be overridden, the advanced options ``SUNDIALS_F77_FUNC_CASE`` and ``SUNDIALS_F77_FUNC_UNDERSCORES`` can be used to
   manually set the name-mangling scheme and bypass trying to infer the scheme.

-  Parts of the main CMakeLists.txt file were moved to new files in the ``src`` and ``example`` directories to make the
   CMake configuration file structure more modular.

### IDA Changes in v3.2.0

Fixed a problem with setting ``sunindextype`` which would occur with some
compilers (e.g. armclang) that did not define ``__STDC_VERSION__``.

Added hybrid MPI/CUDA and MPI/RAJA vectors to allow use of more than one MPI
rank when using a GPU system. The vectors assume one GPU device per MPI rank.

Changed the name of the RAJA ``N_Vector`` library to
``libsundials_nveccudaraja.lib`` from ``libsundials_nvecraja.lib`` to better
reflect that we only support CUDA as a backend for RAJA currently.

Several changes were made to the build system:

* CMake 3.1.3 is now the minimum required CMake version.

* Deprecate the behavior of the :cmakeop:`SUNDIALS_INDEX_TYPE` CMake option and
  added the :cmakeop:`SUNDIALS_INDEX_SIZE` CMake option to select the
  ``sunindextype`` integer size.

* The native CMake FindMPI module is now used to locate an MPI installation.

* If MPI is enabled and MPI compiler wrappers are not set, the build system will
  check if ``CMAKE_<language>_COMPILER`` can compile MPI programs before trying
  to locate and use an MPI installation.

* The previous options for setting MPI compiler wrappers and the executable for
  running MPI programs have been have been depreated. The new options that align
  with those used in native CMake FindMPI module are :cmakeop:`MPI_C_COMPILER`,
  :cmakeop:`MPI_CXX_COMPILER`, :cmakeop:`MPI_Fortran_COMPILER`, and
  :cmakeop:`MPIEXEC_EXECUTABLE`.

* When a Fortran name-mangling scheme is needed (e.g., :cmakeop:`ENABLE_LAPACK`
  is ``ON``) the build system will infer the scheme from the Fortran compiler.
  If a Fortran compiler is not available or the inferred or default scheme needs
  to be overridden, the advanced options :cmakeop:`SUNDIALS_F77_FUNC_CASE` and
  :cmakeop:`SUNDIALS_F77_FUNC_UNDERSCORES` can be used to manually set the
  name-mangling scheme and bypass trying to infer the scheme.

* Parts of the main CMakeLists.txt file were moved to new files in the ``src``
  and ``example`` directories to make the CMake configuration file structure
  more modular.

### IDAS Changes in v2.2.0

Fixed a problem with setting ``sunindextype`` which would occur with some
compilers (e.g. armclang) that did not define ``__STDC_VERSION__``.

Added hybrid MPI/CUDA and MPI/RAJA vectors to allow use of more than one MPI
rank when using a GPU system. The vectors assume one GPU device per MPI rank.

Changed the name of the RAJA ``N_Vector`` library to
``libsundials_nveccudaraja.lib`` from ``libsundials_nvecraja.lib`` to better
reflect that we only support CUDA as a backend for RAJA currently.

Several changes were made to the build system:

* CMake 3.1.3 is now the minimum required CMake version.

* Deprecate the behavior of the :cmakeop:`SUNDIALS_INDEX_TYPE` CMake option and
  added the :cmakeop:`SUNDIALS_INDEX_SIZE` CMake option to select the
  ``sunindextype`` integer size.

* The native CMake FindMPI module is now used to locate an MPI installation.

* If MPI is enabled and MPI compiler wrappers are not set, the build system will
  check if ``CMAKE_<language>_COMPILER`` can compile MPI programs before trying
  to locate and use an MPI installation.

* The previous options for setting MPI compiler wrappers and the executable for
  running MPI programs have been have been depreated. The new options that align
  with those used in native CMake FindMPI module are :cmakeop:`MPI_C_COMPILER`,
  :cmakeop:`MPI_CXX_COMPILER`, :cmakeop:`MPI_Fortran_COMPILER`, and
  :cmakeop:`MPIEXEC_EXECUTABLE`.

* When a Fortran name-mangling scheme is needed (e.g., :cmakeop:`ENABLE_LAPACK`
  is ``ON``) the build system will infer the scheme from the Fortran compiler.
  If a Fortran compiler is not available or the inferred or default scheme needs
  to be overridden, the advanced options :cmakeop:`SUNDIALS_F77_FUNC_CASE` and
  :cmakeop:`SUNDIALS_F77_FUNC_UNDERSCORES` can be used to manually set the
  name-mangling scheme and bypass trying to infer the scheme.

* Parts of the main CMakeLists.txt file were moved to new files in the ``src``
  and ``example`` directories to make the CMake configuration file structure
  more modular.

### KINSOL Changes in v3.2.0

Fixed a problem with setting ``sunindextype`` which would occur with some compilers (e.g. armclang) that
did not define ``__STDC_VERSION__``. Added hybrid MPI/CUDA and MPI/RAJA vectors to allow use of more than
one MPI rank when using a GPU system. The vectors assume one GPU device per MPI rank.  Changed the name
of the RAJA ``N_Vector`` library to ``libsundials_nveccudaraja.lib`` from ``libsundials_nvecraja.lib`` to
better reflect that we only support CUDA as a backend for RAJA currently. Several changes were made to
the build system:

-  CMake 3.1.3 is now the minimum required CMake version.

-  Deprecate the behavior of the ``SUNDIALS_INDEX_TYPE`` CMake option and added the ``SUNDIALS_INDEX_SIZE`` CMake option
   to select the ``sunindextype`` integer size.

-  The native CMake FindMPI module is now used to locate an MPI installation.

-  If MPI is enabled and MPI compiler wrappers are not set, the build system will check if ``CMAKE_<language>_COMPILER``
   can compile MPI programs before trying to locate and use an MPI installation.

-  The previous options for setting MPI compiler wrappers and the executable for running MPI programs have been have
   been depreated. The new options that align with those used in native CMake FindMPI module are ``MPI_C_COMPILER``,
   ``MPI_CXX_COMPILER``, ``MPI_Fortran_COMPILER``, and ``MPIEXEC_EXECUTABLE``.

-  When a Fortran name-mangling scheme is needed (e.g., ``ENABLE_LAPACK`` is ``ON``) the build system will infer the
   scheme from the Fortran compiler. If a Fortran compiler is not available or the inferred or default scheme needs to
   be overridden, the advanced options ``SUNDIALS_F77_FUNC_CASE`` and ``SUNDIALS_F77_FUNC_UNDERSCORES`` can be used to
   manually set the name-mangling scheme and bypass trying to infer the scheme.

-  Parts of the main CMakeLists.txt file were moved to new files in the ``src`` and ``example`` directories to make the
   CMake configuration file structure more modular.


## Changes to SUNDIALS in release 3.1.2

Fixed Windows specific problem where `sunindextype` was not correctly defined
when using 64-bit integers. On Windows `sunindextype` is now defined as the MSVC
basic type `__int64`.

Changed LICENSE install path to `instdir/include/sundials`.

Updated the minimum required version of CMake to 2.8.12 and enabled using rpath
by default to locate shared libraries on OSX.

The misnamed function `CVSpilsSetJacTimesSetupFnBS` in cvodes has been
deprecated and replaced by `CVSpilsSetJacTimesBS`. The deprecated function
`CVSpilsSetJacTimesSetupFnBS` will be removed in the next major release.

Added and updated usage-notes examples from the SUNDIALS website to work with
SUNDIALS 3.x. The new examples are `cvode/cvDisc_dns.c`,
`cvode/cvRoberts_dns_negsol.c`, and `cvodes/cvsRoberts_FSA_dns_Switch.c`.

Added sparse SUNMatrix "Reallocate" routine to allow specification of the
nonzero storage.

Updated the KLU SUNLinearSolver module to set constants for the two
reinitialization types, and fixed a bug in the full reinitialization approach
where the sparse SUNMatrix pointer would go out of scope on some architectures.

Updated the "ScaleAdd" and "ScaleAddI" implementations in the sparse SUNMatrix
module to more optimally handle the case where the target matrix contained
sufficient storage for the sum, but had the wrong sparsity pattern. The sum now
occurs in-place, by performing the sum backwards in the existing storage.
However, it is still more efficient if the user-supplied Jacobian routine
allocates storage for the sum I + gamma J or M + gamma J manually (with zero
entries if needed).

### ARKODE Changes in v2.1.2

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


### CVODE Changes in v3.1.2

The changes in this minor release include the following:

-  Updated the minimum required version of CMake to 2.8.12 and enabled
   using rpath by default to locate shared libraries on OSX.
-  Fixed Windows specific problem where was not correctly defined when
   using 64-bit integers for the SUNDIALS index type. On Windows ``sunindextype`` is
   now defined as the MSVC basic type ``__int64``.
-  Added sparse SUNMatrix "Reallocate" routine to allow specification of
   the nonzero storage.
-  Updated the KLU SUNLinearSolver module to set constants for the two
   reinitialization types, and fixed a bug in the full reinitialization
   approach where the sparse SUNMatrix pointer would go out of scope on
   some architectures.
-  Updated the "ScaleAdd" and "ScaleAddI" implementations in the sparse
   SUNMatrix module to more optimally handle the case where the target
   matrix contained sufficient storage for the sum, but had the wrong
   sparsity pattern. The sum now occurs in-place, by performing the sum
   backwards in the existing storage. However, it is still more
   efficient if the user-supplied Jacobian routine allocates storage for
   the sum :math:`I+\gamma J` manually (with zero entries if needed).
-  Added the following examples from the usage notes page of the
   SUNDIALS website, and updated them to work with SUNDIALS 3.x:

   -  ``cvDisc_dns.c``, which demonstrates using CVODE with discontinuous solutions or RHS.
   -  ``cvRoberts_dns_negsol.c``, which illustrates the use of the RHS function return value to
      control unphysical negative concentrations.

-  Changed the LICENSE install path to `instdir/icnlude/sundials`.

### CVODES Changes in v3.1.2

The changes in this minor release include the following:

-  Updated the minimum required version of CMake to 2.8.12 and enabled using rpath by default to locate shared libraries
   on OSX.

-  Fixed Windows specific problem where ``sunindextype`` was not correctly defined when using 64-bit integers for the
   SUNDIALS index type. On Windows ``sunindextype`` is now defined as the MSVC basic type ``__int64``.

-  Added sparse SUNMatrix "Reallocate" routine to allow specification of the nonzero storage.

-  Updated the KLU SUNLinearSolver module to set constants for the two reinitialization types, and fixed a bug in the
   full reinitialization approach where the sparse SUNMatrix pointer would go out of scope on some architectures.

-  Updated the "ScaleAdd" and "ScaleAddI" implementations in the sparse SUNMatrix module to more optimally handle the
   case where the target matrix contained sufficient storage for the sum, but had the wrong sparsity pattern. The sum
   now occurs in-place, by performing the sum backwards in the existing storage. However, it is still more efficient if
   the user-supplied Jacobian routine allocates storage for the sum :math:`I+\gamma J` manually (with zero entries if
   needed).

-  Added new example, ``cvRoberts_FSA_dns_Switch.c``, which demonstrates switching on/off forward sensitivity
   computations. This example came from the usage notes page of the SUNDIALS website.

-  The misnamed function ``CVSpilsSetJacTimesSetupFnBS`` has been deprecated and replaced by ``CVSpilsSetJacTimesBS``.
   The deprecated function ``CVSpilsSetJacTimesSetupFnBS`` will be removed in the next major release.

-  Changed the LICENSE install path to ``instdir/include/sundials``.

### IDA Changes in v3.1.2

The changes in this minor release include the following:

* Updated the minimum required version of CMake to 2.8.12 and enabled using
  rpath by default to locate shared libraries on OSX.

* Fixed Windows specific problem where ``sunindextype`` was not correctly
  defined when using 64-bit integers for the SUNDIALS index type. On Windows
  ``sunindextype`` is now defined as the MSVC basic type ``__int64``.

* Added sparse SUNMatrix "Reallocate" routine to allow specification of the
  nonzero storage.

* Updated the KLU SUNLinearSolver module to set constants for the two
  reinitialization types, and fixed a bug in the full reinitialization approach
  where the sparse SUNMatrix pointer would go out of scope on some
  architectures.

* Updated the :c:func:`SUNMatScaleAdd` and :c:func:`SUNMatScaleAddI`
  implementations in the sparse SUNMatrix module to more optimally handle the
  case where the target matrix contained sufficient storage for the sum, but had
  the wrong sparsity pattern. The sum now occurs in-place, by performing the sum
  backwards in the existing storage. However, it is still more efficient if the
  user-supplied Jacobian routine allocates storage for the sum
  :math:`I+\gamma J` manually (with zero entries if needed).

* Changed the LICENSE install path to ``instdir/include/sundials``.

### IDAS Changes in v2.1.2

The changes in this minor release include the following:

* Updated the minimum required version of CMake to 2.8.12 and enabled using
  rpath by default to locate shared libraries on OSX.

* Fixed Windows specific problem where ``sunindextype`` was not correctly
  defined when using 64-bit integers for the SUNDIALS index type. On Windows
  ``sunindextype`` is now defined as the MSVC basic type ``__int64``.

* Added sparse SUNMatrix "Reallocate" routine to allow specification of the
  nonzero storage.

* Updated the KLU SUNLinearSolver module to set constants for the two
  reinitialization types, and fixed a bug in the full reinitialization approach
  where the sparse SUNMatrix pointer would go out of scope on some
  architectures.

* Updated the :c:func:`SUNMatScaleAdd` and :c:func:`SUNMatScaleAddI`
  implementations in the sparse SUNMatrix module to more optimally handle the
  case where the target matrix contained sufficient storage for the sum, but had
  the wrong sparsity pattern. The sum now occurs in-place, by performing the sum
  backwards in the existing storage. However, it is still more efficient if the
  user-supplied Jacobian routine allocates storage for the sum
  :math:`I+\gamma J` manually (with zero entries if needed).

* Changed the LICENSE install path to ``instdir/include/sundials``.

### KINSOL Changes in v3.1.2

The changes in this minor release include the following:

-  Updated the minimum required version of CMake to 2.8.12 and enabled using rpath by default to locate shared libraries
   on OSX.

-  Fixed Windows specific problem where ``sunindextype`` was not correctly defined when using 64-bit integers for the
   SUNDIALS index type. On Windows ``sunindextype`` is now defined as the MSVC basic type ``__int64``.

-  Added sparse SUNMatrix "Reallocate" routine to allow specification of the nonzero storage.

-  Updated the KLU SUNLinearSolver module to set constants for the two reinitialization types, and fixed a bug in the
   full reinitialization approach where the sparse SUNMatrix pointer would go out of scope on some architectures.

-  Updated the "ScaleAdd" and "ScaleAddI" implementations in the sparse SUNMatrix module to more optimally handle the
   case where the target matrix contained sufficient storage for the sum, but had the wrong sparsity pattern. The sum
   now occurs in-place, by performing the sum backwards in the existing storage. However, it is still more efficient if
   the user-supplied Jacobian routine allocates storage for the sum :math:`I+\gamma J` manually (with zero entries if
   needed).

-  Changed the LICENSE install path to ``instdir/include/sundials``.


## Changes to SUNDIALS in release 3.1.1

Fixed a minor bug in the CVODE and CVODES `cvSLdet` routine, where a return was
missing in the error check for three inconsistent roots.

Fixed a potential memory leak in the SPGMR and SPFGMR linear solvers: if
"Initialize" was called multiple times then the solver memory was reallocated
(without being freed).

Fixed a minor bug in the `ARKReInit` routine, where a flag was incorrectly set
to indicate that the problem had been resized (instead of just re-initialized).

Fixed C++11 compiler errors/warnings about incompatible use of string literals.

Updated KLU SUNLinearSolver module to use a typedef for the precision-specific
solve function to be used (to avoid compiler warnings).

Added missing typecasts for some (`void*`) pointers to avoid compiler warnings.

Bugfix in `sunmatrix_sparse.c` where `int` was used instead of `sunindextype` in
one location.

Fixed a minor bug in `KINPrintInfo` where a case was missing for
`KIN_REPTD_SYSFUNC_ERR` leading to an undefined info message.

Added missing `#include <stdio.h>` in NVECTOR and SUNMATRIX header files.

Added missing prototypes for `ARKSpilsGetNumMTSetups` in ARKODE and
`IDASpilsGetNumJTSetupEvals` in IDA and IDAS.

Fixed an indexing bug in the CUDA NVECTOR implementation of `N_VWrmsNormMask`
and revised the RAJA NVECTOR implementation of `N_VWrmsNormMask` to work with
mask arrays using values other than zero or one. Replaced `double` with
`realtype` in the RAJA vector test functions.

Fixed compilation issue with GCC 7.3.0 and Fortran programs that do not require
a SUNMatrix or SUNLinearSolver module (e.g. iterative linear solvers, explicit
methods in ARKODE, functional iteration in CVODE, etc.).

### ARKODE Changes in v2.1.1

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


### CVODE Changes in v3.1.1

The changes in this minor release include the following:

-  Fixed a minor bug in the cvSLdet routine, where a return was missing
   in the error check for three inconsistent roots.
-  Fixed a potential memory leak in the SPGMR and SPFGMR linear
   solvers: if "Initialize" was called multiple times then the solver
   memory was reallocated (without being freed).
-  Updated KLU ``SUNLinearSolver`` module to use a for the precision-specific
   solve function to be used (to avoid compiler warnings).
-  Added missing typecasts for some pointers (again, to avoid compiler
   warnings).
-  Bugfix in ``sunmatric_sparse.c`` where we had used instead of in one location.
-  Added missing ``#include <stio.h>`` in ``N_Vector`` and ``SUNMatrix`` header files.
-  Fixed an indexing bug in the CUDA ``N_Vector`` implementation of
   and revised the RAJA ``N_Vector`` implementation of :c:func:`N_VWrmsNormMask` to work with
   mask arrays using values other than zero or one. Replaced ``double`` with ``realtype`` in the
   RAJA vector test functions.
-  Fixed compilation issue with GCC 7.3.0 and Fortran programs that do
   not require a ``SUNMatrix`` or ``SUNLinearSolver`` module (e.g., iterative
   linear solvers or fixed-point iteration).

In addition to the changes above, minor corrections were also made to
the example programs, build system, and user documentation.

### CVODES Changes in v3.1.1

The changes in this minor release include the following:

-  Fixed a minor bug in the cvSLdet routine, where a return was missing in the error check for three inconsistent roots.

-  Fixed a potential memory leak in the SPGMR and SPFGMR linear solvers: if "Initialize" was called multiple
   times then the solver memory was reallocated (without being freed).

-  Updated KLU ``SUNLinearSolver`` module to use a ``typedef`` for the precision-specific solve function to be used (to
   avoid compiler warnings).

-  Added missing typecasts for some ``(void*)`` pointers (again, to avoid compiler warnings).

-  Bugfix in ``sunmatrix_sparse.c`` where we had used ``int`` instead of ``sunindextype`` in one location.

-  Added missing ``#include <stdio.h>`` in ``N_Vector`` and ``SUNMatrix`` header files.

-  Fixed an indexing bug in the CUDA ``N_Vector`` implementation of ``N_VWrmsNormMask`` and revised the
   RAJA ``N_Vector`` implementation of ``N_VWrmsNormMask`` to work with mask arrays using values other than zero
   or one. Replaced ``double`` with ``realtype`` in the RAJA vector test functions.

In addition to the changes above, minor corrections were also made to the
example programs, build system, and user documentation.

### IDA Changes in v3.1.1

The changes in this minor release include the following:

* Fixed a potential memory leak in the :ref:`SUNLINSOL_SPGMR <SUNLinSol.SPGMR>`
  and :ref:`SUNLINSOL_SPFGMR <SUNLinSol.SPFGMR>` linear solvers: if
  "Initialize" was called multiple times then the solver memory was reallocated
  (without being freed).

* Updated KLU ``SUNLinearSolver`` module to use a ``typedef`` for the
  precision-specific solve function to be used (to avoid compiler warnings).

* Added missing typecasts for some ``(void*)`` pointers (again, to avoid
  compiler warnings).

* Bugfix in ``sunmatrix_sparse.c`` where we had used ``int`` instead of
  ``sunindextype`` in one location.

* Added missing ``#include <stdio.h>`` in ``N_Vector`` and ``SUNMatrix`` header
  files.

* Added missing prototype for :c:func:`IDASpilsGetNumJTSetupEvals`.

* Fixed an indexing bug in the CUDA ``N_Vector`` implementation of
  :c:func:`N_VWrmsNormMask` and revised the RAJA ``N_Vector`` implementation of
  :c:func:`N_VWrmsNormMask` to work with mask arrays using values other than
  zero or one. Replaced ``double`` with ``realtype`` in the RAJA vector test
  functions.

* Fixed compilation issue with GCC 7.3.0 and Fortran programs that do not
  require a ``SUNMatrix`` module (e.g., iterative linear solvers).

In addition to the changes above, minor corrections were also made to the
example programs, build system, and user documentation.

### IDAS Changes in v2.1.1

The changes in this minor release include the following:

* Fixed a potential memory leak in the :ref:`SUNLINSOL_SPGMR <SUNLinSol.SPGMR>`
  and :ref:`SUNLINSOL_SPFGMR <SUNLinSol.SPFGMR>` linear solvers: if
  "Initialize" was called multiple times then the solver memory was reallocated
  (without being freed).

* Updated KLU ``SUNLinearSolver`` module to use a ``typedef`` for the
  precision-specific solve function to be used (to avoid compiler warnings).

* Added missing typecasts for some ``(void*)`` pointers (again, to avoid
  compiler warnings).

* Bugfix in ``sunmatrix_sparse.c`` where we had used ``int`` instead of
  ``sunindextype`` in one location.

* Added missing ``#include <stdio.h>`` in ``N_Vector`` and ``SUNMatrix`` header
  files.

* Added missing prototype for :c:func:`IDASpilsGetNumJTSetupEvals`.

* Fixed an indexing bug in the CUDA ``N_Vector`` implementation of
  :c:func:`N_VWrmsNormMask` and revised the RAJA ``N_Vector`` implementation of
  :c:func:`N_VWrmsNormMask` to work with mask arrays using values other than
  zero or one. Replaced ``double`` with ``realtype`` in the RAJA vector test
  functions.

* Fixed compilation issue with GCC 7.3.0 and Fortran programs that do not
  require a ``SUNMatrix`` module (e.g., iterative linear solvers).

In addition to the changes above, minor corrections were also made to the
example programs, build system, and user documentation.

### KINSOL Changes in v3.1.1

The changes in this minor release include the following:

-  Fixed a potential memory leak in the SPGMR and SPFGMR linear solvers: if "Initialize" was called multiple
   times then the solver memory was reallocated (without being freed).

-  Updated KLU SUNLinearSolver module to use a ``typedef`` for the precision-specific solve function to be used (to
   avoid compiler warnings).

-  Added missing typecasts for some ``(void*)`` pointers (again, to avoid compiler warnings).

-  Bugfix in ``sunmatrix_sparse.c`` where we had used ``int`` instead of ``sunindextype`` in one location.

-  Fixed a minor bug in ``KINPrintInfo`` where a case was missing for ``KIN_REPTD_SYSFUNC_ERR`` leading to an undefined
   info message.

-  Added missing ``#include <stdio.h>`` in ``N_Vector`` and ``SUNMatrix`` header files.

-  Fixed an indexing bug in the CUDA ``N_Vector`` implementation of ``N_VWrmsNormMask`` and revised the
   RAJA ``N_Vector`` implementation of ``N_VWrmsNormMask`` to work with mask arrays using values other than zero
   or one. Replaced ``double`` with ``realtype`` in the RAJA vector test functions.

-  Fixed compilation issue with GCC 7.3.0 and Fortran programs that do not require a ``SUNMatrix`` or ``SUNLinearSolver``
   module (e.g., iterative linear solvers or fixed pointer solver).

In addition to the changes above, minor corrections were also made to the example programs, build system, and user
documentation.


## Changes to SUNDIALS in release 3.1.0

Added NVECTOR print functions that write vector data to a specified file (e.g.,
`N_VPrintFile_Serial`).

Added `make test` and `make test_install` options to the build system for
testing SUNDIALS after building with `make` and installing with `make install`
respectively.

Added "Changes in ..." (latest version) to all User Guides.

### ARKODE Changes in v2.1.0

Added NVECTOR print functions that write vector data to a specified
file (e.g. ``N_VPrintFile_Serial``).

Added ``make test`` and ``make test_install`` options to the build
system for testing SUNDIALS after building with ``make`` and
installing with ``make install`` respectively.

### CVODE Changes in v3.1.0

Added ``N_Vector`` print functions that write vector data to a specified
file (e.g., :c:func:`N_VPrintFile_Serial`).

Added ``make test`` and ``make test_install`` options to the build system for testing SUNDIALS after
building with and installing with respectively.

### CVODES Changes in v3.1.0

Added ``N_Vector`` print functions that write vector data to a specified file
(e.g., ``N_VPrintFile_Serial``).

Added ``make test`` and ``make test_install`` options to the build system for
testing SUNDIALS after building with ``make`` and installing with ``make
install`` respectively.

### IDA Changes in v3.1.0

Added ``N_Vector`` print functions that write vector data to a specified file
(e.g., :c:func:`N_VPrintFile_Serial`).

Added ``make test`` and ``make test_install`` options to the build system for
testing SUNDIALS after building with ``make`` and installing with ``make
install`` respectively.

### IDA Changes in v2.1.0

Added ``N_Vector`` print functions that write vector data to a specified file
(e.g., :c:func:`N_VPrintFile_Serial`).

Added ``make test`` and ``make test_install`` options to the build system for
testing SUNDIALS after building with ``make`` and installing with ``make
install`` respectively.

### KINSOL Changes in v3.1.0

Added ``N_Vector`` print functions that write vector data to a specified file (e.g., ``N_VPrintFile_Serial``).

Added ``make test`` and ``make test_install`` options to the build system for testing SUNDIALS after building with
``make`` and installing with ``make install`` respectively.


## Changes to SUNDIALS in release 3.0.0

Added new linear solver and matrix interfaces for all SUNDIALS packages and
updated the existing linear solver and matrix modules. The goal of the redesign
is to provide greater encapsulation and ease interfacing custom linear solvers
with linear solver libraries. Specific changes include:

 * Added generic SUNMATRIX module with three provided implementations:
   dense, banded and sparse.  These replicate previous SUNDIALS Dls and
   Sls matrix structures in a single object-oriented API.

 * Added example problems demonstrating use of generic SUNMATRIX modules.

 * Added generic SUNLINEARSOLVER module with eleven provided
   implementations: dense, banded, LAPACK dense, LAPACK band, KLU,
   SuperLU_MT, SPGMR, SPBCGS, SPTFQMR, SPFGMR, PCG.  These replicate
   previous SUNDIALS generic linear solvers in a single object-oriented
   API.

 * Added example problems demonstrating use of generic SUNLINEARSOLVER
   modules.

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

 * Added Spils interface routines to ARKODE, CVODE, CVODES, IDA and
   IDAS to allow specification of a user-provided "JTSetup" routine.
   This change supports users who wish to set up data structures for
   the user-provided Jacobian-times-vector ("JTimes") routine, and
   where the cost of one JTSetup setup per Newton iteration can be
   amortized between multiple JTimes calls.

Corresponding updates were made to all the example programs.

Two new NVECTOR modules added: for CUDA and RAJA support for GPU systems
(Information on RAJA: <https://software.llnl.gov/RAJA/> )
These vectors are supplied to provide very basic support for running
on GPU architectures.  Users are advised that these vectors both move all data
to the GPU device upon construction, and speedup will only be realized if the
user also conducts the right-hand-side function evaluation on the device.
In addition, these vectors assume the problem fits on one GPU.
For further information about RAJA, users are referred to the web site,
<https://software.llnl.gov/RAJA/.>

Addition of sunindextype option for 32-bit or 64-bit integer data index types
within all SUNDIALS structures

  * sunindextype is defined to be int32_t or int64_t when portable types are
    supported, otherwise it is defined as int or long int.

  * The Fortran interfaces continue to use `long int` for indices, except for
    their sparse matrix interface that now uses the new sunindextype.

  * Includes interfaces to PETSc, hypre, SuperLU_MT, and KLU with either 32-bit
    or 64-bit capabilities depending how the user configures SUNDIALS.

To avoid potential namespace conflicts, the macros defining booleantype
values TRUE and FALSE have been changed to SUNTRUE and SUNFALSE respectively.

Temporary vectors were removed from preconditioner setup and solve
routines for all packages.  It is assumed that all necessary data
for user-provided preconditioner operations will be allocated and
stored in user-provided data structures.

The file include/sundials\_fconfig.h was added.  This file contains
SUNDIALS type information for use in Fortran programs.

Added support for many xSDK-compliant build system keys
(Information on xSDK compliance: <https://xsdk.info/policies/> )
The xSDK is a movement in scientific software to provide a foundation for the
rapid and efficient production of high-quality,
sustainable extreme-scale scientific applications.  More information can
be found at <https://xsdk.info.>

Added functions SUNDIALSGetVersion and SUNDIALSGetVersionNumber to
get SUNDIALS release version information at runtime.

### Build System

Renamed CMake options to enable/disable examples for greater clarity
and added option to enable/disable Fortran 77 examples:

  * Changed `EXAMPLES_ENABLE` to `EXAMPLES_ENABLE_C`
  * Changed `CXX_ENABLE` to `EXAMPLES_ENABLE_CXX`
  * Changed `F90_ENABLE` to `EXAMPLES_ENABLE_F90`
  * Added `EXAMPLES_ENABLE_F77` option

Added separate `BLAS_ENABLE` and `BLAS_LIBRARIES` CMake variables

Fixed minor CMake bugs and included additional error checking during CMake
configuration

Corrections and additions to all User Guides.

Added "Changes in ..." (latest version) section to the introduction to in all
User Guides.

### ARKODE

Added comments to `arkode_butcher.c` regarding which methods should have
coefficients accurate enough for use in quad precision.

Fixed `RCONST` usage in `arkode_butcher.c`.

Fixed bug in `arkInitialSetup` to ensure the mass matrix vector product is
set up before the "msetup" routine is called.

Fixed ARKODE printf-related compiler warnings when building SUNDIALS
with extended precision.

### CVODE and CVODES

In `CVodeFree`, now call `lfree` unconditionally (if non-NULL).

### IDA and IDAS

Added missing prototype for `IDASetMaxBacksIC` in `ida.h` and `idas.h`.

### KINSOL

Corrected KINSOL fcmix name translation for `FKIN_SPFGMR`.

Renamed `KINLocalFn` and `KINCommFn` to `KINBBDLocalFn` and `KINBBDCommFn`
respectively in the BBD preconditioner module for consistency with other
SUNDIALS solvers.

### ARKODE Changes in v2.0.0

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

* Added Spils interface routines to ARKODE, CVODE, CVODES, IDA and
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


### CVODE Changes in v3.0.0

All interfaces to matrix structures and linear solvers have been reworked, and
all example programs have been updated. The goal of the redesign of these
interfaces was to provide more encapsulation and ease in interfacing custom
linear solvers and interoperability with linear solver libraries. Specific
changes include:

-  Added generic SUNMATRIX module with three provided implementations:
   dense, banded and sparse. These replicate previous SUNDIALS Dls and Sls
   matrix structures in a single object-oriented API.

-  Added example problems demonstrating use of generic SUNMATRIX
   modules.

-  Added generic SUNLINEARSOLVER module with eleven provided
   implementations: dense, banded, LAPACK dense, LAPACK band, KLU, SuperLU_MT,
   SPGMR, SPBCGS, SPTFQMR, SPFGMR, PCG. These replicate previous SUNDIALS
   generic linear solvers in a single object-oriented API.

-  Added example problems demonstrating use of generic SUNLINEARSOLVER
   modules.

-  Expanded package-provided direct linear solver (Dls) interfaces and
   scaled, preconditioned, iterative linear solver (Spils) interfaces to utilize
   generic SUNMATRIX and SUNLINEARSOLVER objects.

-  Removed package-specific, linear solver-specific, solver modules
   (e.g. CVDENSE, KINBAND, IDAKLU, ARKSPGMR) since their functionality is
   entirely replicated by the generic Dls/Spils interfaces and
   SUNLINEARSOLVER/SUNMATRIX modules. The exception is CVDIAG, a diagonal
   approximate Jacobian solver available to CVODE and CVODES.

-  Converted all SUNDIALS example problems to utilize new generic
   SUNMATRIX and SUNLINEARSOLVER objects, along with updated Dls and Spils
   linear solver interfaces.

-  Added Spils interface routines to ARKode, CVODE, CVODES, IDA and IDAS
   to allow specification of a user-provided "JTSetup" routine. This change
   supports users who wish to set up data structures for the user-provided
   Jacobian-times-vector ("JTimes") routine, and where the cost of one JTSetup
   setup per Newton iteration can be amortized between multiple JTimes calls.

Two additional ``N_Vector`` implementations were added - one for CUDA and one
for RAJA vectors. These vectors are supplied to provide very basic support for
running on GPU architectures. Users are advised that these vectors both move all
data to the GPU device upon construction, and speedup will only be realized if
the user also conducts the right-hand-side function evaluation on the device. In
addition, these vectors assume the problem fits on one GPU. Further information
about RAJA, users are referred to th web site, https://software.llnl.gov/RAJA/.
These additions are accompanied by additions to various interface functions and
to user documentation.

All indices for data structures were updated to a new ``sunindextype`` that can
be configured to be a 32- or 64-bit integer data index type. ``sunindextype`` is
defined to be ``int32_t`` or ``int64_t`` when portable types are supported,
otherwise it is defined as ``int`` or ``long int``. The Fortran interfaces
continue to use for indices, except for their sparse matrix interface that now
uses the new . This new flexible capability for index types includes interfaces
to PETSc, hypre, SuperLU_MT, and KLU with either 32-bit or 64-bit capabilities
depending how the user configures SUNDIALS.

To avoid potential namespace conflicts, the macros defining ``booleantype``
values ``TRUE`` and ``FALSE`` have been changed to ``SUNTRUE`` and ``SUNFALSE``
respectively.

Temporary vectors were removed from preconditioner setup and solve routines for
all packages. It is assumed that all necessary data for user-provided
preconditioner operations will be allocated and stored in user-provided data
structures.

The file ``include/sundials_fconfig.h`` was added. This file contains SUNDIALS
type information for use in Fortran programs.

Added functions :c:func:`SUNDIALSGetVersion` and
:c:func:`SUNDIALSGetVersionNumber` to get SUNDIALS release version information
at runtime.

The build system was expanded to support many of the xSDK-compliant keys. The
xSDK is a movement in scientific software to provide a foundation for the rapid
and efficient production of high-quality, sustainable extreme-scale scientific
applications. More information can be found at, https://xsdk.info.

In addition, numerous changes were made to the build system. These include the
addition of separate ``BLAS_ENABLE`` and ``BLAS_LIBRARIES`` CMake variables,
additional error checking during CMake configuration, minor bug fixes, and
renaming CMake options to enable/disable examples for greater clarity and an
added option to enable/disable Fortran 77 examples. These changes included
changing ``EXAMPLES_ENABLE`` to ``EXAMPLES_ENABLE_C``, changing ``CXX_ENABLE``
to ``EXAMPLES_ENABLE_CXX``, changing ``F90_ENABLE`` to ``EXAMPLES_ENABLE_F90``,
and adding an ``EXAMPLES_ENABLE_F77`` option.

A bug fix was made in :c:func:`CVodeFree` to call ``lfree`` unconditionally (if
non-NULL).

Corrections and additions were made to the examples, to installation-related
files, and to the user documentation.

### CVODES Changes in v3.0.0

All interfaces to matrix structures and linear solvers have been reworked, and
all example programs have been updated. The goal of the redesign of these
interfaces was to provide more encapsulation and ease in interfacing custom
linear solvers and interoperability with linear solver libraries. Specific
changes include:

-  Added generic SUNMATRIX module with three provided implementations: dense, banded and sparse. These replicate
   previous SUNDIALS Dls and Sls matrix structures in a single object-oriented API.

-  Added example problems demonstrating use of generic SUNMATRIX modules.

-  Added generic SUNLINEARSOLVER module with eleven provided implementations: dense, banded, LAPACK dense, LAPACK band,
   KLU, SuperLU_MT, SPGMR, SPBCGS, SPTFQMR, SPFGMR, PCG. These replicate previous SUNDIALS generic linear solvers
   in a single object-oriented API.

-  Added example problems demonstrating use of generic SUNLINEARSOLVER modules.

-  Expanded package-provided direct linear solver (Dls) interfaces and scaled, preconditioned, iterative linear solver
   (Spils) interfaces to utilize generic SUNMATRIX and SUNLINEARSOLVER objects.

-  Removed package-specific, linear solver-specific, solver modules (e.g. CVDENSE, KINBAND, IDAKLU, ARKSPGMR) since
   their functionality is entirely replicated by the generic Dls/Spils interfaces and SUNLINEARSOLVER/SUNMATRIX modules.
   The exception is CVDIAG, a diagonal approximate Jacobian solver available to CVODE and CVODES.

-  Converted all SUNDIALS example problems to utilize new generic SUNMATRIX and SUNLINEARSOLVER objects, along
   with updated Dls and Spils linear solver interfaces.

-  Added Spils interface routines to ARKode, CVODE, CVODES, IDA and IDAS to allow specification of a user-provided
   "JTSetup" routine. This change supports users who wish to set up data structures for the user-provided
   Jacobian-times-vector ("JTimes") routine, and where the cost of one JTSetup setup per Newton iteration can be
   amortized between multiple JTimes calls.

Two additional ``N_Vector`` implementations were added - one for CUDA and one
for RAJA vectors. These vectors are supplied to provide very basic support for
running on GPU architectures. Users are advised that these vectors both move all
data to the GPU device upon construction, and speedup will only be realized if
the user also conducts the right-hand-side function evaluation on the device. In
addition, these vectors assume the problem fits on one GPU. Further information
about RAJA, users are referred to th web site, https://software.llnl.gov/RAJA/.
These additions are accompanied by additions to various interface functions and
to user documentation.

All indices for data structures were updated to a new ``sunindextype`` that can
be configured to be a 32- or 64-bit integer data index type. ``sunindextype`` is
defined to be ``int32_t`` or ``int64_t`` when portable types are supported,
otherwise it is defined as ``int`` or ``long int``. The Fortran interfaces
continue to use ``long int`` for indices, except for their sparse matrix
interface that now uses the new ``sunindextype``. This new flexible capability
for index types includes interfaces to PETSc, hypre, SuperLU_MT, and KLU with
either 32-bit or 64-bit capabilities depending how the user configures SUNDIALS.

To avoid potential namespace conflicts, the macros defining ``booleantype``
values ``TRUE`` and ``FALSE`` have been changed to ``SUNTRUE`` and ``SUNFALSE``
respectively.

Temporary vectors were removed from preconditioner setup and solve routines for
all packages. It is assumed that all necessary data for user-provided
preconditioner operations will be allocated and stored in user-provided data
structures.

The file ``include/sundials_fconfig.h`` was added. This file contains SUNDIALS
type information for use in Fortran programs.

Added functions ``SUNDIALSGetVersion`` and ``SUNDIALSGetVersionNumber`` to get
SUNDIALS release version information at runtime.

The build system was expanded to support many of the xSDK-compliant keys. The
xSDK is a movement in scientific software to provide a foundation for the rapid
and efficient production of high-quality, sustainable extreme-scale scientific
applications. More information can be found at, https://xsdk.info.

In addition, numerous changes were made to the build system. These include the
addition of separate ``BLAS_ENABLE`` and ``BLAS_LIBRARIES`` CMake variables,
additional error checking during CMake configuration, minor bug fixes, and
renaming CMake options to enable/disable examples for greater clarity and an
added option to enable/disable Fortran 77 examples. These changes included
changing ``EXAMPLES_ENABLE`` to ``EXAMPLES_ENABLE_C``, changing ``CXX_ENABLE``
to ``EXAMPLES_ENABLE_CXX``, changing ``F90_ENABLE`` to ``EXAMPLES_ENABLE_F90``,
and adding an ``EXAMPLES_ENABLE_F77`` option.

A bug fix was made in ``CVodeFree`` to call ``lfree`` unconditionally (if
non-NULL).

Corrections and additions were made to the examples, to installation-related
files, and to the user documentation.

### IDA Changes in v3.0.0

All interfaces to matrix structures and linear solvers have been reworked, and
all example programs have been updated.  The goal of the redesign of these
interfaces was to provide more encapsulation and to ease interfacing of custom
linear solvers and interoperability with linear solver libraries.  Specific
changes include:

* Added generic ``SUNMatrix`` module with three provided implementations: dense,
  banded, and sparse. These replicate previous SUNDIALS Dls and Sls matrix
  structures in a single object-oriented API.

* Added example problems demonstrating use of generic ``SUNMatrix`` modules.

* Added generic ``SUNLinearSolver`` module with eleven provided implementations:
  SUNDIALS native dense, SUNDIALS native banded, LAPACK dense, LAPACK band, KLU,
  SuperLU_MT, SPGMR, SPBCGS, SPTFQMR, SPFGMR, and PCG. These replicate previous
  SUNDIALS generic linear solvers in a single object-oriented API.

* Added example problems demonstrating use of generic ``SUNLinearSolver``
  modules.

* Expanded package-provided direct linear solver (Dls) interfaces and scaled,
  preconditioned, iterative linear solver (Spils) interfaces to utilize generic
  ``SUNMatrix`` and ``SUNLinearSolver`` objects.

* Removed package-specific, linear solver-specific, solver modules
  (e.g. ``CVDENSE``, ``KINBAND``, ``IDAKLU``, ``ARKSPGMR``) since their
  functionality is entirely replicated by the generic Dls/Spils interfaces and
  ``SUNLinearSolver`` and ``SUNMatrix`` modules. The exception is ``CVDIAG``, a
  diagonal approximate Jacobian solver available to CVODE and CVODES.

* Converted all SUNDIALS example problems and files to utilize the new generic
  ``SUNMatrix`` and ``SUNLinearSolver`` objects, along with updated Dls and
  Spils linear solver interfaces.

* Added Spils interface routines to ARKODE, CVODE, CVODES, IDA, and IDAS to
  allow specification of a user-provided "JTSetup" routine.  This change
  supports users who wish to set up data structures for the user-provided
  Jacobian-times-vector ("JTimes") routine, and where the cost of one JTSetup
  setup per Newton iteration can be amortized between multiple JTimes calls.

Two additional ``N_Vector`` implementations were added - one for CUDA and one
for RAJA vectors.  These vectors are supplied to provide very basic support for
running on GPU architectures. Users are advised that these vectors both move all
data to the GPU device upon construction, and speedup will only be realized if
the user also conducts the right-hand-side or residual function evaluation on
the device. In addition, these vectors assume the problem fits on one GPU.
For further information about RAJA, users are referred to the web site,
https://software.llnl.gov/RAJA/.  These additions are accompanied by updates
to various interface functions and to user documentation.

All indices for data structures were updated to a new ``sunindextype`` that can
be configured to be a 32- or 64-bit integer data index type.  ``sunindextype``
is defined to be ``int32_t`` or ``int64_t`` when portable types are supported,
otherwise it is defined as ``int`` or ``long int``.  The Fortran interfaces
continue to use ``long int`` for indices, except for their sparse matrix
interface that now uses the new ``sunindextype``.  This new flexible capability
for index types includes interfaces to PETSc, hypre, SuperLU_MT, and KLU with
either 32-bit or 64-bit capabilities depending how the user configures SUNDIALS.

To avoid potential namespace conflicts, the macros defining ``booleantype``
values ``TRUE`` and ``FALSE`` have been changed to ``SUNTRUE`` and ``SUNFALSE``
respectively.

Temporary vectors were removed from preconditioner setup and solve routines for
all packages. It is assumed that all necessary data for user-provided
preconditioner operations will be allocated and stored in user-provided data
structures.

The file ``include/sundials_fconfig.h`` was added. This file contains SUNDIALS
type information for use in Fortran programs.

The build system was expanded to support many of the xSDK-compliant keys.  The
xSDK is a movement in scientific software to provide a foundation for the rapid
and efficient production of high-quality, sustainable extreme-scale scientific
applications. More information can be found at, https://xsdk.info.

Added functions :c:func:`SUNDIALSGetVersion` and
:c:func:`SUNDIALSGetVersionNumber` to get SUNDIALS release version information
at runtime.

In addition, numerous changes were made to the build system.  These include the
addition of separate ``BLAS_ENABLE`` and ``BLAS_LIBRARIES`` CMake variables,
additional error checking during CMake configuration, minor bug fixes, and
renaming CMake options to enable/disable examples for greater clarity and an
added option to enable/disable Fortran 77 examples.  These changes included
changing ``EXAMPLES_ENABLE`` to :cmakeop:`EXAMPLES_ENABLE_C`, changing
``CXX_ENABLE`` to :cmakeop:`EXAMPLES_ENABLE_CXX`, changing ``F90_ENABLE`` to
:cmakeop:`EXAMPLES_ENABLE_F90`, and adding an :cmakeop:`EXAMPLES_ENABLE_F77`
option.

A bug fix was done to add a missing prototype for :c:func:`IDASetMaxBacksIC` in
``ida.h``.

Corrections and additions were made to the examples, to installation-related
files, and to the user documentation.

### IDAS Changes in v2.0.0

All interfaces to matrix structures and linear solvers have been reworked, and
all example programs have been updated.  The goal of the redesign of these
interfaces was to provide more encapsulation and to ease interfacing of custom
linear solvers and interoperability with linear solver libraries.  Specific
changes include:

* Added generic ``SUNMatrix`` module with three provided implementations: dense,
  banded, and sparse. These replicate previous SUNDIALS Dls and Sls matrix
  structures in a single object-oriented API.

* Added example problems demonstrating use of generic ``SUNMatrix`` modules.

* Added generic ``SUNLinearSolver`` module with eleven provided implementations:
  SUNDIALS native dense, SUNDIALS native banded, LAPACK dense, LAPACK band, KLU,
  SuperLU_MT, SPGMR, SPBCGS, SPTFQMR, SPFGMR, and PCG. These replicate previous
  SUNDIALS generic linear solvers in a single object-oriented API.

* Added example problems demonstrating use of generic ``SUNLinearSolver``
  modules.

* Expanded package-provided direct linear solver (Dls) interfaces and scaled,
  preconditioned, iterative linear solver (Spils) interfaces to utilize generic
  ``SUNMatrix`` and ``SUNLinearSolver`` objects.

* Removed package-specific, linear solver-specific, solver modules
  (e.g. ``CVDENSE``, ``KINBAND``, ``IDAKLU``, ``ARKSPGMR``) since their
  functionality is entirely replicated by the generic Dls/Spils interfaces and
  ``SUNLinearSolver`` and ``SUNMatrix`` modules. The exception is ``CVDIAG``, a
  diagonal approximate Jacobian solver available to CVODE and CVODES.

* Converted all SUNDIALS example problems and files to utilize the new generic
  ``SUNMatrix`` and ``SUNLinearSolver`` objects, along with updated Dls and
  Spils linear solver interfaces.

* Added Spils interface routines to ARKODE, CVODE, CVODES, IDAS, and IDAS to
  allow specification of a user-provided "JTSetup" routine.  This change
  supports users who wish to set up data structures for the user-provided
  Jacobian-times-vector ("JTimes") routine, and where the cost of one JTSetup
  setup per Newton iteration can be amortized between multiple JTimes calls.

Two additional ``N_Vector`` implementations were added - one for CUDA and one
for RAJA vectors.  These vectors are supplied to provide very basic support for
running on GPU architectures. Users are advised that these vectors both move all
data to the GPU device upon construction, and speedup will only be realized if
the user also conducts the right-hand-side or residual function evaluation on
the device. In addition, these vectors assume the problem fits on one GPU.
For further information about RAJA, users are referred to the web site,
https://software.llnl.gov/RAJA/.  These additions are accompanied by updates
to various interface functions and to user documentation.

All indices for data structures were updated to a new ``sunindextype`` that can
be configured to be a 32- or 64-bit integer data index type.  ``sunindextype``
is defined to be ``int32_t`` or ``int64_t`` when portable types are supported,
otherwise it is defined as ``int`` or ``long int``.  The Fortran interfaces
continue to use ``long int`` for indices, except for their sparse matrix
interface that now uses the new ``sunindextype``.  This new flexible capability
for index types includes interfaces to PETSc, hypre, SuperLU_MT, and KLU with
either 32-bit or 64-bit capabilities depending how the user configures SUNDIALS.

To avoid potential namespace conflicts, the macros defining ``booleantype``
values ``TRUE`` and ``FALSE`` have been changed to ``SUNTRUE`` and ``SUNFALSE``
respectively.

Temporary vectors were removed from preconditioner setup and solve routines for
all packages. It is assumed that all necessary data for user-provided
preconditioner operations will be allocated and stored in user-provided data
structures.

The file ``include/sundials_fconfig.h`` was added. This file contains SUNDIALS
type information for use in Fortran programs.

The build system was expanded to support many of the xSDK-compliant keys.  The
xSDK is a movement in scientific software to provide a foundation for the rapid
and efficient production of high-quality, sustainable extreme-scale scientific
applications. More information can be found at, https://xsdk.info.

Added functions :c:func:`SUNDIALSGetVersion` and
:c:func:`SUNDIALSGetVersionNumber` to get SUNDIALS release version information
at runtime.

In addition, numerous changes were made to the build system.  These include the
addition of separate ``BLAS_ENABLE`` and ``BLAS_LIBRARIES`` CMake variables,
additional error checking during CMake configuration, minor bug fixes, and
renaming CMake options to enable/disable examples for greater clarity and an
added option to enable/disable Fortran 77 examples.  These changes included
changing ``EXAMPLES_ENABLE`` to :cmakeop:`EXAMPLES_ENABLE_C`, changing
``CXX_ENABLE`` to :cmakeop:`EXAMPLES_ENABLE_CXX`, changing ``F90_ENABLE`` to
:cmakeop:`EXAMPLES_ENABLE_F90`, and adding an :cmakeop:`EXAMPLES_ENABLE_F77`
option.

A bug fix was done to add a missing prototype for :c:func:`IDASetMaxBacksIC` in
``idas.h``.

Corrections and additions were made to the examples, to installation-related
files, and to the user documentation.

### KINSOL Changes in v3.0.0

All interfaces to matrix structures and linear solvers have been reworked, and all example programs have been updated.
The goal of the redesign of these interfaces was to provide more encapsulation and ease in the interfacing of custom
linear solvers and interoperability with linear solver libraries. Specific changes include:

-  Added generic SUNMATRIX module with three provided implementations: dense, banded and sparse. These replicate
   previous SUNDIALS Dls and Sls matrix structures in a single object-oriented API.

-  Added example problems demonstrating use of generic SUNMATRIX modules.

-  Added generic ``SUNLinearSolver`` module with eleven provided implementations: SUNDIALS native dense,
   SUNDIALS native banded, LAPACK dense, LAPACK band, KLU, SuperLU_MT, SPGMR, SPBCGS, SPTFQMR, SPFGMR, and PCG.
   These replicate previous SUNDIALS generic linear solvers in a single object-oriented API.

-  Added example problems demonstrating use of generic SUNLINEARSOLVER modules.

-  Expanded package-provided direct linear solver (Dls) interfaces and scaled, preconditioned, iterative linear solver
   (Spils) interfaces to utilize generic SUNMATRIX and SUNLINEARSOLVER objects.

-  Removed package-specific, linear solver-specific, solver modules (e.g. CVDENSE, KINBAND, IDAKLU, ARKSPGMR) since
   their functionality is entirely replicated by the generic Dls/Spils interfaces and SUNLINEARSOLVER/SUNMATRIX modules.
   The exception is CVDIAG, a diagonal approximate Jacobian solver available to CVODE and CVODES.

-  Converted all SUNDIALS example problems to utilize new generic SUNMATRIX and SUNLINEARSOLVER objects, along
   with updated Dls and Spils linear solver interfaces.

-  Added Spils interface routines to ARKode, CVODE, CVODES, IDA and IDAS to allow specification of a user-provided
   "JTSetup" routine. This change supports users who wish to set up data structures for the user-provided
   Jacobian-times-vector ("JTimes") routine, and where the cost of one JTSetup setup per Newton iteration can be
   amortized between multiple JTimes calls.

Two additional ``N_Vector`` implementations were added - one for CUDA and one for RAJA vectors. These
vectors are supplied to provide very basic support for running on GPU architectures. Users are advised that these
vectors both move all data to the GPU device upon construction, and speedup will only be realized if the user also
conducts the right-hand-side function evaluation on the device. In addition, these vectors assume the problem fits on
one GPU. Further information about RAJA, users are referred to th web site, https://software.llnl.gov/RAJA/. These
additions are accompanied by additions to various interface functions and to user documentation.

All indices for data structures were updated to a new ``sunindextype`` that can be configured to be a 32- or 64-bit
integer data index type. ``sunindextype`` is defined to be ``int32_t`` or ``int64_t`` when portable types are supported,
otherwise it is defined as ``int`` or ``long int``. The Fortran interfaces continue to use ``long int`` for indices,
except for their sparse matrix interface that now uses the new ``sunindextype``. This new flexible capability for index
types includes interfaces to PETSc, hypre, SuperLU_MT, and KLU with either 32-bit or 64-bit capabilities depending how
the user configures SUNDIALS.

To avoid potential namespace conflicts, the macros defining ``booleantype`` values ``TRUE`` and ``FALSE`` have been
changed to ``SUNTRUE`` and ``SUNFALSE`` respectively.

Temporary vectors were removed from preconditioner setup and solve routines for all packages. It is assumed that all
necessary data for user-provided preconditioner operations will be allocated and stored in user-provided data
structures.

The file ``include/sundials_fconfig.h`` was added. This file contains SUNDIALS type information for use in Fortran
programs.

The build system was expanded to support many of the xSDK-compliant keys. The xSDK is a movement in scientific software
to provide a foundation for the rapid and efficient production of high-quality, sustainable extreme-scale scientific
applications. More information can be found at, https://xsdk.info.

Added functions ``SUNDIALSGetVersion`` and ``SUNDIALSGetVersionNumber`` to get SUNDIALS release version
information at runtime.

In addition, numerous changes were made to the build system. These include the addition of separate ``BLAS_ENABLE`` and
``BLAS_LIBRARIES`` CMake variables, additional error checking during CMake configuration, minor bug fixes, and renaming
CMake options to enable/disable examples for greater clarity and an added option to enable/disable Fortran 77 examples.
These changes included changing ``EXAMPLES_ENABLE`` to ``EXAMPLES_ENABLE_C``, changing ``CXX_ENABLE`` to
``EXAMPLES_ENABLE_CXX``, changing ``F90_ENABLE`` to ``EXAMPLES_ENABLE_F90``, and adding an ``EXAMPLES_ENABLE_F77``
option.

A bug fix was done to correct the fcmix name translation for ``FKIN_SPFGMR``.

Corrections and additions were made to the examples, to installation-related files, and to the user documentation.


## Changes to SUNDIALS in release v2.7.0

- Two new NVECTOR modules added: for _hypre_ ParVector and PETSC.
- In vector API, added new required function, N\_VGetVectorID.
- Upgrades to sparse solver interfaces; now support CSR matrix type with KLU solver.
- In all packages, example codes were changed from using NV\_DATA macro to using N\_VGetArrayPointer\_\* when using the native vectors shipped with SUNDIALS
- In all packages, fixed memory leak in banded preconditioner interface.
- Fixed some examples w.r.t. switch to new macro/function names SUNRexp etc.
- Various minor fixes to installation-related files.
- Corrected name N\_VCloneEmptyVectorArray to N\_VCloneVectorArrayEmpty in all documentation files.
- Updated all packages to return integers from linear solver and preconditioner 'free' functions.
- Removed Matlab interface from distribution as it has not been updated since 2009. We expect to update this interface soon.
- In FKINSOL, FCVODE, and FIDA, added missing Fortran interface routines so that users can supply the sparse Jacobian routine.
- Minor corrections and additions to all User Guides, including removal of references to specific NVECTOR names in usage skeletons.
- Additional example programs added throughout.
- In CVODE
  - in FCVODE, fixed argument order bugs in FCVKLU and FCVSUPERLUMT linear solver interfaces.
- In CVODES
  - changed each \*\*FreeB() to type int; added return(0) to each.
  - in interpolation routines for backward problems, added logic to bypass sensitivity interpolation if input sensitivity argument is NULL.
- In ARKODE
  - updated linear and mass matrix solvers so that 'free' routines return integer instead of void; updated documentation accordingly.
  - fixed initialization of linear solver performance counters.
  - method and embedding for Billington and TRBDF2 explicit Runge-Kutta methods were swapped.
  - fix for user specification of absolute tolerance array along with vector Resize() functionality.
  - fix for user-supplied Butcher tables without embeddings (if fixed time steps or manual adaptivity are employed).
  - multiple documentation updates.
  - added missing ARKSpilsGetNumMtimesEvals() function.
  - implicit predictor algorithms were updated: methods 2 and 3 were improved, a new predictor approach was added, and the default choice was modified.
  - revised handling of integer codes for specifying built-in Butcher tables: a global numbering system is still used, but methods now have #defined names to simplify the user interface.
  - maximum number of Butcher table stages was increased from 8 to 15 to accommodate very high-order methods, and an 8th-order adaptive ERK method was added.
  - added support for the explicit and implicit methods in an additive Runge-Kutta method to utilize different stage times, solution and embedding coefficients, to support new SSP-ARK methods.
  - extended FARKODE interface to include a routine to set scalar/array-valued residual tolerances, to support Fortran applications with non-identity mass-matrices.
- In IDA
  - corrected example idaFoodWeb\_bnd.c in PrintOutput (wrong component printed).
  - added optional input function IDASetMaxBacksIC to limit number of linesearch backtrack operations in IDACalcIC. User guides amended accordingly.
- In IDAS
  - added optional input function IDASetMaxBacksIC to limit number of linesearch backtrack operations in IDACalcIC. User guides amended accordingly.
  - changed each \*\*FreeB() to type int; added return(0) to each.
  - in interpolation routines for backward problems, added logic to bypass sensitivity interpolation if input sensitivity argument is NULL.
- In KINSOL
  - minor bug fix in Picard iteration.
  - minor bug fix in line search to prevent infinite loop when beta condition fails and lambda is below minimum size.

### ARKODE Changes in v1.1.0

We have included numerous bugfixes and enhancements since the
v1.0.2 release.

The bugfixes include:

* For each linear solver, the various solver performance counters are
  now initialized to 0 in both the solver specification function and
  in the solver's ``linit`` function.  This ensures that these solver
  counters are initialized upon linear solver instantiation as well as
  at the beginning of the problem solution.

* The choice of the method vs embedding the Billington and TRBDF2
  explicit Runge--Kutta methods were swapped, since in those the
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

* The ARKODE implicit predictor algorithms were updated: methods 2 and
  3 were improved slightly, a new predictor approach was added, and
  the default choice was modified.

* The underlying sparse matrix structure was enhanced to allow both
  CSR and CSC matrices, with CSR supported by the KLU linear solver
  interface.  ARKODE interfaces to the KLU solver from both C and
  Fortran were updated to enable selection of sparse matrix type, and a
  Fortran-90 CSR example program was added.

* The missing :c:func:`ARKSpilsGetNumMtimesEvals()` function was added
  -- this had been included in the previous documentation but had not
  been implemented.

* The handling of integer codes for specifying built-in ARKODE Butcher
  tables was enhanced.  While a global numbering system is still used,
  methods now have #defined names to simplify the user interface and to
  streamline incorporation of new Butcher tables into ARKODE.

* The maximum number of Butcher table stages was increased from 8 to
  15 to accommodate very high order methods, and an 8th-order adaptive
  ERK method was added.

* Support was added for the explicit and implicit methods in an
  additive Runge--Kutta method to utilize different stage times,
  solution and embedding coefficients, to support new SSP-ARK
  methods.

* The FARKODE interface was extended to include a routine to set
  scalar/array-valued residual tolerances, to support Fortran
  applications with non-identity mass-matrices.

### CVODE Changes in v2.9.0

Two additional ``N_Vector`` implementations were added - one for Hypre
(parallel) ParVector vectors, and one for PETSc vectors. These
additions are accompanied by additions to various interface functions
and to user documentation.

Each ``N_Vector`` module now includes a function, :c:func:`N_VGetVectorID`, that returns the
``N_Vector`` module name.

For each linear solver, the various solver performance counters are now
initialized to 0 in both the solver specification function and in solver
``linit`` function. This ensures that these solver counters are initialized upon
linear solver instantiation as well as at the beginning of the problem
solution.

In FCVODE, corrections were made to three Fortran interface
functions. Missing Fortran interface routines were added so that users
can supply the sparse Jacobian routine when using sparse direct solvers.

A memory leak was fixed in the banded preconditioner interface. In
addition, updates were done to return integers from linear solver and
preconditioner 'free' functions.

The Krylov linear solver Bi-CGstab was enhanced by removing a redundant
dot product. Various additions and corrections were made to the
interfaces to the sparse solvers KLU and SuperLU_MT, including support
for CSR format when using KLU.

New examples were added for use of the OpenMP vector and for use of
sparse direct solvers from Fortran.

Minor corrections and additions were made to the CVODE solver, to
the Fortran interfaces, to the examples, to installation-related files,
and to the user documentation.

### CVODES Changes in v2.9.0

Two additional ``N_Vector`` implementations were added - one for Hypre
(parallel) ParVector vectors, and one for PETSc vectors. These additions are
accompanied by additions to various interface functions and to user
documentation.

Each ``N_Vector`` module now includes a function, ``N_VGetVectorID``, that
returns the ``N_Vector`` module name.

A bug was fixed in the interpolation functions used in solving backward problems
for adjoint sensitivity analysis.

For each linear solver, the various solver performance counters are now
initialized to 0 in both the solver specification function and in solver
``linit`` function. This ensures that these solver counters are initialized upon
linear solver instantiation as well as at the beginning of the problem solution.

A memory leak was fixed in the banded preconditioner interface. In addition,
updates were done to return integers from linear solver and preconditioner
'free' functions.

The Krylov linear solver Bi-CGstab was enhanced by removing a redundant dot
product. Various additions and corrections were made to the interfaces to the
sparse solvers KLU and SuperLU_MT, including support for CSR format when using
KLU.

In interpolation routines for backward problems, added logic to bypass
sensitivity interpolation if input sensitivity argument is NULL.

New examples were added for use of sparse direct solvers within sensitivity
integrations and for use of OpenMP.

Minor corrections and additions were made to the CVODES solver, to the examples,
to installation-related files, and to the user documentation.

### IDA Changes in v2.9.0

Two additional ``N_Vector`` implementations were added - one for Hypre
(parallel) ParVector vectors, and one for PETSc vectors. These additions are
accompanied by additions to various interface functions and to user
documentation.

Each ``N_Vector`` module now includes a function, :c:func:`N_VGetVectorID`, that
returns the ``N_Vector`` module name.

An optional input function was added to set a maximum number of linesearch
backtracks in the initial condition calculation.  Also, corrections were made to
three Fortran interface functions.

For each linear solver, the various solver performance counters are now
initialized to 0 in both the solver specification function and in solver
``linit`` function. This ensures that these solver counters are initialized upon
linear solver instantiation as well as at the beginning of the problem solution.

A memory leak was fixed in the banded preconditioner interface.  In addition,
updates were done to return integers from linear solver and preconditioner
"free" functions.

The Krylov linear solver Bi-CGstab was enhanced by removing a redundant dot
product. Various additions and corrections were made to the interfaces to the
sparse solvers KLU and SuperLU_MT, including support for CSR format when using
KLU.

New examples were added for use of the OpenMP vector.

Minor corrections and additions were made to the IDA solver, to the Fortran
interfaces, to the examples, to installation-related files, and to the user
documentation.

### IDAS Changes in v1.3.0

Two additional ``N_Vector`` implementations were added - one for Hypre
(parallel) ParVector vectors, and one for PETSc vectors. These additions are
accompanied by additions to various interface functions and to user
documentation.

Each ``N_Vector`` module now includes a function, :c:func:`N_VGetVectorID`, that
returns the ``N_Vector`` module name.

An optional input function was added to set a maximum number of linesearch
backtracks in the initial condition calculation.  Also, corrections were made to
three Fortran interface functions.

For each linear solver, the various solver performance counters are now
initialized to 0 in both the solver specification function and in solver
``linit`` function. This ensures that these solver counters are initialized upon
linear solver instantiation as well as at the beginning of the problem solution.

A bug in for-loop indices was fixed in :c:func:`IDAAckpntAllocVectors`. A bug
was fixed in the interpolation functions used in solving backward problems.

A memory leak was fixed in the banded preconditioner interface.  In addition,
updates were done to return integers from linear solver and preconditioner
"free" functions.

The Krylov linear solver Bi-CGstab was enhanced by removing a redundant dot
product. Various additions and corrections were made to the interfaces to the
sparse solvers KLU and SuperLU_MT, including support for CSR format when using
KLU.

New examples were added for use of the OpenMP vector.

Minor corrections and additions were made to the IDAS solver, to the examples,
to installation-related files, and to the user documentation.

### KINSOL Changes in v2.9.0

Two additional ``N_Vector`` implementations were added - one for Hypre (parallel) vectors, and one for PETSc vectors.
These additions are accompanied by additions to various interface functions and to user documentation.

Each ``N_Vector`` module now includes a function, ``N_VGetVectorID``, that returns the ``N_Vector`` module name.

The Picard iteration return was chanegd to always return the newest iterate upon success. A minor bug in the line search
was fixed to prevent an infinite loop when the beta condition fails and lamba is below the minimum size.

For each linear solver, the various solver performance counters are now initialized to 0 in both the solver
specification function and in solver ``linit`` function. This ensures that these solver counters are initialized upon
linear solver instantiation as well as at the beginning of the problem solution.

A memory leak was fixed in the banded preconditioner interface. In addition, updates were done to return integers from
linear solver and preconditioner 'free' functions.

Corrections were made to three Fortran interface functions. The Anderson acceleration scheme was enhanced by use of QR
updating.

The Krylov linear solver Bi-CGstab was enhanced by removing a redundant dot product. Various additions and corrections
were made to the interfaces to the sparse solvers KLU and SuperLU_MT, including support for CSR format when using KLU.

The functions FKINCREATE and FKININIT were added to split the FKINMALLOC routine into two pieces. FKINMALLOC remains for
backward compatibility, but documentation for it has been removed.

A new examples was added for use of the OpenMP vector.

Minor corrections and additions were made to the KINSOL solver, to the Fortran interfaces, to the examples, to
installation-related files, and to the user documentation.


## Changes to SUNDIALS in release v2.6.2

- In IDAS, added missing backward problem support functions: IDALapackDenseB, IDALapackDenseFreeB, IDALapackBandB, IDALapackBandFreeB
- In KINSOL and ARKode, updated Anderson acceleration implementation with QR updating.
- Updated BiCGStab solver to remove redundant dot product call.
- Minor corrections and additions to all User Guides.
- In CVODES and IDAS header files, corrected documentation of backward integration functions, especially the 'which' argument.
- In CVODES, added DVKLUB prototype and corrected CVSuperLUMTB prototype.
- In IDAS, made SuperLUMT call for backward problem consistent with CVODES.
- In CVODES and IDAS, added ReInit and SetOrdering wrappers for backward problems. Fixed potential memory leak in KLU ReInit functions in all solvers.
- In CVODE, IDA, and ARKode, fixed Fortran interfaces to enable calls to \*GetErrWeights, \*GetEstLocalErrors, and \*GetDky within a time step. In ARKode, fixed a bug in one Butcher table.
- In ARKode, fixed error in arkDoErrorTest in recovery after failure.
- In IDAS, fixed for-loop bugs in IDAAckpntAllocVectors Various minor fixes to installation-related files.

### ARKODE Changes in v1.0.2

None?

### CVODE Changes in v2.8.2

None?

### CVODES Changes in v2.8.2

None?

### IDA Changes in v2.8.2

None?

### IDAS Changes in v1.2.2

None?

### KINSOL Changes in v2.8.2

None?

## Changes to SUNDIALS in release v2.6.1

- Fixed loop limit bug in SlsAddMat function.
- In all six solver interfaces to KLU and SuperLUMT, added #include lines, and removed redundant KLU structure allocations.
- Numerous minor documentation improvements
- Minor bug fixes in ARKode

### ARKODE Changes in v1.0.1

None?

### CVODE Changes in v2.8.1

None?

### CVODES Changes in v2.8.1

None?

### IDA Changes in v2.8.1

None?

### IDAS Changes in v1.2.1

None?

### KINSOL Changes in v2.8.1

None?

## Changes to SUNDIALS in release v2.6.0

- Addition of ARKode package of explicit, implicit, and additive Runge-Kutta methods for ODES. This package API is close to CVODE so switching between the two should be straightforward. Thanks go to Daniel Reynolds for the addition of this package.
- Addition of support for two sparse direct solver packages when using the serial vector structure, KLU and SuperLU\_MT. exploits highly sparse systems. SuperLU\_MT supports multithreading in the factorization.
- Addition of openMP and PThreads vector kernels.
- Addition of fixed point and Picard iterative solvers within KINSOL. These are both optionally accelerated with Anderson acceleration.
- Addition of FGMRES support for KINSOL.
- Removal of autotools configuration support. We now exclusively use CMake.
- Numerous bug fixes throughout.

### CVODE Changes in v2.8.0

Two major additions were made to the linear system solvers that are available
for use with the CVODE solver. First, in the serial case, an interface to the
sparse direct solver KLU was added. Second, an interface to SuperLU_MT, the
multi-threaded version of SuperLU, was added as a thread-parallel sparse direct
solver option, to be used with the serial version of the ``N_Vector`` module. As
part of these additions, a sparse matrix (CSC format) structure was added to
CVODE.

Otherwise, only relatively minor modifications were made to the CVODE solver:

In ``cvRootFind``, a minor bug was corrected, where the input array was ignored,
and a line was added to break out of root-search loop if the initial interval
size is below the tolerance ``ttol``.

In ``CVLapackBand``, the line ``smu = MIN(N-1,mu+ml)`` was changed to to correct
an illegal input error for ``DGBTRF/DGBTRS``.

In order to eliminate or minimize the differences between the sources for
private functions in CVODE and CVODES, the names of 48 private functions were
changed from to , and a few other names were also changed.

Two minor bugs were fixed regarding the testing of input on the first call to -
one involving and one involving the initialization of ``*tret``.

In order to avoid possible name conflicts, the mathematical macro and function
names ``MIN``, ``MAX``, ``SQR``, ``RAbs``, ``RSqrt``, ``RExp``, ``RPowerI``, and
were changed to ``SUNMIN``, ``SUNMAX``, ``SUNSQR``, ``SUNRabs``, ``SUNRsqrt``,
``SUNRexp``, ``SUNRpowerI``, and ``SUNRPowerR`` respectively. These names occur
in both the solver and in various example programs.

The example program ``cvAdvDiff_diag_p`` was added to illustrate the use of in
parallel.

In the FCVODE optional input routines ``FCVSETIIN`` and ``FCVSETRIN``, the
optional fourth argument ``key_length`` was removed, with hardcoded key string
lengths passed to all tests.

In all FCVODE examples, integer declarations were revised so that those which
must match a C type ``long int`` are declared ``INTEGER*8``, and a comment was
added about the type match. All other integer declarations are just ``INTEGER``.
Corresponding minor corrections were made to the user guide.

Two new ``N_Vector`` modules have been added for thread-parallel computing
environments - one for OpenMP, denoted ``NVECTOR_OPENMP``, and one for Pthreads,
denoted ``NVECTOR_PTHREADS``.

With this version of SUNDIALS, support and documentation of the Autotools mode
of installation is being dropped, in favor of the CMake mode, which is
considered more widely portable.

### CVODES Changes in v2.8.0

Two major additions were made to the linear system solvers that are available
for use with the CVODES solver. First, in the serial case, an interface to the
sparse direct solver KLU was added. Second, an interface to SuperLU_MT, the
multi-threaded version of SuperLU, was added as a thread-parallel sparse direct
solver option, to be used with the serial version of the ``N_Vector`` module. As
part of these additions, a sparse matrix (CSC format) structure was added to
CVODES.

Otherwise, only relatively minor modifications were made to the CVODES solver:

In ``cvRootfind``, a minor bug was corrected, where the input array ``rootdir``
was ignored, and a line was added to break out of root-search loop if the
initial interval size is below the tolerance ``ttol``.

In ``CVLapackBand``, the line ``smu = MIN(N-1,mu+ml)`` was changed to ``smu = mu
+ ml`` to correct an illegal input error for ``DGBTRF/DGBTRS``.

Some minor changes were made in order to minimize the differences between the
sources for private functions in CVODES and CVODE.

An option was added in the case of Adjoint Sensitivity Analysis with dense or
banded Jacobian: With a call to ``CVDlsSetDenseJacFnBS`` or
``CVDlsSetBandJacFnBS``, the user can specify a user-supplied Jacobian function
of type ``CVDls***JacFnBS``, for the case where the backward problem depends on
the forward sensitivities.

In ``CVodeQuadSensInit``, the line ``cv_mem->cv_fQS_data = ...`` was corrected
(missing ``Q``).

In the User Guide, a paragraph was added in Section 6.2.1 on ``CVodeAdjReInit``,
and a paragraph was added in Section 6.2.9 on ``CVodeGetAdjY``. In the example
``cvsRoberts_ASAi_dns``, the output was revised to include the use of
``CVodeGetAdjY``.

Two minor bugs were fixed regarding the testing of input on the first call to
``CVode`` - one involving ``tstop`` and one involving the initialization of
``*tret``.

For the Adjoint Sensitivity Analysis case in which the backward problem depends
on the forward sensitivities, options have been added to allow for user-supplied
``pset``, ``psolve``, and ``jtimes`` functions.

In order to avoid possible name conflicts, the mathematical macro and function
names ``MIN``, ``MAX``, ``SQR``, ``RAbs``, ``RSqrt``, ``RExp``, ``RPowerI``, and
``RPowerR`` were changed to ``SUNMIN``, ``SUNMAX``, ``SUNSQR``, ``SUNRabs``,
``SUNRsqrt``, ``SUNRexp``, ``SRpowerI``, and ``SUNRpowerR``, respectively. These
names occur in both the solver and example programs.

In the example ``cvsHessian_ASA_FSA``, an error was corrected in the function
``fB2``: ``y2`` in place of ``y3`` in the third term of ``Ith(yBdot,6)``.

Two new ``N_Vector`` modules have been added for thread-parallel computing
environments - one for OpenMP, denoted ``NVECTOR_OPENMP``, and one for Pthreads,
denoted ``NVECTOR_PTHREADS``.

With this version of SUNDIALS, support and documentation of the Autotools mode
of installation is being dropped, in favor of the CMake mode, which is
considered more widely portable.

### IDA Changes in v2.8.0

Two major additions were made to the linear system solvers that are available
for use with the IDA solver. First, in the serial case, an interface to the
sparse direct solver KLU was added.  Second, an interface to SuperLU_MT, the
multi-threaded version of SuperLU, was added as a thread-parallel sparse direct
solver option, to be used with the serial version of the ``N_Vector`` module.
As part of these additions, a sparse matrix (CSC format) structure was added to
IDA.

Otherwise, only relatively minor modifications were made to IDA:

In :c:func:`IDARootfind`, a minor bug was corrected, where the input array
``rootdir`` was ignored, and a line was added to break out of root-search loop
if the initial interval size is below the tolerance ``ttol``.

In ``IDALapackBand``, the line ``smu = MIN(N-1,mu+ml)`` was changed to ``smu =
mu + ml`` to correct an illegal input error for ``DGBTRF/DGBTRS``.

A minor bug was fixed regarding the testing of the input ``tstop`` on the first
call to :c:func:`IDASolve`.

In order to avoid possible name conflicts, the mathematical macro and function
names ``MIN``, ``MAX``, ``SQR``, ``RAbs``, ``RSqrt``, ``RExp``, ``RPowerI``, and
``RPowerR`` were changed to ``SUNMIN``, ``SUNMAX``, ``SUNSQR``, ``SUNRabs``,
``SUNRsqrt``, ``SUNRexp``, ``SRpowerI``, and ``SUNRpowerR``, respectively.
These names occur in both the solver and in various example programs.

In the FIDA optional input routines ``FIDASETIIN``, ``FIDASETRIN``, and
``FIDASETVIN``, the optional fourth argument ``key_length`` was removed, with
hardcoded key string lengths passed to all ``strncmp`` tests.

In all FIDA examples, integer declarations were revised so that those which must
match a C type ``long int`` are declared ``INTEGER*8``, and a comment was added
about the type match. All other integer declarations are just
``INTEGER``. Corresponding minor corrections were made to the user guide.

Two new ``N_Vector`` modules have been added for thread-parallel computing
environments - one for OpenMP, denoted :ref:`NVECTOR_OPENMP <NVectors.OpenMP>`,
and one for Pthreads, denoted :ref:`NVECTOR_PTHREADS <NVectors.Pthreads>`.

With this version of SUNDIALS, support and documentation of the Autotools mode
of installation is being dropped, in favor of the CMake mode, which is
considered more widely portable.

### IDAS Changes in v1.2.0

Two major additions were made to the linear system solvers that are available
for use with the IDAS solver. First, in the serial case, an interface to the
sparse direct solver KLU was added.  Second, an interface to SuperLU_MT, the
multi-threaded version of SuperLU, was added as a thread-parallel sparse direct
solver option, to be used with the serial version of the ``N_Vector`` module.
As part of these additions, a sparse matrix (CSC format) structure was added to
IDAS.

Otherwise, only relatively minor modifications were made to IDAS:

In :c:func:`IDARootfind`, a minor bug was corrected, where the input array
``rootdir`` was ignored, and a line was added to break out of root-search loop
if the initial interval size is below the tolerance ``ttol``.

In ``IDALapackBand``, the line ``smu = MIN(N-1,mu+ml)`` was changed to ``smu =
mu + ml`` to correct an illegal input error for ``DGBTRF/DGBTRS``.

An option was added in the case of Adjoint Sensitivity Analysis with dense or
banded Jacobian: With a call to ``IDADlsSetDenseJacFnBS`` or
``IDADlsSetBandJacFnBS``, the user can specify a user-supplied Jacobian function
of type ``IDADls***JacFnBS``, for the case where the backward problem depends on
the forward sensitivities.

A minor bug was fixed regarding the testing of the input ``tstop`` on the first
call to :c:func:`IDASolve`.

In order to avoid possible name conflicts, the mathematical macro and function
names ``MIN``, ``MAX``, ``SQR``, ``RAbs``, ``RSqrt``, ``RExp``, ``RPowerI``, and
``RPowerR`` were changed to ``SUNMIN``, ``SUNMAX``, ``SUNSQR``, ``SUNRabs``,
``SUNRsqrt``, ``SUNRexp``, ``SRpowerI``, and ``SUNRpowerR``, respectively.
These names occur in both the solver and in various example programs.

In the FIDA optional input routines ``FIDASETIIN``, ``FIDASETRIN``, and
``FIDASETVIN``, the optional fourth argument ``key_length`` was removed, with
hardcoded key string lengths passed to all ``strncmp`` tests.

In all FIDA examples, integer declarations were revised so that those which must
match a C type ``long int`` are declared ``INTEGER*8``, and a comment was added
about the type match. All other integer declarations are just
``INTEGER``. Corresponding minor corrections were made to the user guide.

Two new ``N_Vector`` modules have been added for thread-parallel computing
environments - one for OpenMP, denoted :ref:`NVECTOR_OPENMP <NVectors.OpenMP>`,
and one for Pthreads, denoted :ref:`NVECTOR_PTHREADS <NVectors.Pthreads>`.

With this version of SUNDIALS, support and documentation of the Autotools mode
of installation is being dropped, in favor of the CMake mode, which is
considered more widely portable.

### KINSOL Changes in v2.8.0

Two major additions were made to the globalization strategy options (``KINSol`` argument ``strategy``). One is
fixed-point iteration, and the other is Picard iteration. Both can be accelerated by use of the Anderson acceleration
method. See the relevant paragraphs in Chapter :numref:`KINSOL.Mathematics`.

Three additions were made to the linear system solvers that are available for use with the KINSOL solver. First,
in the serial case, an interface to the sparse direct solver KLU was added. Second, an interface to SuperLU_MT, the
multi-threaded version of SuperLU, was added as a thread-parallel sparse direct solver option, to be used with the
serial version of the ``N_Vector`` module. As part of these additions, a sparse matrix (CSC format) structure was added
to KINSOL. Finally, a variation of GMRES called Flexible GMRES was added.

Otherwise, only relatively minor modifications were made to KINSOL:

In function ``KINStop``, two return values were corrected to make the values of ``uu`` and ``fval`` consistent.

A bug involving initialization of ``mxnewtstep`` was fixed. The error affects the case of repeated user calls to
``KINSol`` with no intervening call to ``KINSetMaxNewtonStep``.

A bug in the increments for difference quotient Jacobian approximations was fixed in function ``kinDlsBandDQJac``.

In ``KINLapackBand``, the line ``smu = MIN(N-1,mu+ml)`` was changed to ``smu = mu + ml`` to correct an illegal input
error for ``DGBTRF/DGBTRS``.

In order to avoid possible name conflicts, the mathematical macro and function names ``MIN``, ``MAX``, ``SQR``,
``RAbs``, ``RSqrt``, ``RExp``, ``RPowerI``, and ``RPowerR`` were changed to ``SUNMIN``, ``SUNMAX``, ``SUNSQR``,
``SUNRabs``, ``SUNRsqrt``, ``SUNRexp``, ``SRpowerI``, and ``SUNRpowerR``, respectively. These names occur in both the
solver and in various example programs.

In the FKINSOL module, an incorrect return value ``ier`` in ``FKINfunc`` was fixed.

In the FKINSOL optional input routines ``FKINSETIIN``, ``FKINSETRIN``, and ``FKINSETVIN``, the optional fourth argument
``key_length`` was removed, with hardcoded key string lengths passed to all ``strncmp`` tests.

In all FKINSOL examples, integer declarations were revised so that those which must match a C type ``long int`` are
declared ``INTEGER*8``, and a comment was added about the type match. All other integer declarations are just
``INTEGER``. Corresponding minor corrections were made to the user guide.

Two new ``N_Vector`` modules have been added for thread-parallel computing environments - one for OpenMP, denoted
``NVECTOR_OPENMP``, and one for Pthreads, denoted ``NVECTOR_PTHREADS``.

With this version of SUNDIALS, support and documentation of the Autotools mode of installation is being dropped,
in favor of the CMake mode, which is considered more widely portable.


## Changes to SUNDIALS in release v2.5.0

- Changes to user interface
  - Problem size and related integers (bandwidth parameters etc.) all have type long int, except in BLAS and LAPACK routines. Function NewIntArray is replaced by a pair NewIntArray/NewLintArray, for int and long int arrays, respectively.

### CVODE Changes in v2.7.0

One significant design change was made with this release: The problem size and
its relatives, bandwidth parameters, related internal indices, pivot arrays, and
the optional output ``lsflag`` have all been changed from type ``int`` to type
``long int``, except for the problem size and bandwidths in user calls to
routines specifying BLAS/LAPACK routines for the dense/band linear solvers. The
function ``NewIntArray`` is replaced by a pair ``NewIntArray`` /
``NewLintArray``, for ``int`` and ``long int`` arrays, respectively.

A large number of minor errors have been fixed. Among these are the following:
In , the logic was changed to avoid a divide by zero. After the solver memory is
created, it is set to zero before being filled. In ``CVSetTqBDF`` each linear
solver interface function, the linear solver memory is freed on an error return,
and the function now includes a line setting to NULL the main memory pointer to
the linear solver memory. In the rootfinding functions ``CVRcheck1``/
``CVRcheck2``, when an exact zero is found, the array ``glo`` of :math:`g`
values at the left endpoint is adjusted, instead of shifting the :math:`t`
location slightly. In the installation files, we modified the treatment of the
macro SUNDIALS_USE_GENERIC_MATH, so that the parameter GENERIC_MATH_LIB is
either defined (with no value) or not defined.

### CVODES Changes in v2.7.0

One significant design change was made with this release: The problem size and
its relatives, bandwidth parameters, related internal indices, pivot arrays, and
the optional output ``lsflag`` have all been changed from type ``int`` to type
``long int``, except for the problem size and bandwidths in user calls to
routines specifying BLAS/LAPACK routines for the dense/band linear solvers. The
function ``NewIntArray`` is replaced by a pair ``NewIntArray`` / ``NewLintArray``,
for ``int`` and ``long int`` arrays, respectively. In a minor change to the user
interface, the type of the index ``which`` in CVODES was changed from ``long
int`` to ``int``.

Errors in the logic for the integration of backward problems were identified and
fixed.

A large number of minor errors have been fixed. Among these are the following:
In ``CVSetTqBDF``, the logic was changed to avoid a divide by zero. After the
solver memory is created, it is set to zero before being filled. In each linear
solver interface function, the linear solver memory is freed on an error return,
and the ``**Free`` function now includes a line setting to NULL the main memory
pointer to the linear solver memory. In the rootfinding functions
``CVRcheck1`` / ``CVRcheck2``, when an exact zero is found, the array ``glo`` of
:math:`g` values at the left endpoint is adjusted, instead of shifting the
:math:`t` location ``tlo`` slightly. In the installation files, we modified the
treatment of the macro SUNDIALS_USE_GENERIC_MATH, so that the parameter
GENERIC_MATH_LIB is either defined (with no value) or not defined.

### IDA Changes in v2.7.0

One significant design change was made with this release: The problem size and
its relatives, bandwidth parameters, related internal indices, pivot arrays, and
the optional output ``lsflag`` have all been changed from type ``int`` to type
``long int``, except for the problem size and bandwidths in user calls to
routines specifying BLAS/LAPACK routines for the dense/band linear solvers. The
function ``NewIntArray`` is replaced by a pair ``NewIntArray`` and
``NewLintArray``, for ``int`` and ``long int`` arrays, respectively.

A large number of minor errors have been fixed. Among these are the following:
After the solver memory is created, it is set to zero before being filled.  To
be consistent with IDAS, IDA uses the function ``IDAGetDky`` for optional output
retrieval.  In each linear solver interface function, the linear solver memory
is freed on an error return, and the ``**Free`` function now includes a line
setting to NULL the main memory pointer to the linear solver memory.  A memory
leak was fixed in two of the ``IDASp***Free`` functions.  In the rootfinding
functions ``IDARcheck1`` and ``IDARcheck2``, when an exact zero is found, the
array ``glo`` of :math:`g` values at the left endpoint is adjusted, instead of
shifting the :math:`t` location ``tlo`` slightly.  In the installation files, we
modified the treatment of the macro SUNDIALS_USE_GENERIC_MATH, so that the
parameter GENERIC_MATH_LIB is either defined (with no value) or not defined.

### IDAS Changes in v1.1.0

One significant design change was made with this release: The problem size and
its relatives, bandwidth parameters, related internal indices, pivot arrays, and
the optional output ``lsflag`` have all been changed from type ``int`` to type
``long int``, except for the problem size and bandwidths in user calls to
routines specifying BLAS/LAPACK routines for the dense/band linear solvers. The
function ``NewIntArray`` is replaced by a pair ``NewIntArray`` and
``NewLintArray``, for ``int`` and ``long int`` arrays, respectively.

Errors in the logic for the integration of backward problems were identified and
fixed. A large number of minor errors have been fixed. Among these are the
following: A missing vector pointer setting was added in
:c:func:`IDASensLineSrch`. In :c:func:`IDACompleteStep`, conditionals around
lines loading a new column of three auxiliary divided difference arrays, for a
possible order increase, were fixed. After the solver memory is created, it is
set to zero before being filled. In each linear solver interface function, the
linear solver memory is freed on an error return, and the ``**Free`` function
now includes a line setting to ``NULL`` the main memory pointer to the linear
solver memory. A memory leak was fixed in two of the ``IDASp***Free`` functions.
In the rootfinding functions ``IDARcheck1`` and ``IDARcheck2``, when an exact
zero is found, the array ``glo`` of ``g`` values at the left endpoint is
adjusted, instead of shifting the ``t`` location ``tlo`` slightly. In the
installation files, we modified the treatment of the macro
``SUNDIALS_USE_GENERIC_MATH``, so that the parameter ``GENERIC_MATH_LIB`` is
either defined (with no value) or not defined.

### KINSOL Changes in v2.7.0

One significant design change was made with this release: The problem size and its relatives, bandwidth parameters,
related internal indices, pivot arrays, and the optional output ``lsflag`` have all been changed from type ``int`` to
type ``long int``, except for the problem size and bandwidths in user calls to routines specifying BLAS/LAPACK routines
for the dense/band linear solvers. The function ``NewIntArray`` is replaced by a pair ``NewIntArray``/``NewLintArray``,
for ``int`` and ``long int`` arrays, respectively.

A large number of errors have been fixed. Three major logic bugs were fixed - involving updating the solution vector,
updating the linesearch parameter, and a missing error return. Three minor errors were fixed - involving setting
``etachoice`` in the Matlab/KINSOL interface, a missing error case in ``KINPrintInfo``, and avoiding an
exponential overflow in the evaluation of ``omega``. In each linear solver interface function, the linear solver memory
is freed on an error return, and the ``**Free`` function now includes a line setting to NULL the main memory pointer to
the linear solver memory. In the installation files, we modified the treatment of the macro SUNDIALS_USE_GENERIC_MATH,
so that the parameter GENERIC_MATH_LIB is either defined (with no value) or not defined.


## Changes to SUNDIALS in release v2.4.0

- New features
  - new linear solver module, based on Blas and Lapack for both dense and banded matrices.
- Changes to user interface
  - reorganization of all linear solver modules into two families (besides the existing family of scaled preconditioned iterative linear solvers, the direct solvers, including the new Lapack-based ones, were also organized into a direct family).
- Changes related to the build system
  - provide CMake-based build option, in addition to that based on autotools.

### CVODE Changes in v2.6.0

Two new features were added in this release: (a) a new linear solver
module, based on BLAS and LAPACK for both dense and banded matrices, and
(b) an option to specify which direction of zero-crossing is to be
monitored while performing rootfinding.

The user interface has been further refined. Some of the API changes
involve: (a) a reorganization of all linear solver modules into two
families (besides the existing family of scaled preconditioned iterative
linear solvers, the direct solvers, including the new LAPACK-based ones,
were also organized into a *direct* family); (b) maintaining a single
pointer to user data, optionally specified through a -type function; and
(c) a general streamlining of the preconditioner modules distributed
with the solver.

### CVODES Changes in v2.6.0

Two new features related to the integration of ODE IVP problems were added in
this release: (a) a new linear solver module, based on BLAS and LAPACK for both
dense and banded matrices, and (b) an option to specify which direction of
zero-crossing is to be monitored while performing rootfinding.

This version also includes several new features related to sensitivity analysis,
among which are: (a) support for integration of quadrature equations depending
on both the states and forward sensitivity (and thus support for forward
sensitivity analysis of quadrature equations), (b) support for simultaneous
integration of multiple backward problems based on the same underlying ODE
(e.g., for use in an *forward-over-adjoint* method for computing second order
derivative information), (c) support for backward integration of ODEs and
quadratures depending on both forward states and sensitivities (e.g., for use in
computing second-order derivative information), and (d) support for
reinitialization of the adjoint module.

The user interface has been further refined. Some of the API changes involve:
(a) a reorganization of all linear solver modules into two families (besides the
existing family of scaled preconditioned iterative linear solvers, the direct
solvers, including the new LAPACK-based ones, were also organized into a
*direct* family); (b) maintaining a single pointer to user data, optionally
specified through a ``Set``-type function; and (c) a general streamlining of the
preconditioner modules distributed with the solver. Moreover, the prototypes of
all functions related to integration of backward problems were modified to
support the simultaneous integration of multiple problems. All backward problems
defined by the user are internally managed through a linked list and identified
in the user interface through a unique identifier.

### IDA Changes in v2.6.0

Two new features were added in this release: (a) a new linear solver module,
based on BLAS and LAPACK for both dense and banded matrices, and (b) option to
specify which direction of zero-crossing is to be monitored while performing
rootfinding.

The user interface has been further refined. Some of the API changes involve:
(a) a reorganization of all linear solver modules into two families (besides the
already present family of scaled preconditioned iterative linear solvers, the
direct solvers, including the new LAPACK-based ones, were also organized into a
*direct* family); (b) maintaining a single pointer to user data, optionally
specified through a ``Set``-type function; (c) a general streamlining of the
band-block-diagonal preconditioner module distributed with the solver.

### KINSOL Changes in v2.6.0

This release introduces a new linear solver module, based on BLAS and LAPACK for both dense and banded matrices.

The user interface has been further refined. Some of the API changes involve: (a) a reorganization of all linear solver
modules into two families (besides the already present family of scaled preconditioned iterative linear solvers, the
direct solvers, including the new LAPACK-based ones, were also organized into a *direct* family); (b) maintaining a
single pointer to user data, optionally specified through a ``Set``-type function; (c) a general streamlining of the
band-block-diagonal preconditioner module distributed with the solver.


## Changes to SUNDIALS in release v2.3.0

- Changes to the user interface
  - modified the functions in the generic dense linear solver (sundials\_dense and sundials\_smalldense) to work for rectangular m by n matrices (m <= n).
  - renamed the factorization and solution functions in the generic dense linear solver to DenseGETRF/denGETRF and DenseGETRS/denGETRS, respectively.
  - renamed the factorization and solution functions in the generic band linear solver to BandGBTRF and BandGBTRS, respectively.
- Changes related to the build system
  - rearranged the entire SUNDIALS source tree
  - all exported header files are now installed in separate subdirectories of the installation include directory
  - header files are included now by specifying the relative path (e.g., #include \<sundials/sundials\_types.h\>)

### CVODE Changes in v2.5.0

The main changes in this release involve a rearrangement of the entire
:ref:`SUNDIALS source tree <CVODE.Organization>`. At the user interface
level, the main impact is in the mechanism of including SUNDIALS header files
which must now include the relative path (e.g. ``#include <cvode/cvode.h>``). Additional changes were made to
the build system: all exported header files are now installed in separate
subdirectories of the instaltion *include* directory.

The functions in the generic dense linear solver (``sundials_dense`` and
``sundials_smalldense``) were modified to work for rectangular :math:`m \times
n` matrices (:math:`m \le n`), while the factorization and solution functions
were renamed to ``DenseGETRF`` / ``denGETRF`` and ``DenseGETRS`` / ``denGETRS``,
respectively. The factorization and solution functions in the generic band
linear solver were renamed ``BandGBTRF`` and ``BandGBTRS``, respectively.

### CVODES Changes in v2.5.0

The main changes in this release involve a rearrangement of the entire SUNDIALS
source tree (see :numref:`CVODES.Organization`). At the user interface level,
the main impact is in the mechanism of including SUNDIALS header files which
must now include the relative path (e.g. ``#include <cvode/cvode.h>``).
Additional changes were made to the build system: all exported header files are
now installed in separate subdirectories of the instaltion *include* directory.

In the adjoint solver module, the following two bugs were fixed: in ``CVodeF``
the solver was sometimes incorrectly taking an additional step before returning
control to the user (in ``CV_NORMAL`` mode) thus leading to a failure in the
interpolated output function; in ``CVodeB``, while searching for the current
check point, the solver was sometimes reaching outside the integration interval
resulting in a segmentation fault.

The functions in the generic dense linear solver (``sundials_dense`` and
``sundials_smalldense``) were modified to work for rectangular :math:`m \times
n` matrices (:math:`m \le n`), while the factorization and solution functions
were renamed to ``DenseGETRF`` / ``denGETRF`` and ``DenseGETRS`` / ``denGETRS``,
respectively. The factorization and solution functions in the generic band
linear solver were renamed ``BandGBTRF`` and ``BandGBTRS``, respectively.

### IDA Changes in v2.5.0

The main changes in this release involve a rearrangement of the entire SUNDIALS
source tree (see :numref:`IDA.Organization`). At the user interface level, the main
impact is in the mechanism of including SUNDIALS header files which must now
include the relative path (e.g. ``#include <cvode/cvode.h>``).  Additional
changes were made to the build system: all exported header files are now
installed in separate subdirectories of the installation *include* directory.

A bug was fixed in the internal difference-quotient dense and banded Jacobian
approximations, related to the estimation of the perturbation (which could have
led to a failure of the linear solver when zero components with sufficiently
small absolute tolerances were present).

The user interface to the consistent initial conditions calculations was
modified.  The :c:func:`IDACalcIC` arguments ``t0``, ``yy0``, and ``yp0`` were
removed and a new function, :c:func:`IDAGetConsistentIC` is provided.

The functions in the generic dense linear solver (``sundials_dense`` and
``sundials_smalldense``) were modified to work for rectangular :math:`m \times
n` matrices (:math:`m \le n`), while the factorization and solution functions
were renamed to ``DenseGETRF / denGETRF`` and ``DenseGETRS / denGETRS``,
respectively.  The factorization and solution functions in the generic band
linear solver were renamed ``BandGBTRF`` and ``BandGBTRS``, respectively.

### KINSOL Changes in v2.5.0

The main changes in this release involve a rearrangement of the entire SUNDIALS source tree (see
:numref:`KINSOL.Organization`). At the user interface level, the main impact is in the mechanism of including
SUNDIALS header files which must now include the relative path (e.g. ``#include <cvode/cvode.h>``). Additional
changes were made to the build system: all exported header files are now installed in separate subdirectories of the
installation *include* directory.

The functions in the generic dense linear solver (``sundials_dense`` and ``sundials_smalldense``) were modified to work
for rectangular :math:`m \times n` matrices (:math:`m \le n`), while the factorization and solution functions were
renamed to ``DenseGETRF``/``denGETRF`` and ``DenseGETRS``/``denGETRS``, respectively. The factorization and solution
functions in the generic band linear solver were renamed ``BandGBTRF`` and ``BandGBTRS``, respectively.


## Changes to SUNDIALS in release v2.2.0

- New features
  - added SPBCG (scaled preconditioned Bi-CGStab) linear solver module
  - added SPTFQMR (scaled preconditioned TFQMR) linear solver module
- Changes related to the build system
  - updated configure script and Makefiles for Fortran examples to avoid C++ compiler errors (now use CC and MPICC to link only if necessary)
  - SUNDIALS shared header files are installed under a sundials subdirectory of the installation include directory
  - the shared object files are now linked into each SUNDIALS library rather than into a separate libsundials\_shared library
- Changes to the user interface
  - added prefix sundials\_ to all shared header files

### CVODE Changes in v2.4.0

CVSPBCG and CVSPTFQMR modules have been added to interface with
the Scaled Preconditioned Bi-CGstab (SPBCG) and Scaled Preconditioned
Transpose-Free Quasi-Minimal Residual (SPTFQMR) linear solver modules,
respectively (for details see :numref:`CVODE.Usage.CC`). Corresponding additions were made
to the Fortran interface module FCVODE. At the same time, function type
names for Scaled Preconditioned Iterative Linear Solvers were added for
the user-supplied Jacobian-times-vector and preconditioner setup and
solve functions.

The deallocation functions now take as arguments the address of the
respective memory block pointer.

To reduce the possibility of conflicts, the names of all header files
have been changed by adding unique prefixes (``cvode_`` and ``sundials_``). When using the
default installation procedure, the header files are exported under
various subdirectories of the target directory. For more details see
:numref:`Installation`.

### CVODES Changes in v2.4.0

CVSPBCG and CVSPTFQMR modules have been added to interface with the Scaled
Preconditioned Bi-CGstab (SPBCG) and Scaled Preconditioned Transpose-Free
Quasi-Minimal Residual (SPTFQMR) linear solver modules, respectively (for
details see Chapter :numref:`CVODES.Usage.SIM`). At the same time,
function type names for Scaled Preconditioned Iterative Linear Solvers were
added for the user-supplied Jacobian-times-vector and preconditioner setup and
solve functions.

A new interpolation method was added to the CVODES adjoint module. The function
``CVadjMalloc`` has an additional argument which can be used to select the
desired interpolation scheme.

The deallocation functions now take as arguments the address of the respective
memory block pointer.

To reduce the possibility of conflicts, the names of all header files have been
changed by adding unique prefixes (``cvodes_`` and ``sundials_``). When using
the default installation procedure, the header files are exported under various
subdirectories of the target ``include`` directory. For more details see
Appendix :numref:`Installation`.

### IDA Changes in v2.4.0

FIDA, a Fortran-C interface module, was added.

IDASPBCG and IDASPTFQMR modules have been added to interface with the Scaled
Preconditioned Bi-CGstab (SPBCG) and Scaled Preconditioned Transpose-Free
Quasi-Minimal Residual (SPTFQMR) linear solver modules, respectively (for
details see :numref:IDA.Usage.CC).  At the same time, function type names for Scaled
Preconditioned Iterative Linear Solvers were added for the user-supplied
Jacobian-times-vector and preconditioner setup and solve functions.

The rootfinding feature was added, whereby the roots of a set of given functions
may be computed during the integration of the DAE system.

A user-callable routine was added to access the estimated local error vector.

The deallocation functions now take as arguments the address of the respective
memory block pointer.

To reduce the possibility of conflicts, the names of all header files have been
changed by adding unique prefixes (``ida_`` and ``sundials_``).  When using the
default installation procedure, the header files are exported under various
subdirectories of the target ``include`` directory. For more details see
Appendix :numref:`Installation`.

### KINSOL Changes in v2.4.0

KINSPBCG, KINSPTFQMR, KINDENSE, and KINBAND modules have been added to interface with the Scaled
Preconditioned Bi-CGStab (SPBCG), Scaled Preconditioned Transpose-Free Quasi-Minimal Residual (SPTFQMR),
DENSE, and BAND linear solver modules, respectively. (For details see Chapter :numref:KINSOL.Usage.CC.)
Corresponding additions were made to the Fortran interface module FKINSOL. At the same time, function type names
for Scaled Preconditioned Iterative Linear Solvers were added for the user-supplied Jacobian-times-vector and
preconditioner setup and solve functions.

Regarding the Fortran interface module FKINSOL, optional inputs are now set using ``FKINSETIIN`` (integer inputs),
``FKINSETRIN`` (real inputs), and ``FKINSETVIN`` (vector inputs). Optional outputs are still obtained from the ``IOUT``
and ``ROUT`` arrays which are owned by the user and passed as arguments to ``FKINMALLOC``.

The KINDENSE and KINBAND linear solver modules include support for nonlinear residual monitoring which can
be used to control Jacobian updating.

To reduce the possibility of conflicts, the names of all header files have been changed by adding unique prefixes
(``kinsol_`` and ``sundials_``). When using the default installation procedure, the header files are exported under
various subdirectories of the target ``include`` directory. For more details see Appendix :numref:`Installation`.


## Changes to SUNDIALS in release v2.1.1

- Changes to the generic NVECTOR module
  - N\_VCloneEmpty was added to the global vector operations table

### CVODE Changes in v2.3.0

None?

### CVODES Changes in v2.3.0

A minor bug was fixed in the interpolation functions of the adjoint CVODES
module.

## Changes to SUNDIALS in release v2.1.0

### CVODE Changes in v2.3.0

The user interface has been further refined. Several functions used for
setting optional inputs were combined into a single one. An optional
user-supplied routine for setting the error weight vector was added.
Additionally, to resolve potential variable scope issues, all SUNDIALS
solvers release user data right after its use. The build systems has
been further improved to make it more robust.

### CVODES Changes in v2.2.0

The user interface has been further refined. Several functions used for setting
optional inputs were combined into a single one. An optional user-supplied
routine for setting the error weight vector was added. Additionally, to resolve
potential variable scope issues, all SUNDIALS solvers release user data right
after its use. The build systems has been further improved to make it more
robust.

### IDA Changes in v2.3.0

The user interface has been further refined. Several functions used for setting
optional inputs were combined into a single one.  An optional user-supplied
routine for setting the error weight vector was added.  Additionally, to resolve
potential variable scope issues, all SUNDIALS solvers release user data right
after its use. The build systems has been further improved to make it more
robust.

### KINSOL Changes in v2.3.0

The user interface has been further refined. Several functions used for setting optional inputs were combined into a
single one. Additionally, to resolve potential variable scope issues, all SUNDIALS solvers release user data right
after its use. The build system has been further improved to make it more robust.


## Changes to SUNDIALS in release v2.0.2

- Changes related to the build system
  - fixed autoconf-related bug to allow configuration with the PGI Fortran compiler
  - modified to use customized detection of the Fortran name mangling scheme (autoconf's AC\_F77\_WRAPPERS routine is problematic on some platforms)

### CVODE Changes in v2.2.2

None?

### CVODES Changes in v2.1.2

A bug was fixed in the ``CVode`` function that was potentially leading to
erroneous behaviour of the rootfinding procedure on the integration first step.

### IDA Changes in v2.2.2

Minor corrections and improvements were made to the build system.  A new chapter
in the User Guide was added - with constants that appear in the user interface.

## Changes to SUNDIALS in release v2.0.1

- Changes related to the build system
  - changed order of compiler directives in header files to avoid compilation errors when using a C++ compiler.
  - changed method of generating sundials\_config.h to avoid potential warnings of redefinition of preprocessor symbols.

### CVODE Changes in v2.2.1

The changes in this minor SUNDIALS release affect only the build system.

### CVODES Changes in v2.1.1

This CVODES release includes bug fixes related to forward sensitivity
computations (possible loss of accuray on a BDF order increase and incorrect
logic in testing user-supplied absolute tolerances). In addition, we have added
the option of activating and deactivating forward sensitivity calculations on
successive CVODES runs without memory allocation/deallocation.

Other changes in this minor SUNDIALS release affect the build system.

### IDA Changes in v2.2.1

The changes in this minor SUNDIALS release affect only the build system.

### KINSOL Changes in v2.2.1

The changes in this minor SUNDIALS release affect only the build system.


## Changes to SUNDIALS in release v2.0.0

- Changes to the generic NVECTOR module
  - removed machEnv, redefined table of vector operations (now contained in the N\_Vector structure itself).
  - all SUNDIALS functions create new N\_Vector variables through cloning, using an N\_Vector passed by the user as a template.
  - a particular NVECTOR implementation is supposed to provide user-callable constructor and destructor functions.
  - removed from structure of vector operations the following functions: N\_VNew, N\_VNew\_S, N\_VFree, N\_VFree\_S, N\_VMake, N\_VDispose, N\_VGetData, N\_VSetData, N\_VConstrProdPos, and N\_VOneMask.
  - added in structure of vector operations the following functions: N\_VClone, N\_VDestroy, N\_VSpace, N\_VGetArrayPointer, N\_VSetArrayPointer, and N\_VWrmsNormMask.
  - Note that nvec\_ser and nvec\_par are now separate modules outside the shared SUNDIALS module.
- Changes to the generic linear solvers
  - in SPGMR, added a dummy N\_Vector argument to be used as a template for cloning.
  - in SPGMR, removed N (problem dimension) from argument list of SpgmrMalloc.
  - iterative.{c,h} replace iterativ.{c,h}
  - modified constant names in iterative.h (preconditioner types are prefixed with 'PREC\_').
  - changed numerical values for MODIFIED\_GS (from 0 to 1) and CLASSICAL\_GS (from 1 to 2).
- Changes to sundialsmath submodule
  - replaced internal routine for estimation of unit roundoff with definition of unit roundoff from float.h
  - modified functions to call appropriate math routines given the precision level specified by the user.
- Changes to sundialstypes submodule
  - removed type 'integertype'.
  - added definitions for 'BIG\_REAL', 'SMALL\_REAL', and 'UNIT\_ROUNDOFF' using values from float.h based on the precision.
  - changed definition of macro RCONST to depend on precision.

### CVODE Changes in v2.2.0

The major changes from the previous version involve a redesign of the user
interface across the entire SUNDIALS suite. We have eliminated the mechanism of
providing optional inputs and extracting optional statistics from the solver
through the `iopt` and `ropt` arrays. Instead, CVODE now provides a set of
routines (with prefix ``CVodeSet``) to change the default values for various
quantities controlling the solver and a set of extraction routines (with prefix
``CVodeGet``) to extract statistics after return from the main solver
routine. Similarly, each linear solver module provides its own set of `Set`- and
`Get`-type routines. For more details see
:numref:`CVODE.Usage.CC.optional_input` and
:numref:`CVODE.Usage.CC.optional_output`.

Additionally, the interfaces to several user-supplied routines (such as those
providing Jacobians and preconditioner information) were simplified by reducing
the number of arguments. The same information that was previously accessible
through such arguments can now be obtained through `Get`-type functions.

The rootfinding feature was added, whereby the roots of a set of given functions
may be computed during the integration of the ODE system.

Installation of CVODE (and all of SUNDIALS) has been completely redesigned and
is now based on configure scripts.

### CVODES Changes in v2.1.0

The major changes from the previous version involve a redesign of the user
interface across the entire SUNDIALS suite. We have eliminated the mechanism of
providing optional inputs and extracting optional statistics from the solver
through the ``iopt`` and ``ropt`` arrays. Instead, CVODES now provides a set of
routines (with prefix ``CVodeSet``) to change the default values for various
quantities controlling the solver and a set of extraction routines (with prefix
``CVodeGet``) to extract statistics after return from the main solver routine.
Similarly, each linear solver module provides its own set of ``Set``- and
``Get``-type routines. For more details see
:numref:`CVODES.Usage.SIM.optional_input` and
:numref:`CVODES.Usage.SIM.optional_output`.

Additionally, the interfaces to several user-supplied routines (such as those
providing Jacobians, preconditioner information, and sensitivity right hand
sides) were simplified by reducing the number of arguments. The same information
that was previously accessible through such arguments can now be obtained
through ``Get``-type functions.

The rootfinding feature was added, whereby the roots of a set of given functions
may be computed during the integration of the ODE system.

Installation of CVODES (and all of SUNDIALS) has been completely redesigned and
is now based on configure scripts.

### IDA Changes in v2.2.0

The major changes from the previous version involve a redesign of the user
interface across the entire SUNDIALS suite. We have eliminated the mechanism of
providing optional inputs and extracting optional statistics from the solver
through the ``iopt`` and ``ropt`` arrays. Instead, IDA now provides a set of
routines (with prefix ``IDASet``) to change the default values for various
quantities controlling the solver and a set of extraction routines (with prefix
``IDAGet``) to extract statistics after return from the main solver routine.
Similarly, each linear solver module provides its own set of ``Set``- and
``Get``-type routines. For more details see :numref:`IDA.Usage.CC.optional_output`.

Additionally, the interfaces to several user-supplied routines (such as those
providing Jacobians and preconditioner information) were simplified by reducing
the number of arguments. The same information that was previously accessible
through such arguments can now be obtained through ``Get``-type functions.

Installation of IDA (and all of SUNDIALS) has been completely redesigned and is
now based on configure scripts.

### KINSOL Changes in v2.2.0

The major changes from the previous version involve a redesign of the user interface across the entire SUNDIALS
suite. We have eliminated the mechanism of providing optional inputs and extracting optional statistics from the solver
through the ``iopt`` and ``ropt`` arrays. Instead, KINSOL now provides a set of routines (with prefix ``KINSet``)
to change the default values for various quantities controlling the solver and a set of extraction routines (with prefix
``KINGet``) to extract statistics after return from the main solver routine. Similarly, each linear solver module
provides its own set of ``Set``- and ``Get``-type routines. For more details see Chapter :numref:KINSOL.Usage.CC.

Additionally, the interfaces to several user-supplied routines (such as those providing Jacobian-vector products and
preconditioner information) were simplified by reducing the number of arguments. The same information that was
previously accessible through such arguments can now be obtained through ``Get``-type functions.

Installation of KINSOL (and all of SUNDIALS) has been completely redesigned and is now based on configure
scripts.



# sundialsTB

sundialsTB is no longer distributed as of sundials v. 2.7.0 as it has not been updated in many years.

## What's new in v2.5.0?

- Bug fixes
  - fixed lines setting etachoice in kimOpts.c
  - in cvm.c and idm.c, fixed size of rootsfound array; added lines to free rootsfound and ckpnt arrays when done using each
- What's new in v2.4.0?
- New Features
  - the Matlab interface to IDAS was extended to provide sensitivity analysis capabilities.
- Changes to user interface
  - the API for adjoint sensitivity analysis (cvodes and idas) was modified to support simultaneous integration of multiple backward problems.

## What's new in v2.3.0?

- New features
  - added Matlab interface to IDA (named idas)
  - on platforms which support configure scripts, installation of sundialsTB can now be enabled while configuring SUNDIALS and installed through make and make install (provided a working MEX compiler is found).
- Bug fixes
  - the installation script install\_STB.m was modified to increase robustness on various platforms (related to path and file names).
- Changes to user interface
  - (cvodes) for improved legibility, some of the keys for forward sensitivity optional inputs were renamed.
  - (cvodes) removed xaxis type option for the internal monitoring function CVodeMonitor.

## What's new in v2.2.0?

- New features
  - modified installation procedure to use a Matlab script
  - added sample Matlab startup file
  - (cvodes) expanded CVodeMonitor
  - (kinsol) added interface to KINSOL's performance monitoring function ('Verbose' option to KINSetOptions)
- Bug fixes
  - (cvodes) fixed bug in interface to quadrature integration which was causing a segmentation violation when monitoring was turned on.
- Changes to user interface
  - updated to reflect changes to the SUNDIALS libraries in v.2.2.0
  - (cvodes) changed the interface for sensitivity analysis (both forward and adjoint) to follow more closely the CVODES calling sequence
  - (cvodes) optional inputs for forward sensitivity analysis are now provided through a separate function, CVodeSensSetOptions
  - removed NVM mex interface
