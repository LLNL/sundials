**New Features**

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
