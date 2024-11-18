**Major Features**

**New Features and Enhancements**

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

Added a time-stepping module to ARKODE for low storage Runge--Kutta methods, 
:ref:`LSRKStep <ARKODE.Usage.LSRKStep>`.  This currently supports five explicit low-storage 
methods: the second-order Runge--Kutta--Chebyshev and Runge--Kutta--Legendre methods, 
and the second- through fourth-order optimal strong stability preserving Runge--Kutta methods.  
All methods include embeddings for temporal adaptivity.

The Trilinos Tpetra NVector interface has been updated to utilize CMake
imported targets added in Trilinos 14 to improve support for different Kokkos
backends with Trilinos. As such, Trilinos 14 or newer is required and the
``Trilinos_INTERFACE_*`` CMake options have been removed.

Example programs using *hypre* have been updated to support v2.20 and newer.

The information level logging output in ARKODE, CVODE(S), and IDA(S) has been
updated to be more uniform across the packages and a new ``tools`` directory
has been added with a Python module, ``suntools``, containing utilities for
parsing logging output. The Python utilities for parsing CSV output have been
relocated from the ``scripts`` directory to the Python module.

The build system has been updated to utilize the CMake LAPACK imported target
which should ease building SUNDIALS with LAPACK libraries that require setting
specific linker flags e.g., MKL.

**Bug Fixes**

Fixed a `bug <https://github.com/LLNL/sundials/issues/581>`__ in the sparse
matrix implementation of :c:func:`SUNMatScaleAddI` which caused out of bounds
writes unless ``indexvals`` were in ascending order for each row/column.

Fixed :c:func:`ARKodeResize` not using the default ``hscale`` when an argument
of ``0`` was provided.

Fixed the loading of ARKStep's default first order explicit method.

Fixed a CMake bug regarding usage of missing "print_warning" macro
that was only triggered when the deprecated ``CUDA_ARCH`` option was used.

Fixed a memory leak that could occur if :c:func:`ARKodeSetDefaults` is called
repeatedly.

Fixed compilation errors when building the Trilinos Teptra NVector with CUDA
support.

Fixed a CMake configuration issue related to aliasing an ``ALIAS`` target when
using ``ENABLE_KLU=ON`` in combination with a static-only build of SuiteSparse.

**Deprecation Notices**

The ARKODE stepper specific functions to retrieve the number of right-hand side
function evaluations have been deprecated. Use :c:func:`ARKodeGetNumRhsEvals`
instead.
