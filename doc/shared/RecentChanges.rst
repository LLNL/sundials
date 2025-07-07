.. For package-specific references use :ref: rather than :numref: so intersphinx
   links to the appropriate place on read the docs

**Major Features**

**New Features and Enhancements**

**Bug Fixes**

The shared library version numbers for the oneMKL dense linear solver and
matrix as well as the PETSc SNES nonlinear solver have been corrected.

Fixed a CMake bug where the MRI H-Tol controller was not included in the ARKODE
Fortran module.

Fixed a bug in the CUDA and HIP implementations of ``SUNMemoryHelper_CopyAsync``
where the execution stream is not extracted correctly from the helper when a
stream is not provided to ``SUNMemoryHelper_CopyAsync``.

Fixed a bug in MRIStep where a segfault would occur when an MRI coupling table
is not explicitly set and an MRI integrator is nested inside another MRI
integrator.

Fixed a bug in MRIStep where MERK methods with unordered stage groups (MERK43
and MERK54) would include stage right-hand side vectors that had not been
computed yet in fast time scale forcing computations. These vectors were scaled
by zero, so in most cases the extraneous computations would not impact results.
However, in cases where these vectors contain ``inf`` or ``nan``, this would
lead to erroneous forcing terms.

**Deprecation Notices**
