.. For package-specific references use :ref: rather than :numref: so intersphinx
   links to the appropriate place on read the docs

**Major Features**

**New Features and Enhancements**

**Bug Fixes**

The shared library version numbers for the oneMKL dense linear solver and
matrix as well as the PETSc SNES nonlinear solver have been corrected.

Fixed a bug in the CUDA and HIP implementations of ``SUNMemoryHelper_CopyAsync``
where the execution stream is not extracted correctly from the helper when a
stream is not provided to ``SUNMemoryHelper_CopyAsync``.

**Deprecation Notices**
