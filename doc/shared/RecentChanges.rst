.. For package-specific references use :ref: rather than :numref: so intersphinx
   links to the appropriate place on read the docs

**Major Features**

**New Features and Enhancements**

The Kokkos N_Vector is now a C++17/Kokkos-4-or-newer same-rank wrapper for
application-owned views (Kokkos 5 requires C++20). Allocation and automatic
host-mirror helpers were removed; use explicit mirrors and ``VectorFor::FromNVector``. Rank 1 uses
``RangePolicy``, ranks 2--6 use native ``MDRangePolicy``, and ranks 7 and
higher are rejected during construction.

Added the public C++ ``sundials::NVectorBase`` interface for implementing
custom NVECTOR types.

Added the utility function, :c:func:`SUNFileFlush` for flushing file
pointers. This is useful when using the Fortran 2003 interfaces.

Added the function :c:func:`SUNContext_HandleError` for reporting an error
through a :c:type:`SUNContext`. It sets the last error and calls the error
handler stack, and can be called from user callbacks and C++ classes.

**Bug Fixes**

Fixed partial-clone cleanup in :c:func:`N_VCloneVectorArray` so all
successfully cloned vectors are destroyed when a later clone fails.

Fixed a segfault that could occur when handling errors without a ``SUNContext``
provided.

:c:func:`SUNLogErrHandlerFn`, :c:func:`SUNAbortErrHandlerFn`, and the package
warning messages no longer interpret ``%`` in an error message as a format
directive. Messages are printed as given.

Fixed a bug in ``FindMAGMA.cmake`` which didn't allow use of MAGMA versions with 
multiple digits in an identifier.

Fixed a bug in the LSRKStep module where the step post-processing function for
the RKL step was called with an incorrect time argument.

**Deprecation Notices**
