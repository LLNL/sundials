.. For package-specific references use :ref: rather than :numref: so intersphinx
   links to the appropriate place on read the docs

**Major Features**

**New Features and Enhancements**

Added the utility function, :c:func:`SUNFileFlush` for flushing file
pointers. This is useful when using the Fortran 2003 interfaces.

Updated :c:func:`ARKodeGetReturnFlagName`, :c:func:`CVodeGetReturnFlagName`,
:c:func:`CVDiagGetReturnFlagName`, :c:func:`IDAGetReturnFlagName`, and
:c:func:`KINGetReturnFlagName` to return ``const char*`` rather than ``char*``.
The returned string should no longer be freed by the user.

**Bug Fixes**

Fixed a segfault that could occur when handling errors without a ``SUNContext``
provided.

Fixed a bug in ``FindMAGMA.cmake`` which didn't allow use of MAGMA versions with 
multiple digits in an identifier.

Fixed a bug in the LSRKStep module where the step post-processing function for
the RKL step was called with an incorrect time argument.

Fixed misleading KLU configuration errors. ``FindKLU.cmake`` now reports a
``KLU_INCLUDE_DIR`` without ``klu.h`` or library paths that do not exist instead
of reporting KLU as found, the 64-bit index check no longer reports a failure to
compile against ``klu.h`` as a ``SuiteSparse_long`` size mismatch, and a failed
check is no longer cached across reconfigures. See GitHub issue `#1016
<https://github.com/LLNL/sundials/issues/1016>`__.

**Deprecation Notices**
