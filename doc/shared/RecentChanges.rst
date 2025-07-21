.. For package-specific references use :ref: rather than :numref: so intersphinx
   links to the appropriate place on read the docs

**Major Features**

**New Features and Enhancements**

**Bug Fixes**

Fixed segfaults in :c:func:`CVodeAdjInit` and :c:func:`IDAAdjInit` when called
after adjoint memory has been freed.

Fixed a CMake bug that would cause the Caliper compile test to fail at configure time.

**Deprecation Notices**
