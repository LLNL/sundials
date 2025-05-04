.. For package-specific references use :ref: rather than :numref: so intersphinx
   links to the appropriate place on read the docs

**Major Features**

**New Features and Enhancements**

The :c:func:`N_VLinearSum` and :c:func:`N_VLinearSumVectorArray` operations are no longer required to
support the use case where the output could be the same as the second input.  While all SUNDIALS
vector implementations still support this use case, the change facilitates user-supplied vectors
and external libraries.

**Bug Fixes**

Fixed a CMake bug that would cause the Caliper compile test to fail at configure time.

**Deprecation Notices**
