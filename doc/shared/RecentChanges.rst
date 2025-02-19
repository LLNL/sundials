**Major Features**

**New Features and Enhancements**

The :c:func:`N_VLinearSum` and :c:func:`N_VLinearSumVectorArray` operations are no longer required to
support the use case where the output could be the same as the second input.  While all SUNDIALS
vector implementations still support this use case, the change facilitates user-supplied vectors
and external libraries.

**Bug Fixes**

Fixed bug in :c:func:`ARKodeSetFixedStep` where it could return ``ARK_SUCCESS``
despite an error occurring.

Fixed the behavior of :cmakeop:`SUNDIALS_ENABLE_ERROR_CHECKS` so additional
runtime error checks are disabled by default with all release build types.
Previously, ``MinSizeRel`` builds enabled additional error checking by default.

Fixed bug in the ARKODE SPRKStep :c:func:`SPRKStepReInit` function and
:c:func:`ARKodeReset` function with SPRKStep that could cause a segmentation
fault when compensated summation is not used.

**Deprecation Notices**
