**Major Features**

**New Features and Enhancements**

**Bug Fixes**

Fixed the behavior of :cmakeop:`SUNDIALS_ENABLE_ERROR_CHECKS` so additional
runtime error checks are disabled by default with all release build types.
Previously, ``MinSizeRel`` builds enabled additional error checking by default.

Fixed bug in the ARKODE SPRKStep :c:func:`SPRKStepReInit` function and
:c:func:`ARKodeReset` function with SPRKStep that could cause a segmentation
fault when compensated summation is not used.

**Deprecation Notices**
