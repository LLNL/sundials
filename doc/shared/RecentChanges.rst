**Major Features**

**New Features and Enhancements**

TODO(SBR): description of order checking improvements

Improved the precision of the coefficients for ``ARKODE_ARK324L2SA_ERK_4_2_3``,
``ARKODE_VERNER_9_5_6``, ``ARKODE_VERNER_10_6_7``, ``ARKODE_VERNER_13_7_8``,
``ARKODE_ARK324L2SA_DIRK_4_2_3``, and ``ARKODE_ESDIRK324L2SA_4_2_3``.

**Bug Fixes**

Fixed bug in the ARKODE SPRKStep :c:func:`SPRKStepReInit` function and
:c:func:`ARKodeReset` function with SPRKStep that could cause a segmentation
fault when compensated summation is not used.

**Deprecation Notices**
