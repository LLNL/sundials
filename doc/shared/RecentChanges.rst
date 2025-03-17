**Major Features**

**New Features and Enhancements**

Added the ``ARKODE_RALSTON_3_1_2`` and ``ARKODE_TSITOURAS_7_4_5`` explicit
Runge-Kutta Butcher tables.

Improved the efficiency of default ARKODE methods with the following changes:

+--------------------+-------------------------------------+--------------------------------------+
| Type               | Old Default                         | New Default                          |
+====================+=====================================+======================================+
| 2nd Order Explicit | ``ARKODE_HEUN_EULER_2_1_2``         | ``ARKODE_RALSTON_3_1_2``             |
+--------------------+-------------------------------------+--------------------------------------+
| 4th Order Explicit | ``ARKODE_ZONNEVELD_5_3_4``          | ``ARKODE_SOFRONIOU_SPALETTA_5_3_4``  |
+--------------------+-------------------------------------+--------------------------------------+
| 5th Order Explicit | ``ARKODE_CASH_KARP_6_4_5``          | ``ARKODE_TSITOURAS_7_4_5``           |
+--------------------+-------------------------------------+--------------------------------------+
| 6th Order Explicit | ``ARKODE_VERNER_8_5_6``             | ``ARKODE_VERNER_9_5_6``              |
+--------------------+-------------------------------------+--------------------------------------+
| 8th Order Explicit | ``ARKODE_FEHLBERG_13_7_8``          | ``ARKODE_VERNER_13_7_8``             |
+--------------------+-------------------------------------+--------------------------------------+
| 2nd Order Implicit | ``ARKODE_SDIRK_2_1_2``              | ``ARKODE_ARK2_DIRK_3_1_2``           |
+--------------------+-------------------------------------+--------------------------------------+
| 3rd Order Implicit | ``ARKODE_ARK324L2SA_DIRK_4_2_3``    | ``ARKODE_ESDIRK325L2SA_5_2_3``       |
+--------------------+-------------------------------------+--------------------------------------+
| 4th Order Implicit | ``ARKODE_SDIRK_5_3_4``              | ``ARKODE_ESDIRK436L2SA_6_3_4``       |
+--------------------+-------------------------------------+--------------------------------------+
| 5th Order Implicit | ``ARKODE_ARK548L2SA_DIRK_8_4_5``    | ``ARKODE_ESDIRK547L2SA2_7_4_5``      |
+--------------------+-------------------------------------+--------------------------------------+
| 4th Order ARK      | ``ARKODE_ARK436L2SA_ERK_6_3_4`` and | ``ARKODE_ARK437L2SA_ERK_7_3_4`` and  |
|                    | ``ARKODE_ARK436L2SA_DIRK_6_3_4``    | ``ARKODE_ARK437L2SA_DIRK_7_3_4``     |
+--------------------+-------------------------------------+--------------------------------------+
| 5th Order ARK      | ``ARKODE_ARK548L2SA_ERK_8_4_5`` and | ``ARKODE_ARK548L2SAb_ERK_8_4_5`` and |
|                    | ``ARKODE_ARK548L2SA_DIRK_8_4_5``    | ``ARKODE_ARK548L2SAb_DIRK_8_4_5``    |
+--------------------+-------------------------------------+--------------------------------------+

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
