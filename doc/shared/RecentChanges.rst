**Major Features**

**New Features and Enhancements**

Improved the precision of the coefficients for ``ARKODE_ARK324L2SA_ERK_4_2_3``,
``ARKODE_VERNER_9_5_6``, ``ARKODE_VERNER_10_6_7``, ``ARKODE_VERNER_13_7_8``,
``ARKODE_ARK324L2SA_DIRK_4_2_3``, and ``ARKODE_ESDIRK324L2SA_4_2_3``.

The Soderlind time step adaptivity controller now starts with an I controller
until there is sufficient history of past time steps and errors.

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

Added support in KINSOL for setting user-supplied functions to compute the
damping factor and depth in fixed-point or Picard iterations. See
:c:func:`KINSetDampingFn` and :c:func:`KINSetDepthFn`, respectively, for more
information.

**Bug Fixes**

Fixed bug in :c:func:`ARKodeResize` which caused it return an error for MRI
methods.

Removed error floors from the :c:type:`SUNAdaptController` implementations
which could unnecessarily limit the time size growth, particularly after the
first step.

Fixed bug in :c:func:`ARKodeSetFixedStep` where it could return ``ARK_SUCCESS``
despite an error occurring.

Fixed the behavior of :cmakeop:`SUNDIALS_ENABLE_ERROR_CHECKS` so additional
runtime error checks are disabled by default with all release build types.
Previously, ``MinSizeRel`` builds enabled additional error checking by default.

Fixed bug in the ARKODE SPRKStep :c:func:`SPRKStepReInit` function and
:c:func:`ARKodeReset` function with SPRKStep that could cause a segmentation
fault when compensated summation is not used.

Fixed a bug in KINSOL where an incorrect damping parameter is applied on the
initial iteration with Anderson acceleration unless :c:func:`KINSetDamping` and
:c:func:`KINSetDampingAA` are both called with the same value when enabling
damping.

Fixed a bug in KINSOL where errors that occurred when computing Anderson
acceleration were not captured.

Added missing return values to :c:func:`KINGetReturnFlagName`.

**Deprecation Notices**

All work space functions, e.g., ``CVodeGetWorkSpace`` and
``ARKodeGetLinWorkSpace``, have been deprecated and will be removed in version
8.0.0.
