**Major Features**

**New Features and Enhancements**

The following changes have been made to the default time step adaptivity
parameters in ARKODE:

+-----------------------+-----------------------+-------------+
| Parameter             | Old Default           | New Default |
+=======================+=======================+=============+
| Controller            | PID (PI for ERKStep)  | I           |
+-----------------------+-----------------------+-------------+
| Safety Factor         | 0.96                  | 0.9         |
+-----------------------+-----------------------+-------------+
| Bias                  | 1.5 (1.2 for ERKStep) | 1.0         |
+-----------------------+-----------------------+-------------+
| Fixed Step Bounds     | [1.0, 1.5]            | [1.0, 1.0]  |
+-----------------------+-----------------------+-------------+
| Adaptivity Adjustment | -1                    | 0           |
+-----------------------+-----------------------+-------------+

The following can be used to restore the old defaults for ERKStep:

.. code-block:: c

  arkode_mem = ...
  SUNAdaptController controller = SUNAdaptController_Soderlind(ctx);
  SUNAdaptController_SetParams_PI(controller, 0.8, -0.31);
  ARKodeSetAdaptController(arkode_mem, controller);
  SUNAdaptController_SetErrorBias(controller, 1.2);
  ARKodeSetSafetyFactor(arkode_mem, 0.96);
  ARKodeSetFixedStepBounds(arkode_mem, 1, 1.5);
  ARKodeSetAdaptivityAdjustment(arkode_mem, -1);
  ...
  SUNAdaptController_Destroy(controller);

The following can be used to restore the old defaults for other ARKODE
integrators:

.. code-block:: c

  arkode_mem = ...
  SUNAdaptController controller = SUNAdaptController_PID(ctx);
  ARKodeSetAdaptController(arkode_mem, controller);
  SUNAdaptController_SetErrorBias(controller, 1.5);
  ARKodeSetSafetyFactor(arkode_mem, 0.96);
  ARKodeSetFixedStepBounds(arkode_mem, 1, 1.5);
  ARKodeSetAdaptivityAdjustment(arkode_mem, -1);
  ...
  SUNAdaptController_Destroy(controller);

Added :c:func:`ARKodeSetAdaptControllerByName` to set a time step adaptivity controller
by a string. There are also three new controllers:
:c:func:`SUNAdaptController_H0211`, :c:func:`SUNAdaptController_H211`, and
:c:func:`SUNAdaptController_H312`.

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

**Deprecation Notices**

All work space functions, e.g., ``CVodeGetWorkSpace`` and
``ARKodeGetLinWorkSpace``, have been deprecated and will be removed in version
8.0.0.
