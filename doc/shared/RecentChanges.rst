.. For package-specific references use :ref: rather than :numref: so intersphinx
   links to the appropriate place on read the docs

**Major Features**

A new discrete adjoint capability for explicit Runge--Kutta methods has been
added to the ARKODE ERKStep and ARKStep stepper modules. This is based on a new
set of shared classes, :c:type:`SUNAdjointStepper` and
:c:type:`SUNAdjointCheckpointScheme`. A new example demonstrating this
capability can be found in
``examples/arkode/C_serial/ark_lotka_volterra_ASA.c``. See the
:ref:`ARKODE.Mathematics.ASA` section of the ARKODE user guide for details.

**New Features and Enhancements**

*ARKODE*

The following changes have been made to the default ERK, DIRK, and ARK methods
in ARKODE to utilize more efficient methods:

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

The old default methods can be loaded using the functions
:c:func:`ERKStepSetTableName` or :c:func:`ERKStepSetTableNum` with ERKStep and
:c:func:`ARKStepSetTableName` or :c:func:`ARKStepSetTableNum` with ARKStep and
passing the desired method name string or constant, respectively. For example,
the following call can be used to load the old default fourth order method with
ERKStep:

.. code-block:: C

   /* Load the old 4th order ERK method using the table name */
   ierr = ERKStepSetTableName(arkode_mem, "ARKODE_ZONNEVELD_5_3_4");

Similarly with ARKStep, the following calls can be used for ERK, DIRK, or ARK
methods, respectively:

.. code-block:: C

   /* Load the old 4th order ERK method by name */
   ierr = ARKStepSetTableName(arkode_mem, "ARKODE_DIRK_NONE",
                              "ARKODE_ZONNEVELD_5_3_4");

   /* Load the old 4th order DIRK method by name */
   ierr = ARKStepSetTableName(arkode_mem, "ARKODE_SDIRK_5_3_4",
                              "ARKODE_ERK_NONE");

   /* Load the old 4th order ARK method by name */
   ierr = ARKStepSetTableName(arkode_mem, "ARKODE_ARK436L2SA_DIRK_6_3_4",
                              "ARKODE_ARK436L2SA_ERK_6_3_4");

Additionally, the following changes have been made to the default time step
adaptivity parameters in ARKODE:

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

The following calls can be used to restore the old defaults for ERKStep:

.. code-block:: c

   SUNAdaptController controller = SUNAdaptController_Soderlind(ctx);
   SUNAdaptController_SetParams_PI(controller, 0.8, -0.31);
   ARKodeSetAdaptController(arkode_mem, controller);
   SUNAdaptController_SetErrorBias(controller, 1.2);
   ARKodeSetSafetyFactor(arkode_mem, 0.96);
   ARKodeSetFixedStepBounds(arkode_mem, 1, 1.5);
   ARKodeSetAdaptivityAdjustment(arkode_mem, -1);

The following calls can be used to restore the old defaults for other ARKODE
integrators:

.. code-block:: c

   SUNAdaptController controller = SUNAdaptController_PID(ctx);
   ARKodeSetAdaptController(arkode_mem, controller);
   SUNAdaptController_SetErrorBias(controller, 1.5);
   ARKodeSetSafetyFactor(arkode_mem, 0.96);
   ARKodeSetFixedStepBounds(arkode_mem, 1, 1.5);
   ARKodeSetAdaptivityAdjustment(arkode_mem, -1);

In both cases above, destroy the controller at the end of the run with
``SUNAdaptController_Destroy(controller);``.

The Soderlind time step adaptivity controller now starts with an I controller
until there is sufficient history of past time steps and errors.

Added :c:func:`ARKodeSetAdaptControllerByName` to set a time step adaptivity controller
with a string. There are also four new controllers:
:c:func:`SUNAdaptController_H0211`, :c:func:`SUNAdaptController_H0321`,
:c:func:`SUNAdaptController_H211`, and :c:func:`SUNAdaptController_H312`.

Added the ``ARKODE_RALSTON_3_1_2`` and ``ARKODE_TSITOURAS_7_4_5`` explicit
Runge-Kutta Butcher tables.

Improved the precision of the coefficients for ``ARKODE_ARK324L2SA_ERK_4_2_3``,
``ARKODE_VERNER_9_5_6``, ``ARKODE_VERNER_10_6_7``, ``ARKODE_VERNER_13_7_8``,
``ARKODE_ARK324L2SA_DIRK_4_2_3``, and ``ARKODE_ESDIRK324L2SA_4_2_3``.

*CVODE / CVODES*

Added support for resizing CVODE and CVODES when solving initial value problems
where the number of equations and unknowns changes over time. Resizing requires
a user supplied history of solution and right-hand side values at the new
problem size, see :c:func:`CVodeResizeHistory` for more information.

*KINSOL*

Added support in KINSOL for setting user-supplied functions to compute the
damping factor and, when using Anderson acceleration, the depth in fixed-point
or Picard iterations. See :c:func:`KINSetDampingFn` and :c:func:`KINSetDepthFn`,
respectively, for more information.

*SUNDIALS Types*

A new type, :c:type:`suncountertype`, was added for the integer type used for
counter variables. It is currently an alias for ``long int``.

**Bug Fixes**

*ARKODE*

Fixed bug in :c:func:`ARKodeResize` which caused it return an error for MRI
methods.

Removed error floors from the :c:type:`SUNAdaptController` implementations
which could unnecessarily limit the time size growth, particularly after the
first step.

Fixed bug in :c:func:`ARKodeSetFixedStep` where it could return ``ARK_SUCCESS``
despite an error occurring.

Fixed bug in the ARKODE SPRKStep :c:func:`SPRKStepReInit` function and
:c:func:`ARKodeReset` function with SPRKStep that could cause a segmentation
fault when compensated summation is not used.

*KINSOL*

Fixed a bug in KINSOL where an incorrect damping parameter is applied on the
initial iteration with Anderson acceleration unless :c:func:`KINSetDamping` and
:c:func:`KINSetDampingAA` are both called with the same value when enabling
damping.

Fixed a bug in KINSOL where errors that occurred when computing Anderson
acceleration were not captured.

Added missing return values to :c:func:`KINGetReturnFlagName`.

*CMake*

Fixed the behavior of :cmakeop:`SUNDIALS_ENABLE_ERROR_CHECKS` so additional
runtime error checks are disabled by default with all release build types.
Previously, ``MinSizeRel`` builds enabled additional error checking by default.

**Deprecation Notices**

All work space functions, e.g., ``CVodeGetWorkSpace`` and
``ARKodeGetLinWorkSpace``, have been deprecated and will be removed in version
8.0.0.
