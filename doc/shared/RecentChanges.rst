**Major Features**

**New Features and Enhancements**

Added Multirate time step adaptivity controllers, based on the recently introduced
`SUNAdaptController` base class, to ARKODE's MRIStep module.

Added functionality to ARKStep and ERKStep to accumulate a temporal error
estimate over multiple time steps.  See the routines :c:func:`ARKStepSetAccumulatedErrorType`,
:c:func:`ARKStepResetAccumulatedError`, :c:func:`ARKStepGetAccumulatedError`,
:c:func:`ERKStepSetAccumulatedErrorType`, :c:func:`ERKStepResetAccumulatedError`,
and :c:func:`ERKStepGetAccumulatedError` for details.

**Bug Fixes**

Fixed the loading of ARKStep's default first order explicit method.

**Deprecation Notices**
