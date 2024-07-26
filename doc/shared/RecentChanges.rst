**Major Features**

**New Features and Enhancements**

Added Multirate time step adaptivity controllers, based on the recently introduced
`SUNAdaptController` base class, to ARKODE's MRIStep module.

Added functionality to ARKODE to accumulate a temporal error
estimate over multiple time steps.  See the routines :c:func:`ARKodeSetAccumulatedErrorType`,
:c:func:`ARKodeResetAccumulatedError`, and :c:func:`ARKodeGetAccumulatedError` for details.

**Bug Fixes**

Fixed the loading of ARKStep's default first order explicit method.

**Deprecation Notices**
