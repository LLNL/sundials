**Major Features**

**New Features and Enhancements**

**Bug Fixes**

Fixed a bug in ARKStep where an extra right-hand side evaluation would occur
each time step when enabling the ``ARKodeSetAutonomous`` option and using an
IMEX method where the DIRK table has an implicit first and is not stiffly
accurate.

**Deprecation Notices**
