**Bug Fixes**

Fixed a bug in ARKStep where an extra right-hand side evaluation would occur
each time step when enabling the :c:func:`ARKodeSetAutonomous` option and using
an IMEX method where the DIRK table has an implicit first stage and is not stiffly
accurate.
