** New Features and Enhancements **

Unit tests were separated from examples. To that end, the following directories 
were moved out of the ``examples/`` directory to the ``test/unit_tests`` directory:
``nvector``, ``sunmatrix``, ``sunlinsol``, and ``sunnonlinsol``.

**Bug Fixes**

Fixed a bug in ARKStep where an extra right-hand side evaluation would occur
each time step when enabling the :c:func:`ARKodeSetAutonomous` option and using
an IMEX method where the DIRK table has an implicit first stage and is not stiffly
accurate.
