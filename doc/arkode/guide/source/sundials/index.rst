.. _SUNDIALS:

=========================================
Using SUNDIALS for C and C++ Applications
=========================================

As discussed in the section, :ref:`Organization`, the six solvers packages
(CVODE(S), IDA(S), ARKODE, KINSOL) that make up SUNDIALS are built upon common
classes/modules for vectors, matrices, and algebraic solvers. In addition, the
six packages all leverage some other common infrastructure, which we discuss
in this section.

.. toctree::
   SUNContext
   Profiling
