..
   Programmer(s): Daniel R. Reynolds @ SMU
   ----------------------------------------------------------------
   Copyright (c) 2017, Southern Methodist University.
   All rights reserved.
   For details, see the LICENSE file.
   ----------------------------------------------------------------

:tocdepth: 3


.. _SUNLinSol.Compatibility:

Compatibility of SUNLinearSolver modules
========================================

We note that not all ``SUNLinearSolver`` types are compatible with all
``SUNMatrix`` and ``N_Vector`` types provided with SUNDIALS.  In Table 
:ref:`SUNLinSol.linsol-matrix` we show the direct linear solvers
available as ``SUNLinearSolver`` modules, and the compatible matrix
implementations.  Recall that Table :ref:`CInterface.solver-vector` shows the
compatibility between all ``SUNLinearSolver`` modules and vector
implementations. 


.. _SUNLinSol.linsol-matrix:

Compatible SUNLinearSolver and SUNMatrix implementations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. cssclass:: table-bordered

================ ===== ====== ====== =============
Linear Solver    Dense Banded Sparse User Supplied
================ ===== ====== ====== =============
Dense            X                   X
Band                   X             X
LapackDense      X                   X
LapackBand             X             X
KLU                           X      X
SuperLU_MT                    X      X
User supplied    X     X      X      X
================ ===== ====== ====== =============


