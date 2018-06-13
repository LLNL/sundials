..
   Programmer(s): Daniel R. Reynolds @ SMU
   ----------------------------------------------------------------
   Copyright (c) 2017, Southern Methodist University.
   All rights reserved.
   For details, see the LICENSE file.
   ----------------------------------------------------------------

:tocdepth: 3


.. _SUNMatrix.Compatibility:

Compatibility of SUNMATRIX types
======================================

We note that not all ``SUNMatrix`` types are compatible with all
``N_Vector`` types provided with SUNDIALS.  This is primarily due to 
the need for compatibility within the ``SUNMatMatvec`` routine;
however, compatibility between ``SUNMatrix`` and ``N_Vector``
implementations is more crucial when considering their interaction
within ``SUNLinearSolver`` objects, as will be described in more detail in
section :ref:`SUNLinSol`.  More specifically, in the Table
:ref:`SUNMatrix.matrix-vector` we show the matrix interfaces available as
``SUNMatrix`` modules, and the compatible vector implementations.


.. _SUNMatrix.matrix-vector:

SUNDIALS matrix interfaces and vector implementations that can be used for each
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. cssclass:: table-bordered

================ ====== ============== ====== ======== ============ ========== ==== ==== ===========
Linear Solver    Serial Parallel (MPI) OpenMP pThreads *hypre* Vec. PETSc Vec. CUDA RAJA User Suppl.
================ ====== ============== ====== ======== ============ ========== ==== ==== ===========
Dense            X                     X      X                                          X
Band             X                     X      X                                          X
Sparse           X                     X      X                                          X
User supplied    X      X              X      X        X            X          X    X    X 
================ ====== ============== ====== ======== ============ ========== ==== ==== ===========
