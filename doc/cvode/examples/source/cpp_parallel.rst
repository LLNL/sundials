..
   Programmer(s): Daniel M. Margolis @ SMU
   ----------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2002-2023, Lawrence Livermore National Security
   and Southern Methodist University.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   ----------------------------------------------------------------

:tocdepth: 3


.. _parallel_cpp:

=====================
Parallel C++ Examples
=====================



.. _cv_heat2D_p:

cv_heat2D_p
======================================================================

Description
------------

Exactly the same as ``cv_heat2D`` in :ref:`serial_cpp` except that here we
are using multiple processors with MPI.


Problem output
---------------

.. literalinclude:: ../../../../examples/cvode/CXX_parallel/cv_heat2D_p_--np_2_2.out
   :language: text


Numerical method
----------------

The same as ``cv_heat2D`` aside from the use of NVECTOR_PARALLEL.
