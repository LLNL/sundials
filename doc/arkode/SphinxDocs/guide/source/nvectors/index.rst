..
   Programmer(s): Daniel R. Reynolds @ SMU
   ----------------------------------------------------------------
   Copyright (c) 2013, Southern Methodist University.
   All rights reserved.
   For details, see the LICENSE file.
   ----------------------------------------------------------------

:tocdepth: 3

.. _NVectors:

=======================
Vector Data Structures
=======================

The SUNDIALS library comes packaged with a variety of NVECTOR
implementations, designed for simulations in serial, shared-memory
parallel, and distributed-memory parallel environments, as well as
interfaces to vector data structures used within external linear
solver libraries.  All native implementations assume that the
process-local data is stored contiguously, and they in turn provide a
variety of standard vector algebra operations that may be performed on
the data. 

In addition, SUNDIALS provides a simple interface for generic vectors
(akin to a C++ *abstract base class*).  All of the major SUNDIALS
solvers (CVODE(s), IDA(s), KINSOL, ARKODE) in turn are constructed to
only depend on these generic vector operations, making them immediately
extensible to new user-defined vector objects.  The only exceptions to
this rule relate to the dense, banded and sparse-direct linear system
solvers, since they rely on particular data storage and access
patterns in the NVECTORS used.



.. toctree::
   :maxdepth: 1

   NVector_Description
   NVector_Operations
   NVector_Serial
   NVector_Parallel
   NVector_OpenMP
   NVector_Pthreads
   NVector_ParHyp
   NVector_PETSc
   NVector_CUDA
   NVector_RAJA
   NVector_Examples
   ARKode_requirements
