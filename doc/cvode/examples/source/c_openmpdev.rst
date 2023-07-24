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

.. _openmpdev_c:

============================================
OpenMP Device Offload C example problems
============================================




.. _cvAdvDiff_kry_ompdev:

cvAdvDiff_kry_ompdev
===============================

Description
------------

This simple example problem is an OpenMP dev implementation of previous
serial C example ``cvAdvDiff_bnd`` except that here we use the SPGMR
linear solver.


Problem output
---------------

.. include:: ../../../../examples/cvode/C_openmpdev/cvAdvDiff_kry_ompdev.out
   :literal:


Numerical method
-----------------

The numerical method is identical to the previous implementation,
except that we now use SUNDIALS' OpenMP-device-offload-enabled vector
kernel module, NVECTOR_OPENMPDEV, and have similarly threaded the
supplied right-hand side residual and banded Jacobian construction
functions.

