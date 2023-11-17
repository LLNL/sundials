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


.. _sycl_cpp:

====================================
SYCL C++ example problems
====================================

.. _cvAdvDiff_kry_sycl:

cvAdvDiff_kry_sycl
==============================

Description
------------

This example problem is a simple demonstration of the method
used to incorporate SYCL into SUNDIALS.  Otherwise, it is exactly
the same problem as ``cvAdvDiff_bnd`` from :ref:`serial_c` except
that here we use the SUNDIALS SPGMR linear solver instead.


Problem output
---------------

.. include:: ../../../../examples/cvode/CPP_sycl/cvAdvDiff_kry_sycl.out
   :literal:


Numerical method
-----------------

As previously mentioned, the problem is exactly the same as ``cvAdvDiff_bnd``
except that here we use SUNLINSOL_SPGMR and NVECTOR_SYCL for demonstration.


