..
   Daniel R. Reynolds @ SMU
   ----------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2002-2024, Lawrence Livermore National Security
   and Southern Methodist University.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   ----------------------------------------------------------------

.. _SUNAdaptController:

#####################################
Time Step Adaptivity Controllers
#####################################

The SUNDIALS library comes packaged with a variety of :c:type:`SUNAdaptController`
implementations, designed to support various forms of error-based time step
adaptivity within SUNDIALS time integrators.  To support applications that may
want to adjust the controller parameters or provide their own implementations,
SUNDIALS defines the :c:type:`SUNAdaptController` base class, along with a
variety of default implementations.

.. toctree::
   :maxdepth: 1

   SUNAdaptController_links.rst
