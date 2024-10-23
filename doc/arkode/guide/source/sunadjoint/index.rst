..
   ----------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2002-2024, Lawrence Livermore National Security
   and Southern Methodist University.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   ----------------------------------------------------------------

.. _SUNAdjoint:

############################
Adjoint Sensitivity Analysis
############################

.. versionadded:: x.y.z

The ``SUNAdjoint`` API consists of a few customizable modules that provide a framework for adjoint
sensitivity analysis (ASA). The API itself does not implement ASA, but it provides a common
interface for ASA capabilities implemented in the SUNDIALS packages. Right now it supports :ref:`the
ASA capabilities in ARKODE <ARKODE.Mathematics.ASA>`, while the ASA capabilities in :ref:`CVODES
<CVODES.Mathematics.ASA>` and :ref:`IDAS <IDAS.Mathematics.ASA>` must be used directly.
*Users should read the package specific sections on ASA capabilities before this section.*

.. toctree::
   :maxdepth: 1

   SUNAdjoint_links.rst
