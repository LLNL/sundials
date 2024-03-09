.. ----------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2002-2024, Lawrence Livermore National Security
   and Southern Methodist University.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End

   intersphinx and :numref: does not work currently see
   https://github.com/sphinx-doc/sphinx/issues/12033
   currently in master branch

   Even if :numref: worked how would it show up? The section might
   not exist or could conflict with other sections
   ----------------------------------------------------------------

**New Features**

**Bug Fixes**

Updated the CMake variable ``HIP_PLATFORM`` to default to ``amd`` as the
previous default, ``hcc``, is no longer recognized in ROCm 5.7.0 or newer. The
new default is also valid in older version of ROCm (at least back to version
4.3.1).

Fixed a bug in the HIP execution policies where ``WAPR_SIZE`` would not be set
with ROCm 6.0.0 or newer.
