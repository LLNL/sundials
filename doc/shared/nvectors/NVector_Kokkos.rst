..
   ----------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2025-2026, Lawrence Livermore National Security,
   University of Maryland Baltimore County, and the SUNDIALS contributors.
   Copyright (c) 2013-2025, Lawrence Livermore National Security
   and Southern Methodist University.
   Copyright (c) 2002-2013, Lawrence Livermore National Security.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   ----------------------------------------------------------------

.. _NVectors.Kokkos:

The NVECTOR_KOKKOS Module
=========================

.. versionchanged:: x.y.z

   NVECTOR_KOKKOS requires Kokkos 4 or newer and C++17. Kokkos 5 and newer
   additionally require C++20. It wraps caller-owned,
   same-rank Kokkos views. The wrapper does not allocate or copy while wrapping.

Usage
-----

Include ``nvector/nvector_kokkos.hpp`` and construct a vector from an ordinary
application-owned view:

.. code-block:: cpp

   Kokkos::View<sunrealtype**, Kokkos::LayoutRight,
                Kokkos::DefaultExecutionSpace::memory_space> data("data", 8, 3);
   using Vector = sundials::kokkos::VectorFor<decltype(data)>;
   Vector y{data, sunctx};

``VectorFor<View>`` and class-template argument deduction accept managed and
explicit unmanaged views. Managed inputs retain a shallow handle to their
allocation; unmanaged inputs are borrowed and retain no owner. The caller must
keep borrowed storage alive until all queued work has completed. Copying a
wrapper is disabled; ``N_VClone`` creates independent managed storage with the
same rank, extents, layout, and supported strides.

The canonical operational view is returned by ``View()``. Use
``VectorFor<View>::FromNVector(v)`` (or ``FromNVectorConst``) to recover the C++
object in callbacks. ``GetLength()`` reports the product of all extents.

Views must be writable ``sunrealtype`` views in an execution-accessible
``LayoutLeft``, ``LayoutRight``, or valid ``LayoutStride`` mapping, with rank
one through six. Rank 1 operations use ``RangePolicy``; ranks 2--6 use native
``MDRangePolicy`` traversal. Ranks 7 and higher are rejected during
construction because Kokkos does not support the required MD policy ranks.
Shapes must match for operations. Host mirrors,
``Kokkos::deep_copy``, and synchronization are application responsibilities;
there are no host-view or allocating-constructor helpers. Raw pointer access is
available only for contiguous rank-1 mappings, and set-pointer operations are
unsupported.

Dense Kokkos matrix and linear-solver integrations intentionally require
rank-1 Kokkos vectors. Higher-rank vectors use the native-rank N_Vector
operations but must use an appropriate solver integration.
