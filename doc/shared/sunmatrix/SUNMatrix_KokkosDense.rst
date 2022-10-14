..
   ----------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2002-2022, Lawrence Livermore National Security
   and Southern Methodist University.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   ----------------------------------------------------------------

.. _SUNMatrix.Kokkos:

The SUNMATRIX_KOKKOSDENSE Module
================================

.. versionadded:: 6.4.0

The SUNMATRIX_KOKKOSDENSE :c:type:`SUNMatrix` implementation provides a data
structure for dense and dense batched (block-diagonal) matrices using Kokkos
:cite:p:`edwards2014kokkos,trott2022kokkos` and KokkosKernels
:cite:p:`trott2021kokkos` to support variety backends including serial, OpenMP,
CUDA, HIP, and SYCL. Since Kokkos is a modern C++ library, the module is also
written in modern C++ (it requires C++14) as a header only library. To utilize
this ``SUNMatrix`` users will need to include
``sunmatrix/sunmatrix_kokkosdense.hpp``. More instructions on building SUNDIALS
with Kokkos and KokkosKernels enabled are given in
:numref:`Installation.CMake.ExternalLibraries`. For instructions on building and
using Kokkos and KokkosKernels, refer to the
`Kokkos <https://kokkos.github.io/kokkos-core-wiki/index.html>`_
and `KokkosKernels <https://github.com/kokkos/kokkos-kernels/wiki>`_.
documentation.

.. _SUNMatrix.Kokkos.Usage:

Using SUNMATRIX_KOKKOSDENSE
----------------------------

The SUNMATRIX_KOKKOSDENSE module is defined by the ``DenseMatrix`` templated
class in the ``sundials::kokkos`` namespace:

.. code-block:: cpp

   template<class ExecSpace = Kokkos::DefaultExecutionSpace,
            class MemSpace = class ExecSpace::memory_space>
   class DenseMatrix : public sundials::impl::BaseMatrix,
                       public sundials::ConvertibleTo<SUNMatrix>

To use the SUNMATRIX_KOKKOSDENSE module, we begin by constructing an instance of
the Kokkos dense matrix e.g.,

.. code-block:: cpp

   // Single matrix using the default execution space
   sundials::kokkos::DenseMatrix<> A{rows, cols, sunctx};

   // Batched (block-diagonal) matrix using the default execution space
   sundials::kokkos::DenseMatrix<> Abatch{blocks, rows, cols, sunctx};

   // Batched (block-diagonal) matrix using the user-selected execution space
   sundials::kokkos::DenseMatrix<my_exec_space> Abatch{blocks, rows, cols, sunctx};

   // Batched (block-diagonal) matrix using the user-selected execution space and
   // execution space instance
   sundials::kokkos::DenseMatrix<my_exec_space> Abatch{blocks, rows, cols,
                                                       my_exec_space_instance,
                                                       sunctx};

Instances of the ``DenseMatrix`` class are implicitly or explicitly (using the
:cpp:func:`~DenseMatrix::Convert` method) convertible to a :c:type:`SUNMatrix`
e.g.,

.. code-block:: cpp

   sundials::kokkos::DenseMatrix<> A{rows, cols, sunctx};
   SUNMatrix B = A;           // implicit conversion to SUNMatrix
   SUNMatrix C = A.Convert(); // explicit conversion to SUNMatrix

No further interaction with a ``DenseMatrix`` is required from this point, and
it is possible to use the :c:type:`SUNMatrix` API to operate on ``B`` or ``C``.

.. warning::

   :c:func:`SUNMatDestroy` should never be called on a ``SUNMatrix`` that was
   created via conversion from a ``sundials::kokkos::DenseMatrix``. Doing so may
   result in a double free.

The underlying ``DenseMatrix`` can be extracted from a ``SUNMatrix`` using
:cpp:func:`GetDenseMat` e.g.,

.. code-block:: cpp

   auto A_dense_mat = GetDenseMat<>(A_sunmat);

The SUNMATRIX_KOKKOSDENSE module is compatible with the NVECTOR_KOKKOS vector
module (see :numref:`NVectors.Kokkos`) and SUNLINEARSOLVER_KOKKOSDENSE linear
solver module (see :numref:`SUNLinSol.Kokkos`).


.. _SUNMatrix.Kokkos.API:

SUNMATRIX_KOKKOSDENSE API
-------------------------

In this section we list the public API of the ``sundials::kokkos::DenseMatrix``
class.

.. cpp:class:: template<class ExecSpace = Kokkos::DefaultExecutionSpace, \
                        class MemSpace = class ExecSpace::memory_space> \
               DenseMatrix : public sundials::impl::BaseMatrix, \
                             public sundials::ConvertibleTo<SUNMatrix>

   .. cpp:function:: DenseMatrix() = default

      Default constructor -- the matrix must be copied or moved to.

   .. cpp:function:: DenseMatrix(sunindextype rows, sunindextype cols, \
                                 SUNContext sunctx)

      Constructs a single DenseMatrix using the default `ExecSpace` space
      instance.

      :param rows: number of matrix rows
      :param cols: number of matrix columns
      :param sunctx: the SUNDIALS simulation context object (:c:type:`SUNContext`)

   .. cpp:function:: DenseMatrix(sunindextype rows, sunindextype cols, \
                                 ExecSpace exec_space, SUNContext sunctx)

      Constructs a single DenseMatrix using the provided `ExecSpace` space
      instance.

      :param rows: number of matrix rows
      :param cols: number of matrix columns
      :param exec_space: a `ExecSpace` instance
      :param sunctx: the SUNDIALS simulation context object (:c:type:`SUNContext`)

   .. cpp:function:: DenseMatrix(sunindextype blocks, sunindextype block_rows, \
                                 sunindextype block_cols, SUNContext sunctx)

      Constructs a batched (block-diagonal) DenseMatrix using the default
      `ExecSpace` space instance.

      :param blocks: number of matrix blocks
      :param block_rows: number of rows in a block
      :param block_cols: number of columns in a block
      :param sunctx: the SUNDIALS simulation context object (:c:type:`SUNContext`)

   .. cpp:function:: DenseMatrix(sunindextype blocks, sunindextype block_rows, \
                                 sunindextype block_cols, ExecSpace exec_space, \
                                 SUNContext sunctx)

      Constructs a batched (block-diagonal) DenseMatrix using the provided
      `ExecSpace` space instance.

      :param blocks: number of matrix blocks
      :param block_rows: number of rows in a block
      :param block_cols: number of columns in a block
      :param exec_space: a `ExecSpace` instance
      :param sunctx: the SUNDIALS simulation context object (:c:type:`SUNContext`)

   .. cpp:function:: DenseMatrix(DenseMatrix&& that_matrix) noexcept

      Move constructor.

   .. cpp:function:: DenseMatrix(const DenseMatrix& that_matrix)

      Copy constructor (performs a deep copy).

   .. cpp:function:: DenseMatrix& operator=(DenseMatrix&& rhs) noexcept

      Move assignment.

   .. cpp:function:: DenseMatrix& operator=(const DenseMatrix& rhs)

      Copy assignment. This is a shallow copy i.e., a new view is not created.

   .. cpp:function:: virtual ~DenseMatrix() = default;

      Default destructor.

   .. cpp:function:: ExecSpace exec_space()

      Get the execution space instance used by the matrix.

   .. cpp:function:: Kokkos::View<sunrealtype***, MemSpace> view()

      Get the underlying Kokkos view with extents
      ``{blocks, block_rows, block_cols}``.

   .. cpp:function:: sunindextype blocks()

      Get the number of blocks i.e., ``extent(0)``.

   .. cpp:function:: sunindextype block_rows()

      Get the number of rows in a block i.e., ``extent(1)``.

   .. cpp:function:: sunindextype block_cols()

      Get the number of columns in a block i.e., ``extent(2)``.

   .. cpp:function:: sunindextype rows()

      Get the number of rows in the block-diagonal matrix i.e.,
      ``extent(0) * extent(1)``.

   .. cpp:function:: sunindextype cols()

      Get the number of columns in the block-diagonal matrix i.e.,
      ``extent(0) * extent(2)``.

   .. cpp:function:: operator SUNMatrix() override

      Implicit conversion to a :c:type:`SUNMatrix`.

   .. cpp:function:: operator SUNMatrix() const override

      Implicit conversion to a :c:type:`SUNMatrix`.

   .. cpp:function:: SUNMatrix Convert() override

      Explicit conversion to a :c:type:`SUNMatrix`.

   .. cpp:function:: SUNMatrix Convert() const override

      Explicit conversion to a :c:type:`SUNMatrix`.

.. cpp:function:: template<class ExecSpace = Kokkos::DefaultExecutionSpace, \
                           class MemSpace = class ExecSpace::memory_space> \
                  inline DenseMatrix<ExecSpace, MemSpace>* GetDenseMat(SUNMatrix A)

   Get the dense matrix wrapped by a SUNMatrix
