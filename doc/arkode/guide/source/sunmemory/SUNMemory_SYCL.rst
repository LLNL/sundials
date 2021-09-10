..
   ----------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2002-2021, Lawrence Livermore National Security
   and Southern Methodist University.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   ----------------------------------------------------------------

:tocdepth: 3

.. _SUNMemory.SYCL:

The SUNMemoryHelper_Sycl Implementation
=======================================

The ``SUNMemoryHelper_Sycl`` module is an implementation of the
``SUNMemoryHelper`` API that interfaces to the SYCL abstraction layer. The
implementation defines the constructor

.. c:function:: SUNMemoryHelper SUNMemoryHelper_Sycl(sycl::queue *Q)

  Allocates and returns a ``SUNMemoryHelper`` object for handling SYCL memory
  using the provided queue. Otherwise it returns ``NULL``.


.. _SUNMemory.SYCL.Operations:

SUNMemoryHelper API Functions
-----------------------------

The implementation provides the following operations defined by the
``SUNMemoryHelper`` API:

.. c:function:: SUNMemory SUNMemoryHelper_Alloc_Sycl(SUNMemoryHelper helper, SUNMemory memptr, size_t mem_size, SUNMemoryType mem_type)

  Allocates a ``SUNMemory`` object whose ``ptr`` field is allocated for
  ``mem_size`` bytes and is of type ``mem_type``. The new object will have
  ownership of ``ptr`` and will be deallocated when ``SUNMemoryHelper_Dealloc``
  is called.

  The ``SUNMemoryType`` supported are

  - ``SUNMEMTYPE_HOST`` -- memory is allocated with a call to ``malloc``
  - ``SUNMEMTYPE_PINNED`` -- memory is allocated with a call to
    ``sycl::malloc_host``
  - ``SUNMEMTYPE_DEVICE`` -- memory is allocated with a call to
    ``sycl::malloc_device``
  - ``SUNMEMTYPE_UVM`` -- memory is allocated with a call to
    ``sycl::malloc_shared``

  **Arguments:**

  - *helper*  -- the ``SUNMemoryHelper`` object
  - *memptr* -- pointer to the allocated ``SUNMemory``
  - *mem_size* -- the size in bytes of the ``ptr``
  - *mem_type* -- the ``SUNMemoryType`` of the ``ptr``

  **Returns:**

    An ``int`` flag indicating success (zero) or failure (non-zero).


.. c:function:: int SUNMemoryHelper_Dealloc_Sycl(SUNMemoryHelper helper, SUNMemory mem)

  Deallocates the ``mem->ptr`` field if it is owned by ``mem``, and then
  deallocates the ``mem`` object.

  **Arguments:**

  - *helper* -- the ``SUNMemoryHelper`` object
  - *mem* -- the ``SUNMemory`` object

  **Returns:**

    An ``int`` flag indicating success (zero) or failure (non-zero).


.. c:function:: int SUNMemoryHelper_Copy_Sycl(SUNMemoryHelper helper, SUNMemory dst, SUNMemory src, size_t mem_size)

  Synchronously copies ``mem_size`` bytes from the the source memory to the
  destination memory.  The copy can be across memory spaces, e.g. host to
  device, or within a memory space, e.g. host to host.  The ``helper``
  object will use the memory types of ``dst`` and ``src`` to determine
  the appropriate transfer type necessary.

  **Arguments:**

  - *helper* -- the ``SUNMemoryHelper`` object
  - *dst* -- the destination memory to copy to
  - *src* -- the source memory to copy from
  - *mem_size* -- the number of bytes to copy

  **Returns:**

    An ``int`` flag indicating success (zero) or failure (non-zero).


.. c:function:: int SUNMemoryHelper_CopyAsync(SUNMemoryHelper helper, SUNMemory dst, SUNMemory src, size_t mem_size, void* ctx)

  Asynchronously copies ``mem_size`` bytes from the the source memory to the
  destination memory.  The copy can be across memory spaces, e.g. host to
  device, or within a memory space, e.g. host to host.  The ``helper`` object
  will use the memory types of ``dst`` and ``src`` to determine the
  appropriate transfer type necessary.

  **Arguments:**

  - *helper* -- the ``SUNMemoryHelper`` object
  - *dst* -- the destination memory to copy to
  - *src* -- the source memory to copy from
  - *mem_size* -- the number of bytes to copy
  - *ctx* -- is unused in this function

  **Returns:**

    An ``int`` flag indicating success (zero) or failure (non-zero).
