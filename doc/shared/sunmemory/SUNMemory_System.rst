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

.. _SUNMemory.Sys:

The SUNMemoryHelper_Sys Implementation
=======================================

The SUNMemoryHelper_Sys module is an implementation of the ``SUNMemoryHelper``
API that interfaces with system memory management. The implementation defines the constructor

.. c:function:: SUNMemoryHelper SUNMemoryHelper_Sys(SUNContext sunctx)

   Allocates and returns a ``SUNMemoryHelper`` object for handling system memory
   if successful. Otherwise, it returns ``NULL``.

.. _SUNMemory.Sys.Operations:

SUNMemoryHelper_Sys API Functions
----------------------------------

The implementation provides the following operations defined by the
``SUNMemoryHelper`` API:

.. c:function:: SUNErrCode SUNMemoryHelper_Alloc_Sys(SUNMemoryHelper helper, \
                                                     SUNMemory* memptr, \
                                                     size_t mem_size, \
                                                     SUNMemoryType mem_type, \
                                                     void* queue)

   Allocates a ``SUNMemory`` object whose ``ptr`` field is allocated for
   ``mem_size`` bytes and is of type ``mem_type``. The new object will have
   ownership of ``ptr`` and will be deallocated when
   :c:func:`SUNMemoryHelper_Dealloc_Sys` is called.

.. c:function:: SUNErrCode SUNMemoryHelper_AllocStrided_Sys(SUNMemoryHelper helper, \
                                                          SUNMemory* memptr, \
                                                          size_t mem_size, \
                                                          size_t stride, \
                                                          SUNMemoryType mem_type, \
                                                          void* queue)

   Allocates a strided ``SUNMemory`` object whose ``ptr`` field is allocated for
   ``mem_size`` bytes with a specified ``stride``. The new object will have
   ownership of ``ptr`` and will be deallocated when
   :c:func:`SUNMemoryHelper_Dealloc_Sys` is called.

.. c:function:: SUNErrCode SUNMemoryHelper_Dealloc_Sys(SUNMemoryHelper helper, \
                                                       SUNMemory mem, void* queue)

   Deallocates the ``mem->ptr`` field if it is owned by ``mem``, and then
   deallocates the ``mem`` object.

.. c:function:: SUNErrCode SUNMemoryHelper_Copy_Sys(SUNMemoryHelper helper, \
                                                    SUNMemory dst, \
                                                    SUNMemory src, \
                                                    size_t memory_size, \
                                                    void* queue)

   Synchronously copies ``memory_size`` bytes from the source memory to the
   destination memory. The copy can be across memory spaces or within a memory space.

.. c:function:: SUNMemoryHelper SUNMemoryHelper_Clone_Sys(SUNMemoryHelper helper)

   Clones the given ``SUNMemoryHelper`` object, creating a new instance with the same
   properties.

.. c:function:: SUNErrCode SUNMemoryHelper_GetAllocStats_Sys(SUNMemoryHelper helper, \
                                                             SUNMemoryType mem_type, \
                                                             unsigned long* num_allocations, \
                                                             unsigned long* num_deallocations, \
                                                             size_t* bytes_allocated, \
                                                             size_t* bytes_high_watermark)

   Returns statistics about memory allocations performed with the helper.

.. c:function:: SUNErrCode SUNMemoryHelper_Destroy_Sys(SUNMemoryHelper helper)

   Destroys the ``SUNMemoryHelper`` object and frees any associated resources.