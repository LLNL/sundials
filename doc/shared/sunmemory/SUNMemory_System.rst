..
   ----------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2002-2025, Lawrence Livermore National Security
   and Southern Methodist University.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   ----------------------------------------------------------------

.. _SUNMemory.Sys:

The SUNMemoryHelper_Sys Implementation
=======================================

The SUNMemoryHelper_Sys module is an implementation of the :c:type:`SUNMemoryHelper`.
API that interfaces with standard library memory management through malloc/free.
The implementation defines the constructor

.. c:function:: SUNMemoryHelper SUNMemoryHelper_Sys(SUNContext sunctx)

   Allocates and returns a :c:type:`SUNMemoryHelper` object for handling system memory
   if successful. Otherwise, it returns ``NULL``.

.. _SUNMemory.Sys.Operations:

SUNMemoryHelper_Sys API Functions
----------------------------------

The implementation provides the following operations defined by the
``SUNMemoryHelper`` API:

* :c:func:`SUNMemoryHelper_Alloc`
* :c:func:`SUNMemoryHelper_AllocStrided`
* :c:func:`SUNMemoryHelper_Dealloc`
* :c:func:`SUNMemoryHelper_Copy`
* :c:func:`SUNMemoryHelper_Clone`
* :c:func:`SUNMemoryHelper_GetAllocStats`
* :c:func:`SUNMemoryHelper_Destroy`
