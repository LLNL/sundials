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

.. _SUNMemory.Description:

The SUNMemoryHelper API
=======================

This API consists of three new sundials types: ``SUNMemoryType``, ``SUNMemory``,
and ``SUNMemoryHelper``, which we now define.

The ``SUNMemory`` structure wraps a pointer to actual data. This structure
is defined as

.. code-block:: c

  typedef struct _SUNMemory
  {
    void*         ptr;
    SUNMemoryType type;
    booleantype   own;
  } *SUNMemory;

The ``SUNMemoryType`` type is an enumeration that defines the four supported
memory types:

.. code-block:: c

  typedef enum
  {
    SUNMEMTYPE_HOST,      /* pageable memory accessible on the host     */
    SUNMEMTYPE_PINNED,    /* page-locked memory accesible on the host   */
    SUNMEMTYPE_DEVICE,    /* memory accessible from the device          */
    SUNMEMTYPE_UVM        /* memory accessible from the host or device  */
  } SUNMemoryType;

Finally, the ``SUNMemoryHelper`` structure is defined as

.. code-block:: c

  struct _SUNMemoryHelper
  {
    void*               content;
    SUNMemoryHelper_Ops ops;
  } *SUNMemoryHelper;

where ``SUNMemoryHelper_Ops`` is defined as

.. code-block:: c

  typedef struct _SUNMemoryHelper_Ops
  {
    /* operations that implementations are required to provide */
    int             (*alloc)(SUNMemoryHelper, SUNMemory* memptr
                             size_t mem_size, SUNMemoryType mem_type);
    int             (*dealloc)(SUNMemoryHelper, SUNMemory mem);
    int             (*copy)(SUNMemoryHelper, SUNMemory dst, SUNMemory src,
                            size_t mem_size);

    /* operations that provide default implementations */
    int             (*copyasync)(SUNMemoryHelper, SUNMemory dst, SUNMemory src,
                                 size_t mem_size, void* ctx);
    SUNMemoryHelper (*clone)(SUNMemoryHelper);
    int             (*destroy)(SUNMemoryHelper);
  } *SUNMemoryHelper_Ops;



.. _SUNMemory.Required:

Implementation defined operations
---------------------------------

The SUNMemory API also defines the following operations which do require
a SUNMemoryHelper instance and **require** the implementation to define
them:

.. c:function:: SUNMemory SUNMemoryHelper_Alloc(SUNMemoryHelper helper, SUNMemory* memptr, size_t mem_size, SUNMemoryType mem_type)

  Allocates a ``SUNMemory`` object whose ``ptr`` field is allocated for
  ``mem_size`` bytes and is of type ``mem_type``. The new object will have
  ownership of ``ptr`` and will be deallocated when ``SUNMemoryHelper_Dealloc``
  is called.

  **Arguments:**

  - *helper*  -- the ``SUNMemoryHelper`` object
  - *memptr* -- pointer to the allocated ``SUNMemory``
  - *mem_size* -- the size in bytes of the ``ptr``
  - *mem_type* -- the ``SUNMemoryType`` of the ``ptr``

  **Returns:**

    An ``int`` flag indicating success (zero) or failure (non-zero).


.. c:function:: int SUNMemoryHelper_Dealloc(SUNMemoryHelper helper, SUNMemory mem)

  Deallocates the ``mem->ptr`` field if it is owned by ``mem``, and then
  deallocates the ``mem`` object.

  **Arguments:**

  - *helper* -- the ``SUNMemoryHelper`` object
  - *mem* -- the ``SUNMemory`` object

  **Returns:**

    An ``int`` flag indicating success (zero) or failure (non-zero).


.. c:function:: int SUNMemoryHelper_Copy(SUNMemoryHelper helper, SUNMemory dst, SUNMemory src, size_t mem_size)

  Synchronously copies ``mem_size`` bytes from the the source memory to the
  destination memory.  The copy can be across memory spaces, e.g. host to
  device, or within a memory space, e.g. host to host.  The ``helper``
  object should use the memory types of ``dst`` and ``src`` to determine
  the appropriate transfer type necessary.

  **Arguments:**

  - *helper* -- the ``SUNMemoryHelper`` object
  - *dst* -- the destination memory to copy to
  - *src* -- the source memory to copy from
  - *mem_size* -- the number of bytes to copy

  **Returns:**

    An ``int`` flag indicating success (zero) or failure (non-zero).



.. _SUNMemory.Utilities:

Utility Functions
-----------------

The SUNMemoryHelper API defines the following functions which do not
require a SUNMemoryHelper instance:

.. c:function:: SUNMemory SUNMemoryHelper_Alias(SUNMemory mem1)

  Returns a ``SUNMemory`` object whose ``ptr`` field points to the same address
  as ``mem1``. The new object *will not* have ownership of ``ptr``, therefore,
  it will not free ``ptr`` when ``SUNMemoryHelper_Dealloc`` is called.

  **Arguments:**

  - *mem1* -- a ``SUNMemory`` object

  **Returns:**

    A ``SUNMemory`` object.


.. c:function:: SUNMemory SUNMemoryHelper_Wrap(void* ptr, SUNMemoryType mem_type)

  Returns a ``SUNMemory`` object whose ``ptr`` field points to the ``ptr``
  argument passed to the function. The new object *will not* have ownership of
  ``ptr``, therefore, it will not free ``ptr`` when ``SUNMemoryHelper_Dealloc``
  is called.

  **Arguments:**

  - *ptr* -- the data pointer to wrap in a ``SUNMemory`` object
  - *mem_type* -- the ``SUNMemoryType`` of the ``ptr``

  **Returns:**

    A ``SUNMemory`` object.


.. c:function:: SUNMemoryHelper SUNMemoryHelper_NewEmpty()

  Returns an empty ``SUNMemoryHelper``. This is useful for building custom
  ``SUNMemoryHelper`` implementations.

  **Returns:**

   A ``SUNMemoryHelper`` object.


.. c:function:: int SUNMemoryHelper_CopyOps(SUNMemoryHelper src, SUNMemoryHelper dst)

  Copies the ``ops`` field of ``src`` to the ``ops`` field of ``dst``.
  This is useful for building custom ``SUNMemoryHelper`` implementations.

  **Arguments:**

  - *src* -- the object to copy from
  - *dst* -- the object to copy to

  **Returns:**

    An ``int`` flag indicating success (zero) or failure (non-zero).


.. _SUNMemory.Overridable:

Implementation overridable operations with defaults
---------------------------------------------------

In addition, the SUNMemoryHelper API defines the following *optionally
overridable* operations which do require a SUNMemoryHelper instance:


.. c:function:: int SUNMemoryHelper_CopyAsync(SUNMemoryHelper helper, SUNMemory dst, SUNMemory src, size_t mem_size, void* ctx)

  Asynchronously copies ``mem_size`` bytes from the the source memory to the
  destination memory.  The copy can be across memory spaces, e.g. host to
  device, or within a memory space, e.g. host to host.  The ``helper`` object
  should use the memory types of ``dst`` and ``src`` to determine the
  appropriate transfer type necessary.  The ``ctx`` argument is used when a
  different execution stream needs to be provided to perform the copy in,
  e.g. with ``CUDA`` this would be a ``cudaStream_t``.

  **Arguments:**

  - *helper* -- the ``SUNMemoryHelper`` object
  - *dst* -- the destination memory to copy to
  - *src* -- the source memory to copy from
  - *mem_size* -- the number of bytes to copy
  - *ctx* -- typically a handle for an object representing an alternate
    execution stream, but it can be any implementation specific data

  **Returns:**

    An ``int`` flag indicating success (zero) or failure (non-zero).

  .. note::

     If this operation is not defined by the implementation, then
     ``SUNMemoryHelper_Copy`` will be used.


.. c:function:: SUNMemoryHelper SUNMemoryHelper_Clone(SUNMemoryHelper helper)

  Clones the ``SUNMemoryHelper`` object itself.

  **Arguments:**

  - *helper* -- the ``SUNMemoryHelper`` object to clone

  **Returns:**

    A ``SUNMemoryHelper`` object.

  .. note::

     If this operation is not defined by the implementation, then the default
     clone will only copy the ``SUNMemoryHelper_Ops`` structure stored in
     ``helper->ops``, and not the ``helper->content`` field.


.. c:function:: int SUNMemoryHelper_Destroy(SUNMemoryHelper helper)

  Destroys (frees) the ``SUNMemoryHelper`` object itself.

  **Arguments:**

  - *helper* -- the ``SUNMemoryHelper`` object to destroy

  **Returns:**

    An ``int`` flag indicating success (zero) or failure (non-zero).

  .. note::

     If this operation is not defined by the implementation, then the default
     destroy will only free the ``helper->ops`` field and the ``helper`` itself.
     The ``helper->content`` field will not be freed.



Implementing a custom SUNMemoryHelper
-------------------------------------

A particular implementation of the SUNMemoryHelper API must:

-  Define and implement the required operations. Note that the names of
   these routines should be unique to that implementation in order to
   permit using more than one SUNMemoryHelper module in the same code.

-  Optionally, specify the *content* field of SUNMemoryHelper.

-  Optionally, define and implement additional user-callable routines
   acting on the newly defined SUNMemoryHelper.

An example of a custom SUNMemoryHelper is given in
examples/utilities/custom_memory_helper.h.
