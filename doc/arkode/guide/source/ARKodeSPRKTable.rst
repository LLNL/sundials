.. ----------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2002-2023, Lawrence Livermore National Security
   and Southern Methodist University.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   ----------------------------------------------------------------

.. _SPRKStorage:

===========================
SPRK Method Table Structure
===========================

.. c:type:: ARKodeSPRKTableMem

   Structure representing the SPRK method that holds the method coefficients.

   .. c:member:: int q

      The method order of accuracy.

   .. c:member:: int stages

      The number of stages.

   .. c:member:: sunrealtype* a

      Array of coefficients that generate the explicit Butcher table.
      ``a[i]`` is the coefficient appearing in column i+1.

      .. math::
         \begin{array}{c|cccc}
         c_1 & 0 & \cdots & 0 & 0 \\
         c_2 & a_1 & 0 & \cdots & \vdots \\
         \vdots & \vdots & \ddots & \ddots & \vdots \\
         c_s & a_1 & \cdots & a_{s-1} & 0 \\
         \hline
         & a_1 & \cdots & a_{s-1} & a_s
         \end{array}.

   .. c:member:: sunrealtype* ahat

      Array of coefficients that generate the diagonally-implicit Butcher table.
      ``ahat[i]`` is the coefficient appearing in column i.

      .. math::
         \begin{array}{c|cccc}
         \hat{c}_1 & \hat{a}_1 & \cdots & 0 & 0 \\
         \hat{c}_2 & \hat{a}_1 & \hat{a}_2 & \cdots & \vdots \\
         \vdots & \vdots & \ddots & \ddots & \vdots \\
         \hat{c}_s & \hat{a}_1 & \hat{a}_2 & \cdots & \hat{a}_{s} \\
         \hline
         & \hat{a}_1 & \hat{a}_2 & \cdots & \hat{a}_{s}
         \end{array}.


.. c:type:: ARKodeSPRKTableMem* ARKodeSPRKTable


ARKodeSPRKTable functions
---------------------------

.. _ARKodeSPRKTable.FunctionsTable:
.. table:: ARKodeSPRKTable functions

   +----------------------------------------------+------------------------------------------------------------+
   | **Function name**                            | **Description**                                            |
   +==============================================+============================================================+
   | :c:func:`ARKodeSPRKTable_Alloc()`            | Allocate an empty storage structure                        |
   +----------------------------------------------+------------------------------------------------------------+
   | :c:func:`ARKodeSPRKTable_Load()`             | Load SPRK method using an identifier                       |
   +----------------------------------------------+------------------------------------------------------------+
   | :c:func:`ARKodeSPRKTable_LoadByName()`       | Load SPRK method using a string version of the identifier  |
   +----------------------------------------------+------------------------------------------------------------+
   | :c:func:`ARKodeSPRKTable_Create()`           | Create a new storage structure                             |
   +----------------------------------------------+------------------------------------------------------------+
   | :c:func:`ARKodeSPRKTable_Copy()`             | Create a copy of a storage structure                       |
   +----------------------------------------------+------------------------------------------------------------+
   | :c:func:`ARKodeSPRKTable_Space()`            | Get the storage structure real and integer workspace size  |
   +----------------------------------------------+------------------------------------------------------------+
   | :c:func:`ARKodeSPRKTable_Free()`             | Deallocate a storage structure                             |
   +----------------------------------------------+------------------------------------------------------------+


.. c:function:: ARKodeSPRKTable ARKodeSPRKTable_Alloc(int stages)

   Allocate memory for an :c:type:`ARKodeSPRKTable` structure with the specified number of stages.

   :param stages: The number of stages.
   :return: :c:type:`ARKodeSPRKTable` structure for the loaded method.

.. c:function:: ARKodeSPRKTable ARKodeSPRKTable_Load(ARKODE_SPRKMethodID id)

   Load the :c:type:`ARKodeSPRKTable` structure for the specified method ID.

   :param id: The ID of the SPRK method, see :ref:`Butcher.sprk`.
   :return: :c:type:`ARKodeSPRKTable` structure for the loaded method.

.. c:function:: ARKodeSPRKTable ARKodeSPRKTable_LoadByName(const char* method)

   Load the :c:type:`ARKodeSPRKTable` structure for the specified method name.

   :param method: The name of the SPRK method, see :ref:`Butcher.sprk`.
   :return: :c:type:`ARKodeSPRKTable` structure for the loaded method.

.. c:function:: ARKodeSPRKTable ARKodeSPRKTable_Copy(ARKodeSPRKTable sprk_table)

   Create a copy of the :c:type:`ARKodeSPRKTable` structure.

   :param sprk_table: The :c:type:`ARKodeSPRKTable` structure to copy.
   :return: Pointer to the copied :c:type:`ARKodeSPRKTable` structure.

.. c:function:: void ARKodeSPRKTable_Write(ARKodeSPRKTable sprk_table, FILE* outfile)

   Write the ARKodeSPRKTable out to the file.

   :param sprk_table: The :c:type:`ARKodeSPRKTable` structure to write.
   :param outfile: The FILE that will be written to.
   :return: void

.. c:function:: void ARKodeSPRKTable_Space(ARKodeSPRKTable sprk_table, sunindextype* liw, sunindextype* lrw)

   Get the workspace sizes required for the :c:type:`ARKodeSPRKTable` structure.

   :param sprk_table: The :c:type:`ARKodeSPRKTable` structure.
   :param liw: Pointer to store the integer workspace size.
   :param lrw: Pointer to store the real workspace size.

.. c:function:: void ARKodeSPRKTable_Free(ARKodeSPRKTable sprk_table)

   Free the memory allocated for the :c:type:`ARKodeSPRKTable` structure.

   :param sprk_table: The :c:type:`ARKodeSPRKTable` structure to free.

.. c:function:: int ARKodeSPRKTable_ToButcher(ARKodeSPRKTable sprk_storage, ARKodeSPRKTable* a_ptr, ARKodeSPRKTable* b_ptr)

   Convert the :c:type:`ARKodeSPRKTable` structure to the Butcher table representation.

   :param sprk_storage: The :c:type:`ARKodeSPRKTable` structure.
   :param a_ptr: Pointer to store the explicit Butcher table.
   :param b_ptr: Pointer to store the diagonally-implicit Butcher table.
