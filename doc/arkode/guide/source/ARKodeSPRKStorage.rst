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

==============================
SPRK Method Storage Structure
==============================

.. _SPRKStorage.id:

Method identifiers
------------------

The following enum values are used to identify different SPRK methods:

.. c:macro:: ARKODE_SPRK_NONE

   Identifier representing no SPRK method, this is solely used to mark the beginning of the enum.

.. c:macro:: ARKODE_SPRK_EULER_1_1

   Identifier for the Symplectic Euler 1st order method with 1 stage.

.. c:macro:: ARKODE_SPRK_LEAPFROG_2_2

   Identifier for the Symplectic Leapfrog 2nd order method with 2 stages.

.. c:macro:: ARKODE_SPRK_PSEUDO_LEAPFROG_2_2

   Identifier for the Symplectic Pseudo Leapfrog 2nd order method with 2 stages.

.. c:macro:: ARKODE_SPRK_RUTH_3_3

   Identifier for the Symplectic Ruth 3rd order method with 3 stages.

.. c:macro:: ARKODE_SPRK_MCLACHLAN_2_2

   Identifier for the Symplectic McLachlan 2nd order method with 2 stages.

.. c:macro:: ARKODE_SPRK_MCLACHLAN_3_3

   Identifier for the Symplectic McLachlan 3rd order method with 3 stages.

.. c:macro:: ARKODE_SPRK_CANDY_ROZMUS_4_4

   Identifier for the Symplectic Candy Rozmus 4th order method with 4 stages.

.. c:macro:: ARKODE_SPRK_MCLACHLAN_4_4

   Identifier for the Symplectic McLachlan 4th order method with 4 stages.

.. c:macro:: ARKODE_SPRK_MCLACHLAN_5_6

   Identifier for the Symplectic McLachlan 5th order method with 6 stages.

.. c:macro:: ARKODE_SPRK_YOSHIDA_6_8

   Identifier for the Symplectic Yoshida 6th order method with 8 stages.

.. c:macro:: ARKODE_SPRK_SUZUKI_UMENO_8_16

   Identifier for the Symplectic Suzuki-Umeno 8th order method with 16 stages.

.. c:macro:: ARKODE_SPRK_SOFRONIOU_10_36

   Identifier for the Symplectic Sofroniou 10th order method with 36 stages.

.. c:type:: ARKodeSPRKTableMem

   Structure representing the SPRK method that holds the method coefficients.

   .. c:member:: int q

      The method order of accuracy.

   .. c:member:: int stages
      
      The number of stages.

   .. c:member:: sunrealtype* a

      Array of coefficients that generate the explicit Butcher table.
      ``a[i]`` is the coefficient appearing in column i+1.

   .. c:member:: sunrealtype* ahat

      Array of coefficients that generate the diagonally-implicit Butcher table.
      ``ahat[i]`` is the coefficient appearing in column i.

.. c:type:: ARKodeSPRKTableMem* ARKodeSPRKTable


ARKodeSPRKTable functions
---------------------------

.. _ARKodeSPRKTable.FunctionsTable:
.. table:: ARKodeSPRKTable functions

   +----------------------------------------------+------------------------------------------------------------+
   | **Function name**                            | **Description**                                            |
   +==============================================+============================================================+
   | :c:func:`ARKodeSPRKTable_Alloc()`          | Allocate an empty storage structure                        |
   +----------------------------------------------+------------------------------------------------------------+
   | :c:func:`ARKodeSPRKTable_Load()`           | Load SPRK method using an identifier                       |
   +----------------------------------------------+------------------------------------------------------------+
   | :c:func:`ARKodeSPRKTable_LoadByName()`     | Load SPRK method using a string version of the identifier  |
   +----------------------------------------------+------------------------------------------------------------+
   | :c:func:`ARKodeSPRKTable_Create()`         | Create a new storage structure                             |
   +----------------------------------------------+------------------------------------------------------------+
   | :c:func:`ARKodeSPRKTable_Copy()`           | Create a copy of a storage structure                       |
   +----------------------------------------------+------------------------------------------------------------+
   | :c:func:`ARKodeSPRKTableMempace()`          | Get the storage structure real and integer workspace size  |
   +----------------------------------------------+------------------------------------------------------------+
   | :c:func:`ARKodeSPRKTable_Free()`           | Deallocate a storage structure                             |
   +----------------------------------------------+------------------------------------------------------------+


.. c:function:: ARKodeSPRKTable ARKodeSPRKTable_Alloc(int stages)

   Allocate memory for an :c:type:`ARKodeSPRKTable`` structure with the specified number of stages.

   :param stages: The number of stages.
   :return: ARKodeSPRKTable structure for the loaded method.

.. c:function:: ARKodeSPRKTable ARKodeSPRKTable_Load(ARKODE_SPRKMethodID id)

   Load the ARKodeSPRKTable structure for the specified method ID.

   :param id: The ID of the SPRK method. One of :ref:`SPRKStorage.id`.
   :return: ARKodeSPRKTable structure for the loaded method.

.. c:function:: ARKodeSPRKTable ARKodeSPRKTable_LoadByName(const char* method)

   Load the ARKodeSPRKTable structure for the specified method name.

   :param method: The name of the SPRK method. Must be one of :ref:`SPRKStorage.id` but as a string.
   :return: ARKodeSPRKTable structure for the loaded method.


.. c:function:: ARKodeSPRKTable ARKodeSPRKTable_Copy(ARKodeSPRKTable B)

   Create a copy of the ARKodeSPRKTable structure.

   :param B: The ARKodeSPRKTable structure to copy.
   :return: Pointer to the copied ARKodeSPRKTable structure.

.. c:function:: void ARKodeSPRKTableMempace(ARKodeSPRKTable B, sunindextype* liw, sunindextype* lrw)

   Get the workspace sizes required for the ARKodeSPRKTable structure.

   :param B: The ARKodeSPRKTable structure.
   :param liw: Pointer to store the integer workspace size.
   :param lrw: Pointer to store the real workspace size.

.. c:function:: void ARKodeSPRKTable_Free(ARKodeSPRKTable B)

   Free the memory allocated for the ARKodeSPRKTable structure.

   :param B: The ARKodeSPRKTable structure to free.

.. c:function:: int ARKodeSPRKTable_ToButcher(ARKodeSPRKTable sprk_storage, ARKodeSPRKTable* a_ptr, ARKodeSPRKTable* b_ptr)

   Convert the ARKodeSPRKTable structure to the Butcher table representation.

   :param sprk_storage: The ARKodeSPRKTable structure.
   :param a_ptr: Pointer to store the explicit Butcher table.
   :param b_ptr: Pointer to store the diagonally-implicit Butcher table.
