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

.. c:macro:: ARKODE_SYMPLECTIC_EULER_1_1

   Identifier for the Symplectic Euler 1st order method with 1 stage.

.. c:macro:: ARKODE_SYMPLECTIC_LEAPFROG_2_2

   Identifier for the Symplectic Leapfrog 2nd order method with 2 stages.

.. c:macro:: ARKODE_SYMPLECTIC_PSEUDO_LEAPFROG_2_2

   Identifier for the Symplectic Pseudo Leapfrog 2nd order method with 2 stages.

.. c:macro:: ARKODE_SYMPLECTIC_RUTH_3_3

   Identifier for the Symplectic Ruth 3rd order method with 3 stages.

.. c:macro:: ARKODE_SYMPLECTIC_MCLACHLAN_2_2

   Identifier for the Symplectic McLachlan 2nd order method with 2 stages.

.. c:macro:: ARKODE_SYMPLECTIC_MCLACHLAN_3_3

   Identifier for the Symplectic McLachlan 3rd order method with 3 stages.

.. c:macro:: ARKODE_SYMPLECTIC_CANDY_ROZMUS_4_4

   Identifier for the Symplectic Candy Rozmus 4th order method with 4 stages.

.. c:macro:: ARKODE_SYMPLECTIC_MCLACHLAN_4_4

   Identifier for the Symplectic McLachlan 4th order method with 4 stages.

.. c:macro:: ARKODE_SYMPLECTIC_MCLACHLAN_5_6

   Identifier for the Symplectic McLachlan 5th order method with 6 stages.

.. c:macro:: ARKODE_SYMPLECTIC_YOSHIDA_6_8

   Identifier for the Symplectic Yoshida 6th order method with 8 stages.

.. c:macro:: ARKODE_SYMPLECTIC_MCLACHLAN_8_16

   Identifier for the Symplectic McLachlan 8th order method with 16 stages.

.. c:macro:: ARKODE_SYMPLECTIC_SOFRONIOU_10_36

   Identifier for the Symplectic Sofroniou 10th order method with 36 stages.

.. c:type:: ARKodeSPRKStorage_s

   Structure representing the SPRK method that holds the method coefficients.

   .. c:member:: int q

      The method order of accuracy.

   .. c:member:: int stages
      
      The number of stages.

   .. c:member:: sunrealtype* a

      Array of coefficients that generate the explicit Butcher table.

   .. c:member:: sunrealtype* b

      Array of coefficients that generate the diagonally-implicit Butcher table.

.. c:type:: ARKodeSPRKStorage_s* ARKodeSPRKStorage


ARKodeSPRKStorage functions
---------------------------

.. _ARKodeSPRKStorage.FunctionsTable:
.. table:: ARKodeButcherTable functions

   +----------------------------------------------+------------------------------------------------------------+
   | **Function name**                            | **Description**                                            |
   +----------------------------------------------+------------------------------------------------------------+
   | :c:func:`ARKodeButcherTable_Alloc()`         | Allocate an empty storage structure                        |
   +----------------------------------------------+------------------------------------------------------------+
   | :c:func:`ARKodeButcherTable_Load()`          | Load SPRK method using an identifier                       |
   +----------------------------------------------+------------------------------------------------------------+
   | :c:func:`ARKodeButcherTable_LoadByName()`    | Load SPRK method using a string version of the identifier  |
   +----------------------------------------------+------------------------------------------------------------+
   | :c:func:`ARKodeButcherTable_Create()`        | Create a new storage structure                             |
   +----------------------------------------------+------------------------------------------------------------+
   | :c:func:`ARKodeButcherTable_Copy()`          | Create a copy of a storage structure                       |
   +----------------------------------------------+------------------------------------------------------------+
   | :c:func:`ARKodeButcherTable_Space()`         | Get the storage structure real and integer workspace size  |
   +----------------------------------------------+------------------------------------------------------------+
   | :c:func:`ARKodeButcherTable_Free()`          | Deallocate a storage structure                             |
   +----------------------------------------------+------------------------------------------------------------+


.. c:function:: ARKodeSPRKStorage ARKodeSPRKMem_Alloc(int stages)

   Allocate memory for the ARKodeSPRKStorage structure.

   :param stages: The number of stages.
   :return: Pointer to the allocated ARKodeSPRKStorage structure.

.. c:function:: ARKodeSPRKStorage ARKodeSPRKMem_Load(ARKODE_SPRKMethodID id)

   Load the ARKodeSPRKStorage structure for the specified method ID.

   :param id: The ID of the SPRK method. One of :ref:`SPRKStorage.id`.
   :return: Pointer to the loaded ARKodeSPRKStorage structure.

.. c:function:: ARKodeSPRKStorage ARKodeSPRKMem_LoadByName(const char* method)

   Load the ARKodeSPRKStorage structure for the specified method name.

   :param method: The name of the SPRK method. Must be one of :ref:`SPRKStorage.id` but as a string.
   :return: Pointer to the loaded ARKodeSPRKStorage structure.


.. c:function:: ARKodeSPRKStorage ARKodeSPRKMem_Copy(ARKodeSPRKStorage B)

   Create a copy of the ARKodeSPRKStorage structure.

   :param B: The original ARKodeSPRKStorage structure.
   :return: Pointer to the copied ARKodeSPRKStorage structure.

.. c:function:: void ARKodeSPRKStorage_Space(ARKodeSPRKStorage B, sunindextype* liw, sunindextype* lrw)

   Get the workspace sizes required for the ARKodeSPRKStorage structure.

   :param B: The ARKodeSPRKStorage structure.
   :param liw: Pointer to store the integer workspace size.
   :param lrw: Pointer to store the real workspace size.

.. c:function:: void ARKodeSPRKMem_Free(ARKodeSPRKStorage B)

   Free the memory allocated for the ARKodeSPRKStorage structure.

   :param B: The ARKodeSPRKStorage structure to free.

.. c:function:: int ARKodeSPRKMem_ToButcher(ARKodeSPRKStorage sprk_storage, ARKodeButcherTable* a_ptr, ARKodeButcherTable* b_ptr)

   Convert the ARKodeSPRKStorage structure to the Butcher table representation.

   :param sprk_storage: The ARKodeSPRKStorage structure.
   :param a_ptr: Pointer to store the explicit Butcher table.
   :param b_ptr: Pointer to store the diagonally-implicit Butcher table.

.. c:function:: ARKodeSPRKStorage ARKodeSymplecticEuler()

   Create the ARKodeSPRKStorage structure for the Symplectic Euler method.

.. c:function:: ARKodeSPRKStorage ARKodeSymplecticLeapfrog2()

   Create the ARKodeSPRKStorage structure for the Symplectic Leapfrog 2-2 method.

.. c:function:: ARKodeSPRKStorage ARKodeSymplecticPseudoLeapfrog2()

   Create the ARKodeSPRKStorage structure for the Symplectic Pseudo Leapfrog 2-2 method.

.. c:function:: ARKodeSPRKStorage ARKodeSymplecticRuth3()

   Create the ARKodeSPRKStorage structure for the Symplectic Ruth 3-3 method.

.. c:function:: ARKodeSPRKStorage ARKodeSymplecticCandyRozmus4()

   Create the ARKodeSPRKStorage structure for the Symplectic Candy Rozmus 4-4 method.

.. c:function:: ARKodeSPRKStorage ARKodeSymplecticMcLachlan2()

   Create the ARKodeSPRKStorage structure for the Symplectic McLachlan 2-2 method.

.. c:function:: ARKodeSPRKStorage ARKodeSymplecticMcLachlan3()

   Create the ARKodeSPRKStorage structure for the Symplectic McLachlan 3-3 method.

.. c:function:: ARKodeSPRKStorage ARKodeSymplecticMcLachlan4()

   Create the ARKodeSPRKStorage structure for the Symplectic McLachlan 4-4 method.

.. c:function:: ARKodeSPRKStorage ARKodeSymplecticMcLachlan5()

   Create the ARKodeSPRKStorage structure for the Symplectic McLachlan 5-6 method.

.. c:function:: ARKodeSPRKStorage ARKodeSymplecticYoshida6()

   Create the ARKodeSPRKStorage structure for the Symplectic Yoshida 6-8 method.

.. c:function:: ARKodeSPRKStorage ARKodeSymplecticMcLachlan8()

   Create the ARKodeSPRKStorage structure for the Symplectic McLachlan 8-16 method.

.. c:function:: ARKodeSPRKStorage ARKodeSymplecticSofroniou10()

   Create the ARKodeSPRKStorage structure for the Symplectic Sofroniou 10-36 method.
