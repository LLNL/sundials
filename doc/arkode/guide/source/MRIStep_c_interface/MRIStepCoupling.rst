..
   Programmer(s): Daniel R. Reynolds @ SMU
   ----------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2002-2020, Lawrence Livermore National Security
   and Southern Methodist University.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   ----------------------------------------------------------------

:tocdepth: 3


.. _MRIStepCoupling:

MRI Coupling Coefficients Data Structure
-------------------------------------------

As mentioned in the section :ref:`MRIStep_CInterface.UserCallable`, a user may supply a custom set of coupling
coefficients between fast and slow time scales (through :c:func:`MRIStepSetCoupling()`).  MRIStep uses a custom
data type, :c:type:`MRIStepCoupling`, to store these coefficients, and provides several related utility
routines.  As described in the Section :ref:`Mathematics.MRIStep`, the coupling from slow to fast time scales
in MRI methods may be encoded by a vector of slow 'stage time' abcissae, :math:`c^S \in \mathbb{R}^{s+1}` and a
set of coupling matrices :math:`\Gamma^{\{k\}}\in\mathbb{R}^{(s+1)\times(s+1)}`.  We therefore define the
:c:type:`MRIStepCoupling` structure to be

.. c:type:: typedef MRIStepCouplingMem *MRIStepCoupling

where ``MRIStepCouplingMem`` is the structure

.. code-block:: c

   struct MRIStepCouplingMem {

     int nmat;
     int stages;
     int q;
     int p;
     realtype ***G;
     realtype *c;

   };

Here,

* ``nmat`` corresponds to the number of :math:`\Gamma^{\{k\}}` matrices used for coupling,
* ``stages`` is the number of entries in ``c``, i.e., :math:`s+1` above,
* ``q`` and ``p`` indicate the orders of accuracy for both the MRI method and the embedding, respectively,
* ``G`` is a three-dimensional array with dimensions ``[nmat][stages][stages]`` containing the set of
  :math:`\Gamma^{\{k\}}` matrices, and
* ``c`` is an array of length ``stages`` containing slow abcissae :math:`c^S` for the MRI method.


.. _MRIStepCoupling.Functions:

MRIStepCoupling functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. cssclass:: table-bordered

===================================== ==============================================================
Function name                         Description
===================================== ==============================================================
:c:func:`MRIStepCoupling_LoadTable()` Loads a pre-defined MRIStepCoupling table
:c:func:`MRIStepCoupling_Alloc()`     Allocate an empty MRIStepCoupling table
:c:func:`MRIStepCoupling_Create()`    Create a new MRIStepCoupling table from coefficients
:c:func:`MRIStepCoupling_MIStoMRI()`  Create a new MRIStepCoupling table from a slow Butcher table
:c:func:`MRIStepCoupling_Copy()`      Create a copy of a MRIStepCoupling table
:c:func:`MRIStepCoupling_Space()`     Get the MRIStepCoupling table real and integer workspace sizes
:c:func:`MRIStepCoupling_Free()`      Deallocate a MRIStepCoupling table
:c:func:`MRIStepCoupling_Write()`     Write the MRIStepCoupling table to an output file
===================================== ==============================================================

.. c:function:: MRIStepCoupling MRIStepCoupling_LoadTable(int imethod)

   Retrieves a specified MRIStepCoupling table. The prototype for this function, as well as the
   integer names for each provided method, are defined in top of the header file
   ``arkode/arkode_mristep.h``.  For further information on the current set of coupling tables
   and their corresponding identifiers, see :ref:`MRIStepCoupling.Tables`.

   **Arguments:**
      * *itable* -- MRIStepCoupling table identifier to load.

   **Return value:**
      * :c:type:`MRIStepCoupling` structure if successful.
      * ``NULL`` pointer if *itable* was invalid or an allocation error occured.

        
.. c:function:: MRIStepCoupling MRIStepCoupling_Alloc(int nmat, int stages)

   Allocates an empty MRIStepCoupling table.

   **Arguments:**
      * *nmat* -- number of :math:`\Gamma^{\{k\}}` matrices in the coupling table.
      * *stages* -- number of stages in the coupling table.

   **Return value:**
      * :c:type:`MRIStepCoupling` structure if successful.
      * ``NULL`` pointer if *stages* was invalid or an allocation error occured.

        
.. c:function:: MRIStepCoupling MRIStepCoupling_Create(int nmat, int stages, int q, int p, realtype *G, realtype *c)

   Allocates an MRIStepCoupling table and fills it with the given values.

   **Arguments:**
      * *nmat* -- number of :math:`\Gamma^{\{k\}}` matrices in the coupling table.
      * *stages* -- number of stages in the MRI method.
      * *q* -- global order of accuracy for the MRI method.
      * *p* -- global order of accuracy for the embedded MRI method.
      * *G* -- array of coefficients defining the :math:`\Gamma^{\{k\}}` matrices. This should be
        stored as a 1D array of size *nmat*stages*stages*, in row-major order.
      * *c* -- array (of length *stages*) of slow abcissae for the MRI method.

   **Return value:**
      * :c:type:`MRIStepCoupling` structure if successful.
      * ``NULL`` pointer if *stages* was invalid or an allocation error occured.

   **Notes:** As embeddings are not currently supported in MRIStep, then *p* should be equal to zero.

   
.. c:function:: MRIStepCoupling MRIStepCoupling_MIStoMRI(ARKodeButcherTable B, int q, int p)

   Creates an MRI coupling table for a traditional MIS method based on the slow Butcher table *B*, following
   the formula shown in :eq:`MIS_to_MRI`.

   **Arguments:**
      * *B* -- the :c:type:`ARKodeButcherTable` for the 'slow' MIS method.
      * *q* -- the overall order of the MIS/MRI method.
      * *p* -- the overall order of the MIS/MRI embedding.

   **Return value:**
      * :c:type:`MRIStepCoupling` structure if successful.
      * ``NULL`` pointer an allocation error occured.

   **Notes:** The :math:`s`-stage slow Butcher table must have an explicit first stage (i.e., :math:`c_1=0` 
   and :math:`A_{1,j}=0` for :math:`1\le j\le s`) and sorted abcissae (i.e., :math:`c_{i} \ge  c_{i-1}` for
   :math:`2\le i\le s`).

   Since an MIS method is at most third order accurate, and even then only if it meets certain compatibility
   criteria (see :eq:`MIS_order3`), the values of *q* and *p* may differ from the method and embedding orders
   of accuracy for the Runge--Kutta encoded in *B*, which is why these arguments should be supplied separately.

   As embeddings are not currently supported in MRIStep, then *p* should be equal to zero.
   
        
.. c:function:: MRIStepCoupling MRIStepCoupling_Copy(MRIStepCoupling C)

   Creates copy of the given coupling table.

   **Arguments:**
      * *C* -- the coupling table to copy.

   **Return value:**
      * :c:type:`MRIStepCoupling` structure if successful.
      * ``NULL`` pointer an allocation error occured.

        
.. c:function:: void MRIStepCoupling_Space(MRIStepCoupling C, sunindextype *liw, sunindextype *lrw)

   Get the real and integer workspace size for an MRIStepCoupling table.

   **Arguments:**
      * *C* -- the coupling table.
      * *lenrw* -- the number of ``realtype`` values in the coupling table workspace.
      * *leniw* -- the number of integer values in the coupling table workspace.

   **Return value:**
      * *ARK_SUCCESS* if successful.
      * *ARK_MEM_NULL* if the Butcher table memory was ``NULL``.

        
.. c:function:: void MRIStepCoupling_Free(MRIStepCoupling C)

   Deallocate the coupling table memory.

   **Arguments:**
      * *C* -- the MRIStepCoupling table.

        
.. c:function:: void MRIStepCoupling_Write(MRIStepCoupling C, FILE *outfile)

   Write the coupling table to the provided file pointer.

   **Arguments:**
      * *C* -- the MRIStepCoupling table.
      * *outfile* -- pointer to use for printing the table.

   **Notes:** The *outfile* argument can be ``stdout`` or ``stderr``, or it
   may point to a specific file created using ``fopen``.





.. _MRIStepCoupling.Tables:

MRIStepCoupling tables
^^^^^^^^^^^^^^^^^^^^^^^^^^^

MRIStep currently includes two classes of MRI coupling tables: those that encode methods that are
explicit at the slow time scale, and those that are diagonally-implicit and solve-decoupled at
the slow time scale.  We list the current identifiers, multirate order of accuracy, and relevant
references for each in the tables below.  For the implicit methods, we also list the number of
implicit solves per step that are required at the slow time scale.


.. cssclass:: table-bordered

Explicit MRI coupling tables:

=========================== ===== ===================
Table name                  Order Reference
=========================== ===== ===================
``MIS_KW3``                 3     [SKAW2009]_
``MRI_GARK_ERK45a``         4     [S2019]_
=========================== ===== ===================


Diagonally-implicit, solve-decoupled MRI coupling tables:

=========================== ===== =============== ===================
Table name                  Order Implicit Solves Reference
=========================== ===== =============== ===================
``MRI_GARK_IRK21a``         2     1               [S2019]_
``MRI_GARK_ESDIRK34a``      4     3               [S2019]_
=========================== ===== =============== ===================


   
