.. ----------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2002-2023, Lawrence Livermore National Security
   and Southern Methodist University.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   ----------------------------------------------------------------

.. _ARKODE.Usage.SPRKStep:

==========================================
Using the SPRKStep time-stepping module
==========================================

This chapter is concerned with the use of the SPRKStep time-stepping module for
the solution of Hamiltonian initial value problems (IVPs) of the form
:eq:`ARKODE_IVP_Hamiltonian` in a C or C++ language setting. The following sections
discuss the header files and the layout of the user's main program, and provide
descriptions of the SPRKStep user-callable functions and user-supplied functions.

The example programs located in the source code ``examples/arkode`` folder, may 
be helpful as templates for new codes. In particular,

* ``examples/arkode/C_serial/ark_harmonic_symplectic.c``
* ``examples/arkode/C_serial/ark_damped_harmonic_symplectic.c``, and
* ``examples/arkode/C_serial/ark_kepler.c``

demonstrate ``SPRKStep`` usage. 

SPRKStep uses the input and output constants from the shared ARKODE infrastructure.
These are defined as needed in this chapter, but for convenience the full list is
provided separately in :numref:`ARKODE.Constants`.

The relevant information on using SPRKStep's C and C++ interfaces is
detailed in the following subsections.

.. toctree::
   :maxdepth: 1

   Skeleton
   User_callable
