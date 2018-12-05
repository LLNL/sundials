..
   Programmer(s): David J. Gardner @ LLNL
   ----------------------------------------------------------------
   Based on ERKStep by Daniel R. Reynolds @ SMU
   ----------------------------------------------------------------
   LLNS/SMU Copyright Start
   Copyright (c) 2018, Southern Methodist University and
   Lawrence Livermore National Security

   This work was performed under the auspices of the U.S. Department
   of Energy by Southern Methodist University and Lawrence Livermore
   National Laboratory under Contract DE-AC52-07NA27344.
   Produced at Southern Methodist University and the Lawrence
   Livermore National Laboratory.

   All rights reserved.
   For details, see the LICENSE file.
   LLNS/SMU Copyright End
   ----------------------------------------------------------------

:tocdepth: 3

.. _MRIStep_CInterface:

==========================================
Using MRIStep for C and C++ Applications
==========================================

This chapter is concerned with the use of the MRIStep time-stepping
module for the solution of two-rate initial value problems (IVPs) in a
C or C++ language setting. The following sections discuss the header
files and the layout of the user's main program, and provide
descriptions of the MRIStep user-callable functions and user-supplied
functions.

The example programs described in the companion document [R2018]_ may
be helpful. Those codes may be used as templates for new codes and are
included in the ARKode package ``examples`` subdirectory.

MRIStep uses the input and output constants from the shared ARKode
infrastructure. These are defined as needed in this chapter, but for
convenience the full list is provided separately in the section
:ref:`Constants`.

The relevant information on using MRIStep's C and C++ interfaces is
detailed in the following sub-sections.

.. toctree::
   :maxdepth: 1

   General
   Skeleton
   User_callable
   User_supplied
