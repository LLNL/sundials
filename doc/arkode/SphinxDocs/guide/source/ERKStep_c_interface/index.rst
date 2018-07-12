..
   Programmer(s): Daniel R. Reynolds @ SMU
   ----------------------------------------------------------------
   Copyright (c) 2013, Southern Methodist University.
   All rights reserved.
   For details, see the LICENSE file.
   ----------------------------------------------------------------

:tocdepth: 3

.. _ERKStep_CInterface:

==========================================
Using ERKStep for C and C++ Applications
==========================================

This chapter is concerned with the use of the ERKStep time-stepping
module for the solution of nonstiff initial value problems (IVPs) in a
C or C++ language setting.  The following sections discuss the header
files and the layout of the user's main program, and provide
descriptions of the ERKStep user-callable functions and user-supplied
functions.

The example programs described in the companion document [R2013]_ may
be helpful. Those codes may be used as templates for new codes and are
included in the ARKode package ``examples`` subdirectory.

ERKStep uses the input and output constants from the shared ARKode
infrastructure. These are defined as needed in this chapter, but for
convenience the full list is provided separately in the section
:ref:`Constants`.

The relevant information on using ERKStep's C and C++ interfaces is
detailed in the following sub-sections:

.. toctree::
   :maxdepth: 1

   General
   Skeleton
   User_callable
   User_supplied
