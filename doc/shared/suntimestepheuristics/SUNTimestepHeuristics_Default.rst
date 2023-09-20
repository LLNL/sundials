..
   Programmer(s): Daniel R. Reynolds @ SMU
   ----------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2002-2023, Lawrence Livermore National Security
   and Southern Methodist University.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   ----------------------------------------------------------------

.. _SUNTimestepHeuristics.Default:

The SUNTimestepHeuristics_Default Module
========================================

The default implementation of the SUNTimestepHeuristics class, SUNTimestepHeuristics_Default,
implements a myriad of heuristics control features used throughout SUNDIALS
integrators.

This is implemented as a derived SUNTimestepHeuristics class, and defines its *content*
field as:

.. code-block:: c

   struct SUNTimestepHeuristicsContent_Default_ {
     realtype     hmax_inv;
     realtype     hmin;
     realtype     etamax;
     realtype     etamx1;
     realtype     etamxf;
     realtype     etamin;
     int          small_nef;
     realtype     etacf;
     SUNExpStabFn expstab;
     void*        estab_data;
     realtype     cfl;
     realtype     safety;
     realtype     growth;
     realtype     lbound;
     realtype     ubound;
     long int     nst_acc;
     long int     nst_exp;
   };

These entries of the *content* field contain parameters to store all options
specified by "set" routines in the base SUNTimestepHeuristics class.  The one
non-obvious of these parameters is ``etamax``, that is changed dynamically
throughout the course of a calculation:

.. math::
   \text{etamax} = \begin{cases}
     \text{etamx1 before the first internal step},\\
     \text{growth following each successful internal step},\\
     \text{1 following either a convergence or error test failure},
   \end{cases}

The header file to be included when using this module is
``suntimestepheuristics/suntimestepheuristics_default.h``.

The SUNTimestepHeuristics_Default class provides implementations of all controller
operations listed in :numref:`SUNTimestepHeuristics.Description.operations`. The
SUNTimestepHeuristics_Default class also provides the following constructor routine:

.. c:function:: SUNTimestepHeuristics SUNTimestepHeuristicsDefault(SUNContext sunctx)

   This constructor function creates and allocates memory for a
   SUNTimestepHeuristics_Default object, and inserts its default parameters.  The only
   argument is the SUNDIALS context object.  Upon successful completion it will
   return a :c:type:`SUNTimestepHeuristics` object; otherwise it will return ``NULL``.
