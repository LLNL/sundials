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

.. _SUNHeuristics.Default:

The SUNHeuristics_Default Module
======================================

The default implementation of the SUNHeuristics class, SUNHeuristics_Default,
implements a myriad of heuristics control features used throughout SUNDIALS
integrators.

This is implemented as a derived SUNHeuristics class, and defines its *content*
field as:

.. code-block:: c

   struct _SUNHeuristicsContent_Default {
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
specified by "set" routines in the base SUNHeuristics class.  The one
non-obvious of these parameters is ``etamax``, that is changed dynamically
throughout the course of a calculation:

.. math::
   \text{etamax} = \begin{cases}
     \text{etamx1 before the first internal step},\\
     \text{growth following each successful internal step},\\
     \text{1 following either a convergence or error test failure},
   \end{cases}

The header file to be included when using this module is
``sunheuristics/sunheuristics_default.h``.


The SUNHeuristics_Default class provides implementations of all controller
operations listed in :numref:`SUNHeuristics.Description.operations`. The
SUNHeuristics_Default class also provides the following constructor routine:


.. c:function:: SUNHeuristics SUNHeuristicsDefault(SUNContext sunctx)

   This constructor function creates and allocates memory for a
   SUNHeuristics_Default object, and inserts its default parameters.  The only
   argument is the SUNDIALS context object.  Upon successful completion it will
   return a :c:type:`SUNHeuristics` object; otherwise it will return ``NULL``.
