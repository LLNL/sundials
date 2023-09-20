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

.. _SUNTimestepHeuristics.Unconstrained:

The SUNTimestepHeuristics_Unconstrained Module
==============================================

We provide an "unconstrained" implementation of the SUNTimestepHeuristics class,
SUNTimestepHeuristics_Unconstrained.  While this object may indeed be called by
integrators, it implements essentially no heuristics constraints on the step
size.

This is implemented as a derived SUNTimestepHeuristics class, and defines its *content*
field as:

.. code-block:: c

   struct SUNTimestepHeuristicsContent_Unconstrained_ { long int nst_acc; };

Note that the only entry in the *content* field is a counter, that keeps track
of the number of calls made to its :c:func:`SUNTimestepHeuristicsConstrainStep` routine.
Since this class implements essentially no heuristic time step control, it
relies on the base class implementation for most of the operations from
:numref:`SUNTimestepHeuristics.Description.operations`, with a few notable exceptions.

* It self-identifies as having :c:type:`SUNTimestepHeuristics_ID` type as
  ``SUN_TIMESTEPHEURISTICS_NULL``.

* At each call to :c:func:`SUNTimestepHeuristicsConstrainStep` it updates the
  ``nst_acc`` counter and sets the "constrained" step to be the same value as
  what was requested.

* Its :c:func:`SUNTimestepHeuristicsConvFail` routine does not adjust the input step
  size, and returns the failure code ``SUNHEURISTICS_CANNOT_DECREASE``.

* Its :c:func:`SUNTimestepHeuristicsReset` routine resets ``nst_acc`` to zero.

* Its :c:func:`SUNTimestepHeuristicsWrite`, :c:func:`SUNTimestepHeuristicsGetNumAccSteps` and
  :c:func:`SUNTimestepHeuristicsSpace` routines perform as expected, given the class'
  storage of the ``nst_acc`` counter.

The header file to be included when using this module is
``suntimestepheuristics/suntimestepheuristics_unconstrained.h``.

In addition to the routines mentioned above, the SUNTimestepHeuristics_Unconstrained
class provides the following constructor routine:

.. c:function:: SUNTimestepHeuristics SUNTimestepHeuristicsUnconstrained(SUNContext sunctx)

   This constructor function creates and allocates memory for a
   SUNTimestepHeuristics_Unconstrained object, and initializes its ``nst_acc`` counter
   to zero.  The only argument is the SUNDIALS context object.  Upon successful
   completion it will return a :c:type:`SUNTimestepHeuristics` object; otherwise it will
   return ``NULL``.
