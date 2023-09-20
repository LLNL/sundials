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

.. _SUNHeuristics.Unconstrained:

The SUNHeuristics_Unconstrained Module
======================================

We provide an "unconstrained" implementation of the SUNHeuristics class,
SUNHeuristics_Unconstrained.  While this object may indeed be called by
integrators, it implements essentially no heuristics constraints on the step
size.

This is implemented as a derived SUNHeuristics class, and defines its *content*
field as:

.. code-block:: c

   struct _SUNHeuristicsContent_Unconstrained { long int nst_acc; };

Note that the only entry in the *content* field is a counter, that keeps track
of the number of calls made to its :c:func:`SUNHeuristicsConstrainStep` routine.
Since this class implements essentially no heuristic time step control, it
relies on the base class implementation for most of the operations from
:numref:`SUNHeuristics.Description.operations`, with a few notable exceptions.

* It self-identifies as having :c:type:`SUNHeuristics_ID` type as
  ``SUNDIALS_HEURISTICS_NULL``.

* At each call to :c:func:`SUNHeuristicsConstrainStep` it updates the
  ``nst_acc`` counter and sets the "constrained" step to be the same value as
  what was requested.

* Its :c:func:`SUNHeuristicsConvFail` routine does not adjust the input step
  size, and returns the failure code ``SUNHEURISTICS_CANNOT_DECREASE``.

* Its :c:func:`SUNHeuristicsReset` routine resets ``nst_acc`` to zero.

* Its :c:func:`SUNHeuristicsWrite`, :c:func:`SUNHeuristicsGetNumAccSteps` and
  :c:func:`SUNHeuristicsSpace` routines perform as expected, given the class'
  storage of the ``nst_acc`` counter.

The header file to be included when using this module is
``sunheuristics/sunheuristics_unconstrained.h``.

In addition to the routines mentioned above, the SUNHeuristics_Unconstrained
class provides the following constructor routine:

.. c:function:: SUNHeuristics SUNHeuristicsUnconstrained(SUNContext sunctx)

   This constructor function creates and allocates memory for a
   SUNHeuristics_Unconstrained object, and initializes its ``nst_acc`` counter
   to zero.  The only argument is the SUNDIALS context object.  Upon successful
   completion it will return a :c:type:`SUNHeuristics` object; otherwise it will
   return ``NULL``.
