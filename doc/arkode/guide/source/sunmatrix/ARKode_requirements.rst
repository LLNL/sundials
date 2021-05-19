..
   Programmer(s): Daniel R. Reynolds @ SMU
   ----------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2002-2021, Lawrence Livermore National Security
   and Southern Methodist University.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   ----------------------------------------------------------------

:tocdepth: 3



.. _SUNMatrix.ARKode:

SUNMATRIX functions required by ARKode
==========================================

In Table :ref:`SUNMatrix.ARKode_Use`, we list the matrix functions in
the ``SUNMatrix`` module used within the ARKode package.  The table
also shows, for each function, which of the code modules uses the
function.  The main ARKode time step modules, ARKStep and ERKStep, do
not call any ``SUNMatrix`` functions directly, so the table columns
are specific to the ARKLS interface and the ARKBANDPRE and ARKBBDPRE
preconditioner modules.   We further note that the ARKLS interface
only utilizes these routines when supplied with a *matrix-based*
linear solver, i.e. the ``SUNMatrix`` object (*J* or *M*) passed to 
:c:func:`ARKStepSetLinearSolver()` or
:c:func:`ARKStepSetMassLinearSolver()` was not ``NULL``.

At this point, we should emphasize that the ARKode user does not need
to know anything about the usage of matrix functions by the ARKode
code modules in order to use ARKode.  The information is presented as
an implementation detail for the interested reader.


.. _SUNMatrix.ARKode_Use:

List of matrix functions usage by ARKode code modules
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. cssclass:: table-bordered

==================  ======  ==========  =========
Routine             ARKLS   ARKBANDPRE  ARKBBDPRE
==================  ======  ==========  =========
SUNMatGetID         X
SUNMatClone         X       
SUNMatDestroy       X       X           X
SUNMatZero          X       X           X
SUNMatCopy          X       X           X
SUNMatScaleAddI     X       X           X
SUNMatScaleAdd      1
SUNMatMatvec        1
SUNMatMatvecSetup   1,2
SUNMatSpace         2       2           2
==================  ======  ==========  =========

1. These matrix functions are only used for problems involving a
   non-identity mass matrix.

2. These matrix functions are optionally used, in that these are only
   called if they are implemented in the ``SUNMatrix`` module that is
   being used (i.e. their function pointers are non-``NULL``).  If not
   supplied, these modules will assume that the matrix requires no
   storage. 


We note that both the ARKBANDPRE and ARKBBDPRE preconditioner modules
are hard-coded to use the SUNDIALS-supplied band ``SUNMatrix`` type,
so the most useful information above for user-supplied ``SUNMatrix``
implementations is the column relating to ARKLS requirements.





