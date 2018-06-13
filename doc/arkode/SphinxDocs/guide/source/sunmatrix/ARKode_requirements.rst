..
   Programmer(s): Daniel R. Reynolds @ SMU
   ----------------------------------------------------------------
   Copyright (c) 2017, Southern Methodist University.
   All rights reserved.
   For details, see the LICENSE file.
   ----------------------------------------------------------------

:tocdepth: 3



.. _SUNMatrix.ARKode:

SUNMATRIX functions required by ARKode
==========================================

In Table :ref:`SUNMatrix.ARKode_Use`, we list the matrix functions in
the ``SUNMatrix`` module used within the ARKode package.  The table
also shows, for each function, which of the code modules uses the
function.  Neither the main ARKode integrator or the ARKSPILS
interface call ``SUNMatrix`` functions directly, so the table columns
are specific to the ARKDLS direct solver interface and the
ARKBANDPRE and ARKBBDPRE preconditioner modules. 

At this point, we should emphasize that the ARKode user does not need
to know anything about the usage of matrix functions by the ARKode
code modules in order to use ARKode.  The information is presented as
an implementation detail for the interested reader.


.. _SUNMatrix.ARKode_Use:

List of matrix functions usage by ARKode code modules
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. cssclass:: table-bordered

==================  ======  ==========  =========
Routine             ARKDLS  ARKBANDPRE  ARKBBDPRE
==================  ======  ==========  =========
SUNMatGetID         X       
SUNMatClone         X       
SUNMatDestroy       X       X           X
SUNMatZero          X       X           X
SUNMatCopy          X       X           X
SUNMatScaleAddI     X       X           X
SUNMatScaleAdd      1
SUNMatMatvec        1
SUNMatSpace         2       2           2
==================  ======  ==========  =========

1. These matrix functions are only used for problems involving a
   non-identity mass matrix.

2. These matrix functions are optionally used, in that these are only
   called if they are implemented in the ``SUNMatrix`` module that is
   being used (i.e. their function pointers are non-``NULL``).






