..
   Programmer(s): Daniel R. Reynolds @ SMU
   ----------------------------------------------------------------
   Copyright (c) 2017, Southern Methodist University.
   All rights reserved.
   For details, see the LICENSE file.
   ----------------------------------------------------------------

:tocdepth: 3



.. _SUNLinSol.ARKode:

SUNLinearSolver functions required by ARKode
==============================================

In the table below, we list the linear solver functions in the
``SUNLinearSolver`` module used within the ARKode package. 
The table also shows, for each function, which of the code modules uses
the function.  In general, the main ARKode integrator considers
three categories of linear solvers, *direct*, *iterative*
and *custom*, with interfaces accessible in the ARKode header
files ``arkode/arkode_direct.h`` (ARKDLS), ``arkode/arkode_spils.h``
(ARKSPILS) and ``arkode/arkode_customls.h`` (ARKCLS), respectively. 
Hence, the the table columns reference the use of ``SUNLinearSolver``
functions by each of these solver interfaces.

As with the ``SUNMatrix`` module, we emphasize that the ARKode user
does not need to know detailed usage of linear solver functions by the
ARKode code modules in order to use ARKode. The information is
presented as an implementation detail for the interested reader.
vector functions in the ``N_Vector``



.. cssclass:: table-bordered


==========================  ======  ========  ======
Routine                     ARKDLS  ARKSPILS  ARKCLS
==========================  ======  ========  ======
SUNLinSolGetType            X       X         O
SUNLinSolSetATimes                  X         O
SUNLinSolSetPreconditioner          X         O
SUNLinSolSetScalingVectors          X         O
SUNLinSolInitialize         X       X         X
SUNLinSolSetup              X       X         X
SUNLinSolSolve              X       X         X
SUNLinSolNumIters                   X         O
SUNLinSolResNorm                    X         O
SUNLinSolResid              
SUNLinSolLastFlag           
SUNLinSolFree               X       X         X
SUNLinSolSpace              O       O         O
==========================  ======  ========  ======


The linear solver functions listed above with a "O" are optionally
used, in that these are only called if they are implemented in the
``SUNLinearSolver`` module that is being used  (i.e. their function
pointers are non-``NULL``).  Also, although ARKode does not call
:c:func:`SUNLinSolLastFlag()` directly, this routine is available for
users to query linear solver issues directly. 


