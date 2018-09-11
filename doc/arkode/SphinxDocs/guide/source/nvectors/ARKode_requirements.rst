..
   Programmer(s): Daniel R. Reynolds @ SMU
   ----------------------------------------------------------------
   Copyright (c) 2013, Southern Methodist University.
   All rights reserved.
   For details, see the LICENSE file.
   ----------------------------------------------------------------

:tocdepth: 3



.. _NVectors.ARKode:

NVECTOR functions required by ARKode
==========================================

In the table below, we list the vector functions in the ``N_Vector``
module that are called within the ARKode package.  The table also
shows, for each function, which ARKode module uses the function.  
The ARKode column shows function usage within the main time-stepping
modules and the shared ARKode infrastructure,  while the remaining
columns show function usage within the ARKode linear solver
interfaces, the ARKBANDPRE and ARKBBDPRE preconditioner modules, and
the FARKODE module.  Here ARKDLS stands for the direct linear solver
interface in ARKode, and ARKSPILS stands for the scaled,
preconditioned, iterative linear solver interface in ARKode.

At this point, we should emphasize that the user does not need to know
anything about ARKode's usage of vector functions in order to use
ARKode.  Instead, this information is provided primarily for users
interested in constructing a custom ``N_Vector`` module.  We note that
a number of ``N_Vector`` functions from the section
:ref:`NVectors.Description` are not listed in the above table.
Therefore a user-supplied ``N_Vector`` module for ARKode could safely
omit these functions from their implementation. 



.. cssclass:: table-bordered

====================  =============  ======  ===========  ==========  =========  =============
Routine               ARKode         ARKDLS  ARKSPILS     ARKBANDPRE  ARKBBDPRE  FARKODE
====================  =============  ======  ===========  ==========  =========  =============
N_VAbs                X                                                          X
N_VAddConst           X                                                          X
N_VClone              X                      X                                   X
N_VCloneEmpty                                                                    X
N_VConst              X              X       X                                   X
N_VDestroy            X                      X                                   X
N_VDiv                X                      X                                   X
N_VDotProd            X\ :sup:`1`            X                                   X\ :sup:`1`
N_VGetArrayPointer                   X                    X           X          X
N_VGetVectorID  
N_VInv                X                                                          X
N_VLinearSum          X              X       X                                   X
N_VMaxNorm            X                                                          X
N_VMin                X                                                          X
N_VProd                                      X   
N_VScale              X              X       X            X           X          X
N_VSetArrayPointer                   X                                           X
N_VSpace              X\ :sup:`2`                                                X\ :sup:`2`
N_VWrmsNorm           X              X       X            X           X          X
N_VLinearCombination  X                      X   
N_VDotProdMulti       X\ :sup:`3`            X\ :sup:`3`
====================  =============  ======  ===========  ==========  =========  =============

1. The :c:func:`N_VDotProd()` function is only used by the main
   ARKode integrator module when the fixed-point nonlinear solver is
   specified; when solving an explicit problem or when using a Newton
   solver with direct linear solver, it need not be
   supplied by the ``N_Vector`` implementation.

2. The :c:func:`N_VSpace()` function is only informational, and need
   not be supplied by the ``N_Vector`` implementation.

3. The optional function :c:func:`N_VDotProdMulti()` is only used when
   the Anderson acceleration fixed point solver is used or when
   classical Gram-Schmidt is used with the SPGMR or SPFGMR iterative
   linear solvers.  The remaining fused vector and vector array
   operations are unused and a user-supplied NVECTOR module for ARKode
   could omit these operations.
