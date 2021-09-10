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



.. _NVectors.ARKode:

NVECTOR functions required by ARKode
==========================================

In the table below, we list the vector functions in the ``N_Vector``
module that are called within the ARKode package.  The table also
shows, for each function, which ARKode module uses the function.
The ARKSTEP and ERKSTEP columns show function usage within the main
time-stepping modules and the shared ARKode infrastructure,  while the
remaining columns show function usage within the ARKLS linear solver
interface, the ARKBANDPRE and ARKBBDPRE preconditioner modules, and
the FARKODE module.

Note that since FARKODE is built on top of ARKode, and therefore
requires the same ``N_Vector`` routines, in the FARKODE column we only
list the routines that the FARKODE interface directly utilizes.

Note that for ARKLS we only list the ``N_Vector`` routines used
directly by ARKLS, each ``SUNLinearSolver`` module may have additional
requirements that are not listed here.  In addition, specific
``SUNNonlinearSolver`` modules attached to ARKode may have additional
``N_Vector`` requirements.  For additional requirements by specific
``SUNLinearSolver`` and ``SUNNonlinearSolver`` modules, please see the
accompanying sections :ref:`SUNLinSol` and :ref:`SUNNonlinSol`.

At this point, we should emphasize that the user does not need to know
anything about ARKode's usage of vector functions in order to use
ARKode.  Instead, this information is provided primarily for users
interested in constructing a custom ``N_Vector`` module.  We note that
a number of ``N_Vector`` functions from the section
:ref:`NVectors.Description` are not listed in the above table.
Therefore a user-supplied ``N_Vector`` module for ARKode could safely
omit these functions from their implementation (although
some may be needed by ``SUNNonlinearSolver`` or ``SUNLinearSolver``
modules).



.. cssclass:: table-bordered

==============================  =======  =======  ===========  ==========  =========  =======
Routine                         ARKSTEP  ERKSTEP  ARKLS        ARKBANDPRE  ARKBBDPRE  FARKODE
==============================  =======  =======  ===========  ==========  =========  =======
N_VGetLength                                      X
N_VAbs                          X        X
N_VAddConst                     X        X
N_VClone                        X        X        X
N_VCloneEmpty                                                                         X
N_VConst                        X        X        X                                   X
N_VDestroy                      X        X        X                                   X
N_VDiv                          X        X
N_VGetArrayPointer                                X\ :sup:`1`  X           X          X
N_VInv                          X        X
N_VLinearSum                    X        X        X
N_VMaxNorm                      X        X
N_VMin                          X        X                                            X
N_VScale                        X        X        X            X           X
N_VSetArrayPointer                                X\ :sup:`1`                         X
N_VSpace\ :sup:`2`              X        X        X            X           X
N_VWrmsNorm                     X        X        X            X           X
N_VLinearCombination\ :sup:`3`  X        X
N_VMinQuotient\ :sup:`5`        X        X
N_VConstrMask\ :sup:`5`         X        X
N_VCompare\ :sup:`5`            X        X
==============================  =======  =======  ===========  ==========  =========  =======

1. This is only required with dense or band matrix-based linear solver modules,
   where the default difference-quotient Jacobian approximation is used.

2. The :c:func:`N_VSpace()` function is only informational, and will
   only be called if provided by the ``N_Vector`` implementation.

3. The :c:func:`N_VLinearCombination()` function is in fact optional;
   if it is not supplied then :c:func:`N_VLinearSum()` will be used instead.

4. The :c:func:`N_VGetLength()` function is only required when an iterative or
   matrix iterative ``SUNLinearSolver`` module is used.

5. The functions :c:func:`N_VMinQuotient`, :c:func:`N_VConstrMask`, and
   :c:func:`N_VCompare` are only used when inequality constraints are enabled
   and may be omitted if this feature is not used.
