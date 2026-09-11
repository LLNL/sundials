..
   Copyright (c) 2002-2026, Lawrence Livermore National Security and
   the SUNDIALS contributors. SPDX-License-Identifier: BSD-3-Clause

.. _CVODE.Examples.Parallel:

Parallel example problems
=========================

A nonstiff example: cvAdvDiff_non_p
-----------------------------------

This MPI program solves

.. math::

   u_t=u_{xx}+0.5u_x,\qquad 0<x<2,

with zero boundary values and :math:`u(0,x)=x(2-x)e^{2x}`. Centered
differences produce a nonstiff ODE system solved with the Adams method and
fixed-point iteration. The distributed vector is partitioned into contiguous
subdomains; each right-hand-side evaluation exchanges neighboring endpoint
values before evaluating local differences.

.. literalinclude:: ../../../../examples/cvode/parallel/cvAdvDiff_non_p.out
   :language: none

A user-preconditioner example: cvDiurnal_kry_p
----------------------------------------------

This is an MPI implementation of :ref:`CVODE.Examples.Serial.Diurnal`.
Processes own rectangular submeshes and exchange ghost-cell data before local
right-hand-side evaluation. The preconditioner retains the local block-diagonal
part of the Newton matrix.

.. literalinclude:: ../../../../examples/cvode/parallel/cvDiurnal_kry_p.out
   :language: none

A CVBBDPRE example: cvDiurnal_kry_bbd_p
---------------------------------------

This variant replaces the user preconditioner with
:ref:`CVBBDPRE <CVODE.Usage.CC.precond.cvbbdpre>`. A local approximation function
and communication callback allow CVBBDPRE to construct difference-quotient
banded blocks independently on each process.

.. literalinclude:: ../../../../examples/cvode/parallel/cvDiurnal_kry_bbd_p.out
   :language: none
