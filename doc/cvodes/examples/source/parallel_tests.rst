..
   Copyright (c) 2002-2026, Lawrence Livermore National Security and
   the SUNDIALS contributors. SPDX-License-Identifier: BSD-3-Clause

.. _CVODES.Examples.ParallelTests:

Parallel tests
==============

These historical scaling tests use a two-species diurnal kinetics
advection-diffusion PDE :cite:p:`SeHi:05,Wit:96`. Central differences are used
except for a biased three-point upwind approximation of horizontal advection.
CVODES applies BDF, GMRES, and a block-diagonal left preconditioner.

The fixed global problem contains 1,280,000 equations. Runs compare state-only
integration with staggered sensitivity analysis for the horizontal and vertical
diffusion coefficients, both without and with sensitivity error control.

.. list-table:: Historical run times in seconds
   :header-rows: 1

   * - Processes
     - States
     - Staggered
     - Staggered with error control
   * - 4
     - 460.31
     - 1414.53
     - 2208.14
   * - 8
     - 211.20
     - 646.59
     - 1064.94
   * - 16
     - 97.16
     - 320.78
     - 417.95
   * - 32
     - 42.78
     - 137.51
     - 210.84
   * - 64
     - 19.50
     - 63.34
     - 83.24
   * - 128
     - 13.78
     - 42.71
     - 55.17
   * - 256
     - 9.87
     - 31.33
     - 47.95

The departure from ideal scaling reflects decreasing preconditioner quality,
increasing communication, and, at low process counts, larger local
factorizations and memory pressure.

.. figure:: figs/cvodes/pvfktTest.png
   :width: 60%
   :align: center

   State-only and staggered-sensitivity scaling results.
