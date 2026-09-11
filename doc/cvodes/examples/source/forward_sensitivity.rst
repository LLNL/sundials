..
   Copyright (c) 2002-2026, Lawrence Livermore National Security and
   the SUNDIALS contributors. SPDX-License-Identifier: BSD-3-Clause

.. _CVODES.Examples.FSA:

Forward sensitivity analysis examples
=====================================

CVODES supports ``CV_SIMULTANEOUS``, ``CV_STAGGERED``, and
``CV_STAGGERED1`` sensitivity methods. Sensitivities may be included in or
excluded from the error test.

A serial nonstiff example: cvsAdvDiff_FSA_non
---------------------------------------------

This example solves

.. math::

   u_t=q_1u_{xx}+q_2u_x,qquad 0\le x\le2,quad0\le t\le5,

with zero boundary values, :math:`u(0,x)=x(2-x)e^{2x}`, and nominal parameters
:math:`q_1=1`, :math:`q_2=0.5`. Central differences on ``MX`` interior points
give

.. math::

   \dot u_i=q_1\frac{u_{i+1}-2u_i+u_{i-1}}{\Delta x^2}
   +q_2\frac{u_{i+1}-u_{i-1}}{2\Delta x}.

CVODES integrates the state and sensitivities
:math:`s^j=\partial u/\partial q_j` with Adams and fixed-point iteration.
``CVodeSensInit1`` selects the sensitivity method; parameter scales and indices
are supplied with ``CVodeSetSensParams``. The example uses internal finite
differences for sensitivity right-hand sides.

.. figure:: figs/cvodes/cvsfwdnonx.png
   :width: 90%
   :align: center

   Solution norm and sensitivities with respect to diffusion and advection.

.. literalinclude:: ../../../../examples/cvodes/serial/cvsAdvDiff_FSA_non_-sensi_sim_t.out
   :language: none

A serial dense example: cvsRoberts_FSA_dns
------------------------------------------

This example computes sensitivities of Robertson kinetics with respect to
:math:`p=(0.04,10^4,3\cdot10^7)`:

.. math::

   \dot y_1&=-p_1y_1+p_2y_2y_3,\\
   \dot y_2&=p_1y_1-p_2y_2y_3-p_3y_2^2,\\
   \dot y_3&=p_3y_2^2.

The program supplies state and sensitivity right-hand sides, a dense Jacobian,
and a custom error-weight function. ``CVodeSensInit1`` initializes three
sensitivity vectors, ``CVodeSetSensParams`` supplies scales, and
``CVodeGetSens`` retrieves results after each integration output.

.. figure:: figs/cvodes/cvsfwddenx.png
   :width: 90%
   :align: center

   :math:`y_1` and its sensitivities to the three reaction rates.

.. literalinclude:: ../../../../examples/cvodes/serial/cvsRoberts_FSA_dns_-sensi_sim_t.out
   :language: none

An MPI example: cvsDiurnal_FSA_kry_p
------------------------------------

This program extends the CVODE diurnal kinetics MPI example to compute
sensitivities with respect to two kinetic-rate parameters. It uses BDF, GMRES,
a user block-diagonal preconditioner, and distributed vectors. State and
sensitivity equations use the same submesh partition and ghost-cell exchanges.

.. figure:: figs/cvodes/cvsfwdkryx_p.png
   :width: 90%
   :align: center

   Diurnal kinetics solution and parameter sensitivities.

.. literalinclude:: ../../../../examples/cvodes/parallel/cvsDiurnal_FSA_kry_p_-sensi_sim_t.out
   :language: none
