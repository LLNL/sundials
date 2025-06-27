.. ----------------------------------------------------------------
   Programmer(s): Mustafa Aggul @ SMU
   ----------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2002-2025, Lawrence Livermore National Security
   and Southern Methodist University.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   ----------------------------------------------------------------

.. _SUNDomEigEst.ARKODE:

ARKODE SUNDomEigEstimator interface
==============================================

In :numref:`SUNDomEigEst.ARKODE.Usage`, we list the SUNDomEigEst module functions used
within the ARKDEE interface.  We emphasize that the ARKODE user does not need to know
detailed usage of dominant eigenvalue estomator functions by the ARKODE code modules
in order to use ARKODE. The information is presented as an implementation detail for
the interested reader.

.. _SUNDomEigEst.ARKODE.Usage:
.. table:: List of SUNDomEigEst functions called by the ARKODE dominant eigenvalue
           estimator interface, depending on the self-identified "id" reported from
           :c:func:`SUNDomEigEstGetID`.  Functions marked with "X" are required;
           functions marked with "O" are only called if they are non-``NULL`` and
           functions marked with "N/A" are not applicable in the ``SUNDomEigEstimator``
           implementation that is being used.
   :align: center

   +----------------------------------------------------+---------------------+---------------------+
   | Routine                                            |   POWER ITERATION   |  ARNOLDI ITERATION  |
   |                                                    |                     |                     |
   +====================================================+=====================+=====================+
   | :c:func:`SUNDomEigEstGetID`                        |          O          |          O          |
   +----------------------------------------------------+---------------------+---------------------+
   | :c:func:`SUNDomEigEstSetATimes`                    |          X          |          X          |
   +----------------------------------------------------+---------------------+---------------------+
   | :c:func:`SUNDomEigEstSetMaxPowerIter`\ :sup:`1`    |          O          |         N/A         |
   +----------------------------------------------------+---------------------+---------------------+
   | :c:func:`SUNDomEigEstSetNumPreProcess`             |          O          |          O          |
   +----------------------------------------------------+---------------------+---------------------+
   | :c:func:`SUNDomEigEstSetTol`\ :sup:`2`             |          O          |         N/A         |
   +----------------------------------------------------+---------------------+---------------------+
   | :c:func:`SUNDomEigEstInitialize`                   |          X          |          X          |
   +----------------------------------------------------+---------------------+---------------------+
   | :c:func:`SUNDomEigEstPreProcess`                   |          O          |          O          |
   +----------------------------------------------------+---------------------+---------------------+
   | :c:func:`SUNDomEigEstComputeHess`\ :sup:`3`        |         N/A         |          X          |
   +----------------------------------------------------+---------------------+---------------------+
   | :c:func:`SUNDomEigEstimate`                        |          X          |          X          |
   +----------------------------------------------------+---------------------+---------------------+
   | :c:func:`SUNDomEigEstNumIters`\ :sup:`4`           |          O          |         N/A         |
   +----------------------------------------------------+---------------------+---------------------+
   | :c:func:`SUNDomEigEstRes`\ :sup:`5`                |          O          |         N/A         |
   +----------------------------------------------------+---------------------+---------------------+
   | :c:func:`SUNDomEigEstFree`\ :sup:`6`               |                     |                     |
   +----------------------------------------------------+---------------------+---------------------+


Notes:

1. :c:func:`SUNDomEigEstSetMaxPowerIter()` might or might not be required depending on
   ``SUNDomEigEstimator`` implementation that is being used. This flag must be left
   ``NULL`` if it is not applicable for an estimator.

2. :c:func:`SUNDomEigEstSetTol()` might or might not be required depending on
   ``SUNDomEigEstimator`` implementation that is being used. This flag must be left
   ``NULL`` if it is not applicable for an estimator.

3. :c:func:`SUNDomEigEstComputeHess()` might or might not be required depending on
   ``SUNDomEigEstimator`` implementation that is being used. This flag must be left
   ``NULL`` if it is not applicable for an estimator.

4. :c:func:`SUNDomEigEstNumIters()` is only used to accumulate overall
   iterative estimator statistics.  If it is not implemented by
   the ``SUNDomEigEstimator`` module, then ARKDEE will consider all
   estimates as requiring zero iterations.

5. Although :c:func:`SUNDomEigEstRes()` is optional, if it is not
   implemented by the ``SUNDomEigEstimator`` then ARKDEE will consider all
   estimates a being *exact*.

6. Although ARKDEE does not call :c:func:`SUNDomEigEstFree()`
   directly, this routine should be available for users to call when
   cleaning up from a simulation.

Since there are a wide range of potential SUNDomEigEst use cases, the following
subsections describe some details of the ARKDEE interface, in the case that
interested users wish to develop custom SUNDomEigEst modules.


.. _SUNDomEigEst.Custom:

Providing a custom SUNDomEigEstimator
-------------------------------------

In certain instances, users may wish to provide a custom SUNDomEigEst
implementation to ARKODE in order to leverage the structure of a problem.  While
the "standard" API for these routines is typically sufficient for most users,
others may need additional ARKODE-specific information on top of what is
provided.  For these purposes, we note the following advanced output functions
available in ARKStep and MRIStep:


**ARKStep advanced outputs**: when solving the Newton nonlinear system of
equations in predictor-corrector form,

.. math::
   \begin{array}{ll}
   G(z_{cor}) \equiv z_{cor} - \gamma f^I\left(t^I_{n,i}, z_{i} \right) - \tilde{a}_i = 0 &\qquad  \text{[$M=I$]},\\
   G(z_{cor}) \equiv M z_{cor} - \gamma f^I\left(t^I_{n,i}, z_{i} \right) - \tilde{a}_i = 0 &\qquad  \text{[$M$ static]},\\
   G(z_{cor}) \equiv M(t^I_{n,i}) (z_{cor} - \tilde{a}_i) - \gamma f^I\left(t^I_{n,i}, z_{i}\right) = 0 &\qquad \text{[$M$ time-dependent]}.
   \end{array}

* :c:func:`ARKStepGetCurrentTime()` -- when called within the computation of a
  step (i.e., within a solve) this returns :math:`t^I_{n,i}`. Otherwise the
  current internal solution time is returned.
* :c:func:`ARKStepGetCurrentState()` -- when called within the computation of a
  step (i.e., within a solve) this returns the current stage vector
  :math:`z_{i} = z_{cor} + z_{pred}`. Otherwise the current internal solution
  is returned.
* :c:func:`ARKStepGetCurrentMassMatrix()` -- returns :math:`M(t)`.


**MRIStep advanced outputs**: when solving the Newton nonlinear system of
equations in predictor-corrector form,

.. math::
   G(z_{cor}) \equiv z_{cor} - \gamma f^I\left(t^S_{n,i}, z_{i}\right) - \tilde{a}_i = 0

* :c:func:`MRIStepGetCurrentTime()` -- when called within the computation of a
  step (i.e., within a solve) this returns :math:`t^S_{n,i}`. Otherwise the
  current internal solution time is returned.
* :c:func:`MRIStepGetCurrentState()` -- when called within the computation of a
  step (i.e., within a solve) this returns the current stage vector
  :math:`z_{i} = z_{cor} + z_{pred}`. Otherwise the current internal solution
  is returned.
