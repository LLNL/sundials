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

.. _SUNDomEigEst.CVODE:

CVODE SUNDomEigEstimator interface
==============================================

In :numref:`SUNDomEigEst.CVODE.Usage`, we list the SUNDomEigEst module functions used
within the CVDEE interface.  We emphasize that the CVODE user does not need to know
detailed usage of dominant eigenvalue estomator functions by the CVODE code modules
in order to use CVODE. The information is presented as an implementation detail for
the interested reader.

.. _SUNDomEigEst.CVODE.Usage:
.. table:: List of SUNDomEigEst functions called by the CVODE dominant eigenvalue
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
   | :c:func:`SUNDomEigEstSetMaxPowerIter`              |          O          |          O          |
   +----------------------------------------------------+---------------------+---------------------+
   | :c:func:`SUNDomEigEstSetNumPreProcess`             |          O          |          O          |
   +----------------------------------------------------+---------------------+---------------------+
   | :c:func:`SUNDomEigEstInitialize`                   |          X          |          X          |
   +----------------------------------------------------+---------------------+---------------------+
   | :c:func:`SUNDomEigEstPreProcess`                   |          O          |          O          |
   +----------------------------------------------------+---------------------+---------------------+
   | :c:func:`SUNDomEigEstComputeHess`\ :sup:`1`        |         N/A         |          X          |
   +----------------------------------------------------+---------------------+---------------------+
   | :c:func:`SUNDomEigEstimate`                        |          X          |          X          |
   +----------------------------------------------------+---------------------+---------------------+
   | :c:func:`SUNDomEigEstNumIters`\ :sup:`2`           |          O          |         N/A         |
   +----------------------------------------------------+---------------------+---------------------+
   | :c:func:`SUNDomEigEstRes`\ :sup:`3`                |          O          |         N/A         |
   +----------------------------------------------------+---------------------+---------------------+
   | :c:func:`SUNDomEigEstFree`\ :sup:`4`               |                     |                     |
   +----------------------------------------------------+---------------------+---------------------+


Notes:

1. :c:func:`SUNDomEigEstComputeHess()` might or might not be required depending on
   ``SUNDomEigEstimator`` implementation that is being used. This flag must be left
   ``NULL`` if it is not applicable for an estimator.

2. :c:func:`SUNDomEigEstNumIters()` is only used to accumulate overall
   iterative estimator statistics.  If it is not implemented by
   the ``SUNDomEigEstimator`` module, then CVDEE will consider all
   estimates as requiring zero iterations.

3. Although :c:func:`SUNDomEigEstRes()` is optional, if it is not
   implemented by the ``SUNDomEigEstimator`` then CVDEE will consider all
   estimates a being *exact*.

4. Although CVDEE does not call :c:func:`SUNDomEigEstFree()`
   directly, this routine should be available for users to call when
   cleaning up from a simulation.

Since there are a wide range of potential SUNDomEigEst use cases, the following
subsections describe some details of the CVDEE interface, in the case that
interested users wish to develop custom SUNDomEigEst modules.
