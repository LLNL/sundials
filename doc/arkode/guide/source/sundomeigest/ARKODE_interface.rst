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
detailed usage of dominant eigenvalue estimator functions by the ARKODE code modules
in order to use ARKODE. The information is presented as an implementation detail for
the interested reader.

.. _SUNDomEigEst.ARKODE.Usage:
.. table:: List of SUNDomEigEst functions called by the ARKODE dominant eigenvalue
           estimator interface.  Functions marked with "X" are required;
           functions marked with "O" are only called if they are non-``NULL`` and
           functions marked with "N/A" are not applicable in the ``SUNDomEigEstimator``
           implementation that is being used.
   :align: center

   +----------------------------------------------------+---------------------+---------------------+
   | Routine                                            |   POWER ITERATION   |  ARNOLDI ITERATION  |
   |                                                    |                     |                     |
   +====================================================+=====================+=====================+
   | :c:func:`SUNDomEigEst_SetATimes`                   |          X          |          X          |
   +----------------------------------------------------+---------------------+---------------------+
   | :c:func:`SUNDomEigEst_SetMaxIters`\ :sup:`1`       |          O          |         N/A         |
   +----------------------------------------------------+---------------------+---------------------+
   | :c:func:`SUNDomEigEst_SetNumPreProcess`            |          O          |          O          |
   +----------------------------------------------------+---------------------+---------------------+
   | :c:func:`SUNDomEigEst_SetTol`\ :sup:`2`            |          O          |         N/A         |
   +----------------------------------------------------+---------------------+---------------------+
   | :c:func:`SUNDomEigEst_Initialize`                  |          X          |          X          |
   +----------------------------------------------------+---------------------+---------------------+
   | :c:func:`SUNDomEigEst_PreProcess`                  |          O          |          O          |
   +----------------------------------------------------+---------------------+---------------------+
   | :c:func:`SUNDomEigEst_ComputeHess`\ :sup:`3`       |         N/A         |          X          |
   +----------------------------------------------------+---------------------+---------------------+
   | :c:func:`SUNDomEig_Estimate`                       |          X          |          X          |
   +----------------------------------------------------+---------------------+---------------------+
   | :c:func:`SUNDomEigEst_GetNumIters`\ :sup:`4`       |          O          |         N/A         |
   +----------------------------------------------------+---------------------+---------------------+
   | :c:func:`SUNDomEigEst_GetRes`\ :sup:`5`            |          O          |         N/A         |
   +----------------------------------------------------+---------------------+---------------------+
   | :c:func:`SUNDomEigEstFree`\ :sup:`6`               |                     |                     |
   +----------------------------------------------------+---------------------+---------------------+


Notes:

1. :c:func:`SUNDomEigEst_SetMaxIters()` might or might not be required depending on
   ``SUNDomEigEstimator`` implementation that is being used. This flag must be left
   ``NULL`` if it is not applicable for an estimator.

2. :c:func:`SUNDomEigEst_SetTol()` might or might not be required depending on
   ``SUNDomEigEstimator`` implementation that is being used. This flag must be left
   ``NULL`` if it is not applicable for an estimator.

3. :c:func:`SUNDomEigEst_ComputeHess()` might or might not be required depending on
   ``SUNDomEigEstimator`` implementation that is being used. This flag must be left
   ``NULL`` if it is not applicable for an estimator.

4. :c:func:`SUNDomEigEst_GetNumIters()` is only used to accumulate overall
   iterative estimator statistics.  If it is not implemented by
   the ``SUNDomEigEstimator`` module, then ARKDEE will consider all
   estimates as requiring zero iterations.

5. Although :c:func:`SUNDomEigEst_GetRes()` is optional, if it is not
   implemented by the ``SUNDomEigEstimator`` then ARKDEE will consider all
   estimates a being *exact*.

6. Although ARKDEE does not call :c:func:`SUNDomEigEstFree()`
   directly, this routine should be available for users to call when
   cleaning up from a simulation.