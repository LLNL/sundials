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

====================================
ARKODE Documentation
====================================

This is the documentation for ARKODE, an adaptive step time
integration package for stiff, nonstiff and mixed stiff/nonstiff
systems of ordinary differential equations (ODEs) using Runge-Kutta
(i.e. one-step, multi-stage) methods.  The ARKODE solver is a
component of the `SUNDIALS
<https://computing.llnl.gov/casc/sundials/main.html>`_ suite of
nonlinear and differential/algebraic equation solvers. It is designed
to have a similar user experience to the `CVODE
<https://computing.llnl.gov/casc/sundials/description/description.html#descr_cvode>`_
solver, including user modes to allow adaptive integration to specified
output times, return after each internal step and root-finding
capabilities, and for calculations in serial, using shared-memory
parallelism (via OpenMP, Pthreads, CUDA, Raja) or distributed-memory
parallelism (via MPI).  The default integration and solver options
should apply to most users, though control over nearly all internal
parameters and time adaptivity algorithms is enabled through optional
interface routines.

ARKODE is written in C, with C++ and Fortran interfaces.

ARKODE is developed by `Southern Methodist University <http://www.smu.edu>`_,
with support by the `US Department of Energy <http://www.doe.gov>`_,
`Office of Science <https://www.energy.gov/science/office-science>`_,
`Office of Advanced Scientific Computing Research <https://www.energy.gov/science/ascr/advanced-scientific-computing-research>`_,
`Scientific Discovery through Advanced Computing (SciDAC) <https://www.scidac.gov/>`_
Program through the `FASTMath <http://www.fastmath-scidac.org/>`_ Institute, under
subcontracts B598130 and B626484 from `Lawrence Livermore National Laboratory
<http://www.llnl.gov>`_, and under DOE award DE-SC0021354.


.. only:: html

   **Authors:** *Daniel R. Reynolds, David J. Gardner, Carol S. Woodward, and Cody J. Balos*

   **Contributors:** The SUNDIALS library has been developed over many years by
   a number of contributors. The current SUNDIALS team consists of Cody J.
   Balos, David J. Gardner, Alan C. Hindmarsh, Daniel R. Reynolds, and Carol S.
   Woodward. We thank Radu Serban for significant and critical past
   contributions.

   Other contributors to SUNDIALS include: James Almgren-Bell, Lawrence E.
   Banks, Peter N. Brown, George Byrne, Rujeko Chinomona, Scott D. Cohen, Aaron
   Collier, Keith E. Grant, Steven L. Lee, Shelby L. Lockhart, John Loffeld,
   Daniel McGreer, Slaven Peles, Cosmin Petra, H. Hunter Schwartz, Jean M.
   Sexton, Dan Shumaker, Steve G. Smith, Allan G. Taylor, Hilari C. Tiedeman,
   Chris White, Ting Yan, and Ulrike M. Yang.

   Table of Contents
   =================

.. toctree::
   :numbered:
   :maxdepth: 1

   Introduction
   Mathematics
   Organization
   sundials/index.rst
   Usage/index.rst
   ARKodeButcherTable
   nvectors/index.rst
   sunmatrix/index.rst
   sunlinsol/index.rst
   sunnonlinsol/index.rst
   sunmemory/index.rst
   Install
   Constants
   Butcher
   History
   References
.. only:: html

   * :ref:`genindex`
