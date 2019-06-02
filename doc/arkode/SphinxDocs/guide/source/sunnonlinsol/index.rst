..
   Programmer(s): Daniel R. Reynolds @ SMU
   ----------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2002-2019, Lawrence Livermore National Security
   and Southern Methodist University.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   ----------------------------------------------------------------

:tocdepth: 3

.. _SUNNonlinSol:

=================================
Nonlinear Solver Data Structures
=================================


.. _SUNNonlinSol.Description:

Description of the SUNNonlinearSolver Module
===============================================

SUNDIALS time integration packages are written in terms of generic nonlinear
solver operations defined by the SUNNonlinSol API and implemented by a
particular SUNNonlinSol module of type ``SUNNonlinearSolver``.
Users can supply their own SUNNonlinSol module, or use one of the modules
provided with SUNDIALS.

The time integrators in SUNDIALS specify a default nonlinear solver module
and as such this chapter is intended for users that wish to use a non-default
nonlinear solver module or would like to provide their own nonlinear solver
implementation. Users interested in using a non-default solver module may skip
the description of the SUNNonlinSol API in section :ref:`SUNNonlinSol.API`
and proceeded to the subsequent sections in this chapter that describe the
SUNNonlinSol modules provided with SUNDIALS.


For users interested in providing their own SUNNonlinSol module, the
following section presents the SUNNonlinSol API and its implementation
beginning with the definition of SUNNonlinSol functions in the
sections :ref:`SUNNonlinSol.CoreFn`, :ref:`SUNNonlinSol.SetFn` and
:ref:`SUNNonlinSol.GetFn`. This is followed by the definition of
functions supplied to a nonlinear solver implementation in the section
:ref:`SUNNonlinSol.SUNSuppliedFn`.  The nonlinear solver return 
codes are given in the section :ref:`SUNNonlinSol.ReturnCodes`. The
``SUNNonlinearSolver`` type and the generic SUNNonlinSol module are defined
in the section :ref:`SUNNonlinSol.Generic`. Finally, the section
:ref:`SUNNonlinSol.Custom` lists the requirements for supplying a custom
SUNNonlinSol module. Users wishing to supply their own SUNNonlinSol module
are encouraged to use the SUNNonlinSol implementations provided with
SUNDIALS as a template for supplying custom nonlinear solver modules.




.. toctree::
   :maxdepth: 1

   SUNNonlinSol_API
   SUNNonlinSol_Newton
   SUNNonlinSol_FixedPoint
