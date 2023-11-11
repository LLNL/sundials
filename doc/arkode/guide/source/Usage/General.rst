.. ----------------------------------------------------------------
   Programmer(s): David J. Gardner @ LLNL
                  Daniel R. Reynolds @ SMU
   ----------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2002-2023, Lawrence Livermore National Security
   and Southern Methodist University.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   ----------------------------------------------------------------

.. _ARKODE.Usage.Headers:

Access to library and header files
==================================

At this point, it is assumed that the installation of ARKODE, following the
procedure described in :numref:`Installation`, has been completed successfully.
In the proceeding text, the directories ``libdir`` and ``incdir`` are the
installation library and include directories, respectively. For a default
installation, these are ``instdir/lib`` and ``instdir/include``, respectively,
where ``instdir`` is the directory where SUNDIALS was installed.

Regardless of where the user's application program resides, its
associated compilation and load commands must make reference to the
appropriate locations for the library and header files required by
ARKODE. ARKODE and the SUNDIALS core symbols are found in
``libdir/libsundials_arkode.lib``. Symbols for additional SUNDIALS
modules, such as vectors and algebraic solvers, are found in

- ``libdir/libsundials_nvec*.lib``,
- ``libdir/libsundials_sunlinsol*.lib``,
- ``libdir/libsundials_sunnonlinsol*.lib``
- ``libdir/libsundials_sunmem*.lib``

The file extension ``.lib`` is typically ``.so`` for shared libraries 
and ``.a`` for static libraries.  

The relevant header files for ARKODE are located in the subdirectories
``incdir/include/arkode``. To use ARKODE the application needs to include 
the header file(s) for the ARKODE time-stepper(s) of choice. 

.. code:: c

  #include <arkode/arkode_erkstep.h>  // ERKStep provides explicit RK methods.
  #include <arkode/arkode_arkstep.h>  // ARKStep provides explicit, implicit, IMEX additive RK methods.
  #include <arkode/arkode_mristep.h>  // MRIStep provides mutlirate RK methods.
  #include <arkode/arkode_sprkstep.h> // SPRKStep provides symplectic partition RK methods.

Each of these define several types and various constants, include function
prototypes, and include the shared ``arkode/arkode.h`` and
``arkode/arkode_ls.h`` header files. No other header files are required to be
included, but there are optional 

If the user includes a non-trivial implicit component to their ODE system in
ARKStep, or if the slow time scale for MRIStep should be treated implicitly,
then each implicit stage will require a nonlinear solver for the resulting
system of algebraic equations -- the default for this is a modified or inexact
Newton iteration, depending on the user's choice of linear solver.  If choosing
to set which nonlinear solver module, or when interacting with a
:numref:`SUNNonlinSol` module directly, the calling program must also include a
:numref:`SUNNonlinSol` header file, of the form ``sunnonlinsol/sunnonlinsol_***.h`` 
where ``***`` is the name of the nonlinear solver module 
(see :numref:`SUNNonlinSol` for more information). 

If using a nonlinear solver that requires the solution of a linear system of the
form :math:`\mathcal{A}x=b` (e.g., the default Newton iteration), then a linear
solver module header file will also be required.  Similarly, if the ODE system
in ARKStep involves a non-identity mass matrix :math:`M \ne I`, then each time
step will require a linear solver for systems of the form :math:`Mx=b`. 
In this case it will be necessary to include the header file for a SUNDIALS linear 
solver, which is of the form ``sunlinsol/sunlinsol_***.h``. 

Other headers may be needed, according to the choice of preconditioner, etc.
For example, if preconditioning for an iterative linear solver were performed
using the ARKBBDPRE module, the header ``arkode/arkode_bbdpre.h`` is needed to
access the preconditioner initialization routines.


