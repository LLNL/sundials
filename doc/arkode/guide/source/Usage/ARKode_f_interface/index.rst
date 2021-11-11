..
   Programmer(s): Daniel R. Reynolds @ SMU and
                  Cody J. Balos @ LLNL
   ----------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2002-2021, Lawrence Livermore National Security
   and Southern Methodist University.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   ----------------------------------------------------------------

.. _Usage.Fortran:

=====================================
Using ARKODE for Fortran Applications
=====================================

Fortran 2003 interfaces to each of the time-stepping modules are provided to
support the use of ARKODE, for the solution of ODE systems, in a mixed Fortran/C
setting. While ARKODE is written in C, it is assumed here that the user's
calling program and user-supplied problem-definining rotuines are written in
Fortran.

.. _Usage.Fortran.Modules:

ARKODE Fortran 2003 Interface Modules
=====================================

The ARKODE Fortran 2003 modules define interfaces to most of the ARKODE C
API using the intrinsic ``iso_c_binding`` module which provides a standardized
mechanism for interoperating with C. ARKODE provides four Fortran 2003 modules:

* ``farkode_arkstep_mod``, ``farkode_erkstep_mod``, ``farkode_mristep_mod`` provide
  interfaces to the ARKStep, ERKStep, and MRIStep time-stepping modules respectively
* ``farkode_mod`` which interfaces to the components of ARKODE which are shared by the
  time-stepping modules

All interfaced functions are named after the corresponding C function, but with a
leading 'F'. For example. the ARKStep function ``ARKStepCreate`` is interfaced as
``FARKStepCreate``. Thus, the steps to use an ARKODE time-stepping module from Fortran
are identical (ignoring language differences) to using it from C/C++.

The Fortran 2003 ARKODE interface modules can be accessed by the ``use`` statement,
i.e. ``use farkode_mod``, and linking to the library ``libsundials_farkode_mod.lib``
in addition to ``libsundials_farkode.lib``. Further information on the location of
installed modules is provided in the Chapter :numref:`Installation`.

The Fortran 2003 interface modules were generated with SWIG Fortran, a
fork of SWIG :cite:p:`Swig-Fortran`. Users who are interested in the SWIG code used
in the generation process should contact the SUNDIALS development team.

.. include:: ../../../../../shared/SundialsF2003.rst
