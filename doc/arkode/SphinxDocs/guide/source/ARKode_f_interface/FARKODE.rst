..
   Programmer(s): Daniel R. Reynolds @ SMU and
                  Cody J. Balos @ LLNL
   ----------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2002-2020, Lawrence Livermore National Security
   and Southern Methodist University.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   ----------------------------------------------------------------

.. _FortranInterface:

=====================================================
FARKODE, an Interface Module for FORTRAN Applications
=====================================================

The FARKODE interface module is a package of C functions which
support the use of the ARKStep time-stepping module for the solution
of ODE systems

.. math::
   M\, \dot{y} = f^E(t,y) + f^I(t,y),

in a mixed Fortran/C setting.  While ARKode is written in C, it is
assumed here that the user's calling program and user-supplied
problem-defining routines are written in Fortran.  We assume only
minimal Fortran capabilities; specifically that the Fortran compiler
support full Fortran77 functionality (although more modern standards
are similarly supported).  This package provides the necessary
interfaces to ARKODE for the majority of supplied serial and parallel
NVECTOR implementations.



.. _FInterface.Portability:

Important note on portability
--------------------------------------

In this package, the names of the interface functions, and the names
of the Fortran user routines called by them, appear as dummy names
which are mapped to actual values by a series of definitions in the
header files.  By default, those mapping definitions depend in turn
on the C macro ``F77_FUNC`` defined in the header file
``sundials_config.h``.  The mapping defined by ``F77_FUNC`` in turn
transforms the C interface names to match the name-mangling approach
used by the supplied Fortran compiler.

By "name-mangling", we mean that due to the case-independent nature of
the Fortran language, Fortran compilers convert all subroutine and
object names to use either all lower-case or all upper-case
characters, and append either zero, one or two underscores as a prefix
or suffix the the name.  For example, the Fortran subroutine
``MyFunction()`` will be changed to one of ``myfunction``,
``MYFUNCTION``, ``myfunction__``, ``MYFUNCTION_``, and so on,
depending on the Fortran compiler used.

SUNDIALS determines this name-mangling scheme at configuration time
(see :ref:`Installation`).



.. _FInterface.DataTypes:

Fortran Data Types
-------------------------

Throughout this documentation, we will refer to data types according
to their usage in C.  The equivalent types to these may vary,
depending on your computer architecture and on how SUNDIALS was
compiled (see :ref:`Installation`).  A Fortran user should first
determine the equivalent types for their architecture and compiler,
and then take care that all arguments passed through this Fortran/C
interface are declared of the appropriate type.


**Integers**: SUNDIALS uses ``int``, ``long int`` and ``sunindextype``
types.  As discussed in :ref:`Installation`, at compilation SUNDIALS
allows the configuration of the 'index' type, that accepts values of
32-bit signed and 64-bit signed. This choice dictates the size of a
SUNDIALS ``sunindextype`` variable.

* ``int`` -- equivalent to an ``INTEGER`` or ``INTEGER*4`` in Fortran

* ``long int`` -- this will depend on the computer architecture:

  * 32-bit architecture -- equivalent to an ``INTEGER`` or ``INTEGER*4`` in Fortran

  * 64-bit architecture -- equivalent to an ``INTEGER*8`` in Fortran

* ``sunindextype`` -- this will depend on the SUNDIALS configuration:

  * 32-bit -- equivalent to an ``INTEGER`` or ``INTEGER*4`` in Fortran

  * 64-bit -- equivalent to an ``INTEGER*8`` in Fortran


**Real numbers**:  As discussed in :ref:`Installation`, at compilation
SUNDIALS allows the configuration option  ``--with-precision``,
that accepts values of ``single``, ``double`` or ``extended`` (the
default is ``double``).  This choice dictates the size of a
``realtype`` variable.  The corresponding Fortran types for these
``realtype`` sizes are:

* ``single`` -- equivalent to a ``REAL`` or ``REAL*4`` in Fortran

* ``double`` -- equivalent to a ``DOUBLE PRECISION`` or ``REAL*8`` in Fortran

* ``extended`` -- equivalent to a ``REAL*16`` in Fortran


We note that when SUNDIALS is compiled with Fortran interfaces
enabled, a file ``sundials/sundials_fconfig.h`` is placed in the
installation's ``include`` directory, containing information about
the Fortran types that correspond to the C types of the configured
SUNDIALS installation.  This file may be "included" by Fortran
routines, as long as the compiler supports the Fortran90 standard (or
higher), as shown in the ARKode example programs ``ark_bruss.f90``,
``ark_bruss1D_FEM_klu.f90`` and ``fark_heat2D.f90``.

Details on the Fortran interface to ARKode are provided in the
following sub-sections:

.. toctree::
   :maxdepth: 1

   Routines
   Usage
   Optional_output
   Rootfinding
   Preconditioning
