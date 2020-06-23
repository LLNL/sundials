..
   Programmer(s): David J. Gardner @ LLNL
   ----------------------------------------------------------------
   Based on ERKStep by Daniel R. Reynolds @ SMU
   ----------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2002-2020, Lawrence Livermore National Security
   and Southern Methodist University.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   ----------------------------------------------------------------

:tocdepth: 3

.. _MRIStep_CInterface.Headers:

Access to library and header files
===========================================

At this point, it is assumed that the installation of ARKode,
following the procedure described in the section :ref:`Installation`,
has been completed successfully.

Regardless of where the user's application program resides, its
associated compilation and load commands must make reference to the
appropriate locations for the library and header files required by
ARKode. The relevant library files are

- ``libdir/libsundials_arkode.lib``,
- ``libdir/libsundials_nvec*.lib``,

where the file extension ``.lib`` is typically ``.so`` for shared
libraries and ``.a`` for static libraries.  The relevant header files
are located in the subdirectories

- ``incdir/include/arkode``
- ``incdir/include/sundials``
- ``incdir/include/nvector``

The directories ``libdir`` and ``incdir`` are the installation library
and include directories, respectively.  For a default installation,
these are ``instdir/lib`` and ``instdir/include``, respectively, where
``instdir`` is the directory where SUNDIALS was installed (see the
section :ref:`Installation` for further details).



.. _MRIStep_CInterface.DataTypes:

Data Types
===========================================

The ``sundials_types.h`` file contains the definition of the variable
type ``realtype``, which is used by the SUNDIALS solvers for all
floating-point data, the definition of the integer type
``sunindextype``, which is used for vector and matrix indices, and
``booleantype``, which is used for certain logic operations within
SUNDIALS.


Floating point types
-----------------------

The type ":index:`realtype`" can be set to
``float``, ``double``, or ``long double``, depending on how SUNDIALS
was installed (with the default being ``double``). The user can change
the precision of the SUNDIALS solvers' floating-point arithmetic at the
configuration stage (see the section :ref:`Installation.CMake.Options`).

Additionally, based on the current precision, ``sundials_types.h``
defines the values :index:`BIG_REAL` to be the largest value
representable as a ``realtype``, :index:`SMALL_REAL` to be the
smallest positive value representable as a ``realtype``, and
:index:`UNIT_ROUNDOFF` to be the smallest realtype number,
:math:`\varepsilon`, such that :math:`1.0 + \varepsilon \ne 1.0`.

Within SUNDIALS, real constants may be set to have the appropriate
precision by way of a macro called :index:`RCONST`.  It is this macro
that needs the ability to branch on the definition ``realtype``.  In
ANSI C, a floating-point constant with no suffix is stored as a
``double``. Placing the suffix "F" at the end of a floating point
constant makes it a ``float``, whereas using the suffix "L" makes it a
``long double``. For example,

.. code-block:: c

   #define A 1.0
   #define B 1.0F
   #define C 1.0L

defines ``A`` to be a ``double`` constant equal to 1.0, ``B`` to be a
``float`` constant equal to 1.0, and ``C`` to be a ``long double`` constant
equal to 1.0.  The macro call ``RCONST(1.0)`` automatically expands to
1.0 if ``realtype`` is ``double``, to ``1.0F`` if ``realtype`` is ``float``, or
to ``1.0L`` if ``realtype`` is ``long double``. SUNDIALS uses the ``RCONST``
macro internally to declare all of its floating-point constants.

A user program which uses the type ``realtype`` and the ``RCONST`` macro
to handle floating-point constants is precision-independent, except for
any calls to precision-specific standard math library functions.
Users can, however, use the types ``double``, ``float``, or ``long
double`` in their code (assuming that this usage is consistent with
the size of ``realtype`` values that are passed to and from SUNDIALS).
Thus, a previously existing piece of ANSI C code can use SUNDIALS
without modifying the code to use ``realtype``, so long as the
SUNDIALS libraries have been compiled using the same precision (for
details see the section :ref:`Installation`).



Integer types used for vector and matrix indices
---------------------------------------------------

The type ``sunindextype`` can be either a 32- or 64-bit *signed* integer.
The default is the portable ``int64_t`` type, and the user can change it
to ``int32_t`` at the configuration stage. The configuration system
will detect if the compiler does not support portable types, and will
replace ``int32_t`` and ``int64_t`` with ``int`` and ``long int``,
respectively, to ensure use of the desired sizes on Linux, Mac OS X, and Windows
platforms. SUNDIALS currently does not support *unsigned* integer types
for vector and matrix indices, although these could be added in the future if there
is sufficient demand.

A user program which uses ``sunindextype`` to handle vector and matrix indices
will work with both index storage types except for any calls to index storage-specific
external libraries. (Our ``C`` and ``C++`` example programs use ``sunindextype``.)
Users can, however, use any one of ``int``, ``long int``, ``int32_t``, ``int64_t`` or
``long long int`` in their code, assuming that this usage is consistent with the typedef
for ``sunindextype`` on their architecture. Thus, a previously existing piece of ANSI
C code can use SUNDIALS without modifying the code to use ``sunindextype``,
so long as the SUNDIALS libraries use the appropriate index storage type (for details
see the section :ref:`Installation`).


Header Files
===========================================

When using MRIStep, the calling program must include several header
files so that various macros and data types can be used. The header
file that is always required is:

- ``arkode/arkode_mristep.h``, the main header file for the MRIStep
  time-stepping module, which defines the several types and various
  constants, includes function prototypes, and includes the shared
  ``arkode/arkode.h`` header file.

Note that ``arkode.h`` includes ``sundials_types.h`` directly, which
defines the types ``realtype``,  ``sunindextype``, and ``booleantype``
and the constants ``SUNFALSE`` and ``SUNTRUE``, so a user program does
not need to include ``sundials_types.h`` directly.

Additionally, the calling program must also include an NVECTOR
implementation header file, of the form ``nvector/nvector_***.h``,
corresponding to the user's preferred data layout and form of
parallelism.  See the section :ref:`NVectors` for details for the
appropriate name.  This file in turn includes the header file
``sundials_nvector.h`` which defines the abstract ``N_Vector`` data
type.

