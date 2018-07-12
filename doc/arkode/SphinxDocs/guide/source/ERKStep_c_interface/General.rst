..
   Programmer(s): Daniel R. Reynolds @ SMU
   ----------------------------------------------------------------
   Copyright (c) 2013, Southern Methodist University.
   All rights reserved.
   For details, see the LICENSE file.
   ----------------------------------------------------------------

:tocdepth: 3

.. _ERKStep_CInterface.Headers:

Access to library and header files
===========================================

As stated previously for the ARKStep module, we assume that the
ERKStep module has been installed with ARKode, following the procedure
described in the section :ref:`Installation`.  Moreover, we will
reference the same header file locations and data types as described
in the Sections :ref:`ARKStep_CInterface.Headers` and
:ref:`ARKStep_CInterface.DataTypes`.

When using ERKStep, the calling program must include several header
files so that various macros and data types can be used. The header
file that is always required is:

- ``arkode/arkode_erkstep.h``, the main header file for the ERKStep
  time-stepping module, which defines the several types and various
  constants, includes function prototypes, and includes the shared
  ``arkode/arkode.h`` header file.

Note that ``arkode.h`` includes ``sundials_types.h`` directly, which
defines the types ``realtype``,  ``sunindextype`` and ``booleantype``
and the constants ``SUNFALSE`` and ``SUNTRUE``, so a user program does
not need to include ``sundials_types.h`` directly.

Additionally, the calling program must also include an NVECTOR
implementation header file, of the form ``nvector/nvector_***.h``,
corresponding to the user's preferred data layout and form of
parallelism.  See the section :ref:`NVectors` for details for the
appropriate name.  This file in turn includes the header file
``sundials_nvector.h`` which defines the abstract ``N_Vector`` data
type.
