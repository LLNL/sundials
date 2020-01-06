..
   Programmer(s): Daniel R. Reynolds @ SMU
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

.. _ARKStep_CInterface.Headers:

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
- ``incdir/include/sunmatrix``
- ``incdir/include/sunlinsol``
- ``incdir/include/sunnonlinsol``

The directories ``libdir`` and ``incdir`` are the installation library
and include directories, respectively.  For a default installation,
these are ``instdir/lib`` and ``instdir/include``, respectively, where
``instdir`` is the directory where SUNDIALS was installed (see the
section :ref:`Installation` for further details).



.. _ARKStep_CInterface.DataTypes:

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

Additionally, SUNDIALS defines several macros for common mathematical
functions *e.g.*}, ``fabs``, ``sqrt``, ``exp``, etc. in ``sundials_math.h``. The
macros are prefixed with ``SUNR`` and expand to the appropriate C function based
on the ``realtype``. For example, the macro ``SUNRabs`` expands to the C
function ``fabs`` when ``realtype`` is ``double``, ``fabsf`` when ``realtype``
is ``float``, and ``fabsl`` when ``realtype`` is ``long double``.

A user program which uses the type ``realtype``, the ``RCONST`` macro, and the
``SUNR`` mathematical function macros is precision-independent except for any
calls to precision-specific library functions. Our example programs use
``realtype``, ``RCONST``, and the ``SUNR`` macros. Users can, however, use the
type ``double``, ``float``, or ``long double`` in their code (assuming that this
usage is consistent with the typedef for ``realtype``) and call the appropriate
math library functions directly. Thus, a previously existing piece of ANSI C
code can use SUNDIALS without modifying the code to use ``realtype``,
``RCONST``, or the ``SUNR`` macros so long as the SUNDIALS libraries use
the correct precision (for details see :ref:`Installation`).




Integer types used for indexing
--------------------------------

The type ``sunindextype`` is used for indexing array entries in SUNDIALS modules
(*e.g.*, vectors lengths and matrix sizes) as well as for storing the total
problem size. During configuration ``sunindextype`` may be selected to be either
a 32- or 64-bit *signed* integer with the default being 64-bit. See the section
:ref:`Installation` for the configuration option to select the desired size of
``sunindextype``. When using a 32-bit integer the total problem size is limited
to :math:`2^{31}-1` and with 64-bit integers the limit is :math:`2^{63}-1`. For
users with problem sizes that exceed the 64-bit limit an advanced configuration
option is available to specify the type used for ``sunindextype``.

A user program which uses ``sunindextype`` to handle indices will work with both
index storage types except for any calls to index storage-specific external
libraries. Our C and C++ example programs use ``sunindextype``. Users can,
however, use any compatible type (*e.g.*,  ``int``, ``long int``, ``int32_t``,
``int64_t`` or ``long long int``) in their code, assuming that this usage is
consistent with the typedef for ``sunindextype`` on their architecture. Thus, a
previously existing piece of ANSI C code can use SUNDIALS without modifying the
code to use ``sunindextype``, so long as the SUNDIALS libraries use the
appropriate index storage type (for details :ref:`Installation`).



Header Files
===========================================

When using ARKStep, the calling program must include several header
files so that various macros and data types can be used. The header
file that is always required is:

- ``arkode/arkode_arkstep.h``, the main header file for the ARKStep
  time-stepping module, which defines the several types and various
  constants, includes function prototypes, and includes the shared
  ``arkode/arkode.h`` and ``arkode/arkode_ls.h`` header files.

Note that ``arkode.h`` includes ``sundials_types.h`` directly, which
defines the types ``realtype``,  ``sunindextype`` and ``booleantype``
and the constants ``SUNFALSE`` and ``SUNTRUE``, so a user program does
not need to include ``sundials_types.h`` directly.

Additionally, the calling program must also include a NVECTOR
implementation header file, of the form ``nvector/nvector_***.h``,
corresponding to the user's preferred data layout and form of
parallelism.  See the section :ref:`NVectors` for details for the
appropriate name.  This file in turn includes the header file
``sundials_nvector.h`` which defines the abstract ``N_Vector`` data
type.

If the user includes a non-trivial implicit component to their
ODE system, then each time step will require a nonlinear solver for
the resulting systems of equations -- the default for this is a
modified Newton iteration.  If using a non-default nonlinear solver
module, or when interacting with a SUNNONLINSOL module directly, the
calling program must also include a SUNNONLINSOL header file, of the
form ``sunnonlinsol/sunnonlinsol_***.h`` where ``***`` is the name of
the nonlinear solver module (see the section :ref:`SUNNonlinSol` for
more information).  This file in turn includes the header file
``sundials_nonlinearsolver.h`` which defines the abstract
``SUNNonlinearSolver`` data type.

If using a nonlinear solver that requires the solution of a linear
system of the form :math:`\mathcal{A}x=b` (e.g., the default Newton
iteration), then a linear solver module header file will also
be required.  Similarly, if the ODE system involves a non-identity
mass matrix :math:`M \ne I`, then each time step will require a linear
solver for systems of the form :math:`Mx=b`.  The header files
corresponding to the SUNDIALS-provided linear solver modules available
for use with ARKode are:

- Direct linear solvers:

  - ``sunlinsol/sunlinsol_dense.h``,
    which is used with the dense linear solver module,
    SUNLINSOL_DENSE;

  - ``sunlinsol/sunlinsol_band.h``,
    which is used with the banded linear solver module,
    SUNLINSOL_BAND;

  - ``sunlinsol/sunlinsol_lapackdense.h``,
    which is used with the LAPACK dense linear solver module,
    SUNLINSOL_LAPACKDENSE;

  - ``sunlinsol/sunlinsol_lapackband.h``,
    which is used with the LAPACK banded linear solver module,
    SUNLINSOL_LAPACKBAND;

  - ``sunlinsol/sunlinsol_klu.h``,
    which is used with the KLU sparse linear solver module,
    SUNLINSOL_KLU;

  - ``sunlinsol/sunlinsol_superlumt.h``,
    which is used with the SuperLU_MT sparse linear solver module,
    SUNLINSOL_SUPERLUMT;

- Iterative linear solvers:

  - ``sunlinsol/sunlinsol_spgmr.h``,
    which is used with the scaled, preconditioned GMRES Krylov linear
    solver module, SUNLINSOL_SPGMR;

  - ``sunlinsol/sunlinsol_spfgmr.h``,
    which is used with the scaled, preconditioned FGMRES Krylov linear
    solver module, SUNLINSOL_SPFGMR;

  - ``sunlinsol/sunlinsol_spbcgs.h``,
    which is used with the scaled, preconditioned Bi-CGStab Krylov
    linear solver module, SUNLINSOL_SPBCGS;

  - ``sunlinsol/sunlinsol_sptfqmr.h``,
    which is used with the scaled, preconditioned TFQMR Krylov linear
    solver module, SUNLINSOL_SPTFQMR;

  - ``sunlinsol/sunlinsol_pcg.h``,
    which is used with the scaled, preconditioned CG Krylov linear
    solver module, SUNLINSOL_PCG;

The header files for the SUNLINSOL_DENSE and SUNLINSOL_LAPACKDENSE
linear solver modules include the file
``sunmatrix/sunmatrix_dense.h``, which defines the SUNMATRIX_DENSE
matrix module, as well as various functions and macros for acting on
such matrices.

The header files for the SUNLINSOL_BAND and SUNLINSOL_LAPACKBAND
linear solver modules include the file ``sunmatrix/sunmatrix_band.h``,
which defines the SUNMATRIX_BAND matrix module, as well as various
functions and macros for acting on such matrices.

The header files for the SUNLINSOL_KLU and SUNLINSOL_SUPERLUMT linear
solver modules include the file ``sunmatrix/sunmatrix_sparse.h``,
which defines the SUNMATRIX_SPARSE matrix module, as well as various
functions and macros for acting on such matrices.

The header files for the Krylov iterative solvers include the file
``sundials/sundials_iterative.h``, which enumerates the
preconditioning type and (for the SPGMR and SPFGMR solvers) the
choices for the Gram-Schmidt orthogonalization process.

Other headers may be needed, according to the choice of
preconditioner, etc.  For example, if preconditioning for an iterative
linear solver were performed using the ARKBBDPRE module, the header
``arkode/arkode_bbdpre.h`` is needed to access the preconditioner
initialization routines.
