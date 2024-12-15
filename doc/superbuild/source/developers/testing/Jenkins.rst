..
   -----------------------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2002-2024, Lawrence Livermore National Security
   and Southern Methodist University.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   -----------------------------------------------------------------------------

.. _Test.Jenkins:

Jenkins CI Testing
==================

With this option enabled tests are run using a Python test runner
``test/testRunner``.


Testing Scripts
---------------

The script ``test_driver.sh`` is provided to setup and run SUNDIALS tests. The
input option ``--testtype`` can be used to select from a predefined set of
test configurations:

* BRANCH -- C90 build tests and a small set of build and run tests.

* PR -- C90 build tests and a larger set of build, run, and install tests using
  the SUNDIALS release tarball.

* RELEASE -- C90 build test and an even larger set of build, run, and install
  tests using the SUNDIALS release tarball.

* CUSTOM -- run a single user-defined test configuration set using additional
  input options.

Run ``test_driver.sh -h`` for more information on the options available.

Note: At this time the testing scripts only run development tests when SUNDIALS
is configured with real type double (either index size can be used).

Testing environment
--------------------

The ``env/setup_env.sh`` script and a machine specific environment script are used
to setup the testing environment. The machine specific environment script used
(listed in the order checked) is:

#. The script specified by the environment variable ``SUNDIALS_ENV_FILE``
#. A user's local environment script: ``<sunrepo>/test/env/env.sh``
#. A machine environment script: ``<sunrepo>/test/env/${HOSTNAME}.sh``
#. A machine environment script: ``<sunrepo>/test/env/${HOST}.sh``
#. The default SUNDIALS environment script: ``<sunrepo>/test/env/default.sh``

Environment scripts must set the following environment variables that are used
when configuring SUNDIALS for testing.

.. code-block:: shell

   CC  = C compiler
   CXX = C++ compiler
   FC  = Fortran compiler

   CFLAGS   = C compiler flags
   CXXFLAGS = C++ compiler flags
   FFLAGS   = Fortran compiler flags

Note the test scripts will append the C standard flag (``-std=c90`` or ``-std=c99``)
and C++ standard flag (``-std=c++11`` or ``-std=c++17``) to the compiler flags
provided by the environment variables.

An environment script may optionally set additional environment variables to
enable or disable third party libraries (TPLs) in the SUNDIALS configuration.
The currently supported TPL environment variables are as follows:

.. code-block:: shell

   # Enable/disable MPI support
   SUNDIALS_MPI = ON or OFF
   MPICC        = MPI C compiler wrapper
   MPICXX       = MPI C++ compiler wrapper
   MPIFC        = MPI Fortran compiler wrapper
   MPIEXEC      = executable for launching MPI runs

   # Enable/disable PThread support
   SUNDIALS_PTHREAD = ON or OFF

   # Enable/disable OpenMP support
   SUNDIALS_OPENMP = ON or OFF

   # Enable/disable OpenMP device offloading support
   SUNDIALS_OPENMPDEV = ON or OFF

   # Enable/disable CUDA support
   SUNDIALS_CUDA = ON or OFF

   # Enable/disable HIP support
   SUNDIALS_HIP = ON or OFF

   # Enable/disable RAJA support
   SUNDIALS_RAJA = ON or OFF
   RAJA_ROOT     = full path to RAJA installation
   RAJA_BACKENDS = RAJA backends

   # Enable/disable SYCL support
   SUNDIALS_SYCL = ON or OFF

   # Enable/disable LAPACK linear solvers
   SUNDIALS_LAPACK  = ON or OFF
   LAPACK_LIBRARIES = full path to LAPACK library

   # Enable/disable KLU linear solver
   SUNDIALS_KLU             = ON or OFF
   SUITE_SPARSE_INCLUDE_DIR = full path to SuiteSparse include directory
   SUITE_SPARSE_LIBRARY_DIR = full path to SuiteSparse library directory

   # Enable/disable SuperLU_MT linear solver
   SUNDIALS_SUPERLU_MT    = ON or OFF
   SUPERLU_MT_INCLUDE_DIR = full path to SuperLU_MT include directory
   SUPERLU_MT_LIBRARY_DIR = full path to SuperLU_MT library directory

   # Enable/disable SuperLU_DIST linear solver
   SUNDIALS_SUPERLU_DIST  = ON or OFF
   SUPERLU_DIST_INCLUDE_DIR = full path to SuperLU_DIST include directory
   SUPERLU_DIST_LIBRARY_DIR = full path to SuperLU_DIST library directory
   SUPERLU_DIST_LIBRARIES   = additional link libraries for SuperLU_DIST

   # Enable/disable MAGMA linear solver
   SUNDIALS_MAGMA = ON or OFF
   MAGMA_ROOT     = full path to MAGMA installation
   MAGMA_BAKCENDS = MAGMA backend

   # Enable/disable hypre support
   SUNDIALS_HYPRE    = ON or OFF
   HYPRE_INCLUDE_DIR = full path to hypre include directory
   HYPRE_LIBRARY_DIR = full path to hypre library directory

   # Enable/disable PETSc support
   SUNDIALS_PETSC = ON or OFF
   PETSC_ROOT     = full path to PETSc installation

   # Enable/disable Trilinos support
   SUNDIALS_TRILINOS = ON or OFF
   TRILINOS_ROOT     = full path to Trilinos installation

   # Enable/disable Trilinos support
   SUNDIALS_XBRAID = ON or OFF
   XBRAID_ROOT     = full path to XBraid installation

To aid in setting the above variables appropriately, ``env/setup_env.sh`` will set
environment variables ``SUNDIALS_PRECISION``, ``SUNDIALS_INDEX_SIZE``, and
``SUNDIALS_LIBRARY_TYPE`` for the real type (``single``, ``double``, or ``extended``),
index size (``32`` or ``64``), and library type (``shared``, ``static`` or ``both``) used
in the test. Additionally, ``SUNDIALS_TPLS`` is set to ``ON`` or ``OFF`` to indicate
if the test will use third-party libraries.

Any additional input options passed to ``test_driver.sh`` will be passed through
to the environment script for additional setup information. For example, the
machine script may accept a compiler spec (e.g., ``gcc@4.9.4``) and/or a build
type (e.g., ``opt`` for an optimized build).

## Using Spack to install TPLs

The TPLs needed for a complete build of SUNDIALS can be easily installed with
spack and the spack environment included in the SUNDIALS repository. Simply
navigate to ``test/spack`` and run ``spack install``. For more information on Spack
environments is the [Spack tutorial](https://spack.readthedocs.io/en/latest/tutorial_environments.html).
