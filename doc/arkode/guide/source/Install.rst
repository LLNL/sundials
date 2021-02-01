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

:tocdepth: 3

.. _Installation:

=============================
ARKode Installation Procedure
=============================

The installation of any SUNDIALS package is accomplished by installing
the SUNDIALS suite as a whole, according to the instructions that
follow.  The same procedure applies whether or not the downloaded
file contains one or all solvers in SUNDIALS.

The SUNDIALS suite (or individual solvers) are distributed as
compressed archives (``.tar.gz``).  The name of the distribution
archive is of the form ``SOLVER-X.Y.Z.tar.gz``, where ``SOLVER`` is
one of: ``sundials``, ``cvode``, ``cvodes``, ``arkode``, ``ida``,
``idas``, or ``kinsol``, and ``X.Y.Z`` represents the version number
(of the SUNDIALS suite or of the individual solver).
To begin the installation, first uncompress and expand the sources, by
issuing

.. code-block:: bash

   % tar -zxf SOLVER-X.Y.Z.tar.gz

This will extract source files under a directory ``SOLVER-X.Y.Z``.

Starting with version 2.6.0 of SUNDIALS, CMake is the only supported
method of installation.  The explanations of the installation
procedure begins with a few common observations:

* The remainder of this chapter will follow these conventions:

  ``SOLVERDIR``
     is the directory ``SOLVER-X.Y.Z`` created above; i.e. the
     directory containing the SUNDIALS sources.

  ``BUILDDIR``
     is the (temporary) directory under which SUNDIALS is built.

  ``INSTDIR``
     is the directory under which the SUNDIALS exported header files
     and libraries will be installed. Typically, header files are
     exported under a directory ``INSTDIR/include`` while libraries
     are installed under ``INSTDIR/lib``, with ``INSTDIR``
     specified at configuration time.

* For SUNDIALS' CMake-based installation, in-source builds are prohibited;
  in other words, the build directory ``BUILDDIR`` can **not** be the
  same as ``SOLVERDIR`` and such an attempt will lead to an error.  This
  prevents "polluting" the source tree and allows efficient builds for
  different configurations and/or options.

* The installation directory ``INSTDIR`` can not be the same as
  the source directory ``SOLVERDIR``.

* By default, only the libraries and header files are exported to the
  installation directory ``INSTDIR``.  If enabled by the user (with the
  appropriate toggle for CMake), the
  examples distributed with SUNDIALS will be built together with
  the solver libraries but the installation step will result in
  exporting (by default in a subdirectory of the installation
  directory) the example sources and sample outputs together with
  automatically generated configuration files that reference the
  *installed* SUNDIALS headers and libraries.  As such, these
  configuration files for the SUNDIALS examples can be used as
  "templates" for your own problems. CMake installs
  ``CMakeLists.txt`` files and also (as an option available only under
  Unix/Linux) ``Makefile`` files. Note this installation approach also
  allows the option of building the SUNDIALS examples without having
  to install them.  (This can be used as a sanity check for the
  freshly built libraries.)

* Even if generation of shared libraries is enabled, only static
  libraries are created for the FCMIX modules.  Because of the use of
  fixed names for the Fortran user-provided subroutines, FCMIX shared
  libraries would result in "undefined symbol" errors at link time.


Further details on the CMake-based installation procedures,
instructions for manual compilation, and a roadmap of the resulting
installed libraries and exported header files, are provided in the
following subsections:

* :ref:`Installation.CMake`
* :ref:`Installation.Results`





.. _Installation.CMake:

CMake-based installation
======================================

CMake-based installation provides a platform-independent build
system. CMake can generate Unix and Linux Makefiles, as well as
KDevelop, Visual Studio, and (Apple) XCode project files from the same
configuration file.  In addition, CMake also provides a GUI front end
and which allows an interactive build and installation process.

The SUNDIALS build process requires CMake version 3.0.2 or
higher and a working C compiler.  On Unix-like operating systems, it
also requires Make (and ``curses``, including its development libraries,
for the GUI front end to CMake, ``ccmake`` or ``cmake-gui``), while on
Windows it requires Visual Studio.  While many Linux distributions
offer CMake, the version included may be out of date.  Many new CMake
features have been added recently, and you should download the latest
version from http://www.cmake.org.  Build instructions for CMake
(only necessary for Unix-like systems) can be found on the CMake website.
Once CMake is installed, Linux/Unix users will be able to use
``ccmake`` or ``cmake-gui`` (depending on the version of CMake),
while Windows users will be able to use ``CMakeSetup``.

As previously noted, when using CMake to configure, build and install
SUNDIALS, it is always required to use a separate build
directory. While in-source builds are possible, they are explicitly
prohibited by the SUNDIALS CMake scripts (one of the reasons being
that, unlike autotools, CMake does not provide a ``make distclean``
procedure and it is therefore difficult to clean-up the source tree
after an in-source build). By ensuring a separate build directory, it
is an easy task for the user to clean-up all traces of the build by
simply removing the build directory. CMake does generate a ``make
clean`` which will remove files generated by the compiler and linker.




.. index:: ccmake

.. _Installation.CMake.Unix:

Configuring, building, and installing on Unix-like systems
----------------------------------------------------------------

The default CMake configuration will build all included solvers and
associated examples and will build static and shared libraries. The
INSTDIR defaults to ``/usr/local`` and can be changed by setting
the ``CMAKE_INSTALL_PREFIX`` variable. Support for FORTRAN and all
other options are disabled.

CMake can be used from the command line with the ``cmake`` command, or
from a ``curses``\ -based GUI by using the ``ccmake`` command, or from
a wxWidgets or QT based GUI by using the ``cmake-gui``
command. Examples for using both text and graphical methods will be
presented.  For the examples shown it is assumed that there is a top
level SUNDIALS directory with appropriate source, build and install
directories:


.. code-block:: bash

   $ mkdir (...)/INSTDIR
   $ mkdir (...)/BUILDDIR
   $ cd (...)/BUILDDIR


.. index:: cmake-gui
.. index:: ccmake


Building with the GUI
^^^^^^^^^^^^^^^^^^^^^^^

Using CMake with the ``ccmake`` GUI follows the general process:

* Select and modify values, run configure (``c`` key)

* New values are denoted with an asterisk

* To set a variable, move the cursor to the variable and press enter

  * If it is a boolean (ON/OFF) it will toggle the value

  * If it is string or file, it will allow editing of the string

  * For file and directories, the ``<tab>`` key can be used to complete

* Repeat until all values are set as desired and the generate option
  is available (``g`` key)

* Some variables (advanced variables) are not visible right away

* To see advanced variables, toggle to advanced mode (``t`` key)

* To search for a variable press ``/`` key, and to repeat the search,
  press the ``n`` key


Using CMake with the ``cmake-gui`` GUI follows a similar process:

* Select and modify values, click ``Configure``

* The first time you click ``Configure``, make sure to pick the
  appropriate generator (the following will ssume generation of Unix
  Makfiles).

* New values are highlighted in red

* To set a variable, click on or move the cursor to the variable and press enter

  * If it is a boolean (``ON/OFF``) it will check/uncheck the box

  * If it is string or file, it will allow editing of the string.
    Additionally, an ellipsis button will appear ``...`` on the far
    right of the entry.  Clicking this button will bring up the file
    or directory selection dialog.

  * For files and directories, the ``<tab>`` key can be used to
    complete

* Repeat until all values are set as desired and click the
  ``Generate`` button

* Some variables (advanced variables) are not visible right away

* To see advanced variables, click the ``advanced`` button



To build the default configuration using the curses GUI, from the
BUILDDIR enter the ``ccmake`` command and point to the ``SOLVERDIR``:

.. code-block:: bash

   $ ccmake (...)/SOLVERDIR

Similarly, to build the default configuration using the wxWidgets GUI,
from the BUILDDIR enter the ``cmake-gui`` command and point to the
``SOLVERDIR``:

.. code-block:: bash

   $ cmake-gui (...)/SOLVERDIR

The default curses configuration screen is shown in
the following figure.

.. _ccmakedefault:

.. figure:: figs/ccmakedefault.png
   :scale: 75 %
   :align: center

   Default configuration screen. Note: Initial screen is empty.
   To get this default configuration, press 'c' repeatedly (accepting
   default values denoted with asterisk) until the 'g' option is
   available.

The default INSTDIR for both SUNDIALS and corresponding examples
can be changed by setting the ``CMAKE_INSTALL_PREFIX`` and
the ``EXAMPLES_INSTALL_PATH`` as shown in the following figure.

.. _ccmakeprefix:

.. figure:: figs/ccmakeprefix.png
   :scale: 75 %
   :align: center

   Changing the INSTDIR for SUNDIALS and corresponding EXAMPLES.


Pressing the ``g`` key or clicking ``generate`` will generate
makefiles including all dependencies and all rules to build SUNDIALS
on this system.  Back at the command prompt, you can now run:

.. code-block:: bash

   $ make

or for a faster parallel build (e.g. using 4 threads), you can run

.. code-block:: bash

   $ make -j 4

To install SUNDIALS in the installation directory specified in the configuration, simply run:

.. code-block:: bash

   $ make install





.. index:: cmake

Building from the command line
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Using CMake from the command line is simply a matter of specifying
CMake variable settings with the ``cmake`` command.  The following
will build the default configuration:

.. code-block:: bash

   $ cmake -DCMAKE_INSTALL_PREFIX=/home/myname/sundials/instdir \
   >  -DEXAMPLES_INSTALL_PATH=/home/myname/sundials/instdir/examples \
   >  ../srcdir
   $ make
   $ make install




.. _Installation.CMake.Options:


Configuration options (Unix/Linux)
-----------------------------------

A complete list of all available options for a CMake-based SUNDIALS
configuration is provide below.  Note that the default values shown
are for a typical configuration on a Linux system and are provided as
illustration only.

:index:`BUILD_ARKODE <BUILD_ARKODE (CMake option)>`
   Build the ARKODE library

   Default: ``ON``

:index:`BUILD_CVODE <BUILD_CVODE (CMake option)>`
   Build the CVODE library

   Default: ``ON``

:index:`BUILD_CVODES <BUILD_CVODES (CMake option)>`
   Build the CVODES library

   Default: ``ON``

:index:`BUILD_IDA <BUILD_IDA (CMake option)>`
   Build the IDA library

   Default: ``ON``

:index:`BUILD_IDAS <BUILD_IDAS (CMake option)>`
   Build the IDAS library

   Default: ``ON``

:index:`BUILD_KINSOL <BUILD_KINSOL (CMake option)>`
   Build the KINSOL library

   Default: ``ON``

:index:`BUILD_SHARED_LIBS <BUILD_SHARED_LIBS (CMake option)>`
   Build shared libraries

   Default: ``ON``

:index:`BUILD_STATIC_LIBS <BUILD_STATIC_LIBS (CMake option)>`
   Build static libraries

   Default: ``ON``

:index:`CMAKE_BUILD_TYPE <CMAKE_BUILD_TYPE (CMake option)>`
   Choose the type of build, options are:
   ``None`` (``CMAKE_C_FLAGS`` used), ``Debug``, ``Release``,
   ``RelWithDebInfo``, and ``MinSizeRel``

   Default:

   .. note:: Specifying a build type will trigger the corresponding
             build type specific compiler flag options below which
             will be appended to the flags set by
             ``CMAKE_<language>_FLAGS``.

:index:`CMAKE_C_COMPILER <CMAKE_C_COMPILER (CMake option)>`
   C compiler

   Default: ``/usr/bin/cc``

:index:`CMAKE_C_FLAGS <CMAKE_C_FLAGS (CMake option)>`
   Flags for C compiler

   Default:

:index:`CMAKE_C_FLAGS_DEBUG <CMAKE_C_FLAGS_DEBUG (CMake option)>`
   Flags used by the C compiler during debug
   builds

   Default: ``-g``

:index:`CMAKE_C_FLAGS_MINSIZEREL <CMAKE_C_FLAGS_MINSIZEREL (CMake option)>`
   Flags used by the C compiler during release minsize builds

   Default: ``-Os -DNDEBUG``

:index:`CMAKE_C_FLAGS_RELEASE <CMAKE_C_FLAGS_RELEASE (CMake option)>`
   Flags used by the C compiler during release
   builds

   Default: ``-O3 -DNDEBUG``

:index:`CMAKE_CXX_COMPILER <CMAKE_CXX_COMPILER (CMake option)>`
   C++ compiler

   Default: ``/usr/bin/c++``

   .. note:: A C++ compiler (and all related options) are only are
             triggered if C++ examples are enabled
             (``EXAMPLES_ENABLE_CXX`` is ON). All SUNDIALS solvers can
             be used from C++ applications by default without setting
             any additional configuration options.

:index:`CMAKE_CXX_FLAGS <CMAKE_CXX_FLAGS (CMake option)>`
   Flags for C++ compiler

   Default:

:index:`CMAKE_CXX_FLAGS_DEBUG <CMAKE_CXX_FLAGS_DEBUG (CMake option)>`
   Flags used by the C++ compiler during debug builds

   Default: ``-g``

:index:`CMAKE_CXX_FLAGS_MINSIZEREL <CMAKE_CXX_FLAGS_MINSIZEREL (CMake option)>`
   Flags used by the C++ compiler during release minsize builds

   Default: ``-Os -DNDEBUG``

:index:`CMAKE_CXX_FLAGS_RELEASE <CMAKE_CXX_FLAGS_RELEASE (CMake option)>`
   Flags used by the C++ compiler during release builds

   Default: ``-O3 -DNDEBUG``

:index:`CMAKE_CXX_STANDARD <CMAKE_CXX_STANDARD (CMake option)>`
   The C++ standard to build C++ parts of SUNDIALS with.

   Default: 11

   Note: Options are 98, 11, 14, 17, 20. This option is only used when a
   C++ compiler is required.

:index:`CMAKE_Fortran_COMPILER <CMAKE_Fortran_COMPILER (CMake option)>`
   Fortran compiler

   Default: ``/usr/bin/gfortran``

   .. note:: Fortran support (and all related options) are triggered only if
             either Fortran-C support is (``FCMIX_ENABLE`` is ON) or
             LAPACK support is enabled (``ENABLE_LAPACK`` is ``ON``).

:index:`CMAKE_Fortran_FLAGS <CMAKE_Fortran_FLAGS (CMake option)>`
   Flags for Fortran compiler

   Default:

:index:`CMAKE_Fortran_FLAGS_DEBUG <CMAKE_Fortran_FLAGS_DEBUG (CMake option)>`
   Flags used by the Fortran compiler during debug builds

   Default: ``-g``

:index:`CMAKE_Fortran_FLAGS_MINSIZEREL <CMAKE_Fortran_FLAGS_MINSIZEREL (CMake option)>`
   Flags used by the Fortran compiler during release minsize builds

   Default: ``-Os``

:index:`CMAKE_Fortran_FLAGS_RELEASE <CMAKE_Fortran_FLAGS_RELEASE (CMake option)>`
   Flags used by the Fortran compiler during release builds

   Default: ``-O3``

:index:`CMAKE_INSTALL_PREFIX <CMAKE_INSTALL_PREFIX (CMake option)>`
   Install path prefix, prepended onto install directories

   Default: ``/usr/local``

   .. note:: The user must have write access to the location specified
             through this option. Exported SUNDIALS header files and libraries
             will be installed under subdirectories ``include`` and ``lib`` of
            ``CMAKE_INSTALL_PREFIX``, respectively.

:index:`ENABLE_CUDA <ENABLE_CUDA (CMake option)>`
   Build the SUNDIALS CUDA modules.

   Default: ``OFF``

:index:`CMAKE_CUDA_ARCHITECTURES <CMAKE_CUDA_ARCHITECTURES (CMake option)>`
   Specifies the CUDA architecture to compile for.

   Default: ``sm_30``

:index:`ENABLE_XBRAID <ENABLE_XBRAID (CMake option)>`
   Enable or disable the ARKStep + XBraid interface.

   Default: ``OFF``

   .. note:: See additional information on building with *XBraid*
             enabled in  :ref:`Installation.CMake.ExternalLibraries`.

:index:`EXAMPLES_ENABLE_C <EXAMPLES_ENABLE_C (CMake option)>`
   Build the SUNDIALS C examples

   Default: ``ON``

:index:`EXAMPLES_ENABLE_CXX <EXAMPLES_ENABLE_CXX (CMake option)>`
   Build the SUNDIALS C++ examples

   Default: ``OFF``

:index:`EXAMPLES_ENABLE_CUDA <EXAMPLES_ENABLE_CUDA (CMake option)>`
   Build the SUNDIALS CUDA examples

   Default: ``OFF``

   .. note:: You need to enable CUDA support to build these examples.

:index:`EXAMPLES_ENABLE_F77 <EXAMPLES_ENABLE_F77 (CMake option)>`
   Build the SUNDIALS Fortran77 examples

   Default: ``ON`` (if ``FCMIX_ENABLE`` is ``ON``)

:index:`EXAMPLES_ENABLE_F90 <EXAMPLES_ENABLE_F90 (CMake option)>`
   Build the SUNDIALS Fortran90 examples

   Default: ``ON`` (if ``BUILD_FORTRAN77_INTERFACE`` is ``ON``)

:index:`EXAMPLES_ENABLE_F2003 <EXAMPLES_ENABLE_F2003 (CMake option)>`
   Build the SUNDIALS Fortran2003 examples

   Default: ``ON`` (if ``BUILD_FORTRAN_MODULE_INTERFACE`` is ``ON``)

:index:`EXAMPLES_INSTALL <EXAMPLES_INSTALL (CMake option)>`
   Install example files

   Default: ``ON``

   .. note:: This option is triggered when any of the SUNDIALS
             example programs are enabled
             (``EXAMPLES_ENABLE_<language>`` is ``ON``). If the user
             requires installation of example programs then the
             sources and sample output files for all SUNDIALS modules
             that are currently enabled will be exported to the
             directory specified by ``EXAMPLES_INSTALL_PATH``. A CMake
             configuration script will also be automatically generated
             and exported to the same directory. Additionally, if the
             configuration is done under a Unix-like system, makefiles
             for the compilation of the example programs (using the
             installed SUNDIALS libraries) will be automatically
             generated and exported to the directory specified by
             ``EXAMPLES_INSTALL_PATH``.

:index:`EXAMPLES_INSTALL_PATH <EXAMPLES_INSTALL_PATH (CMake option)>`
   Output directory for installing example
   files

   Default: ``/usr/local/examples``

   .. note:: The actual default value for this option will be an
             ``examples`` subdirectory created under ``CMAKE_INSTALL_PREFIX``.

:index:`BUILD_FORTRAN77_INTERFACE <BUILD_FORTRAN77_INTERFACE (CMake option)>`
   Enable Fortran77-C interface

   Default: ``OFF``

:index:`BUILD_FORTRAN_MODULE_INTERFACE <BUILD_FORTRAN_MODULE_INTERFACE (CMake option)>`
   Enable Fortran2003 interface

   Default: ``OFF``

:index:`ENABLE_HYPRE <ENABLE_HYPRE (CMake option)>`
   Flag to enable *hypre* support

   Default: ``OFF``

   .. note:: See additional information on building with *hypre*
             enabled in  :ref:`Installation.CMake.ExternalLibraries`.

:index:`HYPRE_INCLUDE_DIR <HYPRE_INCLUDE_DIR (CMake option)>`
   Path to *hypre* header files

   Default: none

:index:`HYPRE_LIBRARY <HYPRE_LIBRARY (CMake option)>`
   Path to *hypre* installed library files

   Default: none

:index:`ENABLE_KLU <F90_ENABLE (CMake option)>`
   Enable KLU support

   Default: ``OFF``

   .. note:: See additional information on building with KLU
             enabled in :ref:`Installation.CMake.ExternalLibraries`.

:index:`KLU_INCLUDE_DIR <KLU_INCLUDE_DIR (CMake option)>`
   Path to SuiteSparse header files

   Default: none

:index:`KLU_LIBRARY_DIR <KLU_LIBRARY_DIR (CMake option)>`
   Path to SuiteSparse installed library files

   Default: none

:index:`ENABLE_LAPACK <ENABLE_LAPACK (CMake option)>`
   Enable LAPACK support

   Default: ``OFF``

   .. note:: Setting this option to ``ON`` will trigger additional CMake
             options. See additional information on building with
             LAPACK enabled in :ref:`Installation.CMake.ExternalLibraries`.

:index:`LAPACK_LIBRARIES <LAPACK_LIBRARIES (CMake option)>`
   LAPACK (and BLAS) libraries

   Default: ``/usr/lib/liblapack.so;/usr/lib/libblas.so``

   .. note:: CMake will search for libraries in your
      ``LD_LIBRARY_PATH`` prior to searching default system
      paths.

:index:`ENABLE_MPI <ENABLE_MPI (CMake option)>`
   Enable MPI support. This will build the parallel nvector
   and the MPI-aware version of the ManyVector library.

   Default: ``OFF``

   .. note:: Setting this option to ``ON`` will trigger several additional
             options related to MPI.

:index:`MPI_C_COMPILER <MPI_C_COMPILER (CMake option)>`
   ``mpicc`` program

   Default:

:index:`MPI_CXX_COMPILER <MPI_CXX_COMPILER (CMake option)>`
   ``mpicxx`` program

   Default:

   .. note:: This option is triggered only if MPI is enabled
             (``ENABLE_MPI`` is ``ON``) and C++ examples are enabled
             (``EXAMPLES_ENABLE_CXX`` is ``ON``). All SUNDIALS
             solvers can be used from C++ MPI applications by default
             without setting any additional configuration options
             other than ``ENABLE_MPI``.

:index:`MPI_Fortran_COMPILER <MPI_Fortran_COMPILER (CMake option)>`
   ``mpif77`` or ``mpif90`` program

   Default:

   .. note:: This option is triggered only if MPI is enabled
             (``ENABLE_MPI`` is ``ON``) and Fortran-C support is
             enabled (``EXAMPLES_ENABLE_F77`` or ``EXAMPLES_ENABLE_F90`` are ``ON``).

:index:`MPIEXEC_EXECUTABLE <MPIEXEC_EXECUTABLE (CMake option)>`
   Specify the executable for running MPI programs

   Default: ``mpirun``

   .. note:: This option is triggered only if MPI is enabled (``ENABLE_MPI`` is ``ON``).

:index:`ENABLE_OPENMP <ENABLE_OPENMP (CMake option)>`
   Enable OpenMP support (build the OpenMP NVector)

   Default: ``OFF``

:index:`ENABLE_PETSC <ENABLE_PETSC (CMake option)>`
   Enable PETSc support

   Default: ``OFF``

   .. note:: See additional information on building with
             PETSc enabled in :ref:`Installation.CMake.ExternalLibraries`.

:index:`PETSC_DIR <PETSC_DIR (CMake option)>`
   Path to PETSc installation

   Default: none

:index:`PETSC_LIBRARIES <PETSC_LIBRARIES (CMake option)>` (advanced option)
   Semi-colon separated list of PETSc link libraries. Unless provided by the
   user, this is autopopulated based on the PETSc installation found in
   ``PETSC_DIR``.

   Default: none

:index:`PETSC_INCLUDES <PETSC_INCLUDES (CMake option)>` (advanced option)
   Semi-colon separated list of PETSc include directroies. Unless provided by
   the user, this is autopopulated based on the PETSc installation found in
   ``PETSC_DIR``.

   Default: none

:index:`ENABLE_PTHREAD <ENABLE_PTHREAD (CMake option)>`
   Enable Pthreads support (build the Pthreads NVector)

   Default: ``OFF``

:index:`ENABLE_RAJA <ENABLE_RAJA (CMake option)>`
   Enable RAJA support.

   Default: OFF

   .. note:: You need to enable CUDA or HIP in order to build the
             RAJA vector module.

:index:`SUNDIALS_RAJA_BACKENDS <SUNDIALS_RAJA_BACKENDS (CMake option)>`
   If building SUNDIALS with RAJA support, this sets the RAJA
   backend to target. Values supported are CUDA and HIP.

   Default: CUDA

:index:`ENABLE_SUPERLUDIST <ENABLE_SUPERLUDIST (CMake option)>`
   Enable SuperLU_DIST support

   Default: ``OFF``

   .. note:: See additional information on building wtih
             SuperLU_DIST enabled in :ref:`Installation.CMake.ExternalLibraries`.

:index:`SUPERLUDIST_INCLUDE_DIR <SUPERLUDIST_INCLUDE_DIR (CMake option)>`
   Path to SuperLU_DIST header files (under a typical SuperLU_DIST
   install, this is typically the SuperLU_DIST ``SRC`` directory)

   Default: none

:index:`SUPERLUDIST_LIBRARY_DIR <SUPERLUDIST_LIBRARY_DIR (CMake option)>`
   Path to SuperLU_DIST installed library files

   Default: none

:index:`SUPERLUDIST_LIBRARIES <SUPERLUDIST_LIBRARIES (CMake option)>`
   Semi-colon separated list of libraries needed for SuperLU_DIST

   Default: none

:index:`SUPERLUDIST_OpenMP <SUPERLUDIST_OpenMP (CMake option)>`
   Enable SUNDIALS support for SuperLU_DIST built with OpenMP

   Default: none

   Note: SuperLU_DIST must be built with OpenMP support for this option to function.
   Additionally the environment variable ``OMP_NUM_THREADS`` must be set to the desired
   number of threads.

:index:`ENABLE_SUPERLUMT <ENABLE_SUPERLUMT (CMake option)>`
   Enable SuperLU_MT support

   Default: ``OFF``

   .. note:: See additional information on building with
             SuperLU_MT enabled in :ref:`Installation.CMake.ExternalLibraries`.

:index:`SUPERLUMT_INCLUDE_DIR <SUPERLUMT_INCLUDE_DIR (CMake option)>`
   Path to SuperLU_MT header files (under a typical SuperLU_MT
   install, this is typically the SuperLU_MT ``SRC`` directory)

   Default: none

:index:`SUPERLUMT_LIBRARY_DIR <SUPERLUMT_LIBRARY_DIR (CMake option)>`
   Path to SuperLU_MT installed library files

   Default: none

:index:`SUPERLUMT_THREAD_TYPE <SUPERLUMT_THREAD_TYPE (CMake option)>`
   Must be set to Pthread or OpenMP, depending on how SuperLU_MT was compiled.

   Default: Pthread

:index:`ENABLE_SYCL <ENABLE_SYCL (CMake option)>`
   Enable SYCL support.

   Default: OFF

   .. note::

      At present the only supported SYCL compiler is the DPC++ (Intel oneAPI)
      compiler. CMake does not currently support autodetection of SYCL compilers
      and ``CMAKE_CXX_COMPILER`` must be set to a valid SYCL compiler i.e.,
      ``dpcpp`` in order to build with SYCL support.

:index:`SUNDIALS_BUILD_WITH_MONITORING <SUNDIALS_BUILD_WITH_MONITORING (CMake option)>`
   Build SUNDIALS with capabilties for fine-grained monitoring of solver progress
   and statistics. This is primarily useful for debugging.

   Default: OFF

   Note: Building with monitoring may result in minor performance degradation
   even if monitoring is not utilized.

:index:`SUNDIALS_F77_FUNC_CASE <SUNDIALS_F77_FUNC_CASE (CMake option)>`
   Specify the case to use in the Fortran name-mangling scheme,
   options are: ``lower`` or ``upper``

   Default:

   Note: The build system will attempt to infer the Fortran
   name-mangling scheme using the Fortran compiler. This option should
   only be used if a Fortran compiler is not available or to override
   the inferred or default (``lower``) scheme if one can not be
   determined. If used, ``SUNDIALS_F77_FUNC_UNDERSCORES`` must also
   be set.

:index:`SUNDIALS_F77_FUNC_UNDERSCORES <SUNDIALS_F77_FUNC_UNDERSCORES (CMake option)>`
   Specify the number of underscores to append in the Fortran
   name-mangling scheme, options are: ``none``, ``one``, or ``two``

   Default:

   Note: The build system will attempt to infer the Fortran
   name-mangling scheme using the Fortran compiler. This option should
   only be used if a Fortran compiler is not available or to override
   the inferred or default (``one``) scheme if one can not be
   determined. If used, ``SUNDIALS_F77_FUNC_CASE`` must also be set.

:index:`SUNDIALS_INDEX_TYPE <SUNDIALS_INDEX_TYPE (CMake option)>` (advanced)
   Integer type used for SUNDIALS indices.  The size must match the size provided for
   the ``SUNDIALS_INDEX_SIZE`` option.

   Default:

   Note: In past SUNDIALS versions, a user could set this option to
   ``INT64_T`` to use 64-bit integers, or ``INT32_T`` to use 32-bit
   integers. Starting in SUNDIALS 3.2.0, these special values are
   deprecated. For SUNDIALS 3.2.0 and up, a user will only need to use
   the ``SUNDIALS_INDEX_SIZE`` option in most cases.

:index:`SUNDIALS_INDEX_SIZE <SUNDIALS_INDEX_SIZE (CMake option)>`
   Integer size (in bits) used for indices in SUNDIALS, options are: ``32`` or ``64``

   Default: ``64``

   Note: The build system tries to find an integer type of appropriate
   size. Candidate 64-bit integer types are (in order of preference):
   ``int64_t``, ``__int64``, ``long long``, and ``long``.  Candidate
   32-bit integers are (in order of preference): ``int32_t``,
   ``int``, and ``long``.  The advanced option,
   ``SUNDIALS_INDEX_TYPE`` can be used to provide a type not listed
   here.

:index:`SUNDIALS_PRECISION <SUNDIALS_PRECISION (CMake option)>`
   Precision used in SUNDIALS, options are: ``double``, ``single`` or
   ``extended``

   Default: ``double``

:index:`SUNDIALS_INSTALL_CMAKEDIR <SUNDIALS_INSTALL_CMAKEDIR (CMake option)>`
  Installation directory for the SUNDIALS cmake files (relative to ``CMAKE_INSTALL_PREFIX``).

  Default: ``CMAKE_INSTALL_PREFIX/cmake/sundials``

:index:`USE_GENERIC_MATH <USE_GENERIC_MATH (CMake option)>`
   Use generic (``stdc``) math libraries

   Default: ``ON``

:index:`XBRAID_DIR <XBRAID_DIR (CMake option)>`
   The root directory of the XBraid installation.

   Default: ``OFF``

:index:`XBRAID_INCLUDES <XBRAID_INCLUDES (CMake option)>`
   Semi-colon separated list of XBraid include directories. Unless provided by
   the user, this is autopopulated based on the XBraid installation found in
   ``XBRAID_DIR``.

   Default: none

:index:`XBRAID_LIBRARIES <XBRAID_LIBRARIES (CMake option)>`
   Semi-colon separated list of XBraid link libraries. Unless provided by
   the user, this is autopopulated based on the XBraid installation found in
   ``XBRAID_DIR``.

   Default: none

:index:`USE_XSDK_DEFAULTS <USE_XSDK_DEFAULTS (xSDK CMake option)>`
   Enable xSDK (see `https://xsdk.info <https://xsdk.info>`_ for more
   information) default configuration settings. This sets ``CMAKE_BUILD_TYPE``
   to ``Debug``, ``SUNDIALS_INDEX_SIZE`` to 32 and ``SUNDIALS_PRECISION`` to
   double.

   Default: ``OFF``


.. _Installation.CMake.Examples:

Configuration examples
-----------------------------------

The following examples will help demonstrate usage of the CMake
configure options.

To configure SUNDIALS using the default C and Fortran compilers,
and default ``mpicc`` and ``mpif77`` parallel compilers,
enable compilation of examples, and install libraries, headers, and
example sources under subdirectories of ``/home/myname/sundials/``, use:

.. code-block:: bash

   % cmake \
   > -DCMAKE_INSTALL_PREFIX=/home/myname/sundials/instdir \
   > -DEXAMPLES_INSTALL_PATH=/home/myname/sundials/instdir/examples \
   > -DENABLE_MPI=ON \
   > -DFCMIX_ENABLE=ON \
   > /home/myname/sundials/srcdir

   % make install


To disable installation of the examples, use:

.. code-block:: bash

   % cmake \
   > -DCMAKE_INSTALL_PREFIX=/home/myname/sundials/instdir \
   > -DEXAMPLES_INSTALL_PATH=/home/myname/sundials/instdir/examples \
   > -DENABLE_MPI=ON \
   > -DFCMIX_ENABLE=ON \
   > -DEXAMPLES_INSTALL=OFF \
   > /home/myname/sundials/srcdir

   % make install




.. _Installation.CMake.ExternalLibraries:

Working with external Libraries
-----------------------------------

The SUNDIALS suite contains many options to enable implementation
flexibility when developing solutions. The following are some notes
addressing specific configurations when using the supported third
party libraries.



.. _Installation.CMake.ExternalLibraries.LAPACK:

Building with LAPACK
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To enable LAPACK, set the ``ENABLE_LAPACK`` option to ``ON``.
If the directory containing the LAPACK library is in the
``LD_LIBRARY_PATH`` environment variable, CMake will set the
``LAPACK_LIBRARIES`` variable accordingly, otherwise CMake will
attempt to find the LAPACK library in standard system locations. To
explicitly tell CMake what library to use, the ``LAPACK_LIBRARIES``
variable can be set to the desired libraries required for LAPACK.


.. code-block:: bash

   % cmake \
   > -DCMAKE_INSTALL_PREFIX=/home/myname/sundials/instdir \
   > -DEXAMPLES_INSTALL_PATH=/home/myname/sundials/instdir/examples \
   > -DENABLE_LAPACK=ON \
   > -DLAPACK_LIBRARIES=/mylapackpath/lib/libblas.so;/mylapackpath/lib/liblapack.so \
   > /home/myname/sundials/srcdir

   % make install

.. note:: If a working Fortran compiler is not available to infer the
          Fortran name-mangling scheme, the options
          ``SUNDIALS_F77_FUNC_CASE`` and
          ``SUNDIALS_F77_FUNC_UNDERSCORES`` *must* be set in order to
          bypass the check for a Fortran compiler and define the
          name-mangling scheme. The defaults for these options in
          earlier versions of SUNDIALS were ``lower`` and ``one``,
          respectively.



.. _Installation.CMake.ExternalLibraries.KLU:

Building with KLU
^^^^^^^^^^^^^^^^^^^^^^^^^^^

The KLU libraries are part of SuiteSparse, a suite of sparse matrix
software, available from the Texas A&M University website:
http://faculty.cse.tamu.edu/davis/suitesparse.html .

SUNDIALS has been tested with SuiteSparse version 5.7.2.  To enable
KLU, set ``ENABLE_KLU`` to ``ON``, set ``KLU_INCLUDE_DIR`` to the
``include`` path of the KLU installation and set ``KLU_LIBRARY_DIR``
to the ``lib`` path of the KLU installation.  The CMake configure will
result in populating the following variables: ``AMD_LIBRARY``,
``AMD_LIBRARY_DIR``,  ``BTF_LIBRARY``, ``BTF_LIBRARY_DIR``,
``COLAMD_LIBRARY``, ``COLAMD_LIBRARY_DIR``, and ``KLU_LIBRARY``.


.. _Installation.CMake.ExternalLibraries.SuperLU_MT:

Building with SuperLU_DIST
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The SuperLU_DIST libraries are available for download from the Lawrence
Berkeley National Laboratory website:
http://crd-legacy.lbl.gov/$\sim$xiaoye/SuperLU/\#superlu_dist.

SUNDIALS has been tested with SuperLU_DIST 6.1.1. To enable
SuperLU_DIST, set  ``ENABLE_SUPERLUDIST`` to ``ON``, set
``SUPERLUDIST_INCLUDE_DIR`` to the ``SRC`` path of the SuperLU_DIST
installation, and set the variable ``SUPERLUMT_LIBRARY_DIR`` to the
``lib`` path of the SuperLU_DIST installation.  At the same time, the
variable ``SUPERLUDIST_LIBRARIES`` must be set to a semi-colon separated list
of other libraries SuperLU_DIST depends on. For example, if SuperLU_DIST
was built with LAPACK, then include the LAPACK library in this list.
If SuperLU_DIST was built with OpenMP support, then you may set
``SUPERLUDIST_OpenMP`` to ``ON`` utilize the OpenMP functionality of
SuperLU_DIST.


Building with SuperLU_MT
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The SuperLU_MT libraries are available for download from the Lawrence
Berkeley National Laboratory website:
http://crd-legacy.lbl.gov/$\sim$xiaoye/SuperLU/\#superlu_mt .

SUNDIALS has been tested with SuperLU_MT version 3.1.  To enable
SuperLU_MT, set  ``ENABLE_SUPERLUMT`` to ``ON``, set
``SUPERLUMT_INCLUDE_DIR`` to the ``SRC`` path of the SuperLU_MT
installation, and set the variable ``SUPERLUMT_LIBRARY_DIR`` to the
``lib`` path of the SuperLU_MT installation. At the same time, the
variable ``SUPERLUMT_LIBRARIES`` must be set to a semi-colon separated
list of other libraries SuperLU_MT depends on. For example, if
SuperLU_MT was build with an external blas library, then include the
full path to the blas library in this list. Additionally, the
variable ``SUPERLUMT_THREAD_TYPE`` must be set to either ``Pthread``
or ``OpenMP``.

Do not mix thread types when building SUNDIALS solvers.
If threading is enabled for SUNDIALS by having either
``ENABLE_OPENMP`` or ``ENABLE_PTHREAD`` set to ``ON`` then SuperLU_MT
should be set to use the same threading type.


.. _Installation.CMake.ExternalLibraries.PETSc:

Building with PETSc
^^^^^^^^^^^^^^^^^^^^^^^^^^^

The PETSc libraries are available for download from the Argonne
National Laboratory website:
http://www.mcs.anl.gov/petsc .

SUNDIALS has been tested with PETSc version 3.10.0 - 3.14.0. To enable PETSc,
set ``ENABLE_PETSC`` to ``ON``, and set ``PETSC_DIR`` to the path of the PETSc
installation. Alternatively, a user can provide a list of inlcude paths in
``PETSC_INCLUDES`` and a list of complete paths to the PETSc libraries in
``PETSC_LIBRARIES``.


.. _Installation.CMake.ExternalLibraries.hypre:

Building with *hypre*
^^^^^^^^^^^^^^^^^^^^^^^^^^^

The *hypre* libraries are available for download from the Lawrence
Livermore National Laboratory website:
`http://computing.llnl.gov/projects/hypre <http://computing.llnl.gov/projects/hypre>`_.
SUNDIALS has been tested with *hypre* version 2.19.0.
To enable *hypre*, set  ``ENABLE_HYPRE`` to ``ON``, set ``HYPRE_INCLUDE_DIR``
to the ``include`` path of the *hypre* installation, and set the variable
``HYPRE_LIBRARY_DIR`` to the ``lib`` path of the *hypre* installation.

Note: SUNDIALS must be configured so that ``SUNDIALS_INDEX_SIZE`` (or
equivalently, ``XSDK_INDEX_SIZE``) equals the precision of
``HYPRE_BigInt`` in the corresponding *hypre* installation.


.. _Installation.CMake.ExternalLibraries.CUDA:

Building with CUDA
^^^^^^^^^^^^^^^^^^^^^^

SUNDIALS CUDA modules and examples have been tested with version 10 and 11
of the CUDA toolkit. To build them, you need to install the Toolkit and compatible
NVIDIA drivers. Both are available for download from the NVIDIA website:
`https://developer.nvidia.com/cuda-downloads
<https://developer.nvidia.com/cuda-downloads>`_. To enable CUDA,
set ``ENABLE_CUDA`` to ``ON``. If CUDA is installed in a nonstandard
location, you may be prompted to set the variable
``CUDA_TOOLKIT_ROOT_DIR`` with your CUDA Toolkit installation
path. To enable CUDA examples, set ``EXAMPLES_ENABLE_CUDA`` to ``ON``.


.. _Installation.CMake.ExternalLibraries.RAJA:

Building with RAJA
^^^^^^^^^^^^^^^^^^^^^

RAJA is a performance portability layer developed by Lawrence
Livermore National Laboratory and can be obtained from
`https://github.com/LLNL/RAJA <https://github.com/LLNL/RAJA>`_.
SUNDIALS RAJA modules and examples have been tested with RAJA
version 0.12.1. Building SUNDIALS RAJA modules requires a CUDA-enabled
RAJA installation. To enable RAJA, set ``ENABLE_CUDA`` and
``ENABLE_RAJA`` to ``ON``. If RAJA is installed in a nonstandard
location you will be prompted to set the variable ``RAJA_DIR`` with
the path to the RAJA CMake configuration file. To enable building the
RAJA examples set ``EXAMPLES_ENABLE_CUDA`` to ``ON``.


.. _Installation.CMake.ExternalLibraries.XBraid:

Building with XBraid
^^^^^^^^^^^^^^^^^^^^

The XBraid library is available for download from the XBraid GitHub:
`https://github.com/XBraid/xbraid <https://github.com/XBraid/xbraid>`_.
SUNDIALS has been tested with XBraid version 3.0.0.
To enable XBraid, set  ``ENABLE_XBRAID`` to ``ON``, set ``XBRAID_DIR``
to the root install location of XBraid or the location of the clone of the
XBraid repository.

Note: At this time the XBraid types ``braid_Int`` and ``braid_Real`` are
hard-coded to ``int`` and ``double`` respectively. As such SUNDIALS must be
configured with ``SUNDIALS_INDEX_SIZE`` set to ``32`` and ``SUNDIALS_PRECISION``
set to ``double``. Additionally, SUNDIALS must be configured with ``ENABLE_MPI``
set to ``ON``.



.. _Installation.CMake.Testing:

Testing the build and installation
---------------------------------------

If SUNDIALS was configured with ``EXAMPLES_ENABLE_<language>`` options
to ``ON``, then a set of regression tests can be run after building
with the ``make`` command by running:

.. code-block:: bash

   % make test

Additionally, if ``EXAMPLES_INSTALL`` was also set to ``ON``, then a
set of smoke tests can be run after installing with the ``make install``
command by running:

.. code-block:: bash

   % make test_install


.. _Installation.CMake.BuildRunExamples:

Building and Running Examples
-------------------------------------

Each of the SUNDIALS solvers is distributed with a set of examples
demonstrating basic usage. To build and install the examples, set at
least of the ``EXAMPLES_ENABLE_<language>`` options to ``ON``, and
set ``EXAMPLES_INSTALL`` to ``ON``. Specify
the installation path for the examples with the variable
``EXAMPLES_INSTALL_PATH``. CMake will generate ``CMakeLists.txt``
configuration files (and ``Makefile`` files if on Linux/Unix) that
reference the *installed* SUNDIALS headers and libraries.

Either the ``CMakeLists.txt`` file or the traditional ``Makefile`` may
be used to build the examples as well as serve as a template for
creating user developed solutions.  To use the supplied ``Makefile``
simply run ``make`` to compile and generate the executables.  To use
CMake from within the installed example directory, run ``cmake`` (or
``ccmake`` or ``cmake-gui`` to use the GUI) followed by ``make`` to
compile the example code.  Note that if CMake is used, it will
overwrite the traditional ``Makefile`` with a new CMake-generated
``Makefile``.

The resulting output from running the examples can be compared with
example output bundled in the SUNDIALS distribution.

NOTE: There will potentially be differences in the output due to
machine architecture, compiler versions, use of third party libraries etc.




.. _Installation.CMake.Windows:

Configuring, building, and installing on Windows
----------------------------------------------------------------

CMake can also be used to build SUNDIALS on Windows. To build SUNDIALS
for use with Visual Studio the following steps should be performed:

1. Unzip the downloaded tar file(s) into a directory. This will be the
   ``SOLVERDIR``

2. Create a separate ``BUILDDIR``

3. Open a Visual Studio Command Prompt and cd to ``BUILDDIR``

4. Run ``cmake-gui ../SOLVERDIR``

   a. Hit Configure

   b. Check/Uncheck solvers to be built

   c. Change ``CMAKE_INSTALL_PREFIX`` to ``INSTDIR``

   d. Set other options as desired

   e. Hit Generate

5. Back in the VS Command Window:

   a. Run ``msbuild ALL_BUILD.vcxproj``

   b. Run ``msbuild INSTALL.vcxproj``


The resulting libraries will be in the ``INSTDIR``.

The SUNDIALS project can also now be opened in Visual Studio.
Double click on the ``ALL_BUILD.vcxproj`` file to open the project.
Build the whole *solution* to create the SUNDIALS libraries.
To use the SUNDIALS libraries in your own projects, you must
set the include directories for your project,
add the SUNDIALS libraries to your project solution,
and set the SUNDIALS libraries as dependencies for your project.




.. _Installation.Results:

Installed libraries and exported header files
====================================================

Using the CMake SUNDIALS build system, the command

.. code-block:: bash

   $ make install

will install the libraries under ``LIBDIR`` and the public header
files under ``INCLUDEDIR``. The values for these directories
are ``INSTDIR/lib`` and ``INSTDIR/include``, respectively.  The
location can be changed by setting the CMake variable
``CMAKE_INSTALL_PREFIX``.  Although all installed libraries reside
under ``LIBDIR/lib``, the public header files are further organized
into subdirectories under ``INCLUDEDIR/include``.

The installed libraries and exported header files are listed for
reference in the :ref:`Table: SUNDIALS libraries and header files
<Installation.Table>`. The file extension ``.LIB`` is typically ``.so``
for shared libraries and ``.a`` for static libraries. Note that, in
this table names are relative to ``LIBDIR`` for libraries and to
``INCLUDEDIR`` for header files.

A typical user program need not explicitly include any of the shared
SUNDIALS header files from under the ``INCLUDEDIR/include/sundials``
directory since they are explicitly included by the appropriate solver
header files (e.g., ``cvode_dense.h`` includes
``sundials_dense.h``). However, it is both legal and safe to do so,
and would be useful, for example, if the functions declared in
``sundials_dense.h`` are to be used in building a preconditioner.


Using SUNDIALS as a Third Party Library in other CMake Projects
---------------------------------------------------------------

The ``make install`` command will also install a `CMake package configuration file
<https://cmake.org/cmake/help/v3.12/manual/cmake-packages.7.html\#package-configuration-file>`_
that other CMake projects can load to get all the information needed to build
against SUNDIALS. In the consuming project's CMake code, the ``find_package``
command may be used to search for the configuration file, which will be
installed to ``instdir/SUNDIALS_INSTALL_CMAKEDIR/SUNDIALSConfig.cmake``
alongside a package version file
``instdir/SUNDIALS_INSTALL_CMAKEDIR/SUNDIALSConfigVersion.cmake``. Together
these files contain all the information the consuming project needs to use
SUNDIALS, including exported CMake targets. The SUNDIALS exported CMake targets
follow the same naming convention as the generated library binaries, e.g. the
exported target for CVODE is ``SUNDIALS::cvode``. The CMake code snipped
below shows how a consuming project might leverage the SUNDIALS package
configuration file to build against SUNDIALS in their own CMake project.

.. code-block:: none

  project(MyProject)

  # Set the variable SUNDIALS_DIR to the SUNDIALS instdir.
  # When using the cmake CLI command, this can be done like so:
  #   cmake -D SUNDIALS_DIR=/path/to/sundials/installation

  find_project(SUNDIALS REQUIRED)

  add_executable(myexec main.c)

  # Link to SUNDIALS libraries through the exported targets.
  # This is just an example, users should link to the targets appropriate
  # for their use case.
  target_link_libraries(myexec PUBLIC SUNDIALS::cvode SUNDIALS::nvecpetsc)



.. _Installation.Table:

.. tabularcolumns:: |\Y{0.3}|\Y{0.1}|\Y{0.6}|

.. Table:: SUNDIALS shared libraries and header files

   +------------------------------+--------------+----------------------------------------------+
   | Shared                       | Headers      | ``sundials/sundials_band.h``                 |
   |                              |              +----------------------------------------------+
   |                              |              | ``sundials/sundials_config.h``               |
   |                              |              +----------------------------------------------+
   |                              |              | ``sundials/sundials_cuda_policies.hpp``      |
   |                              |              +----------------------------------------------+
   |                              |              | ``sundials/sundials_dense.h``                |
   |                              |              +----------------------------------------------+
   |                              |              | ``sundials/sundials_direct.h``               |
   |                              |              +----------------------------------------------+
   |                              |              | ``sundials/sundials_fconfig.h``              |
   |                              |              +----------------------------------------------+
   |                              |              | ``sundials/sundials_fnvector.h``             |
   |                              |              +----------------------------------------------+
   |                              |              | ``sundials/sundials_iterative.h``            |
   |                              |              +----------------------------------------------+
   |                              |              | ``sundials/sundials_linearsolver.h``         |
   |                              |              +----------------------------------------------+
   |                              |              | ``sundials/sundials_nonlinearsolver.h``      |
   |                              |              +----------------------------------------------+
   |                              |              | ``sundials/sundials_matrix.h``               |
   |                              |              +----------------------------------------------+
   |                              |              | ``sundials/sundials_math.h``                 |
   |                              |              +----------------------------------------------+
   |                              |              | ``sundials/sundials_nvector.h``              |
   |                              |              +----------------------------------------------+
   |                              |              | ``sundials/sundials_types.h``                |
   |                              |              +----------------------------------------------+
   |                              |              | ``sundials/sundials_version.h``              |
   |                              |              +----------------------------------------------+
   |                              |              | ``sundials/sundials_xbraid.h``               |
   +------------------------------+--------------+----------------------------------------------+
   |                                                                                            |
   | NVECTOR Modules                                                                            |
   |                                                                                            |
   +------------------------------+--------------+----------------------------------------------+
   | SERIAL                       | Libraries    | ``libsundials_nvecserial.LIB``               |
   |                              |              +----------------------------------------------+
   |                              |              | ``libsundials_fnvecserial.a``                |
   |                              +--------------+----------------------------------------------+
   |                              | Headers      | ``nvector/nvector_serial.h``                 |
   +------------------------------+--------------+----------------------------------------------+
   | PARALLEL                     | Libraries    | ``libsundials_nvecparallel.LIB``             |
   |                              |              +----------------------------------------------+
   |                              |              | ``libsundials_fnvecparallel.a``              |
   |                              +--------------+----------------------------------------------+
   |                              | Headers      | ``nvector/nvector_parallel.h``               |
   +------------------------------+--------------+----------------------------------------------+
   | OPENMP                       | Libraries    | ``libsundials_nvecopenmp.LIB``               |
   |                              |              +----------------------------------------------+
   |                              |              | ``libsundials_fnvecopenmp.a``                |
   |                              +--------------+----------------------------------------------+
   |                              | Headers      | ``nvector/nvector_openmp.h``                 |
   +------------------------------+--------------+----------------------------------------------+
   | PTHREADS                     | Libraries    | ``libsundials_nvecpthreads.LIB``             |
   |                              |              +----------------------------------------------+
   |                              |              | ``libsundials_fnvecpthreads.a``              |
   |                              +--------------+----------------------------------------------+
   |                              | Headers      | ``nvector/nvector_pthreads.h``               |
   +------------------------------+--------------+----------------------------------------------+
   | PARHYP                       | Libraries    | ``libsundials_nvecparhyp.LIB``               |
   |                              +--------------+----------------------------------------------+
   |                              | Headers      | ``nvector/nvector_parhyp.h``                 |
   +------------------------------+--------------+----------------------------------------------+
   | PETSC                        | Libraries    | ``libsundials_nvecpetsc.LIB``                |
   |                              +--------------+----------------------------------------------+
   |                              | Headers      | ``nvector/nvector_petsc.h``                  |
   +------------------------------+--------------+----------------------------------------------+
   | CUDA                         | Libraries    | ``libsundials_nveccuda.LIB``                 |
   |                              +--------------+----------------------------------------------+
   |                              | Headers      | ``nvector/nvector_cuda.h``                   |
   +------------------------------+--------------+----------------------------------------------+
   | HIP                          | Libraries    | ``libsundials_nvechip.LIB``                  |
   |                              +--------------+----------------------------------------------+
   |                              | Headers      | ``nvector/nvector_hip.h``                    |
   +------------------------------+--------------+----------------------------------------------+
   | RAJA                         | Libraries    | ``libsundials_nveccudaraja.LIB``             |
   |                              |              | ``libsundials_nvechipraja.LIB``              |
   |                              +--------------+----------------------------------------------+
   |                              | Headers      | ``nvector/nvector_raja.h``                   |
   +------------------------------+--------------+----------------------------------------------+
   | SYCL                         | Libraries    | ``libsundials_nvecsycl.LIB``                 |
   |                              +--------------+----------------------------------------------+
   |                              | Headers      | ``nvector/nvector_sycl.h``                   |
   +------------------------------+--------------+----------------------------------------------+
   | MANYVECTOR                   | Libraries    | ``libsundials_nvecmanyvector.LIB``           |
   |                              +--------------+----------------------------------------------+
   |                              | Headers      | ``nvector/nvector_manyvector.h``             |
   +------------------------------+--------------+----------------------------------------------+
   | MPIMANYVECTOR                | Libraries    | ``libsundials_nvecmpimanyvector.LIB``        |
   |                              +--------------+----------------------------------------------+
   |                              | Headers      | ``nvector/nvector_mpimanyvector.h``          |
   +------------------------------+--------------+----------------------------------------------+
   | MPIPLUSX                     | Libraries    | ``libsundials_nvecmpiplusx.LIB``             |
   |                              +--------------+----------------------------------------------+
   |                              | Headers      | ``nvector/nvector_mpiplusx.h``               |
   +------------------------------+--------------+----------------------------------------------+
   |                                                                                            |
   | SUNMATRIX Modules                                                                          |
   |                                                                                            |
   +------------------------------+--------------+----------------------------------------------+
   | BAND                         | Libraries    | ``libsundials_sunmatrixband.LIB``            |
   |                              |              +----------------------------------------------+
   |                              |              | ``libsundials_fsunmatrixband.a``             |
   |                              +--------------+----------------------------------------------+
   |                              | Headers      | ``sunmatrix/sunmatrix_band.h``               |
   +------------------------------+--------------+----------------------------------------------+
   | DENSE                        | Libraries    | ``libsundials_sunmatrixdense.LIB``           |
   |                              |              +----------------------------------------------+
   |                              |              | ``libsundials_fsunmatrixdense.a``            |
   |                              +--------------+----------------------------------------------+
   |                              | Headers      | ``sunmatrix/sunmatrix_dense.h``              |
   +------------------------------+--------------+----------------------------------------------+
   | SPARSE                       | Libraries    | ``libsundials_sunmatrixsparse.LIB``          |
   |                              |              +----------------------------------------------+
   |                              |              | ``libsundials_fsunmatrixsparse.a``           |
   |                              +--------------+----------------------------------------------+
   |                              | Headers      | ``sunmatrix/sunmatrix_sparse.h``             |
   +------------------------------+--------------+----------------------------------------------+
   | SLUNRLOC                     | Libraries    | ``libsundials_sunmatrixslunrloc.LIB``        |
   |                              +--------------+----------------------------------------------+
   |                              | Headers      | ``sunmatrix/sunmatrix_slunrloc.h``           |
   +------------------------------+--------------+----------------------------------------------+
   | CUSPARSE                     | Libraries    | ``libsundials_sunmatrixcusparse.LIB``        |
   |                              +--------------+----------------------------------------------+
   |                              | Headers      | ``sunmatrix/sunmatrix_cusparse.h``           |
   +------------------------------+--------------+----------------------------------------------+
   |                                                                                            |
   | SUNLINSOL Modules                                                                          |
   |                                                                                            |
   +------------------------------+--------------+----------------------------------------------+
   | BAND                         | Libraries    | ``libsundials_sunlinsolband.LIB``            |
   |                              |              +----------------------------------------------+
   |                              |              | ``libsundials_fsunlinsolband.a``             |
   |                              +--------------+----------------------------------------------+
   |                              | Headers      | ``sunlinsol/sunlinsol_band.h``               |
   +------------------------------+--------------+----------------------------------------------+
   | DENSE                        | Libraries    | ``libsundials_sunlinsoldense.LIB``           |
   |                              |              +----------------------------------------------+
   |                              |              | ``libsundials_fsunlinsoldense.a``            |
   |                              +--------------+----------------------------------------------+
   |                              | Headers      | ``sunlinsol/sunlinsol_dense.h``              |
   +------------------------------+--------------+----------------------------------------------+
   | KLU                          | Libraries    | ``libsundials_sunlinsolklu.LIB``             |
   |                              |              +----------------------------------------------+
   |                              |              | ``libsundials_fsunlinsolklu.a``              |
   |                              +--------------+----------------------------------------------+
   |                              | Headers      | ``sunlinsol/sunlinsol_klu.h``                |
   +------------------------------+--------------+----------------------------------------------+
   | LAPACKBAND                   | Libraries    | ``libsundials_sunlinsollapackband.LIB``      |
   |                              |              +----------------------------------------------+
   |                              |              | ``libsundials_fsunlinsollapackband.a``       |
   |                              +--------------+----------------------------------------------+
   |                              | Headers      | ``sunlinsol/sunlinsol_lapackband.h``         |
   +------------------------------+--------------+----------------------------------------------+
   | LAPACKDENSE                  | Libraries    | ``libsundials_sunlinsollapackdense.LIB``     |
   |                              |              +----------------------------------------------+
   |                              |              | ``libsundials_fsunlinsollapackdense.a``      |
   |                              +--------------+----------------------------------------------+
   |                              | Headers      | ``sunlinsol/sunlinsol_lapackdense.h``        |
   +------------------------------+--------------+----------------------------------------------+
   | PCG                          | Libraries    | ``libsundials_sunlinsolpcg.LIB``             |
   |                              |              +----------------------------------------------+
   |                              |              | ``libsundials_fsunlinsolpcg.a``              |
   |                              +--------------+----------------------------------------------+
   |                              | Headers      | ``sunlinsol/sunlinsol_pcg.h``                |
   +------------------------------+--------------+----------------------------------------------+
   | SPBCGS                       | Libraries    | ``libsundials_sunlinsolspbcgs.LIB``          |
   |                              |              +----------------------------------------------+
   |                              |              | ``libsundials_fsunlinsolspbcgs.a``           |
   |                              +--------------+----------------------------------------------+
   |                              | Headers      | ``sunlinsol/sunlinsol_spbcgs.h``             |
   +------------------------------+--------------+----------------------------------------------+
   | SPFGMR                       | Libraries    | ``libsundials_sunlinsolspfgmr.LIB``          |
   |                              |              +----------------------------------------------+
   |                              |              | ``libsundials_fsunlinsolspfgmr.a``           |
   |                              +--------------+----------------------------------------------+
   |                              | Headers      | ``sunlinsol/sunlinsol_spfgmr.h``             |
   +------------------------------+--------------+----------------------------------------------+
   | SPGMR                        | Libraries    | ``libsundials_sunlinsolspgmr.LIB``           |
   |                              |              +----------------------------------------------+
   |                              |              | ``libsundials_fsunlinsolspgmr.a``            |
   |                              +--------------+----------------------------------------------+
   |                              | Headers      | ``sunlinsol/sunlinsol_spgmr.h``              |
   +------------------------------+--------------+----------------------------------------------+
   | SPTFQMR                      | Libraries    | ``libsundials_sunlinsolsptfqmr.LIB``         |
   |                              |              +----------------------------------------------+
   |                              |              | ``libsundials_fsunlinsolsptfqmr.a``          |
   |                              +--------------+----------------------------------------------+
   |                              | Headers      | ``sunlinsol/sunlinsol_sptfqmr.h``            |
   +------------------------------+--------------+----------------------------------------------+
   | SUPERLUMT                    | Libraries    | ``libsundials_sunlinsolsuperlumt.LIB``       |
   |                              |              +----------------------------------------------+
   |                              |              | ``libsundials_fsunlinsolsuperlumt.a``        |
   |                              +--------------+----------------------------------------------+
   |                              | Headers      | ``sunlinsol/sunlinsol_superlumt.h``          |
   +------------------------------+--------------+----------------------------------------------+
   | SUPERLUDIST                  | Libraries    | ``libsundials_sunlinsolsuperludist.LIB``     |
   |                              +--------------+----------------------------------------------+
   |                              | Headers      | ``sunlinsol/sunlinsol_superludist.h``        |
   +------------------------------+--------------+----------------------------------------------+
   | CUSOLVERSP_BATCHQR           | Libraries    | ``libsundials_sunlinsolcusolversp.LIB``      |
   |                              +--------------+----------------------------------------------+
   |                              | Headers      | ``sunlinsol/sunlinsol_cusolversp_batchqr.h`` |
   +------------------------------+--------------+----------------------------------------------+
   |                                                                                            |
   | SUNNONLINSOL Modules                                                                       |
   |                                                                                            |
   +------------------------------+--------------+----------------------------------------------+
   | NEWTON                       | Libraries    | ``libsundials_sunnonlinsolnewton.LIB``       |
   |                              |              +----------------------------------------------+
   |                              |              | ``libsundials_fsunnonlinsolnewton.a``        |
   |                              +--------------+----------------------------------------------+
   |                              | Headers      | ``sunnonlinsol/sunnonlinsol_newton.h``       |
   +------------------------------+--------------+----------------------------------------------+
   | FIXEDPOINT                   | Libraries    | ``libsundials_sunnonlinsolfixedpoint.LIB``   |
   |                              |              +----------------------------------------------+
   |                              |              | ``libsundials_fsunnonlinsolfixedpoint.a``    |
   |                              +--------------+----------------------------------------------+
   |                              | Headers      | ``sunnonlinsol/sunnonlinsol_fixedpoint.h``   |
   +------------------------------+--------------+----------------------------------------------+
   | PETSCSNES                    | Libraries    | ``libsundials_sunnonlinsolpetscsnes.LIB``    |
   |                              +--------------+----------------------------------------------+
   |                              | Headers      | ``sunnonlinsol/sunnonlinsol_petscsnes.h``    |
   +------------------------------+--------------+----------------------------------------------+
   |                                                                                            |
   | SUNDIALS Packages                                                                          |
   |                                                                                            |
   +------------------------------+--------------+----------------------------------------------+
   | CVODE                        | Libraries    | ``libsundials_cvode.LIB``                    |
   |                              |              +----------------------------------------------+
   |                              |              | ``libsundials_fcvode.a``                     |
   |                              +--------------+----------------------------------------------+
   |                              | Headers      | ``cvode/cvode.h``                            |
   |                              |              +----------------------------------------------+
   |                              |              | ``cvode/cvode_bandpre.h``                    |
   |                              |              +----------------------------------------------+
   |                              |              | ``cvode/cvode_bbdpre.h``                     |
   |                              |              +----------------------------------------------+
   |                              |              | ``cvode/cvode_diag.h``                       |
   |                              |              +----------------------------------------------+
   |                              |              | ``cvode/cvode_direct.h``                     |
   |                              |              +----------------------------------------------+
   |                              |              | ``cvode/cvode_impl.h``                       |
   |                              |              +----------------------------------------------+
   |                              |              | ``cvode/cvode_ls.h``                         |
   |                              |              +----------------------------------------------+
   |                              |              | ``cvode/cvode_spils.h``                      |
   +------------------------------+--------------+----------------------------------------------+
   | CVODES                       | Libraries    | ``libsundials_cvodes.LIB``                   |
   |                              +--------------+----------------------------------------------+
   |                              | Headers      | ``cvodes/cvodes.h``                          |
   |                              |              +----------------------------------------------+
   |                              |              | ``cvodes/cvodes_bandpre.h``                  |
   |                              |              +----------------------------------------------+
   |                              |              | ``cvodes/cvodes_bbdpre.h``                   |
   |                              |              +----------------------------------------------+
   |                              |              | ``cvodes/cvodes_diag.h``                     |
   |                              |              +----------------------------------------------+
   |                              |              | ``cvodes/cvodes_direct.h``                   |
   |                              |              +----------------------------------------------+
   |                              |              | ``cvodes/cvodes_impl.h``                     |
   |                              |              +----------------------------------------------+
   |                              |              | ``cvodes/cvodes_spils.h``                    |
   +------------------------------+--------------+----------------------------------------------+
   | ARKODE                       | Libraries    | ``libsundials_arkode.LIB``                   |
   |                              |              +----------------------------------------------+
   |                              |              | ``libsundials_farkode.a``                    |
   |                              |              +----------------------------------------------+
   |                              |              | ``libsundials_xbraid.LIB``                   |
   |                              +--------------+----------------------------------------------+
   |                              | Headers      | ``arkode/arkode.h``                          |
   |                              |              +----------------------------------------------+
   |                              |              | ``arkode/arkode_arkstep.h``                  |
   |                              |              +----------------------------------------------+
   |                              |              | ``arkode/arkode_bandpre.h``                  |
   |                              |              +----------------------------------------------+
   |                              |              | ``arkode/arkode_bbdpre.h``                   |
   |                              |              +----------------------------------------------+
   |                              |              | ``arkode/arkode_butcher.h``                  |
   |                              |              +----------------------------------------------+
   |                              |              | ``arkode/arkode_butcher_dirk.h``             |
   |                              |              +----------------------------------------------+
   |                              |              | ``arkode/arkode_butcher_erk.h``              |
   |                              |              +----------------------------------------------+
   |                              |              | ``arkode/arkode_erkstep.h``                  |
   |                              |              +----------------------------------------------+
   |                              |              | ``arkode/arkode_impl.h``                     |
   |                              |              +----------------------------------------------+
   |                              |              | ``arkode/arkode_ls.h``                       |
   |                              |              +----------------------------------------------+
   |                              |              | ``arkode/arkode_xbraid.h``                   |
   +------------------------------+--------------+----------------------------------------------+
   | IDA                          | Libraries    | ``libsundials_ida.LIB``                      |
   |                              |              +----------------------------------------------+
   |                              |              | ``libsundials_fida.a``                       |
   |                              +--------------+----------------------------------------------+
   |                              | Headers      | ``ida/ida.h``                                |
   |                              |              +----------------------------------------------+
   |                              |              | ``ida/ida_bbdpre.h``                         |
   |                              |              +----------------------------------------------+
   |                              |              | ``ida/ida_direct.h``                         |
   |                              |              +----------------------------------------------+
   |                              |              | ``ida/ida_impl.h``                           |
   |                              |              +----------------------------------------------+
   |                              |              | ``ida/ida_ls.h``                             |
   |                              |              +----------------------------------------------+
   |                              |              | ``ida/ida_spils.h``                          |
   +------------------------------+--------------+----------------------------------------------+
   | IDAS                         | Libraries    | ``libsundials_idas.LIB``                     |
   |                              +--------------+----------------------------------------------+
   |                              | Headers      | ``idas/idas.h``                              |
   |                              |              +----------------------------------------------+
   |                              |              | ``idas/idas_bbdpre.h``                       |
   |                              |              +----------------------------------------------+
   |                              |              | ``idas/idas_direct.h``                       |
   |                              |              +----------------------------------------------+
   |                              |              | ``idas/idas_impl.h``                         |
   |                              |              +----------------------------------------------+
   |                              |              | ``idas/idas_spils.h``                        |
   +------------------------------+--------------+----------------------------------------------+
   | KINSOL                       | Libraries    | ``libsundials_kinsol.LIB``                   |
   |                              |              +----------------------------------------------+
   |                              |              | ``libsundials_fkinsol.a``                    |
   |                              +--------------+----------------------------------------------+
   |                              | Headers      | ``kinsol/kinsol.h``                          |
   |                              |              +----------------------------------------------+
   |                              |              | ``kinsol/kinsol_bbdpre.h``                   |
   |                              |              +----------------------------------------------+
   |                              |              | ``kinsol/kinsol_direct.h``                   |
   |                              |              +----------------------------------------------+
   |                              |              | ``kinsol/kinsol_impl.h``                     |
   |                              |              +----------------------------------------------+
   |                              |              | ``kinsol/kinsol_ls.h``                       |
   |                              |              +----------------------------------------------+
   |                              |              | ``kinsol/kinsol_spils.h``                    |
   +------------------------------+--------------+----------------------------------------------+
