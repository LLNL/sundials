.. For package-specific references use :ref: rather than :numref: so intersphinx
   links to the appropriate place on read the docs

**Major Features**

**New Features and Enhancements**

The functions :c:func:`KINSetMAA` and :c:func:`KINSetOrthAA` have been updated
to allow for setting the Anderson acceleration depth and orthogonalization
method after :c:func:`KINInit`. Additionally, :c:func:`KINSetMAA` and
:c:func:`KINSetNumMaxIters` may now be called in any order.

Added the :c:type:`SUNDomEigEstimator` interface for estimating the dominant eigenvalue
value of a system. Two implementations are provided: Power Iteration and Arnoldi 
Iteration. The latter method requires building with LAPACK support enabled.

**Bug Fixes**

The shared library version numbers for the oneMKL dense linear solver and
matrix as well as the PETSc SNES nonlinear solver have been corrected.

Fixed a CMake bug where the MRI H-Tol controller was not included in the ARKODE
Fortran module.

Fixed a bug in the CUDA and HIP implementations of
:c:func:`SUNMemoryHelper_CopyAsync` where the execution stream is not extracted
correctly from the helper when a stream is not provided to
:c:func:`SUNMemoryHelper_CopyAsync`.

Fixed a bug in MRIStep where a segfault would occur when an MRI coupling table
is not explicitly set and an MRI integrator is nested inside another MRI
integrator.

Fixed a bug in MRIStep where MERK methods with unordered stage groups (MERK43
and MERK54) would include stage right-hand side vectors that had not been
computed yet in fast time scale forcing computations. These vectors were scaled
by zero, so in most cases the extraneous computations would not impact results.
However, in cases where these vectors contain ``inf`` or ``nan``, this would
lead to erroneous forcing terms.

Fixed a bug in the ``suntools.logs`` Python module where the ``get_history``
function, when given a ``step_status`` for filtering output from a multirate
method, would only extract values from the fast time scale if the slow time
scale step matched the given status filter. Fixed an additional bug in
``get_history`` with MRI-GARK methods where values would not be extracted from a
fast time scale integration associated with an embedding.

**Deprecation Notices**
