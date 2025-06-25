.. For package-specific references use :ref: rather than :numref: so intersphinx
   links to the appropriate place on read the docs

**Major Features**

**New Features and Enhancements**

A new SUNLinearSolver module, SUNLINEARSOLVER_GINKGOBLOCK, for block/batched systems is included.
A corresponding SUNMatrix module, SUNMATRIX_GINKGOBLOCK, was also added.
These new modules use the `Ginkgo linear solver library <https://ginkgo-project.github.io/>`__.

**Bug Fixes**

The shared library version numbers for the oneMKL dense linear solver and
matrix as well as the PETSc SNES nonlinear solver have been corrected.

**Deprecation Notices**
