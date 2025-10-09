.. For package-specific references use :ref: rather than :numref: so intersphinx
   links to the appropriate place on read the docs

**Major Features**

**New Features and Enhancements**

The functions ``CVodeGetUserDataB`` and ``IDAGetUserDataB`` were added to CVODES
and IDAS, respectively.

**Bug Fixes**

The SPRKStep module now accounts for zero coefficients in the SPRK tables, eliminating
extraneous function evaluations.

A bug preventing a user supplied :c:func:`SUNStepper_ResetCheckpointIndex` function from being
called was fixed. 

**Deprecation Notices**

:c:func:`SUNDIALSFileOpen` and :c:func:`SUNDIALSFileClose` will be removed in the next major release. 
Use :c:func:`SUNFileOpen` and :c:func:`SUNFileClose` instead.

The ``Convert`` methods on the ``sundials::kokkos:Vector``, ``sundials::kokkos::DenseMatrix``,
``sundials::ginkgo::Matrix``, ``sundials::ginkgo::BatchMatrix``, ``sundials::kokkos::DenseLinearSolver``,
``sundials::ginkgo::LinearSolver``, and ``sundials::ginkgo::BatchLinearSolver`` classes have
been deprecated and will be removed in the next major release. The method ``get``, should
be used instead.

The ``CSC_MAT`` and ``CSR_MAT`` macros defined in ``sunmatrix_sparse.h`` will be removed in
the next major release. Use ``SUN_CSC_MAT`` and ``SUN_CSR_MAT`` instead.

The ``N_Vector_S`` typedef to ``N_Vector*`` is deprecated and will be removed in the next major release.
