.. For package-specific references use :ref: rather than :numref: so intersphinx
   links to the appropriate place on read the docs

**Major Features**

**New Features and Enhancements**

The functions :c:func:`IDAGetUserDataB` and :c:func:`CVodeGetUserDataB` were added.

**Bug Fixes**

**Deprecation Notices**

:c:func:`SUNDIALSFileOpen` and :c:func:`SUNDIALSFileClose` will be removed in the next major release. 
Use :c:func:`SUNFileOpen` and :c:func:`SUNFileClose` instead.

The ``Convert`` methods on the ``sundials::kokkos:Vector``, ``sundials::kokkos::DenseMatrix``,
``sundials::ginkgo::Matrix``, ``sundials::ginkgo::BatchMatrix``, ``sundials::kokkos::DenseLinearSolver``,
``sundials::ginkgo::LinearSolver``, and ``sundials::ginkgo::BatchLinearSolver`` classes have
been deprecated and will be removed in the next major release. The method ``get``, should
be used instead.