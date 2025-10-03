// #ifndef _SUNMATRIX_SPARSE_H
//
// #ifdef __cplusplus
// #endif
//
m.attr("SUN_CSC_MAT") = 0;
m.attr("SUN_CSR_MAT") = 1;

auto pyClass_SUNMatrixContent_Sparse =
  nb::class_<_SUNMatrixContent_Sparse>(m, "_SUNMatrixContent_Sparse", "")
    .def(nb::init<>()) // implicit default constructor
  ;

m.def("SUNSparseMatrix", SUNSparseMatrix, nb::arg("M"), nb::arg("N"),
      nb::arg("NNZ"), nb::arg("sparsetype"), nb::arg("sunctx"),
      nb::rv_policy::reference);

m.def("SUNSparseFromDenseMatrix", SUNSparseFromDenseMatrix, nb::arg("A"),
      nb::arg("droptol"), nb::arg("sparsetype"), nb::rv_policy::reference);

m.def("SUNSparseFromBandMatrix", SUNSparseFromBandMatrix, nb::arg("A"),
      nb::arg("droptol"), nb::arg("sparsetype"), nb::rv_policy::reference);

m.def("SUNSparseMatrix_Realloc", SUNSparseMatrix_Realloc, nb::arg("A"));

m.def("SUNSparseMatrix_Reallocate", SUNSparseMatrix_Reallocate, nb::arg("A"),
      nb::arg("NNZ"));

m.def("SUNSparseMatrix_Print", SUNSparseMatrix_Print, nb::arg("A"),
      nb::arg("outfile"));

m.def("SUNSparseMatrix_Rows", SUNSparseMatrix_Rows, nb::arg("A"));

m.def("SUNSparseMatrix_Columns", SUNSparseMatrix_Columns, nb::arg("A"));

m.def("SUNSparseMatrix_NNZ", SUNSparseMatrix_NNZ, nb::arg("A"));

m.def("SUNSparseMatrix_NP", SUNSparseMatrix_NP, nb::arg("A"));

m.def("SUNSparseMatrix_SparseType", SUNSparseMatrix_SparseType, nb::arg("A"));
// #ifdef __cplusplus
//
// #endif
//
// #endif
//
