// #ifndef _SUNMATRIX_DENSE_H
//
// #ifdef __cplusplus
// #endif
//

auto pyClass_SUNMatrixContent_Dense =
  nb::class_<_SUNMatrixContent_Dense>(m, "_SUNMatrixContent_Dense", "")
    .def(nb::init<>()) // implicit default constructor
  ;

m.def("SUNDenseMatrix", SUNDenseMatrix, nb::arg("M"), nb::arg("N"),
      nb::arg("sunctx"), nb::rv_policy::reference);

m.def("SUNDenseMatrix_Print", SUNDenseMatrix_Print, nb::arg("A"),
      nb::arg("outfile"));

m.def("SUNDenseMatrix_Rows", SUNDenseMatrix_Rows, nb::arg("A"));

m.def("SUNDenseMatrix_Columns", SUNDenseMatrix_Columns, nb::arg("A"));

m.def("SUNDenseMatrix_LData", SUNDenseMatrix_LData, nb::arg("A"));
// #ifdef __cplusplus
//
// #endif
//
// #endif
//
