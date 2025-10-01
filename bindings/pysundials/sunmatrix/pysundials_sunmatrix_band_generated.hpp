// #ifndef _SUNMATRIX_BAND_H
//
// #ifdef __cplusplus
// #endif
//

auto pyClass_SUNMatrixContent_Band =
  nb::class_<_SUNMatrixContent_Band>(m, "_SUNMatrixContent_Band", "")
    .def(nb::init<>()) // implicit default constructor
  ;

m.def("SUNBandMatrix", SUNBandMatrix, nb::arg("N"), nb::arg("mu"),
      nb::arg("ml"), nb::arg("sunctx"), nb::rv_policy::reference);

m.def("SUNBandMatrixStorage", SUNBandMatrixStorage, nb::arg("N"), nb::arg("mu"),
      nb::arg("ml"), nb::arg("smu"), nb::arg("sunctx"), nb::rv_policy::reference);

m.def("SUNBandMatrix_Print", SUNBandMatrix_Print, nb::arg("A"),
      nb::arg("outfile"));

m.def("SUNBandMatrix_Rows", SUNBandMatrix_Rows, nb::arg("A"));

m.def("SUNBandMatrix_Columns", SUNBandMatrix_Columns, nb::arg("A"));

m.def("SUNBandMatrix_LowerBandwidth", SUNBandMatrix_LowerBandwidth, nb::arg("A"));

m.def("SUNBandMatrix_UpperBandwidth", SUNBandMatrix_UpperBandwidth, nb::arg("A"));

m.def("SUNBandMatrix_StoredUpperBandwidth", SUNBandMatrix_StoredUpperBandwidth,
      nb::arg("A"));

m.def("SUNBandMatrix_LDim", SUNBandMatrix_LDim, nb::arg("A"));

m.def("SUNBandMatrix_LData", SUNBandMatrix_LData, nb::arg("A"));
// #ifdef __cplusplus
//
// #endif
//
// #endif
//
