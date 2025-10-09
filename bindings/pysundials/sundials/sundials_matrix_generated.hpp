// #ifndef _SUNMATRIX_H
//
// #ifdef __cplusplus
// #endif
//

auto pyEnumSUNMatrix_ID =
  nb::enum_<SUNMatrix_ID>(m, "SUNMatrix_ID", nb::is_arithmetic(), "")
    .value("SUNMATRIX_DENSE", SUNMATRIX_DENSE, "")
    .value("SUNMATRIX_MAGMADENSE", SUNMATRIX_MAGMADENSE, "")
    .value("SUNMATRIX_ONEMKLDENSE", SUNMATRIX_ONEMKLDENSE, "")
    .value("SUNMATRIX_BAND", SUNMATRIX_BAND, "")
    .value("SUNMATRIX_SPARSE", SUNMATRIX_SPARSE, "")
    .value("SUNMATRIX_SLUNRLOC", SUNMATRIX_SLUNRLOC, "")
    .value("SUNMATRIX_CUSPARSE", SUNMATRIX_CUSPARSE, "")
    .value("SUNMATRIX_GINKGO", SUNMATRIX_GINKGO, "")
    .value("SUNMATRIX_GINKGOBATCH", SUNMATRIX_GINKGOBATCH, "")
    .value("SUNMATRIX_KOKKOSDENSE", SUNMATRIX_KOKKOSDENSE, "")
    .value("SUNMATRIX_CUSTOM", SUNMATRIX_CUSTOM, "")
    .export_values();
// #ifndef SWIG
//
// #endif
//

auto pyClass_generic_SUNMatrix_Ops =
  nb::class_<_generic_SUNMatrix_Ops>(m, "_generic_SUNMatrix_Ops", "")
    .def(nb::init<>()) // implicit default constructor
  ;

auto pyClass_generic_SUNMatrix =
  nb::class_<_generic_SUNMatrix>(m, "_generic_SUNMatrix", "")
    .def(nb::init<>()) // implicit default constructor
  ;

m.def("SUNMatGetID", SUNMatGetID, nb::arg("A"));

m.def("SUNMatClone", SUNMatClone, nb::arg("A"), nb::rv_policy::reference);

m.def("SUNMatZero", SUNMatZero, nb::arg("A"));

m.def("SUNMatCopy", SUNMatCopy, nb::arg("A"), nb::arg("B"));

m.def("SUNMatScaleAdd", SUNMatScaleAdd, nb::arg("c"), nb::arg("A"), nb::arg("B"));

m.def("SUNMatScaleAddI", SUNMatScaleAddI, nb::arg("c"), nb::arg("A"));

m.def("SUNMatMatvecSetup", SUNMatMatvecSetup, nb::arg("A"));

m.def("SUNMatMatvec", SUNMatMatvec, nb::arg("A"), nb::arg("x"), nb::arg("y"));

m.def("SUNMatHermitianTransposeVec", SUNMatHermitianTransposeVec, nb::arg("A"),
      nb::arg("x"), nb::arg("y"));
// #ifdef __cplusplus
//
// #endif
//
// #endif
