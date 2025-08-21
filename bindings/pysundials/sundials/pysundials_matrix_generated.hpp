// #ifndef _SUNMATRIX_H
// 
// #ifdef __cplusplus 
// #endif
// 


auto pyClass_generic_SUNMatrix_Ops =
    nb::class_<_generic_SUNMatrix_Ops>
        (m, "_generic_SUNMatrix_Ops", "")
    .def(nb::init<>()) // implicit default constructor
    ;


auto pyClass_generic_SUNMatrix =
    nb::class_<_generic_SUNMatrix>
        (m, "_generic_SUNMatrix", "")
    .def(nb::init<>()) // implicit default constructor
    ;


m.def("SUNMatGetID",
    SUNMatGetID, nb::arg("A"));

m.def("SUNMatClone",
    SUNMatClone, nb::arg("A"));

m.def("SUNMatZero",
    SUNMatZero, nb::arg("A"));

m.def("SUNMatCopy",
    SUNMatCopy, nb::arg("A"), nb::arg("B"));

m.def("SUNMatScaleAdd",
    SUNMatScaleAdd, nb::arg("c"), nb::arg("A"), nb::arg("B"));

m.def("SUNMatScaleAddI",
    SUNMatScaleAddI, nb::arg("c"), nb::arg("A"));

m.def("SUNMatMatvecSetup",
    SUNMatMatvecSetup, nb::arg("A"));

m.def("SUNMatMatvec",
    SUNMatMatvec, nb::arg("A"), nb::arg("x"), nb::arg("y"));

m.def("SUNMatHermitianTransposeVec",
    SUNMatHermitianTransposeVec, nb::arg("A"), nb::arg("x"), nb::arg("y"));
// #ifdef __cplusplus
// 
// #endif
// 
// #endif 
