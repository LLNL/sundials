// #ifndef _SUNMATRIX_H
// 
// #ifdef __cplusplus 
// #endif
// 


auto pyClass_generic_SUNMatrix_Ops =
    nb::class_<_generic_SUNMatrix_Ops>
        (m, "_generic_SUNMatrix_Ops", "Structure containing function pointers to matrix operations")
    .def(nb::init<>()) // implicit default constructor 
    ;


auto pyClass_generic_SUNMatrix =
    nb::class_<_generic_SUNMatrix>
        (m, "_generic_SUNMatrix", " A matrix is a structure with an implementation-dependent\n   'content' field, and a pointer to a structure of matrix\n   operations corresponding to that implementation.")
    .def("__init__", [](_generic_SUNMatrix * self, SUNMatrix_Ops ops = SUNMatrix_Ops(), SUNContext sunctx = SUNContext())
    {
        new (self) _generic_SUNMatrix();  // placement new
        auto r = self;
        r->ops = ops;
        r->sunctx = sunctx;
    },
    nb::arg("ops") = SUNMatrix_Ops(), nb::arg("sunctx") = SUNContext()
    )
    .def_rw("content", &_generic_SUNMatrix::content, "")
    .def_rw("ops", &_generic_SUNMatrix::ops, "")
    .def_rw("sunctx", &_generic_SUNMatrix::sunctx, "")
    ;


m.def("SUNMatNewEmpty",
    SUNMatNewEmpty, nb::arg("sunctx"));

m.def("SUNMatCopyOps",
    SUNMatCopyOps, nb::arg("A"), nb::arg("B"));

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
