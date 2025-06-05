// #ifndef _SUNDIALS_ITERATIVE_H
// 
// #ifdef __cplusplus 
// #endif
// 


auto pyEnumSUN_PREC_ID =
    nb::enum_<SUN_PREC_ID>(m, "SUN_PREC_ID", nb::is_arithmetic(), "")
        .value("SUN_PREC_NONE", SUN_PREC_NONE, "")
        .value("SUN_PREC_LEFT", SUN_PREC_LEFT, "")
        .value("SUN_PREC_RIGHT", SUN_PREC_RIGHT, "")
        .value("SUN_PREC_BOTH", SUN_PREC_BOTH, "");


auto pyEnumSUN_GRAMSCHMIDT_ID =
    nb::enum_<SUN_GRAMSCHMIDT_ID>(m, "SUN_GRAMSCHMIDT_ID", nb::is_arithmetic(), "")
        .value("SUN_MODIFIED_GS", SUN_MODIFIED_GS, "")
        .value("SUN_CLASSICAL_GS", SUN_CLASSICAL_GS, "");


m.def("SUNQRAdd_MGS",
    SUNQRAdd_MGS, nb::arg("Q"), nb::arg("R"), nb::arg("df"), nb::arg("m"), nb::arg("mMax"), nb::arg("QRdata"));

m.def("SUNQRAdd_ICWY",
    SUNQRAdd_ICWY, nb::arg("Q"), nb::arg("R"), nb::arg("df"), nb::arg("m"), nb::arg("mMax"), nb::arg("QRdata"));

m.def("SUNQRAdd_ICWY_SB",
    SUNQRAdd_ICWY_SB, nb::arg("Q"), nb::arg("R"), nb::arg("df"), nb::arg("m"), nb::arg("mMax"), nb::arg("QRdata"));

m.def("SUNQRAdd_CGS2",
    SUNQRAdd_CGS2, nb::arg("Q"), nb::arg("R"), nb::arg("df"), nb::arg("m"), nb::arg("mMax"), nb::arg("QRdata"));

m.def("SUNQRAdd_DCGS2",
    SUNQRAdd_DCGS2, nb::arg("Q"), nb::arg("R"), nb::arg("df"), nb::arg("m"), nb::arg("mMax"), nb::arg("QRdata"));

m.def("SUNQRAdd_DCGS2_SB",
    SUNQRAdd_DCGS2_SB, nb::arg("Q"), nb::arg("R"), nb::arg("df"), nb::arg("m"), nb::arg("mMax"), nb::arg("QRdata"));
// #ifdef __cplusplus
// 
// #endif
// 
// #endif
// 
