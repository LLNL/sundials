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
// #ifdef __cplusplus
// 
// #endif
// 
// #endif
// 
