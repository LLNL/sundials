// #ifndef _SUNDIALS_TYPES_H
// 
// #ifdef __cplusplus 
// #endif
// 
// #ifdef SWIG
// 
// #else
// 
// #endif
// 


auto pyEnumSUNOutputFormat_ =
    nb::enum_<SUNOutputFormat_>(m, "SUNOutputFormat_", nb::is_arithmetic(), "")
        .value("SUN_OUTPUTFORMAT_TABLE", SUN_OUTPUTFORMAT_TABLE, "")
        .value("SUN_OUTPUTFORMAT_CSV", SUN_OUTPUTFORMAT_CSV, "");
// #ifndef SWIG
// 
m.attr("SUN_COMM_NULL") = 0;
// #endif
// 
// #ifdef __cplusplus
// 
// #endif
// 


auto pyEnumSUNDataIOMode_ =
    nb::enum_<SUNDataIOMode_>(m, "SUNDataIOMode_", nb::is_arithmetic(), "")
        .value("SUN_DATAIOMODE_INMEM", SUN_DATAIOMODE_INMEM, "");
// #endif 
