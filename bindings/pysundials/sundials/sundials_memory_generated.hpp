// #ifndef _SUNDIALS_MEMORY_H
//
// #ifdef __cplusplus
// #endif
//

auto pyEnumSUNMemoryType = nb::enum_<SUNMemoryType>(m, "SUNMemoryType",
                                                    nb::is_arithmetic(), "")
                             .value("SUNMEMTYPE_HOST", SUNMEMTYPE_HOST, "")
                             .value("SUNMEMTYPE_PINNED", SUNMEMTYPE_PINNED, "")
                             .value("SUNMEMTYPE_DEVICE", SUNMEMTYPE_DEVICE, "")
                             .value("SUNMEMTYPE_UVM", SUNMEMTYPE_UVM, "")
                             .export_values();
// #ifndef SWIG
//
// #endif
//

auto pyClassSUNMemoryHelper_ =
  nb::class_<SUNMemoryHelper_>(m, "SUNMemoryHelper_", "")
    .def(nb::init<>()) // implicit default constructor
  ;

auto pyClassSUNMemoryHelper_Ops_ =
  nb::class_<SUNMemoryHelper_Ops_>(m, "SUNMemoryHelper_Ops_", "")
    .def(nb::init<>()) // implicit default constructor
  ;

m.def("SUNMemoryHelper_Clone", SUNMemoryHelper_Clone, nb::arg("param_0"),
      nb::rv_policy::reference);

m.def("SUNMemoryHelper_SetDefaultQueue", SUNMemoryHelper_SetDefaultQueue,
      nb::arg("param_0"), nb::arg("queue"));

m.def("SUNMemoryHelper_ImplementsRequiredOps",
      SUNMemoryHelper_ImplementsRequiredOps, nb::arg("param_0"));
// #ifdef __cplusplus
//
// #endif
//
// #endif
//
