// #ifndef _SUNDIALS_MEMORY_H
//
// #ifdef __cplusplus
// #endif
//

auto pyEnumSUNMemoryType_ =
  nb::enum_<SUNMemoryType_>(m, "SUNMemoryType_", nb::is_arithmetic(), "")
    .value("SUN_MEMTYPE_HOST", SUN_MEMTYPE_HOST, "")
    .value("SUN_MEMTYPE_PINNED", SUN_MEMTYPE_PINNED, "")
    .value("SUN_MEMTYPE_DEVICE", SUN_MEMTYPE_DEVICE, "")
    .value("SUN_MEMTYPE_UVM", SUN_MEMTYPE_UVM, "")
    .export_values();

auto pyClassSUNMemoryHelper_ =
  nb::class_<SUNMemoryHelper_>(m, "SUNMemoryHelper_", "")
    .def(nb::init<>()) // implicit default constructor
  ;

auto pyClassSUNMemoryHelper_Ops_ =
  nb::class_<SUNMemoryHelper_Ops_>(m, "SUNMemoryHelper_Ops_", "")
    .def(nb::init<>()) // implicit default constructor
  ;

m.def("SUNMemoryHelper_Clone", SUNMemoryHelper_Clone, nb::arg("param_0"));

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
