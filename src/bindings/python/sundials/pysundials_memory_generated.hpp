// #ifndef _SUNDIALS_MEMORY_H
// 
// #ifdef __cplusplus 
// #endif
// 


auto pyClassSUNMemoryHelper_ =
    nb::class_<SUNMemoryHelper_>
        (m, "SUNMemoryHelper_", "")
    .def("__init__", [](SUNMemoryHelper_ * self, SUNMemoryHelper_Ops ops = SUNMemoryHelper_Ops(), SUNContext sunctx = SUNContext())
    {
        new (self) SUNMemoryHelper_();  // placement new
        auto r = self;
        r->ops = ops;
        r->sunctx = sunctx;
    },
    nb::arg("ops") = SUNMemoryHelper_Ops(), nb::arg("sunctx") = SUNContext()
    )
    .def_rw("content", &SUNMemoryHelper_::content, "")
    .def_rw("queue", &SUNMemoryHelper_::queue, "")
    .def_rw("ops", &SUNMemoryHelper_::ops, "")
    .def_rw("sunctx", &SUNMemoryHelper_::sunctx, "")
    ;


auto pyClassSUNMemoryHelper_Ops_ =
    nb::class_<SUNMemoryHelper_Ops_>
        (m, "SUNMemoryHelper_Ops_", "")
    .def(nb::init<>()) // implicit default constructor 
    ;


m.def("SUNMemoryHelper_Clone",
    SUNMemoryHelper_Clone, nb::arg("param_0"));

m.def("SUNMemoryHelper_SetDefaultQueue",
    SUNMemoryHelper_SetDefaultQueue, nb::arg("param_0"), nb::arg("queue"));

m.def("SUNMemoryHelper_ImplementsRequiredOps",
    SUNMemoryHelper_ImplementsRequiredOps, nb::arg("param_0"));
// #ifdef __cplusplus
// 
// #endif
// 
// #endif
// 
