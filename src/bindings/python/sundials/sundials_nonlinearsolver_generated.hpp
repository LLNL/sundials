// #ifndef _SUNNONLINEARSOLVER_H
// 
// #ifdef __cplusplus 
// #endif
// 


auto pyClass_generic_SUNNonlinearSolver_Ops =
    nb::class_<_generic_SUNNonlinearSolver_Ops>
        (m, "_generic_SUNNonlinearSolver_Ops", "Structure containing function pointers to nonlinear solver operations")
    .def(nb::init<>()) // implicit default constructor 
    ;


auto pyClass_generic_SUNNonlinearSolver =
    nb::class_<_generic_SUNNonlinearSolver>
        (m, "_generic_SUNNonlinearSolver", " A nonlinear solver is a structure with an implementation-dependent 'content'\n   field, and a pointer to a structure of solver nonlinear solver operations\n   corresponding to that implementation.")
    .def("__init__", [](_generic_SUNNonlinearSolver * self, SUNNonlinearSolver_Ops ops = SUNNonlinearSolver_Ops(), SUNContext sunctx = SUNContext())
    {
        new (self) _generic_SUNNonlinearSolver();  // placement new
        auto r = self;
        r->ops = ops;
        r->sunctx = sunctx;
    },
    nb::arg("ops") = SUNNonlinearSolver_Ops(), nb::arg("sunctx") = SUNContext()
    )
    .def_rw("content", &_generic_SUNNonlinearSolver::content, "")
    .def_rw("ops", &_generic_SUNNonlinearSolver::ops, "")
    .def_rw("sunctx", &_generic_SUNNonlinearSolver::sunctx, "")
    ;


m.def("SUNNonlinSolNewEmpty",
    SUNNonlinSolNewEmpty, nb::arg("sunctx"));

m.def("SUNNonlinSolGetType",
    SUNNonlinSolGetType, nb::arg("NLS"));

m.def("SUNNonlinSolInitialize",
    SUNNonlinSolInitialize, nb::arg("NLS"));

m.def("SUNNonlinSolSetup",
    SUNNonlinSolSetup, nb::arg("NLS"), nb::arg("y"), nb::arg("mem"));

m.def("SUNNonlinSolSolve",
    SUNNonlinSolSolve, nb::arg("NLS"), nb::arg("y0"), nb::arg("y"), nb::arg("w"), nb::arg("tol"), nb::arg("callLSetup"), nb::arg("mem"));

m.def("SUNNonlinSolSetSysFn",
    SUNNonlinSolSetSysFn, nb::arg("NLS"), nb::arg("SysFn"));

m.def("SUNNonlinSolSetLSetupFn",
    SUNNonlinSolSetLSetupFn, nb::arg("NLS"), nb::arg("SetupFn"));

m.def("SUNNonlinSolSetLSolveFn",
    SUNNonlinSolSetLSolveFn, nb::arg("NLS"), nb::arg("SolveFn"));

m.def("SUNNonlinSolSetConvTestFn",
    SUNNonlinSolSetConvTestFn, nb::arg("NLS"), nb::arg("CTestFn"), nb::arg("ctest_data"));

m.def("SUNNonlinSolSetMaxIters",
    SUNNonlinSolSetMaxIters, nb::arg("NLS"), nb::arg("maxiters"));

m.def("SUNNonlinSolGetNumIters",
    SUNNonlinSolGetNumIters, nb::arg("NLS"), nb::arg("niters"));

m.def("SUNNonlinSolGetCurIter",
    SUNNonlinSolGetCurIter, nb::arg("NLS"), nb::arg("iter"));

m.def("SUNNonlinSolGetNumConvFails",
    SUNNonlinSolGetNumConvFails, nb::arg("NLS"), nb::arg("nconvfails"));
// #ifdef __cplusplus
// 
// #endif
// 
// #endif
// 
