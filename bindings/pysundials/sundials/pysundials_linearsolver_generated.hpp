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
// #ifndef _SUNLINEARSOLVER_H
// 
// #ifdef __cplusplus 
// #endif
// 


auto pyClass_generic_SUNLinearSolver_Ops =
    nb::class_<_generic_SUNLinearSolver_Ops>
        (m, "_generic_SUNLinearSolver_Ops", "Structure containing function pointers to linear solver operations")
    .def(nb::init<>()) // implicit default constructor 
    ;


auto pyClass_generic_SUNLinearSolver =
    nb::class_<_generic_SUNLinearSolver>
        (m, "_generic_SUNLinearSolver", " A linear solver is a structure with an implementation-dependent\n   'content' field, and a pointer to a structure of linear solver\n   operations corresponding to that implementation.")
    .def("__init__", [](_generic_SUNLinearSolver * self, SUNLinearSolver_Ops ops = SUNLinearSolver_Ops(), SUNContext sunctx = SUNContext())
    {
        new (self) _generic_SUNLinearSolver();  // placement new
        auto r_ctor_ = self;
        r_ctor_->ops = ops;
        r_ctor_->sunctx = sunctx;
    },
    nb::arg("ops") = SUNLinearSolver_Ops(), nb::arg("sunctx") = SUNContext()
    )
    .def_rw("content", &_generic_SUNLinearSolver::content, "")
    .def_rw("ops", &_generic_SUNLinearSolver::ops, "")
    .def_rw("sunctx", &_generic_SUNLinearSolver::sunctx, "")
    ;


m.def("SUNLinSolGetType",
    SUNLinSolGetType, nb::arg("S"));

m.def("SUNLinSolGetID",
    SUNLinSolGetID, nb::arg("S"));

m.def("SUNLinSolSetATimes",
    SUNLinSolSetATimes, nb::arg("S"), nb::arg("A_data"), nb::arg("ATimes"));

m.def("SUNLinSolSetPreconditioner",
    SUNLinSolSetPreconditioner, nb::arg("S"), nb::arg("P_data"), nb::arg("Pset"), nb::arg("Psol"));

m.def("SUNLinSolSetScalingVectors",
    SUNLinSolSetScalingVectors, nb::arg("S"), nb::arg("s1"), nb::arg("s2"));

m.def("SUNLinSolSetZeroGuess",
    SUNLinSolSetZeroGuess, nb::arg("S"), nb::arg("onoff"));

m.def("SUNLinSolInitialize",
    SUNLinSolInitialize, nb::arg("S"));

m.def("SUNLinSolSetup",
    SUNLinSolSetup, nb::arg("S"), nb::arg("A"));

m.def("SUNLinSolSolve",
    SUNLinSolSolve, nb::arg("S"), nb::arg("A"), nb::arg("x"), nb::arg("b"), nb::arg("tol"));

m.def("SUNLinSolNumIters",
    SUNLinSolNumIters, nb::arg("S"));

m.def("SUNLinSolResNorm",
    SUNLinSolResNorm, nb::arg("S"));

m.def("SUNLinSolResid",
    SUNLinSolResid, nb::arg("S"));

m.def("SUNLinSolLastFlag",
    SUNLinSolLastFlag, nb::arg("S"));
// #ifdef __cplusplus
// 
// #endif
// 
// #endif
// 
