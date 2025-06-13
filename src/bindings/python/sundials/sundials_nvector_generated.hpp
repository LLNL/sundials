// #ifndef _NVECTOR_H
//
// #ifdef __cplusplus
// #endif
//

auto pyEnumN_Vector_ID =
  nb::enum_<N_Vector_ID>(m, "N_Vector_ID", nb::is_arithmetic(), "")
    .value("SUNDIALS_NVEC_SERIAL", SUNDIALS_NVEC_SERIAL, "")
    .value("SUNDIALS_NVEC_PARALLEL", SUNDIALS_NVEC_PARALLEL, "")
    .value("SUNDIALS_NVEC_OPENMP", SUNDIALS_NVEC_OPENMP, "")
    .value("SUNDIALS_NVEC_PTHREADS", SUNDIALS_NVEC_PTHREADS, "")
    .value("SUNDIALS_NVEC_PARHYP", SUNDIALS_NVEC_PARHYP, "")
    .value("SUNDIALS_NVEC_PETSC", SUNDIALS_NVEC_PETSC, "")
    .value("SUNDIALS_NVEC_CUDA", SUNDIALS_NVEC_CUDA, "")
    .value("SUNDIALS_NVEC_HIP", SUNDIALS_NVEC_HIP, "")
    .value("SUNDIALS_NVEC_SYCL", SUNDIALS_NVEC_SYCL, "")
    .value("SUNDIALS_NVEC_RAJA", SUNDIALS_NVEC_RAJA, "")
    .value("SUNDIALS_NVEC_KOKKOS", SUNDIALS_NVEC_KOKKOS, "")
    .value("SUNDIALS_NVEC_OPENMPDEV", SUNDIALS_NVEC_OPENMPDEV, "")
    .value("SUNDIALS_NVEC_TRILINOS", SUNDIALS_NVEC_TRILINOS, "")
    .value("SUNDIALS_NVEC_MANYVECTOR", SUNDIALS_NVEC_MANYVECTOR, "")
    .value("SUNDIALS_NVEC_MPIMANYVECTOR", SUNDIALS_NVEC_MPIMANYVECTOR, "")
    .value("SUNDIALS_NVEC_MPIPLUSX", SUNDIALS_NVEC_MPIPLUSX, "")
    .value("SUNDIALS_NVEC_CUSTOM", SUNDIALS_NVEC_CUSTOM, "");

auto pyClass_generic_N_Vector_Ops =
  nb::class_<_generic_N_Vector_Ops>(m,
                                    "_generic_N_Vector_Ops", "Structure containing function pointers to vector operations")
    .def(nb::init<>()) // implicit default constructor
  ;

auto pyClass_generic_N_Vector =
  nb::class_<_generic_N_Vector>(m,
                                "_generic_N_Vector", " A vector is a structure with an implementation-dependent\n   'content' field, and a pointer to a structure of vector\n   operations corresponding to that implementation.")
    .def(
      "__init__",
      [](_generic_N_Vector* self, N_Vector_Ops ops = N_Vector_Ops(),
         SUNContext sunctx = SUNContext())
      {
        new (self) _generic_N_Vector(); // placement new
        auto r    = self;
        r->ops    = ops;
        r->sunctx = sunctx;
      },
      nb::arg("ops") = N_Vector_Ops(), nb::arg("sunctx") = SUNContext())
    .def_rw("content", &_generic_N_Vector::content, "")
    .def_rw("ops", &_generic_N_Vector::ops, "")
    .def_rw("sunctx", &_generic_N_Vector::sunctx, "");

m.def("N_VNewEmpty", N_VNewEmpty, nb::arg("sunctx"),
      "py::return_value_policy::reference", nb::rv_policy::reference);

m.def("N_VCopyOps", N_VCopyOps, nb::arg("w"), nb::arg("v"));

m.def("N_VGetVectorID", N_VGetVectorID, nb::arg("w"));

m.def("N_VClone", N_VClone, nb::arg("w"), "py::return_value_policy::reference",
      nb::rv_policy::reference);

m.def("N_VCloneEmpty", N_VCloneEmpty, nb::arg("w"),
      "py::return_value_policy::reference", nb::rv_policy::reference);

m.def("N_VGetCommunicator", N_VGetCommunicator, nb::arg("v"));

m.def("N_VGetLength", N_VGetLength, nb::arg("v"));

m.def("N_VGetLocalLength", N_VGetLocalLength, nb::arg("v"));

m.def("N_VLinearSum", N_VLinearSum, nb::arg("a"), nb::arg("x"), nb::arg("b"),
      nb::arg("y"), nb::arg("z"));

m.def("N_VConst", N_VConst, nb::arg("c"), nb::arg("z"));

m.def("N_VProd", N_VProd, nb::arg("x"), nb::arg("y"), nb::arg("z"));

m.def("N_VDiv", N_VDiv, nb::arg("x"), nb::arg("y"), nb::arg("z"));

m.def("N_VScale", N_VScale, nb::arg("c"), nb::arg("x"), nb::arg("z"));

m.def("N_VAbs", N_VAbs, nb::arg("x"), nb::arg("z"));

m.def("N_VInv", N_VInv, nb::arg("x"), nb::arg("z"));

m.def("N_VAddConst", N_VAddConst, nb::arg("x"), nb::arg("b"), nb::arg("z"));

m.def("N_VDotProd", N_VDotProd, nb::arg("x"), nb::arg("y"));

m.def("N_VMaxNorm", N_VMaxNorm, nb::arg("x"));

m.def("N_VWrmsNorm", N_VWrmsNorm, nb::arg("x"), nb::arg("w"));

m.def("N_VWrmsNormMask", N_VWrmsNormMask, nb::arg("x"), nb::arg("w"),
      nb::arg("id"));

m.def("N_VMin", N_VMin, nb::arg("x"));

m.def("N_VWL2Norm", N_VWL2Norm, nb::arg("x"), nb::arg("w"));

m.def("N_VL1Norm", N_VL1Norm, nb::arg("x"));

m.def("N_VCompare", N_VCompare, nb::arg("c"), nb::arg("x"), nb::arg("z"));

m.def("N_VInvTest", N_VInvTest, nb::arg("x"), nb::arg("z"));

m.def("N_VConstrMask", N_VConstrMask, nb::arg("c"), nb::arg("x"), nb::arg("m"));

m.def("N_VMinQuotient", N_VMinQuotient, nb::arg("num"), nb::arg("denom"));

m.def("N_VDotProdLocal", N_VDotProdLocal, nb::arg("x"), nb::arg("y"));

m.def("N_VMaxNormLocal", N_VMaxNormLocal, nb::arg("x"));

m.def("N_VMinLocal", N_VMinLocal, nb::arg("x"));

m.def("N_VL1NormLocal", N_VL1NormLocal, nb::arg("x"));

m.def("N_VWSqrSumLocal", N_VWSqrSumLocal, nb::arg("x"), nb::arg("w"));

m.def("N_VWSqrSumMaskLocal", N_VWSqrSumMaskLocal, nb::arg("x"), nb::arg("w"),
      nb::arg("id"));

m.def("N_VInvTestLocal", N_VInvTestLocal, nb::arg("x"), nb::arg("z"));

m.def("N_VConstrMaskLocal", N_VConstrMaskLocal, nb::arg("c"), nb::arg("x"),
      nb::arg("m"));

m.def("N_VMinQuotientLocal", N_VMinQuotientLocal, nb::arg("num"),
      nb::arg("denom"));

m.def("N_VBufSize", N_VBufSize, nb::arg("x"), nb::arg("size"));

m.def("N_VPrint", N_VPrint, nb::arg("v"));

m.def("N_VPrintFile", N_VPrintFile, nb::arg("v"), nb::arg("outfile"));
// #ifdef __cplusplus
//
// #endif
//
// #endif
//
