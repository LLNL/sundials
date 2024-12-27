#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>

#include <sundials/sundials_core.hpp>
#include "sundials_logger_impl.h"
#include "sundials_profiler_impl.h"

namespace nb = nanobind;

void bind_core(nb::module_& m)
{
  // Since Python will automatically garabage collect objects,
  // we need to interface to our C++ RAII views of objects
  // instead of directly to the C objects. Otherwise, a Python
  // user would have to be extremely careful with managing scope
  // and lifetimes of objects (which would yield an non-idomatic API)
  //
  // Since we are using the C++ views, we do not interface the Destroy functions.

  nb::class_<SUNLogger_>(m, "SUNLogger_");
  m.def("SUNLogger_Create",
        [](SUNComm comm, int output_rank)
        {
          SUNLogger logger;
          SUNErrCode err = SUNLogger_Create(comm, output_rank, &logger);
          if (err != SUN_SUCCESS)
          {
            throw std::runtime_error("Failed to create SUNLogger");
          }
          return logger;
        });
  m.def("SUNLogger_CreateFromEnv",
        [](SUNComm comm)
        {
          SUNLogger logger;
          SUNErrCode err = SUNLogger_CreateFromEnv(comm, &logger);
          if (err != SUN_SUCCESS)
          {
            throw std::runtime_error(
              "Failed to create SUNLogger from environment");
          }
          return logger;
        });
  m.def("SUNLogger_SetErrorFilename", &SUNLogger_SetErrorFilename);
  m.def("SUNLogger_SetWarningFilename", &SUNLogger_SetWarningFilename);
  m.def("SUNLogger_SetDebugFilename", &SUNLogger_SetDebugFilename);
  m.def("SUNLogger_SetInfoFilename", &SUNLogger_SetInfoFilename);
  m.def("SUNLogger_QueueMsg",
        [](SUNLogger logger, SUNLogLevel lvl, const char* scope,
           const char* label, const char* msg_txt)
        {
          SUNErrCode err = SUNLogger_QueueMsg(logger, lvl, scope, label, msg_txt);
          if (err != SUN_SUCCESS)
          {
            throw std::runtime_error("Failed to queue message");
          }
        });
  m.def("SUNLogger_Flush", &SUNLogger_Flush);
  m.def("SUNLogger_GetOutputRank", &SUNLogger_GetOutputRank);

  nb::class_<SUNProfiler_>(m, "SUNProfiler_");
  m.def("SUNProfiler_Create",
        [](SUNComm comm, const char* title)
        {
          SUNProfiler profiler;
          SUNErrCode err = SUNProfiler_Create(comm, title, &profiler);
          if (err != SUN_SUCCESS)
          {
            throw std::runtime_error("Failed to create SUNProfiler");
          }
          return profiler;
        });
  m.def("SUNProfiler_Begin", &SUNProfiler_Begin);
  m.def("SUNProfiler_End", &SUNProfiler_End);
  m.def("SUNProfiler_GetTimerResolution", &SUNProfiler_GetTimerResolution);
  m.def("SUNProfiler_GetElapsedTime", &SUNProfiler_GetElapsedTime);
  m.def("SUNProfiler_Print", &SUNProfiler_Print);
  m.def("SUNProfiler_Reset", &SUNProfiler_Reset);

  nb::class_<SUNContext_>(m, "SUNContext_");
  nb::class_<sundials::Context>(m, "SUNContextView")
    .def(nb::init<>())
    .def(nb::init<SUNComm>(), nb::arg("comm"))
    .def("get", nb::overload_cast<>(&sundials::Context::get, nb::const_),
         nb::rv_policy::reference);

  m.def("SUNContext_GetLastError",
        [](SUNContext sunctx) { return SUNContext_GetLastError(sunctx); });
  m.def("SUNContext_PeekLastError",
        [](SUNContext sunctx) { return SUNContext_PeekLastError(sunctx); });
  m.def("SUNContext_PushErrHandler",
        [](SUNContext sunctx, SUNErrHandlerFn err_fn, void* err_user_data)
        { return SUNContext_PushErrHandler(sunctx, err_fn, err_user_data); });
  m.def("SUNContext_PopErrHandler",
        [](SUNContext sunctx) { return SUNContext_PopErrHandler(sunctx); });
  m.def("SUNContext_ClearErrHandlers",
        [](SUNContext sunctx) { return SUNContext_ClearErrHandlers(sunctx); });
  m.def("SUNContext_GetProfiler",
        [](SUNContext sunctx)
        {
          SUNProfiler profiler;
          SUNErrCode err = SUNContext_GetProfiler(sunctx, &profiler);
          if (err != SUN_SUCCESS)
          {
            throw std::runtime_error("Failed to get SUNProfiler");
          }
          return profiler;
        });
  m.def("SUNContext_SetProfiler", [](SUNContext sunctx, SUNProfiler profiler)
        { return SUNContext_SetProfiler(sunctx, profiler); });
  m.def("SUNContext_GetLogger",
        [](SUNContext sunctx)
        {
          SUNLogger logger;
          SUNErrCode err = SUNContext_GetLogger(sunctx, &logger);
          if (err != SUN_SUCCESS)
          {
            throw std::runtime_error("Failed to get SUNLogger");
          }
          return logger;
        });
  m.def("SUNContext_SetLogger", [](SUNContext sunctx, SUNLogger logger)
        { return SUNContext_SetLogger(sunctx, logger); });

  nb::enum_<N_Vector_ID>(m, "N_Vector_ID")
    .value("SUNDIALS_NVEC_SERIAL", SUNDIALS_NVEC_SERIAL)
    .value("SUNDIALS_NVEC_PARALLEL", SUNDIALS_NVEC_PARALLEL)
    .value("SUNDIALS_NVEC_OPENMP", SUNDIALS_NVEC_OPENMP)
    .value("SUNDIALS_NVEC_PTHREADS", SUNDIALS_NVEC_PTHREADS)
    .value("SUNDIALS_NVEC_PARHYP", SUNDIALS_NVEC_PARHYP)
    .value("SUNDIALS_NVEC_PETSC", SUNDIALS_NVEC_PETSC)
    .value("SUNDIALS_NVEC_CUDA", SUNDIALS_NVEC_CUDA)
    .value("SUNDIALS_NVEC_HIP", SUNDIALS_NVEC_HIP)
    .value("SUNDIALS_NVEC_SYCL", SUNDIALS_NVEC_SYCL)
    .value("SUNDIALS_NVEC_RAJA", SUNDIALS_NVEC_RAJA)
    .value("SUNDIALS_NVEC_KOKKOS", SUNDIALS_NVEC_KOKKOS)
    .value("SUNDIALS_NVEC_OPENMPDEV", SUNDIALS_NVEC_OPENMPDEV)
    .value("SUNDIALS_NVEC_TRILINOS", SUNDIALS_NVEC_TRILINOS)
    .value("SUNDIALS_NVEC_MANYVECTOR", SUNDIALS_NVEC_MANYVECTOR)
    .value("SUNDIALS_NVEC_MPIMANYVECTOR", SUNDIALS_NVEC_MPIMANYVECTOR)
    .value("SUNDIALS_NVEC_MPIPLUSX", SUNDIALS_NVEC_MPIPLUSX)
    .value("SUNDIALS_NVEC_CUSTOM", SUNDIALS_NVEC_CUSTOM);

  nb::class_<_generic_N_Vector>(m, "_generic_N_Vector");

  nb::class_<sundials::experimental::NVectorView>(m, "NVectorView")
    .def(nb::init<>())
    .def(nb::init<_generic_N_Vector*>())
    // Option 1: nv.get() must be invoked in Python to convert to N_Vector before calling N_V functions
    .def("get",
         nb::overload_cast<>(&sundials::experimental::NVectorView::get,
                             nb::const_),
         nb::rv_policy::reference)
    // Option 2: nv.GetArrayPointer() must be invoked in Python and we wrap every N_V function as a class method
    .def(
      "GetArrayPointer",
      [](sundials::experimental::NVectorView&& v)
      { return N_VGetArrayPointer(v); },
      nb::rv_policy::reference);

  // I don't think implicit conversion will work unless we make the View classes convertible to the underlying type instead of the pointer type
  // nb::implicitly_convertible<sundials::experimental::NVectorView, _generic_N_Vector*>();

  m.def("N_VGetArrayPointer",
        [](N_Vector v)
        {
          auto ptr = N_VGetArrayPointer(v);
          if (!ptr) { throw std::runtime_error("Failed to get array pointer"); }
          auto owner = nb::find(v);
          size_t shape[1]{static_cast<size_t>(N_VGetLength(v))};
          return nb::ndarray<sunrealtype, nb::numpy, nb::ndim<1>,
                             nb::c_contig>(ptr, 1, shape, owner);
        });
  m.def("N_VGetDeviceArrayPointer",
        [](N_Vector v)
        {
          auto ptr = N_VGetDeviceArrayPointer(v);
          if (!ptr) { throw std::runtime_error("Failed to get array pointer"); }
          auto owner = nb::find(v);
          size_t shape[1]{static_cast<size_t>(N_VGetLength(v))};
          return nb::ndarray<sunrealtype, nb::numpy, nb::ndim<1>,
                             nb::c_contig>(ptr, 1, shape, owner);
        });
  m.def("N_VSetArrayPointer",
        [](nb::ndarray<sunrealtype, nb::numpy, nb::ndim<1>, nb::c_contig> arr,
           N_Vector v)
        {
          if (arr.shape(0) != N_VGetLength(v))
          {
            throw std::runtime_error(
              "Array shape does not match vector length");
          }
          N_VSetArrayPointer(arr.data(), v);
        });

  m.def("N_VNewEmpty", N_VNewEmpty, nb::arg("sunctx"),
        "py::return_value_policy::reference", nb::rv_policy::reference);
  m.def("N_VCopyOps", N_VCopyOps, nb::arg("w"), nb::arg("v"));
  m.def("N_VGetVectorID", N_VGetVectorID, nb::arg("w"));
  m.def("N_VClone", N_VClone, nb::arg("w"));
  m.def("N_VCloneEmpty", N_VCloneEmpty, nb::arg("w"));
  m.def("N_VSpace", N_VSpace, nb::arg("v"), nb::arg("lrw"), nb::arg("liw"));
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
  // m.def("N_VLinearCombination", N_VLinearCombination, nb::arg("nvec"),
  //       nb::arg("c"), nb::arg("X"), nb::arg("z"), "fused vector operations");
  // m.def("N_VScaleAddMulti", N_VScaleAddMulti, nb::arg("nvec"), nb::arg("a"),
  //       nb::arg("x"), nb::arg("Y"), nb::arg("Z"));
  // m.def("N_VDotProdMulti", N_VDotProdMulti, nb::arg("nvec"), nb::arg("x"),
  //       nb::arg("Y"), nb::arg("dotprods"));
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
  // m.def("N_VDotProdMultiLocal", N_VDotProdMultiLocal, nb::arg("nvec"),
  //       nb::arg("x"), nb::arg("Y"), nb::arg("dotprods"));
  // m.def("N_VDotProdMultiAllReduce", N_VDotProdMultiAllReduce,
  //       nb::arg("nvec_total"), nb::arg("x"), nb::arg("sum"));
  m.def("N_VBufSize", N_VBufSize, nb::arg("x"), nb::arg("size"));
  m.def("N_VBufPack", N_VBufPack, nb::arg("x"), nb::arg("buf"));
  m.def("N_VBufUnpack", N_VBufUnpack, nb::arg("x"), nb::arg("buf"));
  m.def("N_VPrint", N_VPrint, nb::arg("v"));
  m.def("N_VPrintFile", N_VPrintFile, nb::arg("v"), nb::arg("outfile"));
}