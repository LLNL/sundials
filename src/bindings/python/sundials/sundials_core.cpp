#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>

#include <sundials/sundials_core.hpp>
#include "sundials_logger_impl.h"
#include "sundials_profiler_impl.h"

namespace nb = nanobind;

void bind_adaptcontroller(nb::module_& m);
void bind_sunadjointcheckpointscheme(nb::module_& m);
void bind_sunadjointstepper(nb::module_& m);
void bind_nvector(nb::module_& m);
void bind_linearsolver(nb::module_& m);

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


  bind_adaptcontroller(m);
  // bind_sunadjointcheckpointscheme(m);
  // bind_sunadjointstepper(m);
  bind_nvector(m);
  bind_linearsolver(m);
}