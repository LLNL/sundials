#include <nanobind/nanobind.h>
#include <nanobind/stl/function.h>

#include <sundials/sundials_core.hpp>

#include <arkode/arkode.h>
#include <arkode/arkode.hpp>
#include <arkode/arkode_ls.h>

#include "sundials_adjointcheckpointscheme_impl.h"

namespace nb = nanobind;

using namespace sundials::experimental;

// Forward declarations of functions defined in other translation units
void bind_arkode_erkstep(nb::module_& m);
void bind_arkode_arkstep(nb::module_& m);

void bind_arkode(nb::module_& m)
{
#include "pysundials_arkode_generated.hpp"

  nb::class_<ARKodeView>(m, "ARKodeView")
    .def_static("Create", &ARKodeView::Create<void*>)
    .def("get", nb::overload_cast<>(&ARKodeView::get, nb::const_),
         nb::rv_policy::reference);

  //
  // arkode_ls.h definitions
  //

  m.def("ARKodeSetLinearSolver", &ARKodeSetLinearSolver, nb::arg("arkode_mem"),
        nb::arg("LS"), nb::arg("A").none());
  m.def("ARKodeSetJacFn", &ARKodeSetJacFn);
  m.def("ARKodeSetJacEvalFrequency", &ARKodeSetJacEvalFrequency);
  m.def("ARKodeSetLinSysFn", &ARKodeSetLinSysFn);
  m.def("ARKodeSetPreconditioner", &ARKodeSetPreconditioner);
  m.def("ARKodeSetLSNormFactor", &ARKodeSetLSNormFactor);
  m.def("ARKodeSetEpsLin", &ARKodeSetEpsLin);
  m.def("ARKodeGetNumLinIters", &ARKodeGetNumLinIters);
  m.def("ARKodeGetNumLinConvFails", &ARKodeGetNumLinConvFails);
  m.def("ARKodeGetNumPrecEvals", &ARKodeGetNumPrecEvals);
  m.def("ARKodeGetNumPrecSolves", &ARKodeGetNumPrecSolves);
  m.def("ARKodeGetNumLinRhsEvals", &ARKodeGetNumLinRhsEvals);
  m.def("ARKodeGetLastLinFlag", &ARKodeGetLastLinFlag);
  m.def("ARKodeGetLinReturnFlagName", &ARKodeGetLinReturnFlagName);

  bind_arkode_erkstep(m);
  bind_arkode_arkstep(m);
}
