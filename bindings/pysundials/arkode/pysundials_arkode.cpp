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

  bind_arkode_erkstep(m);
  bind_arkode_arkstep(m);
}
