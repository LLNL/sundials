#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>

#include <sundials/sundials_core.hpp>

namespace nb = nanobind;

void bind_suncontext(nb::module_& m);
void bind_sunprofiler(nb::module_& m);
void bind_sunlogger(nb::module_& m);
void bind_sunadaptcontroller(nb::module_& m);
void bind_nvector(nb::module_& m);
void bind_sunadjointcheckpointscheme(nb::module_& m);
void bind_sunadjointstepper(nb::module_& m);
void bind_sunlinearsolver(nb::module_& m);
void bind_sunmatrix(nb::module_& m);
void bind_sunnonlinearsolver(nb::module_& m);
void bind_sunstepper(nb::module_& m);

void bind_core(nb::module_& m)
{
  bind_suncontext(m);
  bind_sunprofiler(m);
  
  bind_nvector(m);
  bind_sunlinearsolver(m);
  bind_sunmatrix(m);
  bind_sunnonlinearsolver(m);

  bind_sunadaptcontroller(m);
  bind_sunstepper(m);
  bind_sunlogger(m);

  bind_sunadjointstepper(m);
  // bind_sunadjointcheckpointscheme(m);
}