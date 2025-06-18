#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>

#include <sundials/sundials_core.hpp>

namespace nb = nanobind;

void bind_nvector(nb::module_& m);
void bind_sunadaptcontroller(nb::module_& m);
void bind_sunadjointcheckpointscheme(nb::module_& m);
void bind_sunadjointstepper(nb::module_& m);
void bind_suncontext(nb::module_& m);
void bind_sunlinearsolver(nb::module_& m);
void bind_sunlogger(nb::module_& m);
void bind_sunmatrix(nb::module_& m);
void bind_sunmemory(nb::module_& m);
void bind_sunnonlinearsolver(nb::module_& m);
void bind_sunprofiler(nb::module_& m);
void bind_sunstepper(nb::module_& m);

void bind_core(nb::module_& m)
{
#include "pysundials_types_generated.hpp"
#include "pysundials_errors_generated.hpp"

  bind_nvector(m);
  bind_sunadaptcontroller(m);
  bind_sunadjointcheckpointscheme(m);
  bind_sunadjointstepper(m);
  bind_suncontext(m);
  bind_sunlinearsolver(m);
  bind_sunlogger(m);
  bind_sunmatrix(m);
  bind_sunmemory(m);
  bind_sunnonlinearsolver(m);
  bind_sunprofiler(m);
  bind_sunstepper(m);
}