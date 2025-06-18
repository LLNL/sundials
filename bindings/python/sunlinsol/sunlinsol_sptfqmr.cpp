#include <nanobind/nanobind.h>

#include <sundials/sundials_core.hpp>
#include <sunlinsol/sunlinsol_sptfqmr.h>

namespace nb = nanobind;

void bind_sunlinsol_sptfqmr(nb::module_& m)
{
  m.def("SUNLinSol_SPTFQMR", &SUNLinSol_SPTFQMR, nb::rv_policy::reference);
}