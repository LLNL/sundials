#include <nanobind/nanobind.h>

#include <sundials/sundials_core.hpp>
#include <sunlinsol/sunlinsol_dense.h>

namespace nb = nanobind;

void bind_sunlinsol_dense(nb::module_& m)
{
  m.def("SUNLinSol_Dense", &SUNLinSol_Dense, nb::rv_policy::reference);
}