#include <nanobind/nanobind.h>

#include <sundials/sundials_core.hpp>
#include <sunlinsol/sunlinsol_lapackdense.h>

namespace nb = nanobind;

void bind_sunlinsol_lapackdense(nb::module_& m)
{
  m.def("SUNLinSol_LapackDense", &SUNLinSol_LapackDense, nb::rv_policy::reference);
}