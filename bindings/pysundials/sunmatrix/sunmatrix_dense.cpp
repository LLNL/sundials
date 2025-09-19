#include <nanobind/nanobind.h>
#include <sundials/sundials_core.hpp>
#include <sunmatrix/sunmatrix_dense.h>

namespace nb = nanobind;

void bind_sunmatrix_dense(nb::module_& m)
{
  m.def("SUNDenseMatrix", &SUNDenseMatrix, nb::rv_policy::reference);
}
