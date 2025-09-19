#include <nanobind/nanobind.h>
#include <sundials/sundials_core.hpp>
#include <sunmatrix/sunmatrix_sparse.h>

namespace nb = nanobind;

void bind_sunmatrix_sparse(nb::module_& m)
{
  m.def("SUNSparseMatrix", &SUNSparseMatrix, nb::rv_policy::reference);
}
