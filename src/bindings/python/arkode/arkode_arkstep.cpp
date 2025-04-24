#include <nanobind/nanobind.h>
#include <nanobind/stl/function.h>

#include <arkode/arkode.hpp>
#include <arkode/arkode_arkstep.h>
#include <sundials/sundials_core.hpp>

namespace nb = nanobind;

using ark_rhsfn_type = int(sunrealtype, N_Vector, N_Vector, void*);

struct function_pointers
{
  nb::object rhs_fn;
};

int ark_rhsfn_wrapper(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  auto fn_ptrs = static_cast<function_pointers*>(user_data);
  auto ark_rhs = nb::cast<std::function<ark_rhsfn_type>>(fn_ptrs->rhs_fn);
  return ark_rhs(t, y, ydot, user_data);
}

void bind_arkode_arkstep(nb::module_& m)
{

}
