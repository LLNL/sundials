#include <nanobind/nanobind.h>
#include <nanobind/stl/function.h>

#include <arkode/arkode.hpp>
#include <arkode/arkode_arkstep.h>
#include <sundials/sundials_core.hpp>

namespace nb = nanobind;

using ark_rhsfn_type = int(sunrealtype, N_Vector, N_Vector, void*);

struct function_pointers
{
  nb::object fe, fi;
};

int ark_fe_wrapper(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  auto fn_ptrs = static_cast<function_pointers*>(user_data);
  auto ark_rhs = nb::cast<std::function<ark_rhsfn_type>>(fn_ptrs->fe);
  return ark_rhs(t, y, ydot, user_data);
}

int ark_fi_wrapper(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  auto fn_ptrs = static_cast<function_pointers*>(user_data);
  auto ark_rhs = nb::cast<std::function<ark_rhsfn_type>>(fn_ptrs->fi);
  return ark_rhs(t, y, ydot, user_data);
}

void bind_arkode_arkstep(nb::module_& m)
{
  m.def(
    "ARKStepCreate",
    [](std::function<ark_rhsfn_type> fe, std::function<ark_rhsfn_type> fi,
       sunrealtype t0, N_Vector y0, SUNContext sunctx)
    {
      function_pointers* ptrs = new function_pointers;

      if (fe) { ptrs->fe = nb::borrow(nb::cast(fe)); }
      if (fi) { ptrs->fi = nb::borrow(nb::cast(fi)); }
      void* ark_mem = ARKStepCreate(fe ? ark_fe_wrapper : nullptr,
                                    fi ? ark_fi_wrapper : nullptr, t0, y0,
                                    sunctx);

      ARKodeSetUserData(ark_mem, static_cast<void*>(ptrs));
      return ark_mem;
    },
    // .none() must be added to functions that accept nullptr as a valid argument
    nb::arg("fe").none(), nb::arg("fi").none(), nb::arg("t0"), nb::arg("y0"),
    nb::arg("sunctx"));
}
