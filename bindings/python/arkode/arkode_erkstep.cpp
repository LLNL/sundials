#include <nanobind/nanobind.h>
#include <nanobind/stl/function.h>

#include <arkode/arkode.hpp>
#include <arkode/arkode_erkstep.h>
#include <sundials/sundials_core.hpp>

namespace nb = nanobind;

using erk_rhsfn_type = int(sunrealtype, N_Vector, N_Vector, void*);

struct function_pointers
{
  nb::object rhs_fn;
};

// TODO(CJB): we will need these wrappers for every callback function
// This method relies on being able to sneak in the std::function
// in user_data. Unfortunately, we do not have a dedicated pointer like user_data
// for every single callback function. We will have to change that.

int erk_rhsfn_wrapper(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  auto fn_ptrs = static_cast<function_pointers*>(user_data);
  auto erk_rhs = nb::cast<std::function<erk_rhsfn_type>>(fn_ptrs->rhs_fn);
  return erk_rhs(t, y, ydot, user_data);
}

void bind_arkode_erkstep(nb::module_& m)
{
  m.def("ERKStepCreate",
        [](std::function<erk_rhsfn_type> rhs, sunrealtype t0, N_Vector y0,
           SUNContext sunctx)
        {
          // static nb::object the_static_rhs = ;
          // function_pointers["rhsfn"] = the_static_rhs;
          function_pointers* ptrs =
            new function_pointers; // need to clean this up in ARKodeFree
          ptrs->rhs_fn  = nb::borrow(nb::cast(
            rhs)); // with nb::borrow, we will need to delete the rhs_fn during cleanup
          void* ark_mem = ERKStepCreate(erk_rhsfn_wrapper, t0, y0, sunctx);
          ARKodeSetUserData(ark_mem, static_cast<void*>(ptrs));
          return ark_mem;
        });
  m.def("ERKStepReInit", &ERKStepReInit);
  m.def("ERKStepSetTable", &ERKStepSetTable);
  m.def("ERKStepSetTableNum", &ERKStepSetTableNum);
  m.def("ERKStepSetTableName", &ERKStepSetTableName);
  // m.def("ERKStepGetCurrentButcherTable", &ERKStepGetCurrentButcherTable);
  m.def("ERKStepGetTimestepperStats", &ERKStepGetTimestepperStats);
}
