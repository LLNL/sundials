// #ifndef _ARKSTEP_H
// 
// #ifdef __cplusplus 
// #endif
// 

m.def("ARKStepSetExplicit",
    ARKStepSetExplicit, nb::arg("arkode_mem"));

m.def("ARKStepSetImplicit",
    ARKStepSetImplicit, nb::arg("arkode_mem"));

m.def("ARKStepSetImEx",
    ARKStepSetImEx, nb::arg("arkode_mem"));

m.def("ARKStepSetTables",
    ARKStepSetTables, nb::arg("arkode_mem"), nb::arg("q"), nb::arg("p"), nb::arg("Bi"), nb::arg("Be"));

m.def("ARKStepSetTableNum",
    ARKStepSetTableNum, nb::arg("arkode_mem"), nb::arg("itable"), nb::arg("etable"));

m.def("ARKStepSetTableName",
    ARKStepSetTableName, nb::arg("arkode_mem"), nb::arg("itable"), nb::arg("etable"));

m.def("ARKStepGetTimestepperStats",
    [](void * arkode_mem, long expsteps, long accsteps, long step_attempts, long nfe_evals, long nfi_evals, long nlinsetups, long netfails) -> std::tuple<int, long, long, long, long, long, long, long>
    {
        auto ARKStepGetTimestepperStats_adapt_modifiable_immutable_to_return = [](void * arkode_mem, long expsteps, long accsteps, long step_attempts, long nfe_evals, long nfi_evals, long nlinsetups, long netfails) -> std::tuple<int, long, long, long, long, long, long, long>
        {
            long * expsteps_adapt_modifiable = & expsteps;
            long * accsteps_adapt_modifiable = & accsteps;
            long * step_attempts_adapt_modifiable = & step_attempts;
            long * nfe_evals_adapt_modifiable = & nfe_evals;
            long * nfi_evals_adapt_modifiable = & nfi_evals;
            long * nlinsetups_adapt_modifiable = & nlinsetups;
            long * netfails_adapt_modifiable = & netfails;

            int r = ARKStepGetTimestepperStats(arkode_mem, expsteps_adapt_modifiable, accsteps_adapt_modifiable, step_attempts_adapt_modifiable, nfe_evals_adapt_modifiable, nfi_evals_adapt_modifiable, nlinsetups_adapt_modifiable, netfails_adapt_modifiable);
            return std::make_tuple(r, expsteps, accsteps, step_attempts, nfe_evals, nfi_evals, nlinsetups, netfails);
        };

        return ARKStepGetTimestepperStats_adapt_modifiable_immutable_to_return(arkode_mem, expsteps, accsteps, step_attempts, nfe_evals, nfi_evals, nlinsetups, netfails);
    },     nb::arg("arkode_mem"), nb::arg("expsteps"), nb::arg("accsteps"), nb::arg("step_attempts"), nb::arg("nfe_evals"), nb::arg("nfi_evals"), nb::arg("nlinsetups"), nb::arg("netfails"));
// #ifdef __cplusplus
// 
// #endif
// 
// #endif
// 
