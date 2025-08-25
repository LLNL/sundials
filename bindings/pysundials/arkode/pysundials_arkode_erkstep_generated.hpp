// #ifndef _ERKSTEP_H
// 
// #ifdef __cplusplus 
// #endif
// 

m.def("ERKStepSetTable",
    ERKStepSetTable, nb::arg("arkode_mem"), nb::arg("B"));

m.def("ERKStepSetTableNum",
    ERKStepSetTableNum, nb::arg("arkode_mem"), nb::arg("etable"));

m.def("ERKStepSetTableName",
    ERKStepSetTableName, nb::arg("arkode_mem"), nb::arg("etable"));

m.def("ERKStepGetTimestepperStats",
    [](void * arkode_mem, long expsteps, long accsteps, long step_attempts, long nfevals, long netfails) -> std::tuple<int, long, long, long, long, long>
    {
        auto ERKStepGetTimestepperStats_adapt_modifiable_immutable_to_return = [](void * arkode_mem, long expsteps, long accsteps, long step_attempts, long nfevals, long netfails) -> std::tuple<int, long, long, long, long, long>
        {
            long * expsteps_adapt_modifiable = & expsteps;
            long * accsteps_adapt_modifiable = & accsteps;
            long * step_attempts_adapt_modifiable = & step_attempts;
            long * nfevals_adapt_modifiable = & nfevals;
            long * netfails_adapt_modifiable = & netfails;

            int r = ERKStepGetTimestepperStats(arkode_mem, expsteps_adapt_modifiable, accsteps_adapt_modifiable, step_attempts_adapt_modifiable, nfevals_adapt_modifiable, netfails_adapt_modifiable);
            return std::make_tuple(r, expsteps, accsteps, step_attempts, nfevals, netfails);
        };

        return ERKStepGetTimestepperStats_adapt_modifiable_immutable_to_return(arkode_mem, expsteps, accsteps, step_attempts, nfevals, netfails);
    },     nb::arg("arkode_mem"), nb::arg("expsteps"), nb::arg("accsteps"), nb::arg("step_attempts"), nb::arg("nfevals"), nb::arg("netfails"));
// #ifdef __cplusplus
// 
// #endif
// 
// #endif
// 
