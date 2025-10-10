// #ifndef _SUNADJOINT_CHECKPOINTSCHEME_H
//
// #ifdef __cplusplus
// #endif
//

m.def(
  "SUNAdjointCheckpointScheme_NeedsSaving",
  [](SUNAdjointCheckpointScheme check_scheme, suncountertype step_num,
     suncountertype stage_num,
     sunrealtype t) -> std::tuple<SUNErrCode, sunbooleantype>
  {
    auto SUNAdjointCheckpointScheme_NeedsSaving_adapt_modifiable_immutable_to_return =
      [](SUNAdjointCheckpointScheme check_scheme, suncountertype step_num,
         suncountertype stage_num,
         sunrealtype t) -> std::tuple<SUNErrCode, sunbooleantype>
    {
      sunbooleantype yes_or_no_adapt_modifiable;

      SUNErrCode r =
        SUNAdjointCheckpointScheme_NeedsSaving(check_scheme, step_num, stage_num,
                                               t, &yes_or_no_adapt_modifiable);
      return std::make_tuple(r, yes_or_no_adapt_modifiable);
    };

    return SUNAdjointCheckpointScheme_NeedsSaving_adapt_modifiable_immutable_to_return(check_scheme,
                                                                                       step_num,
                                                                                       stage_num,
                                                                                       t);
  },
  nb::arg("check_scheme"), nb::arg("step_num"), nb::arg("stage_num"),
  nb::arg("t"));

m.def(
  "SUNAdjointCheckpointScheme_LoadVector",
  [](SUNAdjointCheckpointScheme check_scheme, suncountertype step_num,
     suncountertype stage_num,
     sunbooleantype peek) -> std::tuple<SUNErrCode, N_Vector, sunrealtype>
  {
    auto SUNAdjointCheckpointScheme_LoadVector_adapt_modifiable_immutable_to_return =
      [](SUNAdjointCheckpointScheme check_scheme, suncountertype step_num,
         suncountertype stage_num,
         sunbooleantype peek) -> std::tuple<SUNErrCode, N_Vector, sunrealtype>
    {
      N_Vector out_adapt_modifiable;
      sunrealtype tout_adapt_modifiable;

      SUNErrCode r =
        SUNAdjointCheckpointScheme_LoadVector(check_scheme, step_num, stage_num,
                                              peek, &out_adapt_modifiable,
                                              &tout_adapt_modifiable);
      return std::make_tuple(r, out_adapt_modifiable, tout_adapt_modifiable);
    };

    return SUNAdjointCheckpointScheme_LoadVector_adapt_modifiable_immutable_to_return(check_scheme,
                                                                                      step_num,
                                                                                      stage_num,
                                                                                      peek);
  },
  nb::arg("check_scheme"), nb::arg("step_num"), nb::arg("stage_num"),
  nb::arg("peek"), nb::rv_policy::reference);

m.def("SUNAdjointCheckpointScheme_EnableDense",
      SUNAdjointCheckpointScheme_EnableDense, nb::arg("check_scheme"),
      nb::arg("on_or_off"));
// #ifdef __cplusplus
//
// #endif
//
// #endif
