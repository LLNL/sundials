/*---------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ SMU
 *---------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2025, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 *---------------------------------------------------------------
 * This file provides command-line control over optional inputs
 * to ARKODE.
 *--------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sundials/sundials_types.h>
#include "arkode/arkode.h"
#include "arkode/arkode_ls.h"
#include "arkode_impl.h"
#include "sundials_cli.h"

/*---------------------------------------------------------------
  ARKodeSetFromCommandLine:

  Parses the command line to control scalar-valued ARKODE options.
  ---------------------------------------------------------------*/

int ARKodeSetFromCommandLine(void* arkode_mem, const char* arkid, int argc,
                             char* argv[])
{
  ARKodeMem ark_mem;
  if (arkode_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem)arkode_mem;

  /* Set lists of command-line arguments, and the corresponding set routines */
  static const struct sunKeyIntPair int_pairs[] =
    {{"order", ARKodeSetOrder},
     {"interpolant_degree", ARKodeSetInterpolantDegree},
     {"linear", ARKodeSetLinear},
     {"autonomous", ARKodeSetAutonomous},
     {"deduce_implicit_rhs", ARKodeSetDeduceImplicitRhs},
     {"lsetup_frequency", ARKodeSetLSetupFrequency},
     {"predictor_method", ARKodeSetPredictorMethod},
     {"max_nonlin_iters", ARKodeSetMaxNonlinIters},
     {"max_hnil_warns", ARKodeSetMaxHnilWarns},
     {"interpolate_stop_time", ARKodeSetInterpolateStopTime},
     {"max_num_constr_fails", ARKodeSetMaxNumConstrFails},
     {"adaptivity_adjustment", ARKodeSetAdaptivityAdjustment},
     {"small_num_efails", ARKodeSetSmallNumEFails},
     {"max_err_test_fails", ARKodeSetMaxErrTestFails},
     {"max_conv_fails", ARKodeSetMaxConvFails},
     {"linear_solution_scaling", ARKodeSetLinearSolutionScaling}};
  static const int num_int_keys = sizeof(int_pairs) / sizeof(*int_pairs);

  static const struct sunKeyLongPair long_pairs[] = {{"max_num_steps",
                                                      ARKodeSetMaxNumSteps},
                                                     {"jac_eval_frequency",
                                                      ARKodeSetJacEvalFrequency}};
  static const int num_long_keys = sizeof(long_pairs) / sizeof(*long_pairs);

  static const struct sunKeyRealPair real_pairs[] =
    {{"nonlin_crdown", ARKodeSetNonlinCRDown},
     {"nonlin_rdiv", ARKodeSetNonlinRDiv},
     {"delta_gamma_max", ARKodeSetDeltaGammaMax},
     {"nonlin_conv_coef", ARKodeSetNonlinConvCoef},
     {"init_step", ARKodeSetInitStep},
     {"min_step", ARKodeSetMinStep},
     {"max_step", ARKodeSetMaxStep},
     {"stop_time", ARKodeSetStopTime},
     {"fixed_step", ARKodeSetFixedStep},
     {"step_direction", ARKodeSetStepDirection},
     {"cfl_fraction", ARKodeSetCFLFraction},
     {"safety_factor", ARKodeSetSafetyFactor},
     {"error_bias", ARKodeSetErrorBias},
     {"max_growth", ARKodeSetMaxGrowth},
     {"min_reduction", ARKodeSetMinReduction},
     {"max_first_growth", ARKodeSetMaxFirstGrowth},
     {"max_efail_growth", ARKodeSetMaxEFailGrowth},
     {"max_cfail_growth", ARKodeSetMaxCFailGrowth},
     {"eps_lin", ARKodeSetEpsLin},
     {"mass_eps_lin", ARKodeSetMassEpsLin},
     {"ls_norm_factor", ARKodeSetLSNormFactor},
     {"mass_ls_norm_factor", ARKodeSetMassLSNormFactor}};
  static const int num_real_keys = sizeof(real_pairs) / sizeof(*real_pairs);

  static const struct sunKeyTwoRealPair tworeal_pairs[] =
    {{"scalar_tolerances", ARKodeSStolerances},
     {"fixed_step_bounds", ARKodeSetFixedStepBounds}};
  static const int num_tworeal_keys = sizeof(tworeal_pairs) /
                                      sizeof(*tworeal_pairs);

  static const struct sunKeyActionPair action_pairs[] =
    {{"nonlinear", ARKodeSetNonlinear},
     {"clear_stop_time", ARKodeClearStopTime},
     {"no_inactive_root_warn", ARKodeSetNoInactiveRootWarn},
     {"reset_accumulated_error", ARKodeResetAccumulatedError}};
  static const int num_action_keys = sizeof(action_pairs) / sizeof(*action_pairs);

  int idx, j, retval;
  SUNErrCode sunretval;
  sunbooleantype write_parameters = SUNFALSE;
  for (idx = 1; idx < argc; idx++)
  {
    sunbooleantype arg_used = SUNFALSE;

    /* if arkid is supplied, skip command-line arguments that do not begin with arkid;
       else, skip command-line arguments that do not begin with "arkode." */
    size_t offset;
    if (arkid != NULL)
    {
      if (strncmp(argv[idx], arkid, strlen(arkid)) != 0) { continue; }
      offset = strlen(arkid) + 1;
    }
    else
    {
      static const char* prefix = "arkode.";
      if (strncmp(argv[idx], prefix, strlen(prefix)) != 0) { continue; }
      offset = strlen(prefix);
    }

    /* check all "int" command-line options */
    for (j = 0; j < num_int_keys; j++)
    {
      sunretval = sunCheckAndSetIntArg(arkode_mem, &idx, argv, offset,
                                       int_pairs[j].key, int_pairs[j].set,
                                       &arg_used);
      if (sunretval != SUN_SUCCESS)
      {
        retval = ARK_ILL_INPUT;
        arkProcessError(ark_mem, retval, __LINE__, __func__, __FILE__,
                        "error setting command-line argument: %s",
                        int_pairs[j].key);
        return retval;
      }
      if (arg_used) break;
    }
    if (arg_used) continue;

    /* check all long int command-line options */
    for (j = 0; j < num_long_keys; j++)
    {
      sunretval = sunCheckAndSetLongArg(arkode_mem, &idx, argv, offset,
                                        long_pairs[j].key, long_pairs[j].set,
                                        &arg_used);
      if (sunretval != SUN_SUCCESS)
      {
        retval = ARK_ILL_INPUT;
        arkProcessError(ark_mem, retval, __LINE__, __func__, __FILE__,
                        "error setting command-line argument: %s",
                        long_pairs[j].key);
        return retval;
      }
      if (arg_used) break;
    }
    if (arg_used) continue;

    /* check all real command-line options */
    for (j = 0; j < num_real_keys; j++)
    {
      sunretval = sunCheckAndSetRealArg(arkode_mem, &idx, argv, offset,
                                        real_pairs[j].key, real_pairs[j].set,
                                        &arg_used);
      if (sunretval != SUN_SUCCESS)
      {
        retval = ARK_ILL_INPUT;
        arkProcessError(ark_mem, retval, __LINE__, __func__, __FILE__,
                        "error setting command-line argument: %s",
                        real_pairs[j].key);
        return retval;
      }
      if (arg_used) break;
    }
    if (arg_used) continue;

    /* check all pair-of-real command-line options */
    for (j = 0; j < num_tworeal_keys; j++)
    {
      sunretval = sunCheckAndSetTwoRealArg(arkode_mem, &idx, argv, offset,
                                           tworeal_pairs[j].key,
                                           tworeal_pairs[j].set, &arg_used);
      if (sunretval != SUN_SUCCESS)
      {
        retval = ARK_ILL_INPUT;
        arkProcessError(ark_mem, retval, __LINE__, __func__, __FILE__,
                        "error setting command-line argument: %s",
                        tworeal_pairs[j].key);
        return retval;
      }
      if (arg_used) break;
    }
    if (arg_used) continue;

    /* check all action command-line options */
    for (j = 0; j < num_action_keys; j++)
    {
      sunretval = sunCheckAndSetActionArg(arkode_mem, &idx, argv, offset,
                                          action_pairs[j].key,
                                          action_pairs[j].set, &arg_used);
      if (sunretval != SUN_SUCCESS)
      {
        retval = ARK_ILL_INPUT;
        arkProcessError(ark_mem, retval, __LINE__, __func__, __FILE__,
                        "error setting command-line argument: %s",
                        action_pairs[j].key);
        return retval;
      }
      if (arg_used) break;
    }
    if (arg_used) continue;

    /*** handle all remaining command-line options ***/

    if (strcmp(argv[idx] + offset, "interpolant_type") == 0)
    {
      idx++;
      retval = ARK_ILL_INPUT;
      if (strcmp(argv[idx], "ARK_INTERP_HERMITE") == 0)
      {
        retval = ARKodeSetInterpolantType(arkode_mem, ARK_INTERP_HERMITE);
      }
      else if (strcmp(argv[idx], "ARK_INTERP_LAGRANGE") == 0)
      {
        retval = ARKodeSetInterpolantType(arkode_mem, ARK_INTERP_LAGRANGE);
      }
      else if (strcmp(argv[idx], "ARK_INTERP_NONE") == 0)
      {
        retval = ARKodeSetInterpolantType(arkode_mem, ARK_INTERP_NONE);
      }
      if (retval != ARK_SUCCESS)
      {
        arkProcessError(ark_mem, retval, __LINE__, __func__, __FILE__,
                        "error setting command-line argument: %s %s",
                        argv[idx - 1], argv[idx]);
        return retval;
      }
      arg_used = SUNTRUE;
      continue;
    }

    if (strcmp(argv[idx] + offset, "accum_error_type") == 0)
    {
      idx++;
      retval = ARK_ILL_INPUT;
      if (strcmp(argv[idx], "ARK_ACCUMERROR_NONE") == 0)
      {
        retval = ARKodeSetAccumulatedErrorType(arkode_mem, ARK_ACCUMERROR_NONE);
      }
      else if (strcmp(argv[idx], "ARK_ACCUMERROR_MAX") == 0)
      {
        retval = ARKodeSetAccumulatedErrorType(arkode_mem, ARK_ACCUMERROR_MAX);
      }
      else if (strcmp(argv[idx], "ARK_ACCUMERROR_SUM") == 0)
      {
        retval = ARKodeSetAccumulatedErrorType(arkode_mem, ARK_ACCUMERROR_SUM);
      }
      else if (strcmp(argv[idx], "ARK_ACCUMERROR_AVG") == 0)
      {
        retval = ARKodeSetAccumulatedErrorType(arkode_mem, ARK_ACCUMERROR_AVG);
      }
      if (retval != ARK_SUCCESS)
      {
        arkProcessError(ark_mem, retval, __LINE__, __func__, __FILE__,
                        "error setting command-line argument: %s %s",
                        argv[idx - 1], argv[idx]);
        return retval;
      }
      arg_used = SUNTRUE;
      continue;
    }

    if (strcmp(argv[idx] + offset, "write_parameters") == 0)
    {
      write_parameters = SUNTRUE;
      arg_used         = SUNTRUE;
      continue;
    }

    /* Call stepper-specific SetFromCommandLine routine (if supplied) to
       process this command-line argument */
    if (ark_mem->step_setfromcommandline)
    {
      retval = ark_mem->step_setfromcommandline(ark_mem, &idx, argv, offset,
                                                &arg_used);
      if (retval != ARK_SUCCESS) { return retval; }
      if (arg_used) { continue; }
    }

    /* warn for uninterpreted arkid.X arguments */
    arkProcessError(ark_mem, ARK_WARNING, __LINE__, __func__, __FILE__,
                    "WARNING: argument %s was not handled\n", argv[idx]);
  }

  /* Call ARKodeWriteParameters (if requested) now that all
     command-line options have been set -- WARNING: this knows
     nothing about MPI, so it could be redundantly written by all
     processes if requested. */
  if (write_parameters)
  {
    retval = ARKodeWriteParameters(arkode_mem, stdout);
    if (retval != ARK_SUCCESS)
    {
      arkProcessError(ark_mem, retval, __LINE__, __func__, __FILE__,
                      "error writing parameters to stdout");
      return retval;
    }
  }

  return (ARK_SUCCESS);
}

/*===============================================================
  EOF
  ===============================================================*/
