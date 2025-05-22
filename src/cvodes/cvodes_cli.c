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
 * to CVODE.
 *--------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sundials/sundials_types.h>
#include "cvodes/cvodes.h"
#include "cvodes/cvodes_ls.h"
#include "cvodes/cvodes_proj.h"
#include "cvodes_impl.h"
#include "sundials_cli.h"

/*---------------------------------------------------------------
  CVodeSetFromCommandLine:

  Parses the command line to control scalar-valued CVODE options.
  ---------------------------------------------------------------*/

int CVodeSetFromCommandLine(void* cvode_mem, const char* cvid, int argc,
                            char* argv[])
{
  CVodeMem cv_mem;
  if (cvode_mem == NULL)
  {
    cvProcessError(NULL, CV_MEM_NULL, __LINE__, __func__, __FILE__, MSGCV_NO_MEM);
    return (CV_MEM_NULL);
  }
  cv_mem = (CVodeMem)cvode_mem;

  /* Set lists of command-line arguments, and the corresponding set routines */
  static const struct sunKeyIntPair int_pairs[] =
    {{"max_conv_fails", CVodeSetMaxConvFails},
     {"max_err_test_fails", CVodeSetMaxErrTestFails},
     {"max_hnil_warns", CVodeSetMaxHnilWarns},
     {"max_nonlin_iters", CVodeSetMaxNonlinIters},
     {"max_order", CVodeSetMaxOrd},
     {"stab_lim_det", CVodeSetStabLimDet},
     {"interpolate_stop_time", CVodeSetInterpolateStopTime},
     {"num_efails_eta_max_err_fail", CVodeSetNumFailsEtaMaxErrFail},
     {"quad_err_con", CVodeSetQuadErrCon},
     {"sens_err_con", CVodeSetSensErrCon},
     {"sens_max_nonlin_iters", CVodeSetSensMaxNonlinIters},
     {"quad_sens_err_con", CVodeSetQuadSensErrCon},
     {"linear_solution_scaling", CVodeSetLinearSolutionScaling},
     {"proj_err_est", CVodeSetProjErrEst},
     {"max_num_proj_fails", CVodeSetMaxNumProjFails}};
  static const int num_int_keys = sizeof(int_pairs) / sizeof(*int_pairs);

  static const struct sunKeyLongPair long_pairs[] =
    {{"lsetup_frequency", CVodeSetLSetupFrequency},
     {"max_num_steps", CVodeSetMaxNumSteps},
     {"monitor_frequency", CVodeSetMonitorFrequency},
     {"jac_eval_frequency", CVodeSetJacEvalFrequency},
     {"proj_frequency", CVodeSetProjFrequency}};
  static const int num_long_keys = sizeof(long_pairs) / sizeof(*long_pairs);

  static const struct sunKeyRealPair real_pairs[] =
    {{"dgmax_lsetup", CVodeSetDeltaGammaMaxLSetup},
     {"init_step", CVodeSetInitStep},
     {"max_step", CVodeSetMaxStep},
     {"min_step", CVodeSetMinStep},
     {"nonlin_conv_coef", CVodeSetNonlinConvCoef},
     {"stop_time", CVodeSetStopTime},
     {"eta_max_first_step", CVodeSetEtaMaxFirstStep},
     {"eta_max_early_step", CVodeSetEtaMaxEarlyStep},
     {"eta_max", CVodeSetEtaMax},
     {"eta_min", CVodeSetEtaMin},
     {"eta_min_err_fail", CVodeSetEtaMinErrFail},
     {"eta_max_err_fail", CVodeSetEtaMaxErrFail},
     {"eta_conv_fail", CVodeSetEtaConvFail},
     {"delta_gamma_max_bad_jac", CVodeSetDeltaGammaMaxBadJac},
     {"eps_lin", CVodeSetEpsLin},
     {"ls_norm_factor", CVodeSetLSNormFactor},
     {"eps_proj", CVodeSetEpsProj},
     {"proj_fail_eta", CVodeSetProjFailEta}};
  static const int num_real_keys = sizeof(real_pairs) / sizeof(*real_pairs);

  static const struct sunKeyTwoRealPair tworeal_pairs[] =
    {{"eta_fixed_step_bounds", CVodeSetEtaFixedStepBounds},
     {"scalar_tolerances", CVodeSStolerances},
     {"quad_scalar_tolerances", CVodeQuadSStolerances}};
  static const int num_tworeal_keys = sizeof(tworeal_pairs) /
                                      sizeof(*tworeal_pairs);

  static const struct sunKeyTwoIntPair twoint_pairs[] =
    {{"max_order_b", CVodeSetMaxOrdB},
     {"stab_lim_det_b", CVodeSetStabLimDetB},
     {"quad_err_con_b", CVodeSetQuadErrConB},
     {"linear_solution_scaling_b", CVodeSetLinearSolutionScalingB}};
  static const int num_twoint_keys = sizeof(twoint_pairs) / sizeof(*twoint_pairs);

  static const struct sunKeyActionPair action_pairs[] =
    {{"clear_stop_time", CVodeClearStopTime},
     {"adj_no_sensi", CVodeSetAdjNoSensi},
     {"no_inactive_root_warn", CVodeSetNoInactiveRootWarn}};
  static const int num_action_keys = sizeof(action_pairs) / sizeof(*action_pairs);

  static const struct sunKeyIntRealPair int_real_pairs[] =
    {{"sens_dq_method", CVodeSetSensDQMethod},
     {"init_step_b", CVodeSetInitStepB},
     {"min_step_b", CVodeSetMinStepB},
     {"max_step_b", CVodeSetMaxStepB},
     {"eps_lin_b", CVodeSetEpsLinB},
     {"ls_norm_factor_b", CVodeSetLSNormFactorB}};
  static const int num_int_real_keys = sizeof(int_real_pairs) /
                                       sizeof(*int_real_pairs);

  static const struct sunKeyIntRealRealPair int_real_real_pairs[] =
    {{"scalar_tolerances_b", CVodeSStolerancesB},
     {"quad_scalar_tolerances_b", CVodeQuadSStolerancesB}};
  static const int num_int_real_real_keys = sizeof(int_real_real_pairs) /
                                            sizeof(*int_real_real_pairs);

  static const struct sunKeyIntLongPair int_long_pairs[] = {
    {"max_num_steps_b", CVodeSetMaxNumStepsB}};
  static const int num_int_long_keys = sizeof(int_long_pairs) /
                                       sizeof(*int_long_pairs);

  SUNErrCode sunretval;
  int idx, j, retval;
  for (idx = 1; idx < argc; idx++)
  {
    sunbooleantype arg_used = SUNFALSE;

    /* if cvid is supplied, skip command-line arguments that do not begin with cvid;
       else, skip command-line arguments that do not begin with "cvodes." */
    size_t offset;
    if (cvid != NULL)
    {
      if (strncmp(argv[idx], cvid, strlen(cvid)) != 0) { continue; }
      offset = strlen(cvid) + 1;
    }
    else
    {
      static const char* prefix = "cvodes.";
      if (strncmp(argv[idx], prefix, strlen(prefix)) != 0) { continue; }
      offset = strlen(prefix);
    }

    /* check all "int" command-line options */
    for (j = 0; j < num_int_keys; j++)
    {
      sunretval = sunCheckAndSetIntArg(cvode_mem, &idx, argv, offset, int_pairs[j].key,
                                    int_pairs[j].set, &arg_used);
      if (sunretval != CV_SUCCESS)
      {
        retval = CV_ILL_INPUT;
        cvProcessError(cv_mem, retval, __LINE__, __func__, __FILE__,
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
      sunretval = sunCheckAndSetLongArg(cvode_mem, &idx, argv, offset,
                                     long_pairs[j].key, long_pairs[j].set,
                                     &arg_used);
      if (sunretval != CV_SUCCESS)
      {
        retval = CV_ILL_INPUT;
        cvProcessError(cv_mem, retval, __LINE__, __func__, __FILE__,
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
      sunretval = sunCheckAndSetRealArg(cvode_mem, &idx, argv, offset,
                                     real_pairs[j].key, real_pairs[j].set,
                                     &arg_used);
      if (sunretval != CV_SUCCESS)
      {
        retval = CV_ILL_INPUT;
        cvProcessError(cv_mem, retval, __LINE__, __func__, __FILE__,
                       "error setting command-line argument: %s",
                       real_pairs[j].key);
        return retval;
      }
      if (arg_used) break;
    }
    if (arg_used) continue;

    /* check all pair-of-int command-line options */
    for (j = 0; j < num_twoint_keys; j++)
    {
      sunretval = sunCheckAndSetTwoIntArg(cvode_mem, &idx, argv, offset,
                                       twoint_pairs[j].key, twoint_pairs[j].set,
                                       &arg_used);
      if (sunretval != CV_SUCCESS)
      {
        retval = CV_ILL_INPUT;
        cvProcessError(cv_mem, retval, __LINE__, __func__, __FILE__,
                       "error setting command-line argument: %s",
                       twoint_pairs[j].key);
        return retval;
      }
      if (arg_used) break;
    }
    if (arg_used) continue;

    /* check all pair-of-real command-line options */
    for (j = 0; j < num_tworeal_keys; j++)
    {
      sunretval = sunCheckAndSetTwoRealArg(cvode_mem, &idx, argv, offset,
                                        tworeal_pairs[j].key,
                                        tworeal_pairs[j].set, &arg_used);
      if (sunretval != CV_SUCCESS)
      {
        retval = CV_ILL_INPUT;
        cvProcessError(cv_mem, retval, __LINE__, __func__, __FILE__,
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
      sunretval = sunCheckAndSetActionArg(cvode_mem, &idx, argv, offset,
                                       action_pairs[j].key, action_pairs[j].set,
                                       &arg_used);
      if (sunretval != CV_SUCCESS)
      {
        retval = CV_ILL_INPUT;
        cvProcessError(cv_mem, retval, __LINE__, __func__, __FILE__,
                       "error setting command-line argument: %s",
                       action_pairs[j].key);
        return retval;
      }
      if (arg_used) break;
    }
    if (arg_used) continue;

    /* check all int+real command-line options */
    for (j = 0; j < num_int_real_keys; j++)
    {
      sunretval = sunCheckAndSetIntRealArg(cvode_mem, &idx, argv, offset,
                                        int_real_pairs[j].key,
                                        int_real_pairs[j].set, &arg_used);
      if (sunretval != CV_SUCCESS)
      {
        retval = CV_ILL_INPUT;
        cvProcessError(cv_mem, retval, __LINE__, __func__, __FILE__,
                       "error setting command-line argument: %s",
                       int_real_pairs[j].key);
        return retval;
      }
      if (arg_used) break;
    }
    if (arg_used) continue;

    /* check all int+long command-line options */
    for (j = 0; j < num_int_long_keys; j++)
    {
      sunretval = sunCheckAndSetIntLongArg(cvode_mem, &idx, argv, offset,
                                        int_long_pairs[j].key,
                                        int_long_pairs[j].set, &arg_used);
      if (sunretval != CV_SUCCESS)
      {
        retval = CV_ILL_INPUT;
        cvProcessError(cv_mem, retval, __LINE__, __func__, __FILE__,
                       "error setting command-line argument: %s",
                       int_long_pairs[j].key);
        return retval;
      }
      if (arg_used) break;
    }
    if (arg_used) continue;

    /* check all int+real+real command-line options */
    for (j = 0; j < num_int_real_real_keys; j++)
    {
      sunretval = sunCheckAndSetIntRealRealArg(cvode_mem, &idx, argv, offset,
                                            int_real_real_pairs[j].key,
                                            int_real_real_pairs[j].set,
                                            &arg_used);
      if (sunretval != CV_SUCCESS)
      {
        retval = CV_ILL_INPUT;
        cvProcessError(cv_mem, retval, __LINE__, __func__, __FILE__,
                       "error setting command-line argument: %s",
                       int_real_real_pairs[j].key);
        return retval;
      }
      if (arg_used) break;
    }
    if (arg_used) continue;

    /* warn for uninterpreted cvid.X arguments */
    cvProcessError(cv_mem, CV_WARNING, __LINE__, __func__, __FILE__,
                   "WARNING: argument %s was not handled\n", argv[idx]);
  }

  return (CV_SUCCESS);
}

/*===============================================================
  EOF
  ===============================================================*/
