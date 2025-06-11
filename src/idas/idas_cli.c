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
 * to IDA.
 *--------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sundials/sundials_types.h>
#include "idas/idas.h"
#include "idas/idas_ls.h"
#include "idas_impl.h"
#include "sundials_cli.h"

/*---------------------------------------------------------------
  IDASetFromCommandLine:

  Parses the command line to control scalar-valued IDA options.
  ---------------------------------------------------------------*/

int IDASetFromCommandLine(void* ida_mem, const char* idaid, int argc,
                          char* argv[])
{
  IDAMem IDA_mem;
  if (ida_mem == NULL)
  {
    IDAProcessError(NULL, IDA_MEM_NULL, __LINE__, __func__, __FILE__, MSG_NO_MEM);
    return (IDA_MEM_NULL);
  }
  IDA_mem = (IDAMem)ida_mem;

  /* Set lists of command-line arguments, and the corresponding set routines */
  static const struct sunKeyIntPair int_pairs[] =
    {{"max_num_steps_ic", IDASetMaxNumStepsIC},
     {"max_num_jacs_ic", IDASetMaxNumJacsIC},
     {"max_num_iters_ic", IDASetMaxNumItersIC},
     {"linesearch_off_ic", IDASetLineSearchOffIC},
     {"max_backs_ic", IDASetMaxBacksIC},
     {"max_order", IDASetMaxOrd},
     {"max_err_test_fails", IDASetMaxErrTestFails},
     {"suppress_alg", IDASetSuppressAlg},
     {"max_conv_fails", IDASetMaxConvFails},
     {"max_nonlin_iters", IDASetMaxNonlinIters},
     {"quad_err_con", IDASetQuadErrCon},
     {"sens_err_con", IDASetSensErrCon},
     {"sens_max_nonlin_iters", IDASetSensMaxNonlinIters},
     {"quad_sens_err_con", IDASetQuadSensErrCon},
     {"linear_solution_scaling", IDASetLinearSolutionScaling}};
  static const int num_int_keys = sizeof(int_pairs) / sizeof(*int_pairs);

  static const struct sunKeyLongPair long_pairs[] = {
    {"max_num_steps", IDASetMaxNumSteps}};
  static const int num_long_keys = sizeof(long_pairs) / sizeof(*long_pairs);

  static const struct sunKeyRealPair real_pairs[] =
    {{"nonlin_conv_coef_ic", IDASetNonlinConvCoefIC},
     {"step_tolerance_ic", IDASetStepToleranceIC},
     {"delta_cj_lsetup", IDASetDeltaCjLSetup},
     {"init_step", IDASetInitStep},
     {"max_step", IDASetMaxStep},
     {"min_step", IDASetMinStep},
     {"stop_time", IDASetStopTime},
     {"eta_min", IDASetEtaMin},
     {"eta_max", IDASetEtaMax},
     {"eta_low", IDASetEtaLow},
     {"eta_min_err_fail", IDASetEtaMinErrFail},
     {"eta_conv_fail", IDASetEtaConvFail},
     {"nonlin_conv_coef", IDASetNonlinConvCoef},
     {"eps_lin", IDASetEpsLin},
     {"ls_norm_factor", IDASetLSNormFactor},
     {"increment_factor", IDASetIncrementFactor}};
  static const int num_real_keys = sizeof(real_pairs) / sizeof(*real_pairs);

  static const struct sunKeyTwoRealPair tworeal_pairs[] =
    {{"eta_fixed_step_bounds", IDASetEtaFixedStepBounds},
     {"scalar_tolerances", IDASStolerances},
     {"quad_scalar_tolerances", IDAQuadSStolerances}};
  static const int num_tworeal_keys = sizeof(tworeal_pairs) /
                                      sizeof(*tworeal_pairs);

  static const struct sunKeyTwoIntPair twoint_pairs[] =
    {{"max_order_b", IDASetMaxOrdB},
     {"suppress_alg_b", IDASetSuppressAlgB},
     {"quad_err_con_b", IDASetQuadErrConB},
     {"linear_solution_scaling_b", IDASetLinearSolutionScalingB}};
  static const int num_twoint_keys = sizeof(twoint_pairs) / sizeof(*twoint_pairs);

  static const struct sunKeyActionPair action_pairs[] =
    {{"clear_stop_time", IDAClearStopTime},
     {"no_inactive_root_warn", IDASetNoInactiveRootWarn},
     {"sens_toggle_off", IDASensToggleOff},
     {"adj_no_sensi", IDAAdjSetNoSensi}};
  static const int num_action_keys = sizeof(action_pairs) / sizeof(*action_pairs);

  static const struct sunKeyIntRealPair int_real_pairs[] =
    {{"sens_dq_method", IDASetSensDQMethod},
     {"init_step_b", IDASetInitStepB},
     {"max_step_b", IDASetMaxStepB},
     {"eps_lin_b", IDASetEpsLinB},
     {"ls_norm_factor_b", IDASetLSNormFactorB},
     {"increment_factor_b", IDASetIncrementFactorB}};
  static const int num_int_real_keys = sizeof(int_real_pairs) /
                                       sizeof(*int_real_pairs);

  static const struct sunKeyIntRealRealPair int_real_real_pairs[] =
    {{"scalar_tolerances_b", IDASStolerancesB},
     {"quad_scalar_tolerances_b", IDAQuadSStolerancesB}};
  static const int num_int_real_real_keys = sizeof(int_real_real_pairs) /
                                            sizeof(*int_real_real_pairs);

  static const struct sunKeyIntLongPair int_long_pairs[] = {
    {"max_num_steps_b", IDASetMaxNumStepsB}};
  static const int num_int_long_keys = sizeof(int_long_pairs) /
                                       sizeof(*int_long_pairs);

  SUNErrCode sunretval;
  int idx, j, retval;
  for (idx = 1; idx < argc; idx++)
  {
    sunbooleantype arg_used = SUNFALSE;

    /* if idaid is supplied, skip command-line arguments that do not begin with idaid;
       else, skip command-line arguments that do not begin with "idas." */
    size_t offset;
    if (idaid != NULL)
    {
      if (strncmp(argv[idx], idaid, strlen(idaid)) != 0) { continue; }
      offset = strlen(idaid) + 1;
    }
    else
    {
      static const char* prefix = "idas.";
      if (strncmp(argv[idx], prefix, strlen(prefix)) != 0) { continue; }
      offset = strlen(prefix);
    }

    /* check all "int" command-line options */
    for (j = 0; j < num_int_keys; j++)
    {
      sunretval = sunCheckAndSetIntArg(ida_mem, &idx, argv, offset,
                                       int_pairs[j].key, int_pairs[j].set,
                                       &arg_used);
      if (sunretval != IDA_SUCCESS)
      {
        retval = IDA_ILL_INPUT;
        IDAProcessError(IDA_mem, retval, __LINE__, __func__, __FILE__,
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
      sunretval = sunCheckAndSetLongArg(ida_mem, &idx, argv, offset,
                                        long_pairs[j].key, long_pairs[j].set,
                                        &arg_used);
      if (sunretval != IDA_SUCCESS)
      {
        retval = IDA_ILL_INPUT;
        IDAProcessError(IDA_mem, retval, __LINE__, __func__, __FILE__,
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
      sunretval = sunCheckAndSetRealArg(ida_mem, &idx, argv, offset,
                                        real_pairs[j].key, real_pairs[j].set,
                                        &arg_used);
      if (sunretval != IDA_SUCCESS)
      {
        retval = IDA_ILL_INPUT;
        IDAProcessError(IDA_mem, retval, __LINE__, __func__, __FILE__,
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
      sunretval = sunCheckAndSetTwoIntArg(ida_mem, &idx, argv, offset,
                                          twoint_pairs[j].key,
                                          twoint_pairs[j].set, &arg_used);
      if (sunretval != IDA_SUCCESS)
      {
        retval = IDA_ILL_INPUT;
        IDAProcessError(IDA_mem, retval, __LINE__, __func__, __FILE__,
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
      sunretval = sunCheckAndSetTwoRealArg(ida_mem, &idx, argv, offset,
                                           tworeal_pairs[j].key,
                                           tworeal_pairs[j].set, &arg_used);
      if (sunretval != IDA_SUCCESS)
      {
        retval = IDA_ILL_INPUT;
        IDAProcessError(IDA_mem, retval, __LINE__, __func__, __FILE__,
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
      sunretval = sunCheckAndSetActionArg(ida_mem, &idx, argv, offset,
                                          action_pairs[j].key,
                                          action_pairs[j].set, &arg_used);
      if (sunretval != IDA_SUCCESS)
      {
        retval = IDA_ILL_INPUT;
        IDAProcessError(IDA_mem, retval, __LINE__, __func__, __FILE__,
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
      sunretval = sunCheckAndSetIntRealArg(ida_mem, &idx, argv, offset,
                                           int_real_pairs[j].key,
                                           int_real_pairs[j].set, &arg_used);
      if (sunretval != IDA_SUCCESS)
      {
        retval = IDA_ILL_INPUT;
        IDAProcessError(IDA_mem, retval, __LINE__, __func__, __FILE__,
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
      sunretval = sunCheckAndSetIntLongArg(ida_mem, &idx, argv, offset,
                                           int_long_pairs[j].key,
                                           int_long_pairs[j].set, &arg_used);
      if (sunretval != IDA_SUCCESS)
      {
        retval = IDA_ILL_INPUT;
        IDAProcessError(IDA_mem, retval, __LINE__, __func__, __FILE__,
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
      sunretval = sunCheckAndSetIntRealRealArg(ida_mem, &idx, argv, offset,
                                               int_real_real_pairs[j].key,
                                               int_real_real_pairs[j].set,
                                               &arg_used);
      if (sunretval != IDA_SUCCESS)
      {
        retval = IDA_ILL_INPUT;
        IDAProcessError(IDA_mem, retval, __LINE__, __func__, __FILE__,
                        "error setting command-line argument: %s",
                        int_real_real_pairs[j].key);
        return retval;
      }
      if (arg_used) break;
    }
    if (arg_used) continue;

    /* warn for uninterpreted idaid.X arguments */
    IDAProcessError(IDA_mem, IDA_WARNING, __LINE__, __func__, __FILE__,
                    "WARNING: argument %s was not handled\n", argv[idx]);
  }

  return (IDA_SUCCESS);
}

/*===============================================================
  EOF
  ===============================================================*/
