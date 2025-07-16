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
#include "ida/ida.h"
#include "ida/ida_ls.h"
#include "ida_impl.h"
#include "sundials_cli.h"

static int idaSetFromCommandLine(void* ida_mem, const char* idaid,
                                 int argc, char* argv[]);

/*---------------------------------------------------------------
  IDASetOptions:

  Sets IDA options using strings.
  ---------------------------------------------------------------*/

int IDASetOptions(void* ida_mem, const char* idaid,
                  const char* file_name, int argc, char* argv[])
{
  if (file_name != NULL && strlen(file_name) > 0)
  {
      int retval = IDA_ILL_INPUT;
      IDAProcessError(ida_mem, retval, __LINE__, __func__, __FILE__,
                      "file-based options are not currently supported.");
      return retval;
  }

  if (argc > 0 && argv != NULL)
  {
    int retval = idaSetFromCommandLine(ida_mem, idaid, argc, argv);
    if (retval != IDA_SUCCESS) { return retval; }
  }

  return IDA_SUCCESS;
}

/* -----------------------------------------------------------------
 * Function to control IDA options from the command line
 */

static int idaSetFromCommandLine(void* ida_mem, const char* idaid,
                                 int argc, char* argv[])
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
     {"line_search_off_ic", IDASetLineSearchOffIC},
     {"max_backs_ic", IDASetMaxBacksIC},
     {"max_order", IDASetMaxOrd},
     {"max_err_test_fails", IDASetMaxErrTestFails},
     {"suppress_alg", IDASetSuppressAlg},
     {"max_conv_fails", IDASetMaxConvFails},
     {"max_nonlin_iters", IDASetMaxNonlinIters},
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
     {"scalar_tolerances", IDASStolerances}};
  static const int num_tworeal_keys = sizeof(tworeal_pairs) /
                                      sizeof(*tworeal_pairs);

  static const struct sunKeyActionPair action_pairs[] =
    {{"clear_stop_time", IDAClearStopTime},
     {"no_inactive_root_warn", IDASetNoInactiveRootWarn}};
  static const int num_action_keys = sizeof(action_pairs) / sizeof(*action_pairs);

  SUNErrCode sunretval;
  int idx, j, retval;
  for (idx = 1; idx < argc; idx++)
  {
    sunbooleantype arg_used = SUNFALSE;

    /* if idaid is supplied, skip command-line arguments that do not begin with idaid;
       else, skip command-line arguments that do not begin with "ida." */
    size_t offset;
    if (idaid != NULL && strlen(idaid) > 0)
    {
      if (strncmp(argv[idx], idaid, strlen(idaid)) != 0) { continue; }
      offset = strlen(idaid) + 1;
    }
    else
    {
      static const char* prefix = "ida.";
      if (strncmp(argv[idx], prefix, strlen(prefix)) != 0) { continue; }
      offset = strlen(prefix);
    }

    /* check all "int" command-line options */
    sunretval = sunCheckAndSetIntArgs(ida_mem, &idx, argv, offset, int_pairs,
                                      num_int_keys, &arg_used, &j);
    if (sunretval != SUN_SUCCESS)
    {
      retval = IDA_ILL_INPUT;
      IDAProcessError(IDA_mem, retval, __LINE__, __func__, __FILE__,
                      "error setting key: %s",
                      int_pairs[j].key);
      return retval;
    }
    if (arg_used) continue;

    /* check all long int command-line options */
    sunretval = sunCheckAndSetLongArgs(ida_mem, &idx, argv, offset, long_pairs,
                                       num_long_keys, &arg_used, &j);
    if (sunretval != SUN_SUCCESS)
    {
      retval = IDA_ILL_INPUT;
      IDAProcessError(IDA_mem, retval, __LINE__, __func__, __FILE__,
                      "error setting key: %s",
                      long_pairs[j].key);
      return retval;
    }
    if (arg_used) continue;

    /* check all real command-line options */
    sunretval = sunCheckAndSetRealArgs(ida_mem, &idx, argv, offset, real_pairs,
                                       num_real_keys, &arg_used, &j);
    if (sunretval != SUN_SUCCESS)
    {
      retval = IDA_ILL_INPUT;
      IDAProcessError(IDA_mem, retval, __LINE__, __func__, __FILE__,
                      "error setting key: %s",
                      real_pairs[j].key);
      return retval;
    }
    if (arg_used) continue;

    /* check all pair-of-real command-line options */
    sunretval = sunCheckAndSetTwoRealArgs(ida_mem, &idx, argv, offset,
                                          tworeal_pairs, num_tworeal_keys,
                                          &arg_used, &j);
    if (sunretval != SUN_SUCCESS)
    {
      retval = IDA_ILL_INPUT;
      IDAProcessError(IDA_mem, retval, __LINE__, __func__, __FILE__,
                      "error setting key: %s",
                      tworeal_pairs[j].key);
      return retval;
    }
    if (arg_used) continue;

    /* check all action command-line options */
    sunretval = sunCheckAndSetActionArgs(ida_mem, &idx, argv, offset,
                                         action_pairs, num_action_keys,
                                         &arg_used, &j);
    if (sunretval != SUN_SUCCESS)
    {
      retval = IDA_ILL_INPUT;
      IDAProcessError(IDA_mem, retval, __LINE__, __func__, __FILE__,
                      "error setting key: %s",
                      action_pairs[j].key);
      return retval;
    }
    if (arg_used) continue;

    /* warn for uninterpreted idaid.X arguments */
    IDAProcessError(IDA_mem, IDA_WARNING, __LINE__, __func__, __FILE__,
                    "WARNING: key %s was not handled\n", argv[idx]);
  }

  return (IDA_SUCCESS);
}

/*===============================================================
  EOF
  ===============================================================*/
