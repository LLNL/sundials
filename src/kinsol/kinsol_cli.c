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
 * to KINSOL.
 *--------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sundials/sundials_types.h>
#include "kinsol/kinsol.h"
#include "kinsol/kinsol_ls.h"
#include "kinsol_impl.h"
#include "sundials_cli.h"

/*---------------------------------------------------------------
  KINSetFromCommandLine:

  Parses the command line to control scalar-valued KINSOL options.
  ---------------------------------------------------------------*/

int KINSetFromCommandLine(void* kinmem, const char* kinid, int argc, char* argv[])
{
  KINMem kin_mem;
  if (kinmem == NULL)
  {
    KINProcessError(NULL, KIN_MEM_NULL, __LINE__, __func__, __FILE__, MSG_NO_MEM);
    return (KIN_MEM_NULL);
  }
  kin_mem = (KINMem)kinmem;

  /* Set lists of command-line arguments, and the corresponding set routines */
  static const struct sunKeyIntPair int_pairs[] =
    {{"orth_aa", KINSetOrthAA},
     {"return_newest", KINSetReturnNewest},
     {"no_init_setup", KINSetNoInitSetup},
     {"no_res_mon", KINSetNoResMon},
     {"eta_form", KINSetEtaForm},
     {"no_min_eps", KINSetNoMinEps}};
  static const int num_int_keys = sizeof(int_pairs) / sizeof(*int_pairs);

  static const struct sunKeyLongPair long_pairs[] =
    {{"m_aa", KINSetMAA},
     {"delay_aa", KINSetDelayAA},
     {"num_max_iters", KINSetNumMaxIters},
     {"max_setup_calls", KINSetMaxSetupCalls},
     {"max_sub_setup_calls", KINSetMaxSubSetupCalls},
     {"max_beta_fails", KINSetMaxBetaFails}};
  static const int num_long_keys = sizeof(long_pairs) / sizeof(*long_pairs);

  static const struct sunKeyRealPair real_pairs[] =
    {{"damping", KINSetDamping},
     {"damping_aa", KINSetDampingAA},
     {"eta_const_value", KINSetEtaConstValue},
     {"res_mon_const_value", KINSetResMonConstValue},
     {"max_newton_step", KINSetMaxNewtonStep},
     {"rel_err_func", KINSetRelErrFunc},
     {"func_norm_tol", KINSetFuncNormTol},
     {"scaled_step_tol", KINSetScaledStepTol}};
  static const int num_real_keys = sizeof(real_pairs) / sizeof(*real_pairs);

  static const struct sunKeyTwoRealPair tworeal_pairs[] =
    {{"eta_params", KINSetEtaParams}, {"res_mon_params", KINSetResMonParams}};
  static const int num_tworeal_keys = sizeof(tworeal_pairs) /
                                      sizeof(*tworeal_pairs);

  SUNErrCode sunretval;
  int idx, j, retval;
  for (idx = 1; idx < argc; idx++)
  {
    sunbooleantype arg_used = SUNFALSE;

    /* if kinid is supplied, skip command-line arguments that do not begin with kinid;
       else, skip command-line arguments that do not begin with "kinsol." */
    size_t offset;
    if (kinid != NULL)
    {
      if (strncmp(argv[idx], kinid, strlen(kinid)) != 0) { continue; }
      offset = strlen(kinid) + 1;
    }
    else
    {
      static const char* prefix = "kinsol.";
      if (strncmp(argv[idx], prefix, strlen(prefix)) != 0) { continue; }
      offset = strlen(prefix);
    }

    /* check all "int" command-line options */
    sunretval = sunCheckAndSetIntArgs(kinmem, &idx, argv, offset,
                                     int_pairs, num_int_keys,
                                     &arg_used, &j);
    if (sunretval != SUN_SUCCESS)
    {
      retval = KIN_ILL_INPUT;
      KINProcessError(kin_mem, retval, __LINE__, __func__, __FILE__,
                      "error setting command-line argument: %s",
                      int_pairs[j].key);
      return retval;
    }
    if (arg_used) continue;

    /* check all long int command-line options */
    sunretval = sunCheckAndSetLongArgs(kinmem, &idx, argv, offset,
                                       long_pairs, num_long_keys,
                                       &arg_used, &j);
    if (sunretval != SUN_SUCCESS)
    {
      retval = KIN_ILL_INPUT;
      KINProcessError(kin_mem, retval, __LINE__, __func__, __FILE__,
                      "error setting command-line argument: %s",
                      long_pairs[j].key);
      return retval;
    }
    if (arg_used) continue;

    /* check all real command-line options */
    sunretval = sunCheckAndSetRealArgs(kinmem, &idx, argv, offset,
                                       real_pairs, num_real_keys,
                                       &arg_used, &j);
    if (sunretval != SUN_SUCCESS)
    {
      retval = KIN_ILL_INPUT;
      KINProcessError(kin_mem, retval, __LINE__, __func__, __FILE__,
                      "error setting command-line argument: %s",
                      real_pairs[j].key);
      return retval;
    }
    if (arg_used) continue;

    /* check all pair-of-real command-line options */
    sunretval = sunCheckAndSetTwoRealArgs(kinmem, &idx, argv, offset,
                                          tworeal_pairs, num_tworeal_keys,
                                          &arg_used, &j);
    if (sunretval != SUN_SUCCESS)
    {
      retval = KIN_ILL_INPUT;
      KINProcessError(kin_mem, retval, __LINE__, __func__, __FILE__,
                      "error setting command-line argument: %s",
                      tworeal_pairs[j].key);
      return retval;
    }
    if (arg_used) continue;

    /* warn for uninterpreted kinid.X arguments */
    KINProcessError(kin_mem, KIN_WARNING, __LINE__, __func__, __FILE__,
                    "WARNING: argument %s was not handled\n", argv[idx]);
  }

  return (KIN_SUCCESS);
}

/*===============================================================
  EOF
  ===============================================================*/
