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
#include "arkode_impl.h"

/*===============================================================
  Command-line input utility routines -- create more of these: TwoInt, TwoReal, String, TwoString.  Remove direct calls to arkProcessError, and instead have the calling routine return a slightly more generic error message that omits the extra arguments.  Propagate both this set of routines and the function types to the SUNDIALS level, so that they can be used across all integrators.
  ===============================================================*/

int arkCheckAndSetIntArg(ARKodeMem ark_mem, int* i, char* argv[],
                         const size_t offset, const char* argtest,
                         arkIntSetFn fname, sunbooleantype* arg_used)
{
  *arg_used = SUNFALSE;
  if (strcmp(argv[*i] + offset, argtest) == 0)
  {
    (*i) += 1;
    int iarg   = atoi(argv[*i]);
    int retval = fname((void*)ark_mem, iarg);
    if (retval != ARK_SUCCESS)
    {
      arkProcessError(ark_mem, retval, __LINE__, __func__, __FILE__,
                      "error setting command-line argument: %s %s",
                      argv[(*i) - 1], argv[*i]);
      return retval;
    }
    *arg_used = SUNTRUE;
  }
  return ARK_SUCCESS;
}

int arkCheckAndSetLongArg(ARKodeMem ark_mem, int* i, char* argv[],
                          const size_t offset, const char* argtest,
                          arkLongSetFn fname, sunbooleantype* arg_used)
{
  *arg_used = SUNFALSE;
  if (strcmp(argv[*i] + offset, argtest) == 0)
  {
    (*i) += 1;
    long int iarg = atol(argv[*i]);
    int retval    = fname((void*)ark_mem, iarg);
    if (retval != ARK_SUCCESS)
    {
      arkProcessError(ark_mem, retval, __LINE__, __func__, __FILE__,
                      "error setting command-line argument: %s %s",
                      argv[(*i) - 1], argv[*i]);
      return retval;
    }
    *arg_used = SUNTRUE;
  }
  return ARK_SUCCESS;
}

int arkCheckAndSetRealArg(ARKodeMem ark_mem, int* i, char* argv[],
                          const size_t offset, const char* argtest,
                          arkRealSetFn fname, sunbooleantype* arg_used)
{
  *arg_used = SUNFALSE;
  if (strcmp(argv[*i] + offset, argtest) == 0)
  {
    (*i) += 1;
    sunrealtype rarg = atof(argv[*i]);
    int retval       = fname((void*)ark_mem, rarg);
    if (retval != ARK_SUCCESS)
    {
      arkProcessError(ark_mem, retval, __LINE__, __func__, __FILE__,
                      "error setting command-line argument: %s %s",
                      argv[(*i) - 1], argv[*i]);
      return retval;
    }
    *arg_used = SUNTRUE;
  }
  return ARK_SUCCESS;
}

int arkCheckAndSetActionArg(ARKodeMem ark_mem, int* i, char* argv[],
                            const size_t offset, const char* argtest,
                            arkActionSetFn fname, sunbooleantype* arg_used)
{
  *arg_used = SUNFALSE;
  if (strcmp(argv[*i] + offset, argtest) == 0)
  {
    int retval = fname((void*)ark_mem);
    if (retval != ARK_SUCCESS)
    {
      arkProcessError(ark_mem, retval, __LINE__, __func__, __FILE__,
                      "error setting command-line argument: %s", argv[*i]);
      return retval;
    }
    *arg_used = SUNTRUE;
  }
  return ARK_SUCCESS;
}


/*---------------------------------------------------------------
  ARKodeSetFromCommandLine:

  Parses the command line to control scalar-valued ARKODE options.

  (this leverages a multiple typedefs and static utility routines)
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

  /* Set list of integer command-line arguments, and the corresponding set routine */
  static struct arkKeyIntPair int_pairs[] =
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
     {"max_conv_fails", ARKodeSetMaxConvFails}};
  static const int num_int_keys = sizeof(int_pairs) / sizeof(*int_pairs);

  static struct arkKeyLongPair long_pairs[] = {
    {"max_num_steps", ARKodeSetMaxNumSteps}};
  static const int num_long_keys = sizeof(long_pairs) / sizeof(*long_pairs);

  static struct arkKeyRealPair real_pairs[] =
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
     {"max_cfail_growth", ARKodeSetMaxCFailGrowth}};
  static const int num_real_keys = sizeof(real_pairs) / sizeof(*real_pairs);

  static struct arkKeyActionPair action_pairs[] =
    {{"nonlinear", ARKodeSetNonlinear},
     {"clear_stop_time", ARKodeClearStopTime},
     {"no_inactive_root_warn", ARKodeSetNoInactiveRootWarn},
     {"reset_accumulated_error", ARKodeResetAccumulatedError}};
  static const int num_action_keys = sizeof(action_pairs) / sizeof(*action_pairs);

  int i, j, retval;
  sunbooleantype write_parameters = SUNFALSE;
  for (i = 1; i < argc; i++)
  {
    sunbooleantype arg_used = SUNFALSE;

    /* if arkid is supplied, skip command-line arguments that do not begin with arkid;
       else, skip command-line arguments that do not begin with "arkode." */
    size_t offset;
    if (strlen(arkid) > 0) {
      if (strncmp(argv[i], arkid, strlen(arkid)) != 0) { continue; }
      offset = strlen(arkid) + 1;
    }
    else
    {
      static const char* prefix = "arkode.";
      if (strncmp(argv[i], prefix, strlen(prefix)) != 0) { continue; }
      offset = strlen(prefix);
    }

    /* check all "int" command-line options */
    for (j = 0; j < num_int_keys; j++)
    {
      retval = arkCheckAndSetIntArg(ark_mem, &i, argv, offset, int_pairs[j].key,
                                    int_pairs[j].set, &arg_used);
      if (retval != ARK_SUCCESS) { return retval; }
      if (arg_used) break;
    }
    if (arg_used) continue;

    /* check all long int command-line options */
    for (j = 0; j < num_long_keys; j++)
    {
      retval = arkCheckAndSetLongArg(ark_mem, &i, argv, offset, long_pairs[j].key,
                                     long_pairs[j].set, &arg_used);
      if (retval != ARK_SUCCESS) { return retval; }
      if (arg_used) break;
    }
    if (arg_used) continue;

    /* check all real command-line options */
    for (j = 0; j < num_real_keys; j++)
    {
      retval = arkCheckAndSetRealArg(ark_mem, &i, argv, offset, real_pairs[j].key,
                                     real_pairs[j].set, &arg_used);
      if (retval != ARK_SUCCESS) { return retval; }
      if (arg_used) break;
    }
    if (arg_used) continue;

    /* check all action command-line options */
    for (j = 0; j < num_action_keys; j++)
    {
      retval = arkCheckAndSetActionArg(ark_mem, &i, argv, offset,
                                       action_pairs[j].key, action_pairs[j].set,
        &arg_used);
      if (retval != ARK_SUCCESS) { return retval; }
      if (arg_used) break;
    }
    if (arg_used) continue;

    /*** handle all remaining command-line options ***/

    if (strcmp(argv[i] + offset, "interpolant_type") == 0)
    {
      i++;
      retval = ARK_ILL_INPUT;
      if (strcmp(argv[i], "ARK_INTERP_HERMITE") == 0)
      {
        retval = ARKodeSetInterpolantType(arkode_mem, ARK_INTERP_HERMITE);
      }
      else if (strcmp(argv[i], "ARK_INTERP_LAGRANGE") == 0)
      {
        retval = ARKodeSetInterpolantType(arkode_mem, ARK_INTERP_LAGRANGE);
      }
      else if (strcmp(argv[i], "ARK_INTERP_NONE") == 0)
      {
        retval = ARKodeSetInterpolantType(arkode_mem, ARK_INTERP_NONE);
      }
      if (retval != ARK_SUCCESS)
      {
        arkProcessError(ark_mem, retval, __LINE__, __func__, __FILE__,
                        "error setting command-line argument: %s %s",
                        argv[i - 1], argv[i]);
        return retval;
      }
      arg_used = SUNTRUE;
      continue;
    }

    if (strcmp(argv[i] + offset, "scalar_tolerances") == 0)
    {
      i++;
      sunrealtype rtol = atof(argv[i]);
      i++;
      sunrealtype atol = atof(argv[i]);
      retval           = ARKodeSStolerances(arkode_mem, rtol, atol);
      if (retval != ARK_SUCCESS)
      {
        arkProcessError(ark_mem, retval, __LINE__, __func__, __FILE__,
                        "error setting command-line argument: %s %s %s",
                        argv[i - 2], argv[i - 1], argv[i]);
        return retval;
      }
      arg_used = SUNTRUE;
      continue;
    }

    if (strcmp(argv[i] + offset, "fixed_step_bounds") == 0)
    {
      i++;
      sunrealtype lb = atof(argv[i]);
      i++;
      sunrealtype ub = atof(argv[i]);
      retval         = ARKodeSetFixedStepBounds(arkode_mem, lb, ub);
      if (retval != ARK_SUCCESS)
      {
        arkProcessError(ark_mem, retval, __LINE__, __func__, __FILE__,
                        "error setting command-line argument: %s %s %s",
                        argv[i - 2], argv[i - 1], argv[i]);
        return retval;
      }
      arg_used = SUNTRUE;
      continue;
    }

    if (strcmp(argv[i] + offset, "accum_error_type") == 0)
    {
      i++;
      retval = ARK_ILL_INPUT;
      if (strcmp(argv[i], "ARK_ACCUMERROR_NONE") == 0)
      {
        retval = ARKodeSetAccumulatedErrorType(arkode_mem, ARK_ACCUMERROR_NONE);
      }
      else if (strcmp(argv[i], "ARK_ACCUMERROR_MAX") == 0)
      {
        retval = ARKodeSetAccumulatedErrorType(arkode_mem, ARK_ACCUMERROR_MAX);
      }
      else if (strcmp(argv[i], "ARK_ACCUMERROR_SUM") == 0)
      {
        retval = ARKodeSetAccumulatedErrorType(arkode_mem, ARK_ACCUMERROR_SUM);
      }
      else if (strcmp(argv[i], "ARK_ACCUMERROR_AVG") == 0)
      {
        retval = ARKodeSetAccumulatedErrorType(arkode_mem, ARK_ACCUMERROR_AVG);
      }
      if (retval != ARK_SUCCESS)
      {
        arkProcessError(ark_mem, retval, __LINE__, __func__, __FILE__,
                        "error setting command-line argument: %s %s",
                        argv[i - 1], argv[i]);
        return retval;
      }
      arg_used = SUNTRUE;
      continue;
    }

    if (strcmp(argv[i] + offset, "writeparameters") == 0)
    {
      write_parameters = SUNTRUE;
      arg_used = SUNTRUE;
      continue;
    }

    /* Call stepper-specific SetFromCommandLine routine (if supplied) to
       process this command-line argument */
    if (ark_mem->step_setfromcommandline)
    {
      retval = ark_mem->step_setfromcommandline(ark_mem, &i, argv,
                                                offset, &arg_used);
      if (retval != ARK_SUCCESS)
      {
        arkProcessError(ark_mem, retval, __LINE__, __func__, __FILE__,
                        "error setting command-line argument: %s", argv[i]);
        return retval;
      }
      if (arg_used) { continue; }
    }

    /* warn for uninterpreted arkid.X arguments */
    arkProcessError(ark_mem, ARK_WARNING, __LINE__, __func__, __FILE__,
                    "WARNING: argument %s was not handled\n", argv[i]);
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
