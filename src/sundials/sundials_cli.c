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

#include "sundials_cli.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sundials/sundials_types.h>
#include <sundials/sundials_math.h>

/*===============================================================
  Command-line input utility routines
  ===============================================================*/

SUNErrCode sunCheckAndSetIntArg(void* mem, int* i, char* argv[], const size_t offset,
                         const char* argtest, sunIntSetFn fname,
                         sunbooleantype* arg_used)
{
  *arg_used = SUNFALSE;
  if (strcmp(argv[*i] + offset, argtest) == 0)
  {
    (*i) += 1;
    int iarg   = atoi(argv[*i]);
    int retval = fname(mem, iarg);
    if (retval != SUN_SUCCESS) { return retval; }
    *arg_used = SUNTRUE;
  }
  return SUN_SUCCESS;
}

SUNErrCode sunCheckAndSetTwoIntArg(void* mem, int* i, char* argv[],
                            const size_t offset, const char* argtest,
                            sunTwoIntSetFn fname, sunbooleantype* arg_used)
{
  *arg_used = SUNFALSE;
  if (strcmp(argv[*i] + offset, argtest) == 0)
  {
    (*i) += 1;
    int iarg1 = atoi(argv[*i]);
    (*i) += 1;
    int iarg2  = atoi(argv[*i]);
    int retval = fname(mem, iarg1, iarg2);
    if (retval != SUN_SUCCESS) { return retval; }
    *arg_used = SUNTRUE;
  }
  return SUN_SUCCESS;
}

SUNErrCode sunCheckAndSetLongArg(void* mem, int* i, char* argv[], const size_t offset,
                          const char* argtest, sunLongSetFn fname,
                          sunbooleantype* arg_used)
{
  *arg_used = SUNFALSE;
  if (strcmp(argv[*i] + offset, argtest) == 0)
  {
    (*i) += 1;
    long int iarg = atol(argv[*i]);
    int retval    = fname(mem, iarg);
    if (retval != SUN_SUCCESS) { return retval; }
    *arg_used = SUNTRUE;
  }
  return SUN_SUCCESS;
}

SUNErrCode sunCheckAndSetIntRealArg(void* mem, int* i, char* argv[],
                             const size_t offset, const char* argtest,
                             sunIntRealSetFn fname, sunbooleantype* arg_used)
{
  *arg_used = SUNFALSE;
  if (strcmp(argv[*i] + offset, argtest) == 0)
  {
    (*i) += 1;
    int iarg = atoi(argv[*i]);
    (*i) += 1;
    sunrealtype rarg = SUNStrToReal(argv[*i]);
    int retval       = fname(mem, iarg, rarg);
    if (retval != SUN_SUCCESS) { return retval; }
    *arg_used = SUNTRUE;
  }
  return SUN_SUCCESS;
}

SUNErrCode sunCheckAndSetIntRealRealArg(void* mem, int* i, char* argv[],
                                 const size_t offset, const char* argtest,
                                 sunIntRealRealSetFn fname,
                                 sunbooleantype* arg_used)
{
  *arg_used = SUNFALSE;
  if (strcmp(argv[*i] + offset, argtest) == 0)
  {
    (*i) += 1;
    int iarg = atoi(argv[*i]);
    (*i) += 1;
    sunrealtype rarg1 = SUNStrToReal(argv[*i]);
    (*i) += 1;
    sunrealtype rarg2 = SUNStrToReal(argv[*i]);
    int retval        = fname(mem, iarg, rarg1, rarg2);
    if (retval != SUN_SUCCESS) { return retval; }
    *arg_used = SUNTRUE;
  }
  return SUN_SUCCESS;
}

SUNErrCode sunCheckAndSetIntLongArg(void* mem, int* i, char* argv[],
                             const size_t offset, const char* argtest,
                             sunIntLongSetFn fname, sunbooleantype* arg_used)
{
  *arg_used = SUNFALSE;
  if (strcmp(argv[*i] + offset, argtest) == 0)
  {
    (*i) += 1;
    int iarg = atoi(argv[*i]);
    (*i) += 1;
    long int large = atol(argv[*i]);
    int retval     = fname(mem, iarg, large);
    if (retval != SUN_SUCCESS) { return retval; }
    *arg_used = SUNTRUE;
  }
  return SUN_SUCCESS;
}

SUNErrCode sunCheckAndSetRealArg(void* mem, int* i, char* argv[], const size_t offset,
                          const char* argtest, sunRealSetFn fname,
                          sunbooleantype* arg_used)
{
  *arg_used = SUNFALSE;
  if (strcmp(argv[*i] + offset, argtest) == 0)
  {
    (*i) += 1;
    sunrealtype rarg = SUNStrToReal(argv[*i]);
    int retval       = fname(mem, rarg);
    if (retval != SUN_SUCCESS) { return retval; }
    *arg_used = SUNTRUE;
  }
  return SUN_SUCCESS;
}

SUNErrCode sunCheckAndSetTwoRealArg(void* mem, int* i, char* argv[],
                             const size_t offset, const char* argtest,
                             sunTwoRealSetFn fname, sunbooleantype* arg_used)
{
  *arg_used = SUNFALSE;
  if (strcmp(argv[*i] + offset, argtest) == 0)
  {
    (*i) += 1;
    sunrealtype rarg1 = SUNStrToReal(argv[*i]);
    (*i) += 1;
    sunrealtype rarg2 = SUNStrToReal(argv[*i]);
    int retval        = fname(mem, rarg1, rarg2);
    if (retval != SUN_SUCCESS) { return retval; }
    *arg_used = SUNTRUE;
  }
  return SUN_SUCCESS;
}

SUNErrCode sunCheckAndSetCharArg(void* mem, int* i, char* argv[], const size_t offset,
                          const char* argtest, sunCharSetFn fname,
                          sunbooleantype* arg_used)
{
  *arg_used = SUNFALSE;
  if (strcmp(argv[*i] + offset, argtest) == 0)
  {
    (*i) += 1;
    int retval = fname(mem, argv[*i]);
    if (retval != SUN_SUCCESS) { return retval; }
    *arg_used = SUNTRUE;
  }
  return SUN_SUCCESS;
}

SUNErrCode sunCheckAndSetTwoCharArg(void* mem, int* i, char* argv[],
                             const size_t offset, const char* argtest,
                             sunTwoCharSetFn fname, sunbooleantype* arg_used)
{
  *arg_used = SUNFALSE;
  if (strcmp(argv[*i] + offset, argtest) == 0)
  {
    int retval = fname(mem, argv[*i + 1], argv[*i + 2]);
    (*i) += 2;
    if (retval != SUN_SUCCESS) { return retval; }
    *arg_used = SUNTRUE;
  }
  return SUN_SUCCESS;
}

SUNErrCode sunCheckAndSetActionArg(void* mem, int* i, char* argv[],
                            const size_t offset, const char* argtest,
                            sunActionSetFn fname, sunbooleantype* arg_used)
{
  *arg_used = SUNFALSE;
  if (strcmp(argv[*i] + offset, argtest) == 0)
  {
    int retval = fname(mem);
    if (retval != SUN_SUCCESS) { return retval; }
    *arg_used = SUNTRUE;
  }
  return SUN_SUCCESS;
}
