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
#include <sundials/sundials_math.h>
#include <sundials/sundials_types.h>

/*===============================================================
  Command-line input utility routines
  ===============================================================*/

SUNErrCode sunCheckAndSetIntArg(void* mem, int* argidx, char* argv[],
                                const size_t offset, const char* argtest,
                                sunIntSetFn fname, sunbooleantype* arg_used)
{
  *arg_used = SUNFALSE;
  if (strcmp(argv[*argidx] + offset, argtest) == 0)
  {
    (*argidx) += 1;
    int iarg   = atoi(argv[*argidx]);
    int retval = fname(mem, iarg);
    if (retval != SUN_SUCCESS) { return retval; }
    *arg_used = SUNTRUE;
  }
  return SUN_SUCCESS;
}

SUNErrCode sunCheckAndSetIntArgs(void* mem, int* argidx,
                                 char* argv[], const size_t offset,
                                 const struct sunKeyIntPair* testpairs,
                                 int numpairs, sunbooleantype* arg_used,
                                 int *failedarg)
{
  for (int j = 0; j < numpairs; j++)
  {
    int retval = sunCheckAndSetIntArg(mem, argidx, argv, offset,
                                      testpairs[j].key, testpairs[j].set,
                                      arg_used);
    if (retval != SUN_SUCCESS)
    {
      *failedarg = j;
      return retval;
    }
    if (*arg_used) { return SUN_SUCCESS; }
  }
  return SUN_SUCCESS;
}

SUNErrCode sunCheckAndSetTwoIntArg(void* mem, int* argidx, char* argv[],
                                   const size_t offset, const char* argtest,
                                   sunTwoIntSetFn fname, sunbooleantype* arg_used)
{
  *arg_used = SUNFALSE;
  if (strcmp(argv[*argidx] + offset, argtest) == 0)
  {
    (*argidx) += 1;
    int iarg1 = atoi(argv[*argidx]);
    (*argidx) += 1;
    int iarg2  = atoi(argv[*argidx]);
    int retval = fname(mem, iarg1, iarg2);
    if (retval != SUN_SUCCESS) { return retval; }
    *arg_used = SUNTRUE;
  }
  return SUN_SUCCESS;
}

SUNErrCode sunCheckAndSetTwoIntArgs(void* mem, int* argidx,
                                    char* argv[], const size_t offset,
                                    const struct sunKeyTwoIntPair* testpairs,
                                    int numpairs, sunbooleantype* arg_used,
                                    int *failedarg)
{
  for (int j = 0; j < numpairs; j++)
  {
    int retval = sunCheckAndSetTwoIntArg(mem, argidx, argv, offset,
                                         testpairs[j].key, testpairs[j].set,
                                         arg_used);
    if (retval != SUN_SUCCESS)
    {
      *failedarg = j;
      return retval;
    }
    if (*arg_used) { return SUN_SUCCESS; }
  }
  return SUN_SUCCESS;
}

SUNErrCode sunCheckAndSetLongArg(void* mem, int* argidx, char* argv[],
                                 const size_t offset, const char* argtest,
                                 sunLongSetFn fname, sunbooleantype* arg_used)
{
  *arg_used = SUNFALSE;
  if (strcmp(argv[*argidx] + offset, argtest) == 0)
  {
    (*argidx) += 1;
    long int iarg = atol(argv[*argidx]);
    int retval    = fname(mem, iarg);
    if (retval != SUN_SUCCESS) { return retval; }
    *arg_used = SUNTRUE;
  }
  return SUN_SUCCESS;
}

SUNErrCode sunCheckAndSetLongArgs(void* mem, int* argidx,
                                  char* argv[], const size_t offset,
                                  const struct sunKeyLongPair* testpairs,
                                  int numpairs, sunbooleantype* arg_used,
                                  int *failedarg)
{
  for (int j = 0; j < numpairs; j++)
  {
    int retval = sunCheckAndSetLongArg(mem, argidx, argv, offset,
                                       testpairs[j].key, testpairs[j].set,
                                       arg_used);
    if (retval != SUN_SUCCESS)
    {
      *failedarg = j;
      return retval;
    }
    if (*arg_used) { return SUN_SUCCESS; }
  }
  return SUN_SUCCESS;
}

SUNErrCode sunCheckAndSetIntRealArg(void* mem, int* argidx, char* argv[],
                                    const size_t offset, const char* argtest,
                                    sunIntRealSetFn fname,
                                    sunbooleantype* arg_used)
{
  *arg_used = SUNFALSE;
  if (strcmp(argv[*argidx] + offset, argtest) == 0)
  {
    (*argidx) += 1;
    int iarg = atoi(argv[*argidx]);
    (*argidx) += 1;
    sunrealtype rarg = SUNStrToReal(argv[*argidx]);
    int retval       = fname(mem, iarg, rarg);
    if (retval != SUN_SUCCESS) { return retval; }
    *arg_used = SUNTRUE;
  }
  return SUN_SUCCESS;
}

SUNErrCode sunCheckAndSetIntRealArgs(void* mem, int* argidx,
                                     char* argv[], const size_t offset,
                                     const struct sunKeyIntRealPair* testpairs,
                                     int numpairs, sunbooleantype* arg_used,
                                     int *failedarg)
{
  for (int j = 0; j < numpairs; j++)
  {
    int retval = sunCheckAndSetIntRealArg(mem, argidx, argv, offset,
                                          testpairs[j].key, testpairs[j].set,
                                          arg_used);
    if (retval != SUN_SUCCESS)
    {
      *failedarg = j;
      return retval;
    }
    if (*arg_used) { return SUN_SUCCESS; }
  }
  return SUN_SUCCESS;
}

SUNErrCode sunCheckAndSetIntRealRealArg(void* mem, int* argidx, char* argv[],
                                        const size_t offset, const char* argtest,
                                        sunIntRealRealSetFn fname,
                                        sunbooleantype* arg_used)
{
  *arg_used = SUNFALSE;
  if (strcmp(argv[*argidx] + offset, argtest) == 0)
  {
    (*argidx) += 1;
    int iarg = atoi(argv[*argidx]);
    (*argidx) += 1;
    sunrealtype rarg1 = SUNStrToReal(argv[*argidx]);
    (*argidx) += 1;
    sunrealtype rarg2 = SUNStrToReal(argv[*argidx]);
    int retval        = fname(mem, iarg, rarg1, rarg2);
    if (retval != SUN_SUCCESS) { return retval; }
    *arg_used = SUNTRUE;
  }
  return SUN_SUCCESS;
}

SUNErrCode sunCheckAndSetIntRealRealArgs(void* mem, int* argidx,
                                         char* argv[], const size_t offset,
                                         const struct sunKeyIntRealRealPair* testpairs,
                                         int numpairs, sunbooleantype* arg_used,
                                         int *failedarg)
{
  for (int j = 0; j < numpairs; j++)
  {
    int retval = sunCheckAndSetIntRealRealArg(mem, argidx, argv, offset,
                                              testpairs[j].key, testpairs[j].set,
                                              arg_used);
    if (retval != SUN_SUCCESS)
    {
      *failedarg = j;
      return retval;
    }
    if (*arg_used) { return SUN_SUCCESS; }
  }
  return SUN_SUCCESS;
}

SUNErrCode sunCheckAndSetIntLongArg(void* mem, int* argidx, char* argv[],
                                    const size_t offset, const char* argtest,
                                    sunIntLongSetFn fname,
                                    sunbooleantype* arg_used)
{
  *arg_used = SUNFALSE;
  if (strcmp(argv[*argidx] + offset, argtest) == 0)
  {
    (*argidx) += 1;
    int iarg = atoi(argv[*argidx]);
    (*argidx) += 1;
    long int large = atol(argv[*argidx]);
    int retval     = fname(mem, iarg, large);
    if (retval != SUN_SUCCESS) { return retval; }
    *arg_used = SUNTRUE;
  }
  return SUN_SUCCESS;
}

SUNErrCode sunCheckAndSetIntLongArgs(void* mem, int* argidx,
                                     char* argv[], const size_t offset,
                                     const struct sunKeyIntLongPair* testpairs,
                                     int numpairs, sunbooleantype* arg_used,
                                     int *failedarg)
{
  for (int j = 0; j < numpairs; j++)
  {
    int retval = sunCheckAndSetIntLongArg(mem, argidx, argv, offset,
                                          testpairs[j].key, testpairs[j].set,
                                          arg_used);
    if (retval != SUN_SUCCESS)
    {
      *failedarg = j;
      return retval;
    }
    if (*arg_used) { return SUN_SUCCESS; }
  }
  return SUN_SUCCESS;
}

SUNErrCode sunCheckAndSetRealArg(void* mem, int* argidx, char* argv[],
                                 const size_t offset, const char* argtest,
                                 sunRealSetFn fname, sunbooleantype* arg_used)
{
  *arg_used = SUNFALSE;
  if (strcmp(argv[*argidx] + offset, argtest) == 0)
  {
    (*argidx) += 1;
    sunrealtype rarg = SUNStrToReal(argv[*argidx]);
    int retval       = fname(mem, rarg);
    if (retval != SUN_SUCCESS) { return retval; }
    *arg_used = SUNTRUE;
  }
  return SUN_SUCCESS;
}

SUNErrCode sunCheckAndSetRealArgs(void* mem, int* argidx,
                                  char* argv[], const size_t offset,
                                  const struct sunKeyRealPair* testpairs,
                                  int numpairs, sunbooleantype* arg_used,
                                  int *failedarg)
{
  for (int j = 0; j < numpairs; j++)
  {
    int retval = sunCheckAndSetRealArg(mem, argidx, argv, offset,
                                       testpairs[j].key, testpairs[j].set,
                                       arg_used);
    if (retval != SUN_SUCCESS)
    {
      *failedarg = j;
      return retval;
    }
    if (*arg_used) { return SUN_SUCCESS; }
  }
  return SUN_SUCCESS;
}

SUNErrCode sunCheckAndSetTwoRealArg(void* mem, int* argidx, char* argv[],
                                    const size_t offset, const char* argtest,
                                    sunTwoRealSetFn fname,
                                    sunbooleantype* arg_used)
{
  *arg_used = SUNFALSE;
  if (strcmp(argv[*argidx] + offset, argtest) == 0)
  {
    (*argidx) += 1;
    sunrealtype rarg1 = SUNStrToReal(argv[*argidx]);
    (*argidx) += 1;
    sunrealtype rarg2 = SUNStrToReal(argv[*argidx]);
    int retval        = fname(mem, rarg1, rarg2);
    if (retval != SUN_SUCCESS) { return retval; }
    *arg_used = SUNTRUE;
  }
  return SUN_SUCCESS;
}

SUNErrCode sunCheckAndSetTwoRealArgs(void* mem, int* argidx,
                                     char* argv[], const size_t offset,
                                     const struct sunKeyTwoRealPair* testpairs,
                                     int numpairs, sunbooleantype* arg_used,
                                     int *failedarg)
{
  for (int j = 0; j < numpairs; j++)
  {
    int retval = sunCheckAndSetTwoRealArg(mem, argidx, argv, offset,
                                          testpairs[j].key, testpairs[j].set,
                                          arg_used);
    if (retval != SUN_SUCCESS)
    {
      *failedarg = j;
      return retval;
    }
    if (*arg_used) { return SUN_SUCCESS; }
  }
  return SUN_SUCCESS;
}

SUNErrCode sunCheckAndSetCharArg(void* mem, int* argidx, char* argv[],
                                 const size_t offset, const char* argtest,
                                 sunCharSetFn fname, sunbooleantype* arg_used)
{
  *arg_used = SUNFALSE;
  if (strcmp(argv[*argidx] + offset, argtest) == 0)
  {
    (*argidx) += 1;
    int retval = fname(mem, argv[*argidx]);
    if (retval != SUN_SUCCESS) { return retval; }
    *arg_used = SUNTRUE;
  }
  return SUN_SUCCESS;
}

SUNErrCode sunCheckAndSetCharArgs(void* mem, int* argidx,
                                  char* argv[], const size_t offset,
                                  const struct sunKeyCharPair* testpairs,
                                  int numpairs, sunbooleantype* arg_used,
                                  int *failedarg)
{
  for (int j = 0; j < numpairs; j++)
  {
    int retval = sunCheckAndSetCharArg(mem, argidx, argv, offset,
                                       testpairs[j].key, testpairs[j].set,
                                       arg_used);
    if (retval != SUN_SUCCESS)
    {
      *failedarg = j;
      return retval;
    }
    if (*arg_used) { return SUN_SUCCESS; }
  }
  return SUN_SUCCESS;
}

SUNErrCode sunCheckAndSetTwoCharArg(void* mem, int* argidx, char* argv[],
                                    const size_t offset, const char* argtest,
                                    sunTwoCharSetFn fname,
                                    sunbooleantype* arg_used)
{
  *arg_used = SUNFALSE;
  if (strcmp(argv[*argidx] + offset, argtest) == 0)
  {
    int retval = fname(mem, argv[*argidx + 1], argv[*argidx + 2]);
    (*argidx) += 2;
    if (retval != SUN_SUCCESS) { return retval; }
    *arg_used = SUNTRUE;
  }
  return SUN_SUCCESS;
}

SUNErrCode sunCheckAndSetTwoCharArgs(void* mem, int* argidx,
                                     char* argv[], const size_t offset,
                                     const struct sunKeyTwoCharPair* testpairs,
                                     int numpairs, sunbooleantype* arg_used,
                                     int *failedarg)
{
  for (int j = 0; j < numpairs; j++)
  {
    int retval = sunCheckAndSetTwoCharArg(mem, argidx, argv, offset,
                                          testpairs[j].key, testpairs[j].set,
                                          arg_used);
    if (retval != SUN_SUCCESS)
    {
      *failedarg = j;
      return retval;
    }
    if (*arg_used) { return SUN_SUCCESS; }
  }
  return SUN_SUCCESS;
}

SUNErrCode sunCheckAndSetActionArg(void* mem, int* argidx, char* argv[],
                                   const size_t offset, const char* argtest,
                                   sunActionSetFn fname, sunbooleantype* arg_used)
{
  *arg_used = SUNFALSE;
  if (strcmp(argv[*argidx] + offset, argtest) == 0)
  {
    int retval = fname(mem);
    if (retval != SUN_SUCCESS) { return retval; }
    *arg_used = SUNTRUE;
  }
  return SUN_SUCCESS;
}

SUNErrCode sunCheckAndSetActionArgs(void* mem, int* argidx,
                                    char* argv[], const size_t offset,
                                    const struct sunKeyActionPair* testpairs,
                                    int numpairs, sunbooleantype* arg_used,
                                    int *failedarg)
{
  for (int j = 0; j < numpairs; j++)
  {
    int retval = sunCheckAndSetActionArg(mem, argidx, argv, offset,
                                         testpairs[j].key, testpairs[j].set,
                                         arg_used);
    if (retval != SUN_SUCCESS)
    {
      *failedarg = j;
      return retval;
    }
    if (*arg_used) { return SUN_SUCCESS; }
  }
  return SUN_SUCCESS;
}
