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
#include "sundials_cli.h"

/*===============================================================
  Command-line input utility routines
  ===============================================================*/

int sunCheckAndSetIntArg(void* mem, int* i, char* argv[], const size_t offset,
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

int sunCheckAndSetTwoIntArg(void* mem, int* i, char* argv[],
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

int sunCheckAndSetLongArg(void* mem, int* i, char* argv[], const size_t offset,
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

int sunCheckAndSetRealArg(void* mem, int* i, char* argv[], const size_t offset,
                          const char* argtest, sunRealSetFn fname,
                          sunbooleantype* arg_used)
{
  *arg_used = SUNFALSE;
  if (strcmp(argv[*i] + offset, argtest) == 0)
  {
    (*i) += 1;
    sunrealtype rarg = atof(argv[*i]);
    int retval       = fname(mem, rarg);
    if (retval != SUN_SUCCESS) { return retval; }
    *arg_used = SUNTRUE;
  }
  return SUN_SUCCESS;
}

int sunCheckAndSetTwoRealArg(void* mem, int* i, char* argv[],
                             const size_t offset, const char* argtest,
                             sunTwoRealSetFn fname, sunbooleantype* arg_used)
{
  *arg_used = SUNFALSE;
  if (strcmp(argv[*i] + offset, argtest) == 0)
  {
    (*i) += 1;
    sunrealtype rarg1 = atof(argv[*i]);
    (*i) += 1;
    sunrealtype rarg2 = atof(argv[*i]);
    int retval        = fname(mem, rarg1, rarg2);
    if (retval != SUN_SUCCESS) { return retval; }
    *arg_used = SUNTRUE;
  }
  return SUN_SUCCESS;
}

int sunCheckAndSetCharArg(void* mem, int* i, char* argv[], const size_t offset,
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

int sunCheckAndSetTwoCharArg(void* mem, int* i, char* argv[],
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

int sunCheckAndSetActionArg(void* mem, int* i, char* argv[],
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
