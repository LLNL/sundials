/* -----------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2025, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * SUNDIALS command-line utility definitions.
 * -----------------------------------------------------------------*/

#ifndef _SUNDIALS_CLI_H
#define _SUNDIALS_CLI_H

#include <stdio.h>
#include "sundials/sundials_errors.h"
#include "sundials/sundials_types.h"

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

/* utilities for integer "set" routines */
typedef int (*sunIntSetFn)(void*, int);

struct sunKeyIntPair
{
  const char* key;
  sunIntSetFn set;
};

SUNDIALS_EXPORT int sunCheckAndSetIntArg(void* mem, int* i, char* argv[],
                                         const size_t offset,
                                         const char* argtest, sunIntSetFn fname,
                                         sunbooleantype* arg_used);

/* utilities for pair-of-integer "set" routines */
typedef int (*sunTwoIntSetFn)(void*, int, int);

struct sunKeyTwoIntPair
{
  const char* key;
  sunTwoIntSetFn set;
};

SUNDIALS_EXPORT int sunCheckAndSetTwoIntArg(void* mem, int* i, char* argv[],
                                            const size_t offset,
                                            const char* argtest,
                                            sunTwoIntSetFn fname,
                                            sunbooleantype* arg_used);

/* utilities for long int "set" routines */
typedef int (*sunLongSetFn)(void*, long int);

struct sunKeyLongPair
{
  const char* key;
  sunLongSetFn set;
};

SUNDIALS_EXPORT int sunCheckAndSetLongArg(void* mem, int* i, char* argv[],
                                          const size_t offset,
                                          const char* argtest, sunLongSetFn fname,
                                          sunbooleantype* arg_used);

/* utilities for sunrealtype "set" routines */
typedef int (*sunRealSetFn)(void*, sunrealtype);

struct sunKeyRealPair
{
  const char* key;
  sunRealSetFn set;
};

SUNDIALS_EXPORT int sunCheckAndSetRealArg(void* mem, int* i, char* argv[],
                                          const size_t offset,
                                          const char* argtest, sunRealSetFn fname,
                                          sunbooleantype* arg_used);

/* utilities for pair-of-sunrealtype "set" routines */
typedef int (*sunTwoRealSetFn)(void*, sunrealtype, sunrealtype);

struct sunKeyTwoRealPair
{
  const char* key;
  sunTwoRealSetFn set;
};

SUNDIALS_EXPORT int sunCheckAndSetTwoRealArg(void* mem, int* i, char* argv[],
                                             const size_t offset,
                                             const char* argtest,
                                             sunTwoRealSetFn fname,
                                             sunbooleantype* arg_used);

/* utilities for char* "set" routines */
typedef int (*sunCharSetFn)(void*, const char*);

struct sunKeyCharPair
{
  const char* key;
  sunCharSetFn set;
};

SUNDIALS_EXPORT int sunCheckAndSetCharArg(void* mem, int* i, char* argv[],
                                          const size_t offset,
                                          const char* argtest, sunCharSetFn fname,
                                          sunbooleantype* arg_used);

/* utilities for pair-of-char* "set" routines */
typedef int (*sunTwoCharSetFn)(void*, const char*, const char*);

struct sunKeyTwoCharPair
{
  const char* key;
  sunTwoCharSetFn set;
};

SUNDIALS_EXPORT int sunCheckAndSetTwoCharArg(void* mem, int* i, char* argv[],
                                             const size_t offset,
                                             const char* argtest,
                                             sunTwoCharSetFn fname,
                                             sunbooleantype* arg_used);

/* utilities for action "set" routines */
typedef int (*sunActionSetFn)(void*);

struct sunKeyActionPair
{
  const char* key;
  sunActionSetFn set;
};

SUNDIALS_EXPORT int sunCheckAndSetActionArg(void* mem, int* i, char* argv[],
                                            const size_t offset,
                                            const char* argtest,
                                            sunActionSetFn fname,
                                            sunbooleantype* arg_used);

#ifdef __cplusplus
}
#endif

#endif
