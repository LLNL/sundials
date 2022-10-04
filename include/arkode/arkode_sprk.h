/* -----------------------------------------------------------------
 * Programmer(s): Cody J. Balos @ LLNL
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2022, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * -----------------------------------------------------------------*/

#ifndef _ARKODE_SPRKMEM_H
#define _ARKODE_SPRKMEM_H

#include <sundials/sundials_types.h>
#include <arkode/arkode_butcher.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

struct ARKodeSprkMem_s {

  int q;           /* method order of accuracy        */
  int stages;      /* number of stages                */
  sunrealtype* b;  /* diagonally implicit cofficients */
  sunrealtype* B;  /* explicit table coefficients     */

};

typedef _SUNDIALS_STRUCT_ ARKodeSprkMem_s *ARKodeSprkMem;

/* Utility routines to allocate/free/output SPRK structures */
SUNDIALS_EXPORT ARKodeSprkMem ARKodeSprkMem_Alloc(int stages);
SUNDIALS_EXPORT ARKodeSprkMem ARKodeSprkMem_Copy(ARKodeSprkMem B);
SUNDIALS_EXPORT void ARKodeSprkMem_Space(ARKodeSprkMem B, sunindextype *liw, sunindextype *lrw);
SUNDIALS_EXPORT void ARKodeSprkMem_Free(ARKodeSprkMem B);
SUNDIALS_EXPORT int ARKodeSprkMem_ToButcher(ARKodeSprkMem sprk_mem, ARKodeButcherTable* b_ptr, ARKodeButcherTable* B_ptr);

typedef enum {
  ARKODE_SPRK_NONE = -1, /* ensure enum is signed int */
  ARKODE_MIN_SPRK_NUM = 0,
  ARKODE_SYMPLECTIC_EULER = ARKODE_MIN_SPRK_NUM,
  ARKODE_MAX_SPRK_NUM = ARKODE_SYMPLECTIC_EULER
} ARKODE_SPRKMethodID;

#ifdef __cplusplus
}
#endif

#endif
