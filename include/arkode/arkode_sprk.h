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

#ifndef _ARKODE_SPRKSTORAGE_H
#define _ARKODE_SPRKSTORAGE_H

#include <arkode/arkode_butcher.h>
#include <sundials/sundials_types.h>

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

typedef enum
{
  ARKODE_SPRK_NONE            = -1, /* ensure enum is signed int */
  ARKODE_MIN_SPRK_NUM         = 0,
  ARKODE_SPRK_EULER_1_1 = ARKODE_MIN_SPRK_NUM,
  ARKODE_SPRK_LEAPFROG_2_2,
  ARKODE_SPRK_PSEUDO_LEAPFROG_2_2,
  ARKODE_SPRK_RUTH_3_3,
  ARKODE_SPRK_MCLACHLAN_2_2,
  ARKODE_SPRK_MCLACHLAN_3_3,
  ARKODE_SPRK_CANDY_ROZMUS_4_4,
  ARKODE_SPRK_MCLACHLAN_4_4,
  ARKODE_SPRK_MCLACHLAN_5_6,
  ARKODE_SPRK_YOSHIDA_6_8,
  ARKODE_SPRK_SUZUKI_UMENO_8_16,
  ARKODE_SPRK_SOFRONIOU_10_36,
  ARKODE_MAX_SPRK_NUM = ARKODE_SPRK_SOFRONIOU_10_36
} ARKODE_SPRKMethodID;

struct ARKodeSPRKStorage_s
{
  /* method order of accuracy */
  int q;
  /* number of stages */
  int stages;
  /* the a_i coefficients generate the explicit Butcher table */
  sunrealtype* a;
  /* the ahat_i coefficients generate the diagonally-implicit Butcher table */
  sunrealtype* ahat;
};

typedef _SUNDIALS_STRUCT_ ARKodeSPRKStorage_s* ARKodeSPRKStorage;

/* Utility routines to allocate/free/output SPRK structures */
SUNDIALS_EXPORT ARKodeSPRKStorage ARKodeSPRKStorage_Alloc(int stages);
SUNDIALS_EXPORT ARKodeSPRKStorage ARKodeSPRKStorage_Load(ARKODE_SPRKMethodID id);
SUNDIALS_EXPORT ARKodeSPRKStorage ARKodeSPRKStorage_LoadByName(const char* method);
SUNDIALS_EXPORT ARKodeSPRKStorage
ARKodeSPRKStorage_Copy(ARKodeSPRKStorage that_sprk_storage);
SUNDIALS_EXPORT void ARKodeSPRKStorage_Space(ARKodeSPRKStorage sprk_storage,
                                             sunindextype* liw,
                                             sunindextype* lrw);
SUNDIALS_EXPORT void ARKodeSPRKStorage_Free(ARKodeSPRKStorage sprk_storage);
SUNDIALS_EXPORT int ARKodeSPRKStorage_ToButcher(ARKodeSPRKStorage sprk_storage,
                                                ARKodeButcherTable* erk_ptr,
                                                ARKodeButcherTable* dirk_ptr);

/* Different methods */

#ifdef __cplusplus
}
#endif

#endif
