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
  ARKODE_SYMPLECTIC_EULER_1_1 = ARKODE_MIN_SPRK_NUM,
  ARKODE_SYMPLECTIC_LEAPFROG_2_2,
  ARKODE_SYMPLECTIC_PSEUDO_LEAPFROG_2_2,
  ARKODE_SYMPLECTIC_RUTH_3_3,
  ARKODE_SYMPLECTIC_MCLACHLAN_2_2,
  ARKODE_SYMPLECTIC_MCLACHLAN_3_3,
  ARKODE_SYMPLECTIC_CANDY_ROZMUS_4_4,
  ARKODE_SYMPLECTIC_MCLACHLAN_4_4,
  ARKODE_SYMPLECTIC_MCLACHLAN_5_6,
  ARKODE_SYMPLECTIC_YOSHIDA_6_8,
  ARKODE_SYMPLECTIC_MCLACHLAN_8_16,
  ARKODE_SYMPLECTIC_SOFRONIOU_10_36,
  ARKODE_MAX_SPRK_NUM = ARKODE_SYMPLECTIC_SOFRONIOU_10_36
} ARKODE_SPRKMethodID;

struct ARKodeSPRKStorage_s
{
  int q;          /* method order of accuracy         */
  int stages;     /* number of stages                 */
  sunrealtype* a; /* coefficients that generate the explicit Butcher table */
  sunrealtype* b; /* coefficients that generate the diagonally-implicit Butcher
                     table */

  /* the a_i coefficients generate the explicit Butcher table */
  /* the b_i coefficients generate the diagonally-implicit Butcher table */
};

typedef _SUNDIALS_STRUCT_ ARKodeSPRKStorage_s* ARKodeSPRKStorage;

/* Utility routines to allocate/free/output SPRK structures */
SUNDIALS_EXPORT ARKodeSPRKStorage ARKodeSPRKStorage_Alloc(int stages);
SUNDIALS_EXPORT ARKodeSPRKStorage ARKodeSPRKStorage_Load(ARKODE_SPRKMethodID id);
SUNDIALS_EXPORT ARKodeSPRKStorage ARKodeSPRKStorage_LoadByName(const char* method);
SUNDIALS_EXPORT ARKodeSPRKStorage
ARKodeSPRKStorage_Copy(ARKodeSPRKStorage that_sprk_storage);
SUNDIALS_EXPORT void ARKodeSPRKStorage_Space(ARKodeSPRKStorage B,
                                             sunindextype* liw,
                                             sunindextype* lrw);
SUNDIALS_EXPORT void ARKodeSPRKStorage_Free(ARKodeSPRKStorage sprk_storage);
SUNDIALS_EXPORT int ARKodeSPRKStorage_ToButcher(ARKodeSPRKStorage sprk_storage,
                                                ARKodeButcherTable* a_ptr,
                                                ARKodeButcherTable* b_ptr);

/* Different methods */


#ifdef __cplusplus
}
#endif

#endif
