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

typedef enum {
  ARKODE_SPRK_NONE = -1, /* ensure enum is signed int */
  ARKODE_MIN_SPRK_NUM = 0,
  ARKODE_SYMPLECTIC_EULER_1 = ARKODE_MIN_SPRK_NUM,
  ARKODE_SYMPLECTIC_LEAPFROG_2,
  ARKODE_SYMPLECTIC_PSEUDO_LEAPFROG_2,
  ARKODE_SYMPLECTIC_RUTH_3,
  ARKODE_SYMPLECTIC_MCLACHLAN_2,
  ARKODE_SYMPLECTIC_MCLACHLAN_3,
  ARKODE_SYMPLECTIC_CANDY_ROZMUS_4,
  ARKODE_SYMPLECTIC_MCLACHLAN_4,
  ARKODE_SYMPLECTIC_MCLACHLAN_5,
  ARKODE_SYMPLECTIC_YOSHIDA_6,
  ARKODE_SYMPLECTIC_MCLACHLAN_8,
  ARKODE_SYMPLECTIC_SOFRONIOU_10,
  ARKODE_MAX_SPRK_NUM = ARKODE_SYMPLECTIC_SOFRONIOU_10
} ARKODE_SPRKMethodID;

struct ARKodeSPRKMem_s {

  int q;           /* method order of accuracy         */
  int stages;      /* number of stages                 */
  sunrealtype* a;  /* coefficients multiplying q'     */
  sunrealtype* b;  /* coefficients multiplying p'     */

  /* the a_i coefficients generate the explicit Butcher table */
  /* the b_i coefficients generate the diagonally-implicit Butcher table */

};

typedef _SUNDIALS_STRUCT_ ARKodeSPRKMem_s *ARKodeSPRKMem;

/* Utility routines to allocate/free/output SPRK structures */
SUNDIALS_EXPORT ARKodeSPRKMem ARKodeSPRKMem_Alloc(int stages);
SUNDIALS_EXPORT ARKodeSPRKMem ARKodeSPRKMem_Load(ARKODE_SPRKMethodID id);
SUNDIALS_EXPORT ARKodeSPRKMem ARKodeSPRKMem_Copy(ARKodeSPRKMem B);
SUNDIALS_EXPORT void ARKodeSPRKMem_Space(ARKodeSPRKMem B, sunindextype *liw, sunindextype *lrw);
SUNDIALS_EXPORT void ARKodeSPRKMem_Free(ARKodeSPRKMem B);
SUNDIALS_EXPORT int ARKodeSPRKMem_ToButcher(ARKodeSPRKMem sprk_mem, ARKodeButcherTable* b_ptr, ARKodeButcherTable* B_ptr);

/* Different methods */

ARKodeSPRKMem ARKodeSymplecticEuler();
ARKodeSPRKMem ARKodeSymplecticLeapfrog2();
ARKodeSPRKMem ARKodeSymplecticPseudoLeapfrog2();
ARKodeSPRKMem ARKodeSymplecticRuth3();
ARKodeSPRKMem ARKodeSymplecticCandyRozmus4();
ARKodeSPRKMem ARKodeSymplecticMcLachlan2();
ARKodeSPRKMem ARKodeSymplecticMcLachlan3();
ARKodeSPRKMem ARKodeSymplecticMcLachlan4();
ARKodeSPRKMem ARKodeSymplecticMcLachlan5();
ARKodeSPRKMem ARKodeSymplecticYoshida6();
ARKodeSPRKMem ARKodeSymplecticMcLachlan8();
ARKodeSPRKMem ARKodeSymplecticSofroniou10();

#ifdef __cplusplus
}
#endif

#endif
