/* -----------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ SMU
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2023, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * This is the header file for ARKode's built-in DIRK Butcher tables.
 * -----------------------------------------------------------------*/

#ifndef _ARKODE_DIRK_TABLES_H
#define _ARKODE_DIRK_TABLES_H

#include <arkode/arkode_butcher.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

typedef enum {
  ARKODE_DIRK_NONE = -1, /* ensure enum is signed int */
  ARKODE_MIN_DIRK_NUM = 100,
  ARKODE_SDIRK_2_1_2 = ARKODE_MIN_DIRK_NUM,
  ARKODE_BILLINGTON_3_3_2,
  ARKODE_TRBDF2_3_3_2,
  ARKODE_KVAERNO_4_2_3,
  ARKODE_ARK324L2SA_DIRK_4_2_3,
  ARKODE_CASH_5_2_4,
  ARKODE_CASH_5_3_4,
  ARKODE_SDIRK_5_3_4,
  ARKODE_KVAERNO_5_3_4,
  ARKODE_ARK436L2SA_DIRK_6_3_4,
  ARKODE_KVAERNO_7_4_5,
  ARKODE_ARK548L2SA_DIRK_8_4_5,
  ARKODE_ARK437L2SA_DIRK_7_3_4,
  ARKODE_ARK548L2SAb_DIRK_8_4_5,
  ARKODE_ESDIRK324L2SA_4_2_3,
  ARKODE_ESDIRK325L2SA_5_2_3,
  ARKODE_ESDIRK32I5L2SA_5_2_3,
  ARKODE_ESDIRK436L2SA_6_3_4,
  ARKODE_ESDIRK43I6L2SA_6_3_4,
  ARKODE_QESDIRK436L2SA_6_3_4,
  ARKODE_ESDIRK437L2SA_7_3_4,
  ARKODE_ESDIRK547L2SA_7_4_5,
  ARKODE_ESDIRK547L2SA2_7_4_5,
  ARKODE_ARK2_DIRK_3_1_2,
  ARKODE_MAX_DIRK_NUM = ARKODE_ARK2_DIRK_3_1_2
} ARKODE_DIRKTableID;

/* Accessor routine to load built-in DIRK table */
SUNDIALS_EXPORT ARKodeButcherTable ARKodeButcherTable_LoadDIRK(ARKODE_DIRKTableID imethod);

/* Accessor routine to load built-in DIRK table */
SUNDIALS_EXPORT ARKodeButcherTable ARKodeButcherTable_LoadDIRKByName(const char *imethod);

#ifdef __cplusplus
}
#endif

#endif
