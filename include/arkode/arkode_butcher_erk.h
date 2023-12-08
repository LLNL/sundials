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
 * This is the header file for ARKode's built-in ERK Butcher tables.
 * -----------------------------------------------------------------*/

#ifndef _ARKODE_ERK_TABLES_H
#define _ARKODE_ERK_TABLES_H

#include <arkode/arkode_butcher.h>

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

typedef enum
{
  ARKODE_ERK_NONE         = -1, /* ensure enum is signed int */
  ARKODE_MIN_ERK_NUM      = 0,
  ARKODE_HEUN_EULER_2_1_2 = ARKODE_MIN_ERK_NUM,
  ARKODE_BOGACKI_SHAMPINE_4_2_3,
  ARKODE_ARK324L2SA_ERK_4_2_3,
  ARKODE_ZONNEVELD_5_3_4,
  ARKODE_ARK436L2SA_ERK_6_3_4,
  ARKODE_SAYFY_ABURUB_6_3_4,
  ARKODE_CASH_KARP_6_4_5,
  ARKODE_FEHLBERG_6_4_5,
  ARKODE_DORMAND_PRINCE_7_4_5,
  ARKODE_ARK548L2SA_ERK_8_4_5,
  ARKODE_VERNER_8_5_6,
  ARKODE_FEHLBERG_13_7_8,
  ARKODE_KNOTH_WOLKE_3_3,
  ARKODE_ARK437L2SA_ERK_7_3_4,
  ARKODE_ARK548L2SAb_ERK_8_4_5,
  ARKODE_ARK2_ERK_3_1_2,
  ARKODE_SOFRONIOU_SPALETTA_5_3_4,
  ARKODE_SHU_OSHER_3_2_3,
  ARKODE_VERNER_9_5_6,
  ARKODE_VERNER_10_6_7,
  ARKODE_VERNER_13_7_8,
  ARKODE_VERNER_16_8_9,
  ARKODE_MAX_ERK_NUM = ARKODE_VERNER_16_8_9
} ARKODE_ERKTableID;

/* Accessor routine to load built-in ERK table */
SUNDIALS_EXPORT ARKodeButcherTable
ARKodeButcherTable_LoadERK(ARKODE_ERKTableID emethod);

SUNDIALS_EXPORT ARKodeButcherTable
ARKodeButcherTable_LoadERKByName(const char* emethod);

#ifdef __cplusplus
}
#endif

#endif
