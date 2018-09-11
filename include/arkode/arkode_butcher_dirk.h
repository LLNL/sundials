/*---------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ SMU
 *---------------------------------------------------------------
 * LLNS/SMU Copyright Start
 * Copyright (c) 2018, Southern Methodist University and 
 * Lawrence Livermore National Security
 *
 * This work was performed under the auspices of the U.S. Department 
 * of Energy by Southern Methodist University and Lawrence Livermore 
 * National Laboratory under Contract DE-AC52-07NA27344.
 * Produced at Southern Methodist University and the Lawrence 
 * Livermore National Laboratory.
 *
 * All rights reserved.
 * For details, see the LICENSE file.
 * LLNS/SMU Copyright End
 *---------------------------------------------------------------
 * Header file with Butcher table IDs for built-in DIRK methods.
 *--------------------------------------------------------------*/

#ifndef _ARKODE_DIRK_TABLES_H
#define _ARKODE_DIRK_TABLES_H

#include <arkode/arkode_butcher.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/* Butcher table accessor IDs */
#define SDIRK_2_1_2             12
#define BILLINGTON_3_3_2        13
#define TRBDF2_3_3_2            14
#define KVAERNO_4_2_3           15
#define ARK324L2SA_DIRK_4_2_3   16
#define CASH_5_2_4              17
#define CASH_5_3_4              18
#define SDIRK_5_3_4             19
#define KVAERNO_5_3_4           20
#define ARK436L2SA_DIRK_6_3_4   21
#define KVAERNO_7_4_5           22
#define ARK548L2SA_DIRK_8_4_5   23

/* Utility #defines to ensure valid input IDs for DIRK tables */
#define MIN_DIRK_NUM            12
#define MAX_DIRK_NUM            23

/* Accessor routine to load built-in DIRK table */  
ARKodeButcherTable ARKodeLoadButcherTable_DIRK(int imethod);

  
#ifdef __cplusplus
}
#endif

#endif
