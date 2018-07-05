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
 * Header file with Butcher table IDs for built-in ERK methods.
 *--------------------------------------------------------------*/

#ifndef _ARKODE_ERK_TABLES_H
#define _ARKODE_ERK_TABLES_H

#include <arkode/arkode_butcher.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/* Butcher table accessor IDs -- Take care when adding to this list 
   to ensure that the IDs for DIRK methods are updated as well, 
   since ARKStep module requires ERK->DIRK ordering so that 
   MIN_ERK_NUM, MAX_ERK_NUM, MIN_DIRK_NUM and MAX_DIRK_NUM are valid */

#define HEUN_EULER_2_1_2         0
#define BOGACKI_SHAMPINE_4_2_3   1
#define ARK324L2SA_ERK_4_2_3     2
#define ZONNEVELD_5_3_4          3
#define ARK436L2SA_ERK_6_3_4     4
#define SAYFY_ABURUB_6_3_4       5
#define CASH_KARP_6_4_5          6
#define FEHLBERG_6_4_5           7
#define DORMAND_PRINCE_7_4_5     8
#define ARK548L2SA_ERK_8_4_5     9
#define VERNER_8_5_6            10
#define FEHLBERG_13_7_8         11

/* Utility #defines to ensure valid input IDs for ERK tables */
#define MIN_ERK_NUM              0
#define MAX_ERK_NUM             11

/* Accessor routine to load built-in ERK table */  
ARKodeButcherTable ARKodeLoadButcherTable_ERK(int imethod);

  
#ifdef __cplusplus
}
#endif

#endif
