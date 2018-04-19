/*---------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ SMU
 *---------------------------------------------------------------
 * LLNS/SMU Copyright Start
 * Copyright (c) 2017, Southern Methodist University and 
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
 * Header file for Butcher table structures in ARKode.
 *--------------------------------------------------------------*/

#ifndef _ARKODE_BUTCHER_H
#define _ARKODE_BUTCHER_H

#include <sundials/sundials_types.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/* Butcher table accessors -- explicit */
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

#define DEFAULT_ERK_2           HEUN_EULER_2_1_2
#define DEFAULT_ERK_3           BOGACKI_SHAMPINE_4_2_3
#define DEFAULT_ERK_4           ZONNEVELD_5_3_4
#define DEFAULT_ERK_5           CASH_KARP_6_4_5
#define DEFAULT_ERK_6           VERNER_8_5_6
#define DEFAULT_ERK_8           FEHLBERG_13_7_8

#define MIN_ERK_NUM              0
#define MAX_ERK_NUM             11

/* Butcher table accessors -- implicit */
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

#define DEFAULT_DIRK_2          SDIRK_2_1_2
#define DEFAULT_DIRK_3          ARK324L2SA_DIRK_4_2_3
#define DEFAULT_DIRK_4          SDIRK_5_3_4
#define DEFAULT_DIRK_5          ARK548L2SA_DIRK_8_4_5

#define MIN_DIRK_NUM            12
#define MAX_DIRK_NUM            23

/* Butcher table accessors -- ImEx */
#define DEFAULT_ARK_ETABLE_3    ARK324L2SA_ERK_4_2_3
#define DEFAULT_ARK_ETABLE_4    ARK436L2SA_ERK_6_3_4
#define DEFAULT_ARK_ETABLE_5    ARK548L2SA_ERK_8_4_5
#define DEFAULT_ARK_ITABLE_3    ARK324L2SA_DIRK_4_2_3
#define DEFAULT_ARK_ITABLE_4    ARK436L2SA_DIRK_6_3_4
#define DEFAULT_ARK_ITABLE_5    ARK548L2SA_DIRK_8_4_5



/*---------------------------------------------------------------
 Types : struct ARKodeButcherTableMem, ARKodeButcherTable
---------------------------------------------------------------*/
typedef struct ARKodeButcherTableMem {

  int q;           /* method order of accuracy       */
  int p;           /* embedding order of accuracy    */
  int stages;      /* number of stages               */
  realtype **A;    /* Butcher table coefficients     */
  realtype *c;     /* canopy node coefficients       */
  realtype *b;     /* root node coefficients         */
  realtype *d;     /* embedding coefficients         */

} *ARKodeButcherTable;


/* Utility routines to allocate/free/output Butcher table structures */
ARKodeButcherTable AllocButcherTable(int stages, booleantype embedded);
void ButcherTableSpace(ARKodeButcherTable B, sunindextype *liw, sunindextype *lrw);
void FreeButcherTable(ARKodeButcherTable B);
void WriteButcherTable(ARKodeButcherTable B, FILE *outfile);

/* Utility routine to fill a pre-defined Butcher table structure */  
ARKodeButcherTable ARKodeLoadButcherTable(int imethod);

#ifdef __cplusplus
}
#endif

#endif
