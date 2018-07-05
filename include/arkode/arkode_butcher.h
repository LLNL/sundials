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

#ifdef __cplusplus
}
#endif

#endif
