/* -----------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ SMU
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2024, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * This is the header file for ARKode Butcher table structures.
 * -----------------------------------------------------------------*/

#ifndef _ARKODE_BUTCHER_H
#define _ARKODE_BUTCHER_H

#include <stdio.h>
#include <sundials/sundials_types.h>

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

/*---------------------------------------------------------------
  Types : struct ARKodeButcherTableMem, ARKodeButcherTable
  ---------------------------------------------------------------*/
struct ARKodeButcherTableMem
{
  int q;           /* method order of accuracy       */
  int p;           /* embedding order of accuracy    */
  int stages;      /* number of stages               */
  sunrealtype** A; /* Butcher table coefficients     */
  sunrealtype* c;  /* canopy node coefficients       */
  sunrealtype* b;  /* root node coefficients         */
  sunrealtype* d;  /* embedding coefficients         */
};

typedef _SUNDIALS_STRUCT_ ARKodeButcherTableMem* ARKodeButcherTable;

/* Utility routines to allocate/free/output Butcher table structures */
SUNDIALS_EXPORT ARKodeButcherTable
ARKodeButcherTable_Alloc(int stages, sunbooleantype embedded);
SUNDIALS_EXPORT ARKodeButcherTable ARKodeButcherTable_Create(int s, int q, int p,
                                                             sunrealtype* c,
                                                             sunrealtype* A,
                                                             sunrealtype* b,
                                                             sunrealtype* d);
SUNDIALS_EXPORT ARKodeButcherTable ARKodeButcherTable_Copy(ARKodeButcherTable B);
SUNDIALS_EXPORT void ARKodeButcherTable_Space(ARKodeButcherTable B,
                                              sunindextype* liw,
                                              sunindextype* lrw);
SUNDIALS_EXPORT void ARKodeButcherTable_Free(ARKodeButcherTable B);
SUNDIALS_EXPORT void ARKodeButcherTable_Write(ARKodeButcherTable B,
                                              FILE* outfile);
SUNDIALS_EXPORT sunbooleantype
ARKodeButcherTable_IsStifflyAccurate(ARKodeButcherTable B);
SUNDIALS_EXPORT int ARKodeButcherTable_CheckOrder(ARKodeButcherTable B, int* q,
                                                  int* p, FILE* outfile);
SUNDIALS_EXPORT int ARKodeButcherTable_CheckARKOrder(ARKodeButcherTable B1,
                                                     ARKodeButcherTable B2,
                                                     int* q, int* p,
                                                     FILE* outfile);

#ifdef __cplusplus
}
#endif

#endif
