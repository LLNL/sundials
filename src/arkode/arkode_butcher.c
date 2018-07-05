/*---------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ SMU
 *---------------------------------------------------------------
 * LLNS/SMU Copyright Start
 * Copyright (c) 2015, Southern Methodist University and 
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
 * This is the implementation file for Butcher table structure
 * for the ARKode infrastructure.
 *--------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include "arkode_impl.h"
#include <sundials/sundials_math.h>

#if defined(SUNDIALS_EXTENDED_PRECISION)
#define RSYM ".32Lg"
#else
#define RSYM ".16g"
#endif


/*---------------------------------------------------------------
  Utility routine to allocate a table structure
  ---------------------------------------------------------------*/
ARKodeButcherTable AllocButcherTable(int stages, booleantype embedded) {

  int i;
  ARKodeButcherTable B;

  /* Check for legal 'stages' value */
  if (stages < 1)  return(NULL);

  /* Allocate Butcher table structure */
  B = (ARKodeButcherTable) malloc(sizeof(struct ARKodeButcherTableMem));

  /* set stages into B structure */
  B->stages = stages;
  
  /* Allocate fields within Butcher table structure */
  B->A = (realtype **) calloc( stages, sizeof(realtype*) );
  for (i=0; i<stages; i++)
    B->A[i] = (realtype *) calloc( stages, sizeof(realtype) );
  B->b = (realtype *) calloc( stages, sizeof(realtype) );
  B->c = (realtype *) calloc( stages, sizeof(realtype) );
  if (embedded)
    B->d = (realtype *) calloc( stages, sizeof(realtype) );

  /* initialize order parameters */
  B->q = 0;
  B->p = 0;

  return(B);
}


/*---------------------------------------------------------------
  Utility routine to query the size of a Butcher table structure
  ---------------------------------------------------------------*/
void ButcherTableSpace(ARKodeButcherTable B, sunindextype *liw,
                       sunindextype *lrw)
{
  /* initialize outputs and return if B is not allocated */
  *liw = 0;  *lrw = 0;
  if (B == NULL)  return;

  /* fill outputs based on B */
  *liw = 3;
  if (B->d != NULL) {
    *lrw = B->stages * (B->stages + 3);
  } else {
    *lrw = B->stages * (B->stages + 2);
  }
}


/*---------------------------------------------------------------
  Utility routine to free a Butcher table structure
  ---------------------------------------------------------------*/
void FreeButcherTable(ARKodeButcherTable B) {

  int i;
  
  /* Free each field within Butcher table structure, and then 
     free structure itself */
  if (B != NULL) {
    if (B->d != NULL)  free(B->d);
    if (B->c != NULL)  free(B->c);
    if (B->b != NULL)  free(B->b);
    if (B->A != NULL) {
      for (i=0; i<B->stages; i++)
        if (B->A[i] != NULL)  free(B->A[i]);
      free(B->A);
    }

    free(B);
  }
}


/*---------------------------------------------------------------
  Utility routine to print a Butcher table structure
  ---------------------------------------------------------------*/
void WriteButcherTable(ARKodeButcherTable B, FILE *outfile) {
  int i, j;
  if (B == NULL) return;
  fprintf(outfile, "  A = \n");
  for (i=0; i<B->stages; i++) {
    fprintf(outfile, "      ");
    for (j=0; j<B->stages; j++)
      fprintf(outfile, "%"RSYM"  ", B->A[i][j]);
    fprintf(outfile, "\n");
  }
  fprintf(outfile, "  c = ");
  for (i=0; i<B->stages; i++) 
    fprintf(outfile, "%"RSYM"  ", B->c[i]);
  fprintf(outfile, "\n");
  fprintf(outfile, "  b = ");
  for (i=0; i<B->stages; i++) 
    fprintf(outfile, "%"RSYM"  ", B->b[i]);
  fprintf(outfile, "\n");
  fprintf(outfile, "  d = ");
  for (i=0; i<B->stages; i++) 
    fprintf(outfile, "%"RSYM"  ", B->d[i]);
  fprintf(outfile, "\n");
}


/*---------------------------------------------------------------
  EOF
  ---------------------------------------------------------------*/
