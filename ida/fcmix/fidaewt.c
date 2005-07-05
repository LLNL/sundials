/*
 * -----------------------------------------------------------------
 * $Revision: 1.3 $
 * $Date: 2005-07-05 22:35:14 $
 * ----------------------------------------------------------------- 
 * Programmer(s): Aaron Collier @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2005, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see sundials/ida/LICENSE.
 * -----------------------------------------------------------------
 * Fortran/C interface routines for IDA, for the case of a 
 * user-supplied error weight calculation routine.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include "fida.h"           /* actual function names, prototypes and global
			       variables                                      */
#include "nvector.h"        /* definitions of type N_Vector and vector macros */
#include "sundialstypes.h"  /* definition of type realtype                    */

/*************************************************/

/* Prototype of user-supplied Fortran routine */

#ifdef __cplusplus  /* wrapper to enable C++ usage (IDAEwtFn) */
extern "C" {
#endif

extern void FIDA_EWT(realtype*, realtype*, int*);

#ifdef __cplusplus
}
#endif

/*************************************************/

/* 
 * User-callable function to interface to IDASetEwtFn.
 */

void FIDA_EWTSET(int *flag, int *ier)
{
  *ier = 0;

  if (*flag != 0) *ier = IDASetEwtFn(IDA_idamem, (IDAEwtFn) FIDAEwtSet, NULL);

  return;
}

/*************************************************/

/* 
 * C function to interface between IDA and a Fortran subroutine FIDAVEWT.
 */

int FIDAEwtSet(N_Vector y, N_Vector ewt, void *e_data)
{
  int ier;
  realtype *y_data, *ewt_data;

  /* Initialize all pointers to NULL */
  y_data = ewt_data = NULL;

  /* NOTE: The user-supplied routine should set ier to an
     appropriate value, but we preset the value to zero
     (meaning SUCCESS) so the user need only reset the
     value if an error occurred */
  ier = 0;

  y_data   = N_VGetArrayPointer(y);
  ewt_data = N_VGetArrayPointer(ewt);

  /* Call user-supplied routine */
  FIDA_EWT(y_data, ewt_data, &ier);

  return(ier);
}
