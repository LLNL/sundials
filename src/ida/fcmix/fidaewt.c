/*
 * -----------------------------------------------------------------
 * $Revision: 1.3 $
 * $Date: 2007-04-30 19:29:00 $
 * ----------------------------------------------------------------- 
 * Programmer(s): Aaron Collier @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2005, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * -----------------------------------------------------------------
 * Fortran/C interface routines for IDA, for the case of a 
 * user-supplied error weight calculation routine.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include "fida.h"     /* actual function names, prototypes and global vars.*/
#include "ida_impl.h" /* definition of IDAMem type                         */

/*************************************************/

/* Prototype of user-supplied Fortran routine */

#ifdef __cplusplus  /* wrapper to enable C++ usage (IDAEwtFn) */
extern "C" {
#endif

  extern void FIDA_EWT(realtype*, realtype*,   /* Y, EWT */ 
                       long int*, realtype*,   /* IPAR, RPAR */
                       int*);                  /* IER */

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

  if (*flag != 0) {
    *ier = IDAWFtolerances(IDA_idamem,  FIDAEwtSet);
  }

  return;
}

/*************************************************/

/* 
 * C function to interface between IDA and a Fortran subroutine FIDAVEWT.
 */

int FIDAEwtSet(N_Vector y, N_Vector ewt, void *user_data)
{
  int ier;
  realtype *y_data, *ewt_data;
  FIDAUserData IDA_userdata;

  /* Initialize all pointers to NULL */
  y_data = ewt_data = NULL;

  /* NOTE: The user-supplied routine should set ier to an
     appropriate value, but we preset the value to zero
     (meaning SUCCESS) so the user need only reset the
     value if an error occurred */
  ier = 0;

  y_data   = N_VGetArrayPointer(y);
  ewt_data = N_VGetArrayPointer(ewt);

  IDA_userdata = (FIDAUserData) user_data;

  /* Call user-supplied routine */
  FIDA_EWT(y_data, ewt_data, IDA_userdata->ipar, IDA_userdata->rpar, &ier);

  return(ier);
}
