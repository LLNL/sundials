/*
 * -----------------------------------------------------------------
 * $Revision: 1.5 $
 * $Date: 2006-01-24 22:17:29 $
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

#include "fida.h"             /* actual function names, prototypes and global
			         variables                                      */
#include "ida_impl.h"         /* definition of IDAMem type                      */
#include "sundials_nvector.h" /* definitions of type N_Vector and vector macros */
#include "sundials_types.h"   /* definition of type realtype                    */

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
  IDAMem ida_mem;

  *ier = 0;

  if (*flag != 0) {
    ida_mem = (IDAMem) IDA_idamem;
    *ier = IDASetEwtFn(IDA_idamem, (IDAEwtFn) FIDAEwtSet, ida_mem->ida_rdata);
  }

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

  IDA_userdata = (FIDAUserData) e_data;

  /* Call user-supplied routine */
  FIDA_EWT(y_data, ewt_data, IDA_userdata->ipar, IDA_userdata->rpar, &ier);

  return(ier);
}
