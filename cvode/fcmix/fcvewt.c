/*
 * -----------------------------------------------------------------
 * $Revision: 1.1.2.2 $
 * $Date: 2005-04-06 23:32:53 $
 * ----------------------------------------------------------------- 
 * Programmer: Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see sundials/cvode/LICENSE.
 * -----------------------------------------------------------------
 * Fortran/C interface routines for CVODE, for the case of a 
 * user-supplied error weight calculation routine.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include "fcvode.h"         /* actual fn. names, prototypes and global vars.  */
#include "nvector.h"        /* definitions of type N_Vector and vector macros */
#include "sundialstypes.h"  /* definition of type realtype                    */

/***************************************************************************/

/* Prototype of the Fortran routine */

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif
  extern void FCV_EWT(realtype*, realtype*, int*);
#ifdef __cplusplus
}
#endif

/***************************************************************************/

void FCV_EWTSET(int *flag, int *ier)
{
  if (*flag == 1)
    *ier = CVodeSetEwtFn(CV_cvodemem, FCVEwtSet, NULL);
}

/***************************************************************************/

/* 
 * C function to interface between CVODE and a Fortran subroutine FCVEWT.
 */

int FCVEwtSet(N_Vector y, N_Vector ewt, void *e_data)
{
  int ier = 0;
  realtype *ydata, *ewtdata;

  ydata  = N_VGetArrayPointer(y);
  ewtdata = N_VGetArrayPointer(ewt);

  FCV_EWT(ydata, ewtdata, &ier);

  return(ier);
}
