/*
 * -----------------------------------------------------------------
 * $Revision: 1.2 $
 * $Date: 2007-04-23 23:37:19 $
 * ----------------------------------------------------------------- 
 * Programmer: Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * -----------------------------------------------------------------
 * Fortran/C interface routines for CVODE, for the case of a 
 * user-supplied error weight calculation routine.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include "fcvode.h"           /* actual fn. names, prototypes and global vars.  */
#include "cvode_impl.h"       /* definition of CVodeMem type                    */

/***************************************************************************/

/* Prototype of the Fortran routine */

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif
  extern void FCV_EWT(realtype*, realtype*,  /* Y, EWT */ 
                      long int*, realtype*,  /* IPAR, RPAR */
                      int*);                 /* IER */
#ifdef __cplusplus
}
#endif

/***************************************************************************/

/* 
 * User-callable function to interface to CVodeSetEwtFn.
 */

void FCV_EWTSET(int *flag, int *ier)
{
  CVodeMem cv_mem;

  if (*flag != 0) {
    cv_mem = (CVodeMem) CV_cvodemem;
    *ier = CVodeWFtolerances(CV_cvodemem, FCVEwtSet);
  }
}

/***************************************************************************/

/* 
 * C function to interface between CVODE and a Fortran subroutine FCVEWT.
 */

int FCVEwtSet(N_Vector y, N_Vector ewt, void *f_data)
{
  int ier = 0;
  realtype *ydata, *ewtdata;
  FCVUserData CV_userdata;

  ydata  = N_VGetArrayPointer(y);
  ewtdata = N_VGetArrayPointer(ewt);

  CV_userdata = (FCVUserData) f_data;

  FCV_EWT(ydata, ewtdata, CV_userdata->ipar, CV_userdata->rpar, &ier);

  return(ier);
}
