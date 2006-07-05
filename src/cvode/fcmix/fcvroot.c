/*
 * -----------------------------------------------------------------
 * $Revision: 1.1 $
 * $Date: 2006-07-05 15:32:33 $
 * ----------------------------------------------------------------- 
 * Programmer(s): Aaron Collier @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * -----------------------------------------------------------------
 * The FCVROOT module contains the routines necessary to use
 * the rootfinding feature of the CVODE module and to interface
 * with the user-supplied Fortran subroutine.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include "fcvode.h"     /* actual fn. names, prototypes and global variables */
#include "fcvroot.h"    /* prototypes of interfaces to CVODE                 */
#include "cvode_impl.h" /* definition of CVodeMem type                       */

/***************************************************************************/

/* Prototype of the Fortran routine */

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif
  extern void FCV_ROOTFN(realtype *, realtype*, realtype*,  /* T, Y, G    */
                         long int*, realtype*,              /* IPAR, RPAR */
                         int *ier);                         /* IER        */
#ifdef __cplusplus
}
#endif

/***************************************************************************/

void FCV_ROOTINIT(int *nrtfn, int *ier)
{
  CVodeMem cv_mem;
  
  cv_mem = (CVodeMem) CV_cvodemem;
  *ier = CVodeRootInit(CV_cvodemem, *nrtfn, (CVRootFn) FCVrootfunc, cv_mem->cv_f_data);
  CV_nrtfn = *nrtfn;

  return;
}

/***************************************************************************/

void FCV_ROOTINFO(int *nrtfn, int *info, int *ier)
{
  *ier = CVodeGetRootInfo(CV_cvodemem, info);
  return; 
}

/***************************************************************************/

void FCV_ROOTFREE(void)
{
  CVodeRootInit(CV_cvodemem, 0, NULL, NULL);

  return;
}

/***************************************************************************/

int FCVrootfunc(realtype t, N_Vector y, realtype *gout, void *g_data)
{
  int ier;
  realtype *ydata;
  FCVUserData CV_userdata;

  ydata = N_VGetArrayPointer(y);

  CV_userdata = (FCVUserData) g_data;

  FCV_ROOTFN(&t, ydata, gout, CV_userdata->ipar, CV_userdata->rpar, &ier);

  return(ier);
}

