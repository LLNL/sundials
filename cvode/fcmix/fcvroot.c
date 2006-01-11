/*
 * -----------------------------------------------------------------
 * $Revision: 1.17 $
 * $Date: 2006-01-11 21:13:45 $
 * ----------------------------------------------------------------- 
 * Programmer(s): Aaron Collier @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see sundials/cvode/LICENSE.
 * -----------------------------------------------------------------
 * The FCVROOT module contains the routines necessary to use
 * the rootfinding feature of the CVODE module and to interface
 * with the user-supplied Fortran subroutine.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include "cvode.h"            /* CVODE constants and prototypes                    */
#include "fcvode.h"           /* actual fn. names, prototypes and global variables */
#include "fcvroot.h"          /* prototypes of interfaces to CVODE                 */
#include "cvode_impl.h"       /* definition of CVodeMem type                       */
#include "sundials_nvector.h" /* definition of type N_Vector                       */
#include "sundials_types.h"   /* definition of SUNDIALS type realtype              */

/***************************************************************************/

/* Prototype of the Fortran routine */

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif
  extern void FCV_ROOTFN(realtype *, realtype*, realtype*,  /* T, Y, G */
                         long int*, realtype*);             /* IPAR, RPAR */
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

void FCVrootfunc(realtype t, N_Vector y, realtype *gout, void *g_data)
{
  realtype *ydata;
  FCVUserData CV_userdata;

  ydata = N_VGetArrayPointer(y);

  CV_userdata = (FCVUserData) g_data;

  FCV_ROOTFN(&t, ydata, gout, CV_userdata->ipar, CV_userdata->rpar);

  return;
}

