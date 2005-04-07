/*
 * -----------------------------------------------------------------
 * $Revision: 1.14 $
 * $Date: 2005-04-07 23:28:22 $
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

#include "cvode.h"          /* CVODE constants and prototypes        */
#include "fcvode.h"         /* actual function names, prototypes and
			       global variables                      */
#include "fcvroot.h"        /* prototypes of interfaces to CVROOT    */
#include "nvector.h"        /* definition of type N_Vector           */
#include "sundialstypes.h"  /* definition of SUNDIALS type realtype  */

/***************************************************************************/

/* Prototype of the Fortran routine */

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif
  extern void FCV_ROOTFN(realtype *, realtype*, realtype*);
#ifdef __cplusplus
}
#endif

/***************************************************************************/

void FCV_ROOTINIT(int *nrtfn, int *ier)
{
  *ier = CVodeRootInit(CV_cvodemem, *nrtfn, (CVRootFn) FCVrootfunc, NULL);
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

  ydata = N_VGetArrayPointer(y);

  FCV_ROOTFN(&t, ydata, gout);

  return;
}

