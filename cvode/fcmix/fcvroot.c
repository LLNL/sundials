/*
 * -----------------------------------------------------------------
 * $Revision: 1.11 $
 * $Date: 2004-10-21 20:55:05 $
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

void FCV_ROOTFN(realtype *, realtype*, realtype*);

/***************************************************************************/

void FCV_ROOTINIT(int *nrtfn, int *ier)
{
  *ier = CVodeRootInit(CV_cvodemem, (CVRootFn) FCVrootfunc, *nrtfn);
  CV_nrtfn = *nrtfn;

  return;
}

/***************************************************************************/

void FCV_ROOTINFO(int *nrtfn, int *info, int *ier)
{
  int *rootsfound;
  int i;

  *ier = CVodeGetRootInfo(CV_cvodemem, &rootsfound);

  for (i = 0; i < *nrtfn; i++) info[i] = rootsfound[i];

  return; 
}

/***************************************************************************/

void FCV_ROOTFREE(void)
{
  CVodeRootInit(CV_cvodemem, NULL, 0);

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
