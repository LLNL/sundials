/*
 * -----------------------------------------------------------------
 * $Revision: 1.7 $
 * $Date: 2004-09-22 00:15:44 $
 * ----------------------------------------------------------------- 
 * Programmer(s): Aaron Collier @ LLNL
 * -----------------------------------------------------------------
 * The FCVROOT module contains the routines necessary to use
 * the rootfinding feature of the CVODE module and to interface
 * with the user-supplied Fortran subroutine.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>
#include "sundialstypes.h"  /* definition of SUNDIALS type realtype */
#include "nvector.h"        /* definition of type N_Vector */
#include "cvode.h"          /* CVODE constants and prototypes */
#include "fcvode.h"         /* actual function names, prototypes and
			       global variables */
#include "fcvroot.h"        /* prototypes of interfaces to CVROOT */

/***************************************************************************/

/* Prototype of the Fortran routine */

void FCV_ROOTFN(realtype *, realtype*, realtype*);

/***************************************************************************/

void FCV_ROOTINIT(int *nrtfn, int *ier)
{
  *ier = CVodeRootInit(CV_cvodemem, (RootFn) FCVrootfunc, *nrtfn);
  CV_nrtfn = *nrtfn;

  return;
}

/***************************************************************************/

void FCV_ROOTINFO(int *nrtfn, int **info, int *ier)
{
  int *rootsfound, i;

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
