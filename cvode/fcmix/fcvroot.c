/*
 * -----------------------------------------------------------------
 * $Revision: 1.3 $
 * $Date: 2004-06-09 18:10:06 $
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

void FCV_ROOTFREE(void)
{
  CVodeRootInit(CV_cvodemem, NULL, 0);

  return;
}

/***************************************************************************/

void FCVrootfunc(realtype t, N_Vector y, realtype *gout, void *g_data)
{
  realtype *ydata;

  ydata = N_VGetData(y);

  FCV_ROOTFN(&t, ydata, gout);

  return;
}
