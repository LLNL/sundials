/*******************************************************************
 * File          : fkinpreco.c                                     *
 * Programmers   : Allan G Taylor, Alan C. Hindmarsh, and          *
 *                 Radu Serban @ LLNL                              *
 * Version of    : 07 February 2004                                *
 *-----------------------------------------------------------------*
 * This C function FKINPSet is to interface between KINSOL and the *
 * Fortran user-supplied preconditioner setup routine.             *
 * Note the use of generic name FK_PSET below.                     *
 *******************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "sundialstypes.h" /* definitions of type realtype               */
#include "nvector.h"       /* definitions of type N_Vector               */
#include "kinsol.h"        /* KINSOL constants and prototypes            */
#include "kinspgmr.h"      
#include "fkinsol.h"       /* prototypes of interfaces, global variables */

/*********************************************************************/

/* Prototype of the user-supplied Fortran routine */
void FK_PSET(realtype*, realtype*, realtype*, realtype*, 
             realtype*, realtype*, int*);


/***************************************************************************/

void FKIN_SPGMRSETPSET(int *flag, int *ier)
{
  if (*flag == 0) KINSpgmrSetPrecSetupFn(KIN_mem, NULL);
  else            KINSpgmrSetPrecSetupFn(KIN_mem, FKINPSet);
}

/*********************************************************************/

/* C function FKINPSet to interface between KINSpgmr and FKPSET, 
   the user-supplied Fortran preconditioner setup routine. */

int FKINPSet(N_Vector uu, N_Vector uscale,
             N_Vector fval, N_Vector fscale,
             void *P_data,
             N_Vector vtemp1, N_Vector vtemp2)
{
  realtype *udata,*uscaledata, *fdata, *fscaledata, *vtemp1data, *vtemp2data;
  int retcode;
  
  udata        = N_VGetData(uu);
  uscaledata   = N_VGetData(uscale);
  fdata        = N_VGetData(fval);
  fscaledata   = N_VGetData(fscale);
  vtemp1data   = N_VGetData(vtemp1);
  vtemp2data   = N_VGetData(vtemp2);
  
  FK_PSET(udata, uscaledata, fdata, fscaledata, vtemp1data, vtemp2data, &retcode);

 /* Note: there is no need to use N_VSetData since we are not getting back any
    information that should go into an N_Vector */

 return(retcode);
}
