/*
 * ----------------------------------------------------------------
 * $Revision: 1.9 $
 * $Date: 2004-05-03 21:24:50 $
 * ----------------------------------------------------------------
 * Programmer(s): Allan Taylor, Alan Hindmarsh and
 *                Radu Serban @ LLNL
 * ----------------------------------------------------------------
 * This routine interfaces between KINSOL and the user-supplied
 * Fortran routine FK_PSOL.
 * ----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>
#include "sundialstypes.h"  /* definition of type realtype */
#include "nvector.h"        /* definition of type N_Vector */
#include "kinsol.h"         /* KINSOL constants and prototypes */
#include "kinspgmr.h"
#include "fkinsol.h"        /* prototypes of interfaces and global variables */

/*
 * ----------------------------------------------------------------
 * prototype of the user-supplied fortran routine
 * ----------------------------------------------------------------
 */

void FK_PSOL(realtype*, realtype*, realtype*, realtype*, 
             realtype*, realtype*, int*);

/*
 * ----------------------------------------------------------------
 * Function : FKIN_SPGMRSETPSOL
 * ----------------------------------------------------------------
 */

void FKIN_SPGMRSETPSOL(int *flag, int *ier)
{
  if ((*flag) == 0) KINSpgmrSetPrecSolveFn(KIN_mem, NULL);
  else KINSpgmrSetPrecSolveFn(KIN_mem, FKINPSol);
}

/*
 * ----------------------------------------------------------------
 * Function : FKINPSol
 * ----------------------------------------------------------------
 * C function FKINPSol is used to interface between FK_PSOL and
 * the user-supplied Fortran preconditioner solve routine.
 * ----------------------------------------------------------------
 */

int FKINPSol(N_Vector uu, N_Vector uscale, 
             N_Vector fval, N_Vector fscale, 
             N_Vector vv, void *P_data,
             N_Vector ftem)
{
  realtype *udata,*uscaledata, *fdata,*fscaledata, *vvdata, *ftemdata;
  int retcode;

  udata      = N_VGetData(uu);
  uscaledata = N_VGetData(uscale);
  fdata      = N_VGetData(fval);
  fscaledata = N_VGetData(fscale);
  vvdata     = N_VGetData(vv);
  ftemdata   = N_VGetData(ftem);
  
  FK_PSOL(udata, uscaledata, fdata, fscaledata, vvdata, ftemdata, &retcode);
 
  N_VSetData(vvdata, vv);

  return(retcode);
}
