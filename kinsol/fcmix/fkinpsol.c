/*
 * ----------------------------------------------------------------
 * $Revision: 1.12 $
 * $Date: 2004-07-27 23:52:49 $
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
  realtype *udata, *uscaledata, *fdata, *fscaledata, *vvdata, *ftemdata;
  int retcode;

  udata      = N_VGetArrayPointer(uu);
  uscaledata = N_VGetArrayPointer(uscale);
  fdata      = N_VGetArrayPointer(fval);
  fscaledata = N_VGetArrayPointer(fscale);
  vvdata     = N_VGetArrayPointer(vv);
  ftemdata   = N_VGetArrayPointer(ftem);
  
  FK_PSOL(udata, uscaledata, fdata, fscaledata, vvdata, ftemdata, &retcode);
 
  N_VSetArrayPointer(vvdata, vv);

  return(retcode);
}
