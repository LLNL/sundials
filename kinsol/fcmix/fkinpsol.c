/*
 * -----------------------------------------------------------------
 * $Revision: 1.15 $
 * $Date: 2004-12-06 20:22:19 $
 * -----------------------------------------------------------------
 * Programmer(s): Allan Taylor, Alan Hindmarsh and
 *                Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see sundials/kinsol/LICENSE.
 * -----------------------------------------------------------------
 * This routine interfaces between KINSOL and the user-supplied
 * Fortran routine FK_PSOL.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include "fkinsol.h"        /* prototypes of interfaces and global variables */
#include "kinsol.h"         /* KINSOL constants and prototypes               */
#include "kinspbcg.h"       /* prototypes of KINSPBCG interface routines     */
#include "kinspgmr.h"       /* prototypes of KINSPGMR interface routines     */
#include "nvector.h"        /* definition of type N_Vector                   */
#include "sundialstypes.h"  /* definition of type realtype                   */

/*
 * ----------------------------------------------------------------
 * prototype of the user-supplied fortran routine
 * ----------------------------------------------------------------
 */

extern void FK_PSOL(realtype*, realtype*, realtype*, realtype*, 
		    realtype*, realtype*, int*);

/*
 * ----------------------------------------------------------------
 * Function : FKIN_SPBCGSETPSOL
 * ----------------------------------------------------------------
 */

void FKIN_SPBCGSETPSOL(int *flag, int *ier)
{
  if ((*flag) == 0) KINSpbcgSetPrecSolveFn(KIN_mem, NULL);
  else KINSpbcgSetPrecSolveFn(KIN_mem, FKINPSol);
}

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

  return(retcode);
}
