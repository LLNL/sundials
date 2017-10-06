/*
 * -----------------------------------------------------------------
 * $Revision: 4075 $
 * $Date: 2014-04-24 10:46:58 -0700 (Thu, 24 Apr 2014) $
 * -----------------------------------------------------------------
 * Programmer(s): Allan Taylor, Alan Hindmarsh and
 *                Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * LLNS Copyright Start
 * Copyright (c) 2014, Lawrence Livermore National Security
 * This work was performed under the auspices of the U.S. Department 
 * of Energy by Lawrence Livermore National Laboratory in part under 
 * Contract W-7405-Eng-48 and in part under Contract DE-AC52-07NA27344.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * LLNS Copyright End
 * -----------------------------------------------------------------
 * This file contains the interfaces between KINSOL and the
 * user-supplied Fortran routines FK_PSET and FK_PSOL.
 *
 * The C function FKINPSet is used to interface between KINSOL and
 * the Fortran user-supplied preconditioner setup routine.
 *
 * The C function FKINPSol is used to interface between KINSOL and
 * the Fortran user-supplied preconditioner solve routine.
 *
 * Note: The use of the generic names FK_PSET and FK_PSOL below.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include "fkinsol.h"
#include "kinsol_impl.h"

#include <kinsol/kinsol_spils.h>

/*
 * ----------------------------------------------------------------
 * prototype of the user-supplied fortran routine
 * ----------------------------------------------------------------
 */

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

extern void FK_PSET(realtype*, realtype*, realtype*, realtype*, 
		    realtype*, realtype*, int*);
extern void FK_PSOL(realtype*, realtype*, realtype*, realtype*, 
		    realtype*, realtype*, int*);

#ifdef __cplusplus
}
#endif

/*
 * ----------------------------------------------------------------
 * Function : FKIN_SPILSSETPREC
 * ----------------------------------------------------------------
 */

void FKIN_SPILSSETPREC(int *flag, int *ier)
{
  if ((*flag) == 0) {
    *ier = KINSpilsSetPreconditioner(KIN_kinmem, NULL, NULL);
  } else {
    *ier = KINSpilsSetPreconditioner(KIN_kinmem, FKINPSet, FKINPSol);
  }

  return;
}

/*
 * ----------------------------------------------------------------
 * Function : FKINPSet
 * ----------------------------------------------------------------
 * C function FKINPSet is used to interface between FK_PSET and
 * the user-supplied Fortran preconditioner setup routine.
 * ----------------------------------------------------------------
 */

int FKINPSet(N_Vector uu, N_Vector uscale,
             N_Vector fval, N_Vector fscale,
             void *user_data,
             N_Vector vtemp1, N_Vector vtemp2)
{
  realtype *udata, *uscaledata, *fdata, *fscaledata, *vtemp1data, *vtemp2data;
  int retcode;

  udata = uscaledata = fdata = fscaledata = vtemp1data = vtemp2data = NULL;

  udata      = N_VGetArrayPointer(uu);
  uscaledata = N_VGetArrayPointer(uscale);
  fdata      = N_VGetArrayPointer(fval);
  fscaledata = N_VGetArrayPointer(fscale);
  vtemp1data = N_VGetArrayPointer(vtemp1);
  vtemp2data = N_VGetArrayPointer(vtemp2);

  FK_PSET(udata, uscaledata, fdata, fscaledata, vtemp1data, vtemp2data, &retcode);

 /* Note: There is no need to use N_VSetArrayPointer since we are not getting back
    any information that should go into an N_Vector */

 return(retcode);
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
             N_Vector vv, void *user_data,
             N_Vector ftem)
{
  realtype *udata, *uscaledata, *fdata, *fscaledata, *vvdata, *ftemdata;
  int retcode;

  udata = uscaledata = fdata = fscaledata = vvdata = ftemdata = NULL;

  udata      = N_VGetArrayPointer(uu);
  uscaledata = N_VGetArrayPointer(uscale);
  fdata      = N_VGetArrayPointer(fval);
  fscaledata = N_VGetArrayPointer(fscale);
  vvdata     = N_VGetArrayPointer(vv);
  ftemdata   = N_VGetArrayPointer(ftem);

  FK_PSOL(udata, uscaledata, fdata, fscaledata, vvdata, ftemdata, &retcode);

  return(retcode);
}
