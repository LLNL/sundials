/*
 * -----------------------------------------------------------------
 * $Revision: 1.19 $
 * $Date: 2004-11-30 21:20:57 $
 * -----------------------------------------------------------------
 * Programmer(s): Allan Taylor, Alan Hindmarsh, Radu Serban, and
 *                Aaron Collier @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see sundials/kinsol/LICENSE.
 * -----------------------------------------------------------------
 * This module contains the routines necessary to interface with
 * the KINBBDPRE module and user-supplied Fortran routines. Generic
 * names are used (e.g. FK_COMMFN). The routines here call the
 * generically named routines and provide a standard interface to
 * the C code of the KINBBDPRE package.
 * ----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include "fkinbbd.h"        /* prototypes of interfaces to KINBBDPRE        */
#include "fkinsol.h"        /* prototypes of standard interfaces and global
                               variables                                    */
#include "kinbbdpre.h"      /* prototypes of KINBBDPRE functions and macros */
#include "kinsol.h"         /* KINSOL constants and prototypes              */
#include "kinspgmr.h"       /* prototypes of KINSPGMR interface routines    */
#include "nvector.h"        /* definition of type N_Vector                  */
#include "sundialstypes.h"  /* definition of type realtype                  */

/*
 * ----------------------------------------------------------------
 * private constants
 * ----------------------------------------------------------------
 */

#define ZERO RCONST(0.0)

/*
 * ----------------------------------------------------------------
 * prototypes of the user-supplied fortran routines
 * ----------------------------------------------------------------
 */

extern void FK_LOCFN(long int*, realtype*, realtype*);
extern void FK_COMMFN(long int*, realtype*);

/*
 * ----------------------------------------------------------------
 * Function : FKIN_BBDINIT
 * ----------------------------------------------------------------
 */

void FKIN_BBDINIT(long int *nlocal, long int *mu, long int *ml, int *ier)
{
  KBBD_Data = KINBBDPrecAlloc(KIN_mem, *nlocal, *mu, *ml, ZERO, FKINgloc, FKINgcomm);
  if (KBBD_Data == NULL) { *ier = -1; return; }
}

/*
 * ----------------------------------------------------------------
 * Function : FKIN_BBDSPGMR
 * ----------------------------------------------------------------
 */

void FKIN_BBDSPGMR(int *maxl, int *maxlrst, int *ier)
{
  *ier = KINBBDSpgmr(KIN_mem, *maxl, KBBD_Data);
  if (*ier != KINSPGMR_SUCCESS) return;

  *ier = KINSpgmrSetMaxRestarts(KIN_mem, *maxlrst);
  if (*ier != KINSPGMR_SUCCESS) return;
}

/*
 * ----------------------------------------------------------------
 * Function : FKINgloc
 * ----------------------------------------------------------------
 * C function FKINgloc is the interface between the KINBBDPRE
 * module and the Fortran subroutine FK_LOCFN.
 * ----------------------------------------------------------------
 */

void FKINgloc(long int Nloc, N_Vector uu, N_Vector gval, void *f_data)
{
  realtype *uloc, *gloc;

  uloc = N_VGetArrayPointer(uu);
  gloc = N_VGetArrayPointer(gval);

  FK_LOCFN(&Nloc, uloc, gloc);

  N_VSetArrayPointer(gloc, gval);
}

/*
 * ----------------------------------------------------------------
 * Function : FKINgcomm
 * ----------------------------------------------------------------
 * C function FKINgcomm is the interface between the KINBBDPRE
 * module and the Fortran subroutine FK_COMMFN.
 * ----------------------------------------------------------------
 */

void FKINgcomm(long int Nloc, N_Vector uu, void *f_data)
{
  realtype *uloc;

  uloc = N_VGetArrayPointer(uu);
  
  FK_COMMFN(&Nloc, uloc);
}

/*
 * ----------------------------------------------------------------
 * Function : FKIN_BBDOPT
 * ----------------------------------------------------------------
 * C function FKIN_BBDOPT is used to access optional outputs stored
 * in the KBBD_Data data structure.
 * ----------------------------------------------------------------
 */

void FKIN_BBDOPT(long int *lenrpw, long int *lenipw, long int *nge)
{
  KINBBDPrecGetWorkSpace(KBBD_Data, lenrpw, lenipw);
  KINBBDPrecGetNumGfnEvals(KBBD_Data, nge);
}

/*
 * ----------------------------------------------------------------
 * Function : FKIN_BBDFREE
 * ----------------------------------------------------------------
 * C function FKIN_BBDFREE is used to interface with KINBBDPrecFree
 * in order to free memory created by KINBBDPrecAlloc.
 * ----------------------------------------------------------------
 */

void FKIN_BBDFREE(void)
{
  KINBBDPrecFree(KBBD_Data);
}
