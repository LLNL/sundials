/****************************************************************************
 * File         : fkinbbd.c                                                 *
 * Programmers  : Allan G Taylor, Alan C. Hindmarsh, and Radu Serban @ LLNL * 
 * Version of   : 07 February 2004                                          *
 ****************************************************************************
 *                                                                          *
 * This module contains the routines necessary to interface with the        *
 * KINBBDPRE module and user-supplied Fortran routines. Generic names are   *
 * used (e.g. K_COMMFN). The routines here call the generically named       *
 * routines and provide a standard interface to the C code of the           *
 * KINBBDPRE package.                                                       *
 *                                                                          *
 ***************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "sundialstypes.h" /* definitions of type realtype                        */
#include "nvector.h"       /* definitions of type N_Vector                        */
#include "kinsol.h"        /* KINSOL constants and prototypes                     */
#include "fkinsol.h"       /* prototypes of standard interfaces, global variables */
#include "fkinbbd.h"       /* prototypes of interfaces to KINBBDPRE               */
#include "kinspgmr.h"      /* prototypes of KINSPGMR interface routines           */
#include "kinbbdpre.h"     /* prototypes of KINBBDPRE functions, macros           */

/***************************************************************************/

/* Prototypes of the user-supplied Fortran routines */

void FK_LOCFN(long int*, realtype*, realtype*);
void FK_COMMFN(long int*, realtype*);

/***************************************************************************/

void FKIN_BBDINIT(long int *nlocal, long int *mu, long int *ml, int *ier)
{
  KBBD_Data = KBBDPrecAlloc(KIN_mem, *nlocal, *mu, *ml, 0.0, FKINgloc, FKINgcomm);
  if (KBBD_Data == NULL) { *ier = -1; return; }
}

/***************************************************************************/

void FKIN_BBDSPGMR(int *maxl, int *maxlrst, int *ier)
{
  *ier = KBBDSpgmr(KIN_mem, *maxl, KBBD_Data);
  if (*ier != 0) return;

  *ier = KINSpgmrSetMaxRestarts(KIN_mem, *maxlrst);
  if (*ier != 0) return;
}

/***************************************************************************/

/* C function FKINgloc to interface between KINBBDPRE module and a Fortran 
   subroutine FKLOCFN. */

void FKINgloc(long int Nloc, N_Vector uu, N_Vector gval, void *f_data)
{
  realtype *uloc, *gloc;

  uloc = N_VGetData(uu);
  gloc = N_VGetData(gval);

  FK_LOCFN(&Nloc, uloc, gloc);

  N_VSetData(gloc, gval);

}

/***************************************************************************/

/* C function FKINgcomm to interface between KINBBDPRE module and a Fortran 
   subroutine FKCOMMFN. */


void FKINgcomm(long int Nloc, N_Vector uu, void *f_data)
{
  realtype *uloc;

  uloc = N_VGetData(uu);
  
  FK_COMMFN(&Nloc, uloc);

}

/***************************************************************************/

/* C function FKINBBDOPT to access optional outputs from KBBD_Data */

void FKIN_BBDOPT(long int *lenrpw, long int *lenipw, long int *nge)
{
  KBBDPrecGetIntWorkSpace(KBBD_Data, lenipw);
  KBBDPrecGetRealWorkSpace(KBBD_Data, lenrpw);
  KBBDPrecGetNumGfnEvals(KBBD_Data, nge);
}


/***************************************************************************/

/* C function FKINBBDFREE to interface to KBBDFree, to free memory 
   created by KBBDAlloc */

void FKIN_BBDFREE()
{
  KBBDPrecFree(KBBD_Data);
}



