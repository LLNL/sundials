/****************************************************************************
 *                                                                          *
 * File         : fkinbbd.c                                                 *
 * Programmers  : Allan G Taylor, Alan C. Hindmarsh, and Radu Serban @ LLNL * 
 * Version of   : 5 August 2003                                             *
 *                                                                          *
 ****************************************************************************
 *                                                                          *
 * This module contains the routines necessary to interface with the        *
 * KINBBDPRE module and user-supplied Fortran routines. Generic names are   *
 * used (e.g. K_COMMFN, see fcmixpar.h). The routines here call the         *
 * generically named routines and provide a standard interface to the C code*
 * of the KINBBDPRE package.                                                *
 *                                                                          *
 ***************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "sundialstypes.h" /* definitions of types realtype and integertype   */
#include "nvector.h"   /* definitions of type N_Vector                        */
#include "kinsol.h"    /* KINSOL constants and prototypes                     */
#include "fcmixpar.h"  /* definition of global F2C_machEnv variable           */
#include "fkinsol.h"   /* prototypes of standard interfaces, global variables */
#include "fkinbbd.h"   /* prototypes of interfaces to KINBBDPRE               */
#include "kinspgmr.h"  /* prototypes of KINSPGMR interface routines           */
#include "kinbbdpre.h" /* prototypes of KINBBDPRE functions, macros           */

/***************************************************************************/

/* Prototypes of the user-supplied Fortran routines */

void K_LOCFN(integertype*, realtype*, realtype*);
void K_COMMFN(integertype*, realtype*);

/***************************************************************************/

void F_KINBBDINIT(integertype *nlocal, int *maxl, int *maxlrst,
                  integertype *mu, integertype *ml, int *ier)
{

  KBBD_Data = KBBDPrecAlloc(KIN_mem, *nlocal, *mu, *ml, 0.0, KINgloc, KINgcomm);
  if (KBBD_Data == NULL) { *ier = -1; return; }

  *ier = KBBDSpgmr(KIN_mem, *maxl, KBBD_Data);
  if (*ier != 0) return;

  *ier = KINSpgmrSetMaxRestarts(KIN_mem, *maxlrst);
  if (*ier != 0) return;
}

/***************************************************************************/

/* C function KINgloc to interface between KINBBDPRE module and a Fortran 
   subroutine KLOCFN. */

void KINgloc(integertype Nloc, N_Vector uu, N_Vector gval, void *f_data)
{
  realtype *uloc, *gloc;

  uloc = N_VGetData(uu);
  gloc = N_VGetData(gval);

  K_LOCFN(&Nloc, uloc, gloc);

  N_VSetData(gloc, gval);

}

/***************************************************************************/

/* C function KINgcomm to interface between KINBBDPRE module and a Fortran 
   subroutine KCOMMFN. */


void KINgcomm(integertype Nloc, N_Vector uu, void *f_data)
{
  realtype *uloc;

  uloc = N_VGetData(uu);
  
  K_COMMFN(&Nloc, uloc);

}

/***************************************************************************/

/* C function FKINBBDOPT to access optional outputs from KBBD_Data */

void F_KINBBDOPT(integertype *lenrpw, integertype *lenipw, int *nge)
{
  KBBDPrecGetIntWorkSpace(KBBD_Data, lenipw);
  KBBDPrecGetRealWorkSpace(KBBD_Data, lenrpw);
  KBBDPrecGetNumGfnEvals(KBBD_Data, nge);
}


/***************************************************************************/

/* C function FKINBBDFREE to interface to KBBDFree, to free memory 
   created by KBBDAlloc */

void F_KINBBDFREE()
{
  KBBDPrecFree(KBBD_Data);
}



