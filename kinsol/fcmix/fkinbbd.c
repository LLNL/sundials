/****************************************************************************
 *                                                                          *
 * File         : fkinbbd.c                                                 *
 * Programmers  : Allan G Taylor, Alan C. Hindmarsh, and Radu Serban @ LLNL * 
 * Version of   : 27 June 2002                                              *
 *                                                                          *
 ****************************************************************************
 *                                                                          *
 * This module contains the routines necessary to interface with the        *
 * KINBBDPRE module and user-supplied Fortran routines. Generic names are   *
 * used (e.g. K_COMMFN, see fcmixpar.h). The routines here call the         *
 * generically named routines and provide a standard interface to the C code*
 * of the KINBBDPRE package.  The routine F_KINBBDINIT0 has a counterpart   *
 * F_KINBBDINIT1 in a separate file: fkinbbdinit1.c .                       *
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

/* Prototypes of the Fortran routines */

void K_LOCFN(integertype*, realtype*, realtype*);
void K_COMMFN(integertype*, realtype*);

/***************************************************************************/

void F_KINBBDINIT0(integertype *nlocal, int *maxl, int *maxlrst, int *msbpre,
                   integertype *mu, integertype *ml, int *ier)
{
  /* First call KBBDAlloc to initialize KINBBDPRE module:
     *mu, *ml      are the half-bandwidths for the preconditioner blocks
     0.0           is the value for dq_rel_uu -- 0.0 triggers the default
     KINgloc       is a pointer to the KINLocalFn function
     KINgcomm      is a pointer to the KINCommFn function
     NULL          is the pointer to f_data.  NULL is used here since handling
                   of local data in the func routines is done only in Fortran */

  KBBD_Data = KBBDAlloc(*nlocal, *mu, *ml, 0.0, KINgloc, KINgcomm, NULL, KIN_kmem);
  if (KBBD_Data == NULL) { *ier = -1; return; }

  /* Call KINSpgmr to specify the SPGMR linear solver:
     KIN_kmem    is the pointer to the KINSOL memory block
     *maxl       is the maximum Krylov dimension
     *maxlrst    is the max number of linear solver restarts
     *msbpre     is the max number of steps w/o calling the precond setup fcn
     KBBDPrecon  is a pointer to the preconditioner setup routine
     KBBDPSol    is a pointer to the preconditioner solve routine
     KBBD_Data   is the pointer to P_data                                  */

  KINSpgmr(KIN_kmem, *maxl, *maxlrst, *msbpre,
           KBBDPrecon, KBBDPSol, NULL, KBBD_Data);

  *ier = 0;
}

/***************************************************************************/

/* C function KINgloc to interface between KINBBDPRE module and a Fortran 
   subroutine KLOCFN. */

void KINgloc(integertype Nloc, N_Vector uu, N_Vector gval, void *f_data)
{
  realtype *uloc, *gloc;

  uloc = N_VGetData(uu);
  gloc = N_VGetData(gval) ;

  K_LOCFN(&Nloc, uloc, gloc);

  N_VSetData(gloc, gval);

}

/***************************************************************************/

/* C function KINgcomm to interface between KINBBDPRE module and a Fortran 
   subroutine KCOMMFN. */


void KINgcomm(integertype Nloc, realtype *uloc, void *f_data)
{

  K_COMMFN(&Nloc, uloc);

}

/***************************************************************************/

/* C function FKINBBDOPT to access optional outputs from KBBD_Data */

void F_KINBBDOPT(integertype *lenrpw, integertype *lenipw, integertype *nge)
{
  KBBDData pdata;
  pdata = (KBBDData)(KBBD_Data);
  *lenrpw = KBBD_RPWSIZE(pdata);
  *lenipw = KBBD_IPWSIZE(pdata);
  *nge = KBBD_NGE(pdata);
}


/***************************************************************************/

/* C function FKINBBDFREE to interface to KBBDFree, to free memory 
   created by KBBDAlloc */

void F_KINBBDFREE()
{
  KBBDFree(KBBD_Data);
}



