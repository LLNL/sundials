/****************************************************************************
 *                                                                          *
 * File          : fkinbbd.c                                                *
 * Programmers   : Allan G Taylor and Alan C. Hindmarsh @ LLNL              * 
 * Version of    : 18 January 2001                                          *
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
#include "llnltyps.h"  /* definitions of types real and integer               */
#include "nvector.h"   /* definitions of type N_Vector and vector macros      */
#include "kinsol.h"    /* KINSOL constants and prototypes                     */
#include "fkinsol.h"   /* prototypes of standard interfaces, global variables */
#include "fkinbbd.h"   /* prototypes of interfaces to KINBBDPRE               */
#include "kinbbdpre.h" /* prototypes of KINBBDPRE functions, macros           */

/***************************************************************************/

void F_KINBBDINIT0 (int *maxl, int *maxlrst, int *msbpre,
		    integer *mu, integer *ml, int *ier)
{
  integer Nlocal;

  /* First call KBBDAlloc to initialize KINBBDPRE module:
     *Nlocal       is the local vector size
     *mu, *ml      are the half-bandwidths for the preconditioner blocks
     0.0           is the value for dq_rel_uu -- 0.0 triggers the default
     KINgloc       is a pointer to the KINLocalFn function
     KINgcomm      is a pointer to the KINCommFn function
     NULL          is the pointer to f_data.  NULL is used here since handling
                   of local data in the func routines is done only in Fortran */

  Nlocal = ((machEnvType)KIN_machEnv)->local_vec_length;
  
  KBBD_Data = KBBDAlloc (Nlocal, *mu, *ml, 0.0, KINgloc, KINgcomm, NULL, 
                         KIN_kmem, KIN_machEnv);
  if (KBBD_Data == NULL) { *ier = -1; return; }

  /* Call KINSpgmr to specify the SPGMR linear solver:
     KIN_kmem    is the pointer to the KINSOL memory block
     *maxl       is the maximum Krylov dimension
     *maxlrst    is the max number of linear solver restarts
     *msbpre     is the max number of steps w/o calling the precond setup fcn
     KBBDPrecon  is a pointer to the preconditioner setup routine
     KBBDPSol    is a pointer to the preconditioner solve routine
     KBBD_Data   is the pointer to P_data                                  */

  KINSpgmr (KIN_kmem, *maxl, *maxlrst, *msbpre,
           KBBDPrecon, KBBDPSol, NULL, KBBD_Data);

  *ier = 0;
}

/***************************************************************************/

/* C function KINgloc to interface between KINBBDPRE module and a Fortran 
   subroutine KLOCFN. */

void KINgloc(integer Nloc, N_Vector uu, N_Vector gval, void *f_data)
{
  real *uloc, *gloc;

  uloc = N_VDATA(uu);
  gloc = N_VDATA(gval) ;
  K_LOCFN(&Nloc, uloc, gloc);

}

/***************************************************************************/

/* C function KINgcomm to interface between KINBBDPRE module and a Fortran 
   subroutine KCOMMFN. */


void KINgcomm(integer Nloc, real *uloc, void *f_data)
{

  K_COMMFN(&Nloc, uloc);

}

/***************************************************************************/

/* C function FKINBBDOPT to access optional outputs from KBBD_Data */

void F_KINBBDOPT(integer *lenrpw, integer *lenipw, integer *nge)
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



