/********************************************************************
 *                                                                  *
 * File        : fkinbbdinit1.c                                     *
 * Programmers : Allan G Taylor, Alan C. Hindmarsh, and Radu Serban *
 *               @ LLNL                                             *
 * Version of  : 29 July 2002                                       *
 *------------------------------------------------------------------*
 * This is the implementation file for the routine fkinbbdinit1,    *
 * used as part of the KINSOL Fortran/C interface package.          *
 *                                                                  *
 * The routine fkinbbdinit1 does the initialization necessary for   *
 * the KINBBDPRE package. In particular, this version assumes that  *
 * the user has supplied an ATimes routine FATIMES in Fortran       *
 * that is called by KINUATimes in the C code. This module uses     *
 * generic names for Fortran names (e.g. K_COMMFN, ... )            *
 *******************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "sundialstypes.h" /* definitions of types realtype and integertype   */
#include "nvector.h"   /* definitions of type N_Vector                        */
#include "kinsol.h"    /* KINSOL constants and prototypes                     */
#include "fcmixpar.h"  /* definition of global F2C_machEnv variable           */
#include "fkinsol.h"   /* prototypes of standard interfaces, global variables */
#include "fkinbbd.h"   /* prototypes of interfaces to KINBBDPRE               */
#include "kinspgmr.h"  /* prototypes of KINSPGMR interface routines           */
#include "kinbbdpre.h" /* prototypes of KINBBDPRE                             */

/***************************************************************************/

void F_KINBBDINIT1 (integertype *nlocal, int *maxl, int *maxlrst, int *msbpre,
                    integertype *mu, integertype *ml, int *ier)
{

  /* First call KBBDAlloc to initialize KINBBDPRE module:
     *nlocal       is the local vector size
     *mu, *ml      are the half-bandwidths for the preconditioner blocks
     0.0           is the value for dq_rel_uu -- 0.0 triggers the default
     KINgloc       is a pointer to the KINLocalFn function
     KINgcomm      is a pointer to the KINCommFn function
     NULL          is the pointer to f_data.   NULL is used here since handling
                   of local data in the func routines is done only in Fortran */

  KBBD_Data = KBBDAlloc (*nlocal, *mu, *ml, 0.0, KINgloc, KINgcomm, NULL,
                         KIN_kmem);
  if (KBBD_Data == NULL) { *ier = -1; return; }
  
  /* Call KINSpgmr to specify the SPGMR linear solver:
     KIN_kmem is the pointer to the KINSOL memory block
     *maxl       is the maximum Krylov dimension
     *maxlrst    is the max number of linear solver restarts
     *msbpre     is the max number of steps w/o calling the precond setup fcn
     KBBDPrecon  is a pointer to the preconditioner setup routine
     KBBDPSol    is a pointer to the preconditioner solve routine
     KINUAtimes  is a pointer to the standard interface routine to the 
                 user-supplied Atimes routine FATIMES
     KBBD_Data   is the pointer to P_data                                  */
  
  KINSpgmr (KIN_kmem, *maxl, *maxlrst, *msbpre,
            KBBDPrecon, KBBDPSol, KINUAtimes, KBBD_Data);
  
  *ier = 0;
}
