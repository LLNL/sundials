/****************************************************************************
 * File         : fcvbbdin1.c                                               *
 * Programmers  : Alan C. Hindmarsh and Radu Serban @ LLNL                  * 
 * Version of   : 27 March 2002                                             *
 *                                                                          *
 ****************************************************************************
 *                                                                          *
 * Fortran/C and C/Fortran interface routine for use of CVODE with CVBBDPRE *
 * preconditioner module, for the case of a user-supplied Jtimes routine.   *
 *                                                                          *
 ***************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "llnltyps.h" /* definitions of types real and integer           */
#include "nvector.h"  /* definition of type N_Vector                     */
#include "cvode.h"    /* CVODE constants and prototypes                  */
#include "fcmixpar.h" /* definition of global F2C_machEnv variable       */
#include "fcvode.h"   /* actual function names, prototypes, global vars. */
#include "fcvbbd.h"   /* prototypes of interfaces to PVBBDPRE            */
#include "cvspgmr.h"  /* prototypes of CVSPGMR interface routines        */
#include "cvbbdpre.h" /* prototypes of PVBBDPRE functions, macros        */

/***************************************************************************/

void FCV_BBDIN1(integer *Nloc, integer *mudq, integer *mldq, 
                integer *mu, integer *ml,
                real* dqrely, int *pretype, int *gstype, int *maxl,
                real *delt, int *ier)
{

  /* First call PVBBDAlloc to initialize PVBBDPRE module:
     *Nloc       is the local vector size
     *mudq,*mldq are the half-bandwidths for computing preconditioner blocks
     *mu, *ml    are the half-bandwidths of the retained preconditioner blocks
     *dqrely     is the difference quotient relative increment factor
     CVgloc      is a pointer to the CVLocalFn function
     CVcfn       is a pointer to the CVCommFn function
     NULL        is the pointer to f_data                             */

  CVBBD_Data = CVBBDAlloc(*Nloc, *mudq, *mldq, *mu, *ml, *dqrely, 
                          CVgloc, CVcfn, NULL, CV_cvodemem);
  if (CVBBD_Data == NULL) { *ier = -1; return; }

  /* Call CVSpgmr to specify the SPGMR linear solver:
     CV_cvodemem is the pointer to the CVODE memory block
     *pretype    is the preconditioner type
     *gstype     is the Gram-Schmidt process type
     *maxl       is the maximum Krylov dimension
     *delt       is the linear convergence tolerance factor
     CVBBDPrecon is a pointer to the preconditioner setup routine
     CVBBDPSol   is a pointer to the preconditioner solve routine
     CVBBD_Data  is the pointer to P_data
     CVjtimes    is a pointer to the Jtimes routine
     NULL        is the pointer to jac_data                               */

  *ier = CVSpgmr(CV_cvodemem, *pretype, *gstype, *maxl, *delt, 
                 CVBBDPrecon, CVBBDPSol, CVBBD_Data, CVJtimes, NULL);

}

/***************************************************************************/

void FCV_REINBBD1(integer *Nloc, integer *mudq, integer *mldq, 
                  integer *mu, integer *ml,
                  real* dqrely, int *pretype, int *gstype, int *maxl,
                  real *delt, int *ier)
{
  int flag;

  /* First call CVReInitBBD to re-initialize CVBBDPRE module:
     CVBBD_Data  is the pointer to P_data
     *Nloc       is the local vector size
     *mudq,*mldq are the half-bandwidths for computing preconditioner blocks
     *mu, *ml    are the half-bandwidths of the retained preconditioner blocks
     *dqrely     is the difference quotient relative increment factor
     CVgloc      is a pointer to the CVLocalFn function
     CVcfn       is a pointer to the CVCommFn function
     NULL        is the pointer to f_data                             */

  flag = CVReInitBBD(CVBBD_Data, *Nloc, *mudq, *mldq, *mu, *ml,
                     *dqrely, CVgloc, CVcfn, NULL);

  /* Call CVReInitSpgmr to re-initialize the SPGMR linear solver:
     CV_cvodemem is the pointer to the CVODE memory block
     *pretype    is the preconditioner type
     *gstype     is the Gram-Schmidt process type
     *maxl       is the maximum Krylov dimension
     *delt       is the linear convergence tolerance factor
     CVBBDPrecon is a pointer to the preconditioner setup routine
     CVBBDPSol   is a pointer to the preconditioner solve routine
     CVBBD_Data  is the pointer to P_data
     CVjtimes    is a pointer to the Jtimes routine
     NULL        is the pointer to jac_data                               */

  *ier = CVReInitSpgmr(CV_cvodemem, *pretype, *gstype, *maxl, *delt, 
                       CVBBDPrecon, CVBBDPSol, CVBBD_Data, CVJtimes, NULL);

}
