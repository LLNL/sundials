/****************************************************************************
 * File         : fcvbbd.c                                                  *
 * Programmers  : Alan C. Hindmarsh and Radu Serban @ LLNL                  * 
 * Version of   : 26 June 2002                                              *
 *                                                                          *
 ****************************************************************************
 *                                                                          *
 * This module contains the routines necessary to interface with the        *
 * CVBBDPRE module and user-supplied Fortran routines.                      *
 * The routines here call the generically named routines and provide a      *
 * standard interface to the C code of the CVBBDPRE package.                *
 * The routine FCV_BBDIN0 has a counterpart FCV_BBDIN1 in a separate        *
 * file: fcvbbdin1.c.                                                       *
 *                                                                          *
 ***************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "sundialstypes.h" /* definitions of types realtype and integertype   */
#include "nvector.h"       /* definitions of type N_Vector                    */
#include "cvode.h"         /* CVODE constants and prototypes                  */
#include "fcmixpar.h"      /* definition of global F2C_machEnv variable       */
#include "fcvode.h"        /* actual function names, prototypes, global vars. */
#include "fcvbbd.h"        /* prototypes of interfaces to CVBBDPRE            */
#include "cvspgmr.h"       /* prototypes of CVSPGMR interface routines        */
#include "cvbbdpre.h"      /* prototypes of CVBBDPRE functions, macros        */

/***************************************************************************/

/* Prototypes of the Fortran routines */

void FCV_GLOCFN(integertype*, realtype*, realtype*, realtype*);
void FCV_COMMFN(integertype*, realtype*, realtype*);

/***************************************************************************/

void FCV_BBDIN0(integertype *Nloc, integertype *mudq, integertype *mldq, 
                integertype *mu, integertype *ml, realtype* dqrely,
                int *pretype, int *gstype, int *maxl, realtype *delt, int *ier)
{

  /* First call CVBBDAlloc to initialize CVBBDPRE module:
     *Nloc       is the local vector size
     *mudq,*mldq are the half-bandwidths for computing preconditioner blocks
     *mu, *ml    are the half-bandwidths of the retained preconditioner blocks
     *dqrely     is the difference quotient relative increment factor
     CVgloc      is a pointer to the CVLocalFn function
     CVcfn       is a pointer to the CVCommFn function
     NULL        is the pointer to f_data                             */

  CVBBD_Data = CVBBDAlloc (*Nloc, *mudq, *mldq, *mu, *ml, *dqrely, 
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
     NULL        is a pointer to the Jtimes routine
     NULL        is the pointer to jac_data                               */

  *ier = CVSpgmr (CV_cvodemem, *pretype, *gstype, *maxl, *delt, 
                  CVBBDPrecon, CVBBDPSol, CVBBD_Data, NULL, NULL);

}

/***************************************************************************/

void FCV_REINBBD0(integertype *Nloc, integertype *mudq, integertype *mldq, 
                  integertype *mu, integertype *ml,
                  realtype* dqrely, int *pretype, int *gstype, int *maxl,
                  realtype *delt, int *ier)
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
     NULL        is a pointer to the Jtimes routine
     NULL        is the pointer to jac_data                               */

  *ier = CVReInitSpgmr(CV_cvodemem, *pretype, *gstype, *maxl, *delt, 
                       CVBBDPrecon, CVBBDPSol, CVBBD_Data, NULL, NULL);

}

/***************************************************************************/

/* C function CVgloc to interface between CVBBDPRE module and a Fortran 
   subroutine CVLOCFN. */

void CVgloc(integertype Nloc, realtype t, realtype *yloc, realtype *gloc,
            void *f_data)
{

  FCV_GLOCFN(&Nloc, &t, yloc, gloc);

}

/***************************************************************************/

/* C function CVcfn to interface between CVBBDPRE module and a Fortran 
   subroutine CVCOMMF. */


void CVcfn(integertype Nloc, realtype t, N_Vector y, void *f_data)
{
  realtype *yloc;

  yloc = N_VGetData(y);

  FCV_COMMFN(&Nloc, &t, yloc);

  /* Note: there is no need to use N_VSetData here because the data in
     the N_Vector y is input only to FCV_COMMFN.  */

}

/***************************************************************************/

/* C function FCVBBDOPT to access optional outputs from CVBBD_Data */

void FCV_BBDOPT(integertype *lenrpw, integertype *lenipw, integertype *nge)
{
  CVBBDData pdata;
  pdata = (CVBBDData)(CVBBD_Data);
  *lenrpw = CVBBD_RPWSIZE(pdata);
  *lenipw = CVBBD_IPWSIZE(pdata);
  *nge = CVBBD_NGE(pdata);
}


/***************************************************************************/

/* C function FCVBBDF to interface to CVBBDFree, to free memory 
   created by CVBBDAlloc */

void FCV_BBDF ()
{
  CVBBDFree(CVBBD_Data);
}
