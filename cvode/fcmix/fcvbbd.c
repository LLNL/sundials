/****************************************************************************
 * File         : fcvbbd.c                                                  *
 * Programmers  : Alan C. Hindmarsh and Radu Serban @ LLNL                  * 
 * Version of   : 31 July 2003                                              *
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

void FCV_BBDIN(integertype *Nloc, integertype *mudq, integertype *mldq, 
               integertype *mu, integertype *ml, realtype* dqrely,
               int *pretype, int *gstype, int *maxl, realtype *delt, int *ier)
{

  /* First call CVBBDPrecAlloc to initialize CVBBDPRE module:
     *Nloc       is the local vector size
     *mudq,*mldq are the half-bandwidths for computing preconditioner blocks
     *mu, *ml    are the half-bandwidths of the retained preconditioner blocks
     *dqrely     is the difference quotient relative increment factor
     CVgloc      is a pointer to the CVLocalFn function
     CVcfn       is a pointer to the CVCommFn function */

  CVBBD_Data = CVBBDPrecAlloc(CV_cvodemem, *Nloc, *mudq, *mldq, *mu, *ml, 
                              *dqrely, CVgloc, CVcfn);
  if (CVBBD_Data == NULL) { *ier = -1; return; }

  /* Call CVBBDSpgmr to specify the SPGMR linear solver:
     CV_cvodemem is the pointer to the CVODE memory block
     *pretype    is the preconditioner type
     *gstype     is the Gram-Schmidt process type
     *maxl       is the maximum Krylov dimension
     *delt       is the linear convergence tolerance factor */

  *ier = CVBBDSpgmr(CV_cvodemem, *pretype, *maxl, CVBBD_Data);
  if (*ier != 0) return;

  *ier = CVSpgmrSetGSType(CV_cvodemem, *gstype);
  if (*ier != 0) return;

  *ier = CVSpgmrSetDelt(CV_cvodemem, *delt);
  if (*ier != 0) return;

}

/***************************************************************************/

void FCV_REINBBD(integertype *Nloc, integertype *mudq, integertype *mldq, 
                 realtype* dqrely, int *pretype, int *gstype,
                 realtype *delt, int *ier)
{
  int flag;

  /* First call CVReInitBBD to re-initialize CVBBDPRE module:
     CVBBD_Data  is the pointer to P_data
     *mudq,*mldq are the half-bandwidths for computing preconditioner blocks
     *dqrely     is the difference quotient relative increment factor
     CVgloc      is a pointer to the CVLocalFn function
     CVcfn       is a pointer to the CVCommFn function */

  flag = CVBBDPrecReInit(CVBBD_Data, *mudq, *mldq,
                         *dqrely, CVgloc, CVcfn);

  /* Call CVSetSpgmr* to re-initialize the SPGMR linear solver:
     *pretype    is the preconditioner type
     *gstype     is the Gram-Schmidt process type
     *delt       is the linear convergence tolerance factor */
   
  *ier = CVSpgmrResetPrecType(CV_cvodemem, *pretype);
  if (*ier != 0) return;

  *ier = CVSpgmrSetGSType(CV_cvodemem, *gstype);
  if (*ier != 0) return;

  *ier = CVSpgmrSetDelt(CV_cvodemem, *delt);
  if (*ier != 0) return;

}

/***************************************************************************/

/* C function CVgloc to interface between CVBBDPRE module and a Fortran 
   subroutine CVLOCFN. */

void CVgloc(integertype Nloc, realtype t, N_Vector yloc, N_Vector gloc,
            void *f_data)
{
  realtype *yloc_data, *gloc_data;
  
  yloc_data = N_VGetData(yloc);
  gloc_data = N_VGetData(gloc);

  FCV_GLOCFN(&Nloc, &t, yloc_data, gloc_data);

  N_VSetData(gloc_data, gloc);

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

void FCV_BBDOPT(integertype *lenrpw, integertype *lenipw, int *nge)
{
  CVBBDPrecGetIntWorkSpace(CVBBD_Data, lenipw);
  CVBBDPrecGetRealWorkSpace(CVBBD_Data, lenrpw);
  CVBBDPrecGetNumGfnEvals(CVBBD_Data, nge);

}


/***************************************************************************/

/* C function FCVBBDF to interface to CVBBDPrecFree, to free memory 
   created by CVBBDPrecAlloc */

void FCV_BBDF ()
{
  CVBBDPrecFree(CVBBD_Data);
}
