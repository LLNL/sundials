/******************************************************************
 * File          : fcvode.c                                       *
 * Programmers   : Alan C. Hindmarsh and Radu Serban @ LLNL       *
 * Version of    : 27 March 2002                                  *
 *----------------------------------------------------------------*
 * This is the implementation file for the Fortran interface to   *
 * the CVODE package. See fcvode.h for usage.                     *
 * NOTE: some routines are necessarily stored elsewhere to avoid  *
 * linking problems. See also, therefore, fcvpreco.c, fcvpsol.c   *
 * fcvjtimes.c, and the five fcvspgmr**.c where ** = 01, 10, 11   *
 * 20, 21, for all the options available                          *
 ******************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "llnltyps.h" /* definitions of types real and integer             */
#include "nvector.h"  /* definitions of type N_Vector and vector macros    */
#include "cvode.h"    /* CVODE constants and prototypes                    */
#include "cvdiag.h"   /* prototypes for CVDIAG interface routines          */
#include "cvspgmr.h"  /* prototypes for CVSPGMR interface routines         */
#include "fcmixpar.h" /* global F2C_machEnv variable                       */
#include "fcvode.h"   /* actual function names, prototypes, global vars.   */

/***************************************************************************/

/* Prototypes of the Fortran routines */
void FCV_FUN(integer*, real*, real*, real*);

/**************************************************************************/
void FCV_MALLOC(integer *neq, real *t0, real *y0, 
                integer *meth, integer *itmeth, integer *iatol, 
                real *rtol, real *atol,
                integer *optin, long int *iopt, real *ropt, int *ier)
{
  int lmm, iter, itol;
  N_Vector atolvec;
  void *atolptr;

  CV_yvec = N_VMake(*neq, y0, F2C_machEnv);
  lmm = (*meth == 1) ? ADAMS : BDF;
  iter = (*itmeth == 1) ? FUNCTIONAL : NEWTON;
  if (*iatol == 1)
    { itol = SS; atolptr = atol; }
  else
    { atolvec = N_VMake(*neq, atol, F2C_machEnv);
      itol = SV; atolptr = atolvec; }

  /* Call CVodeMalloc to initialize CVODE: 
     *neq    is the problem size
     CVf     is the user's right-hand side function in y'=f(t,y)
     *t0     is the initial time
     CV_yvec is the initial dependent variable vector
     lmm     is the method specifier
     iter    is the iteration method specifier
     itol    specifies tolerance type
     rtol    is a pointer to the scalar relative tolerance
     atolptr is the absolute tolerance pointer (to scalar or vector)
     NULL    is the pointer to f_data
     NULL    is the pointer to the error message file
     *optin  is the flag indicating optional inputs in iopt and ropt
     iopt    is an array used to communicate optional integer input and output
     ropt    is an array used to communicate optional real input and output
     F2C_machEnv is the pointer to the machine environment block

     A pointer to CVODE problem memory is returned and stored in CV_cvodemem. */

  CV_cvodemem = CVodeMalloc(*neq, CVf, *t0, CV_yvec, lmm, iter, 
                            itol, rtol, atolptr, NULL, NULL, 
                            *optin, iopt, ropt, F2C_machEnv);

  *ier = (CV_cvodemem == NULL) ? -1 : 0 ;
}

/***************************************************************************/

void FCV_REINIT(real *t0, real *y0, integer *meth, integer *itmeth,
                integer *iatol, real *rtol, real *atol, integer *optin,
                long int *iopt, real *ropt, int *ier)
{
  int lmm, iter, itol;
  N_Vector atolvec;
  void *atolptr;
  integer neq;

  neq = ((CVodeMem)CV_cvodemem)->cv_N;

  N_VSetData(y0, CV_yvec);
  lmm = (*meth == 1) ? ADAMS : BDF;
  iter = (*itmeth == 1) ? FUNCTIONAL : NEWTON;
  if (*iatol == 1)
    { itol = SS; atolptr = atol; }
  else
    { atolvec = N_VMake(neq, atol, F2C_machEnv);
      itol = SV; atolptr = atolvec; }

  /* Call CVReInit to re-initialize CVODE: 
     CV_cvodemem is the pointer to the CVODE memory block
     CVf     is the user's right-hand side function in y'=f(t,y)
     *t0     is the initial time
     CV_yvec is the initial dependent variable vector
     lmm     is the method specifier
     iter    is the iteration method specifier
     itol    specifies tolerance type
     rtol    is a pointer to the scalar relative tolerance
     atolptr is the absolute tolerance pointer (to scalar or vector)
     NULL    is the pointer to f_data
     NULL    is the pointer to the error message file
     *optin  is the flag indicating optional inputs in iopt and ropt
     iopt    is an array used to communicate optional integer input and output
     ropt    is an array used to communicate optional real input and output
     F2C_machEnv is the pointer to the machine environment block

     A pointer to CVODE problem memory is returned and stored in CV_cvodemem. */

  *ier = CVReInit(CV_cvodemem, CVf, *t0, CV_yvec, lmm, iter, 
                  itol, rtol, atolptr, NULL, NULL,
                  *optin, iopt, ropt, F2C_machEnv);
}

/***************************************************************************/

void FCV_DIAG(int *ier)
{
  /* Call CVDiag:
     CV_cvodemem is the pointer to the CVODE memory block  */

  *ier = CVDiag(CV_cvodemem);
}

/***************************************************************************/

void FCV_SPGMR00(int *gstype, int *maxl, real *delt, int *ier)
{
  /* Call CVSpgmr to specify the SPGMR linear solver:
     CV_cvodemem is the pointer to the CVODE memory block
     0           is the preconditioner type (none)
     *gstype     is the Gram-Schmidt process type
     *maxl       is the maximum Krylov dimension
     *delt       is the linear convergence tolerance factor
     NULL        is a pointer to the preconditioner setup routine
     NULL        is a pointer to the preconditioner solve routine
     NULL        is the pointer to P_data
     NULL        is a pointer to the Jtimes routine
     NULL        is the pointer to J_data                             */

  *ier = CVSpgmr(CV_cvodemem, 0, *gstype, *maxl, *delt, NULL, NULL,
                 NULL, NULL, NULL);
}

/***************************************************************************/

void FCV_REINSPGMR00(int *gstype, int *maxl, real *delt, int *ier)
{
  /* Call CVReInitSpgmr to specify the SPGMR linear solver:
     CV_cvodemem is the pointer to the CVODE memory block
     0           is the preconditioner type (none)
     *gstype     is the Gram-Schmidt process type
     *maxl       is the maximum Krylov dimension
     *delt       is the linear convergence tolerance factor
     NULL        is a pointer to the preconditioner setup routine
     NULL        is a pointer to the preconditioner solve routine
     NULL        is the pointer to P_data
     NULL        is a pointer to the Jtimes routine
     NULL        is the pointer to J_data                             */

  *ier = CVReInitSpgmr(CV_cvodemem, 0, *gstype, *maxl, *delt, NULL, NULL,
                       NULL, NULL, NULL);
}

/***************************************************************************/

void FCV_CVODE(real *tout, real *t, real *y, int *itask, int *ier)
{
  /* Call CVode:
     CV_cvodemem is the pointer to the CVODE memory block
     *tout       is the t value where output is desired
     CV_yvec     is the N_Vector containing the solution on return
     t           is the returned independent variable value
     *itask      is the task indicator (NORMAL or ONE_STEP) */

  *ier = CVode(CV_cvodemem, *tout, CV_yvec, t, *itask);

  y = N_VGetData(CV_yvec);

}

/***************************************************************************/

void FCV_DKY (real *t, int *k, real *dky, int *ier)
{
  /* Call CVodeDky:
     CV_cvodemem is the pointer to the CVODE memory block
     *t          is the t value where output is desired
     *k          is the derivative order
     CV_yvec     is the N_Vector containing the solution derivative on return */

  *ier = CVodeDky (CV_cvodemem, *t, *k, CV_yvec);

  dky = N_VGetData(CV_yvec);

}

/***************************************************************************/

void FCV_FREE ()
{
  /* Call CVodeFree:
     CV_cvodemem is the pointer to the CVODE memory block */

  CVodeFree (CV_cvodemem);
}

/***************************************************************************/

/* C function CVf to interface between CVODE and a Fortran subroutine CVFUN.
   Addresses of N, t, y, and ydot are passed to CVFUN, using the
   routine N_VGetData from the NVECTOR module.
   Auxiliary data is assumed to be communicated by Common. */

void CVf(integer N, real t, N_Vector y, N_Vector ydot, void *f_data)
{
  real *ydata, *dydata;

  ydata = N_VGetData(y);
  dydata = N_VGetData(ydot);

  FCV_FUN (&N, &t, ydata, dydata);

  N_VSetData(dydata, ydot);

}
