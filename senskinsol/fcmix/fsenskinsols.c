/*******************************************************************
 *                                                                 *
 * File          : fsenskinsols.c                                  *
 * Programmers   : Allan G Taylor, Alan C. Hindmarsh, and          *
 *                 Keith Grant @ LLNL                              *
 * Version of    : 27 Oct 2000  Modified for Sens_KINSol           *
 *               : 25 Mar 1999  Removed machEnv from KINSol args   *
 *-----------------------------------------------------------------*
 * This is the implementation file for the Fortran interface to    *
 * the SensKINSOL package. See fsenskinsol.h for usage.            *
 * NOTE: some routines are necessarily stored elsewhere to avoid   *
 * linking problems. See also, therefore, fkinpreco.c, fkinpsol.c  *
 * fkinuatimes.c, and the five fkinspgmr**.c where ** = 01, 10, 11 *
 * 20, 21, for all the options available.                          *
 *******************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "llnltyps.h"      /* definitions of types real and integer      */
#include "nvector.h"       /* definitions and macros for type N_Vector   */
#include "kinsol.h"        /* KINSOL constants and prototypes            */
#include "sens_kinsol.h"   /* SensKINSOL constants and prototypes        */
#include "sens_kinspgmr.h" /* SensKINSol linear solver definitions       */
#include "fkinsol.h"       /* FORTRAN prototypes etc. for KINSol         */
#include "fsenskinsol.h"   /* FORTRAN prototypes etc. for SensKINSolbles */

/***************************************************************************/

void F_SSENSKINMALLOC(integer *neq, integer *np, integer *diffrhs,
         real *diffscl, real *p, real *pbar, integer *ier)
{

  /***********************************************************************
   * Call SensKINMalloc to initialize memory for SensKINSOL:
   *
   *  *neq        is the problem size
   *
   *  *np         is the number of parameters for sensitivity analysis
   *
   *  *diffrhs    is the difference option for the sensitivity RHS
   *
   *  *diffscl    is the sensitivity difference increment scaling factor
   *
   *  p           is the array of sensitivity parameters
   *
   *  pbar        is the array of nominal magnitudes of the parameters
   *
   *  NULL        is the pointer to the error message file
   *
   *  NULL        is the NULL pointer to the (parallel only) machine
   *                 environment block
   *
   * A pointer to Sens_KINSOL problem memory is returned and stored in
   * the global variable SENSKIN_smem.
   *
   ***********************************************************************/

  SENSKIN_smem = SensKINMalloc(*neq, *np, *diffrhs, *diffscl, p, pbar,
         NULL, NULL);

  *ier = (SENSKIN_smem == NULL) ? -1 : 0 ;
}


/***************************************************************************/

void F_SENSKINSPGMR00( integer *maxl, integer *maxlrst, integer *msbpre)
{
  /***********************************************************************
   * Call SensKINSpgmr to specify the SPGMR linear solver:
   *
   * This is the base case: no preconditioning routines and no
   * user ATimes routine.
   *
   *
   *  SENSKIN_smem   is the pointer to the Sens_KINSOL memory block
   *
   *  *maxl          is the maximum Krylov dimension
   *
   *  *maxlrst       is the max number of lin solver restarts
   *
   *  *msbpre        is the maximum number of steps calling the
   *                    preconditioning solver w/o calling the
   *                    preconditioner setup routine
   *
   *  NULL           is a pointer to the preconditioner setup interface
   *                    routine
   *
   *  NULL           is a pointer to the preconditioner solve interface
   *                    routine
   *
   *  NULL           is a pointer to the user ATimes interface routine
   *
   *  NULL           is a pointer to the P_data memory structure
   *
   ***********************************************************************/

   SensKINSpgmr (SENSKIN_smem, *maxl, *maxlrst, *msbpre, NULL, NULL,
      NULL, NULL);
}

/***************************************************************************/

void F_SENSKINSOL(integer *Neq, real *uu, integer *globalstrategy,
	       real *uscale , real *fscale, real *fnormtol, real *scsteptol,
          real *constraints, integer *inopt, long int *iopt, real *ropt,
          integer *ier)

{
  N_Vector uuvec, uscalevec, fscalevec, constraintsvec;
  boole binopt;
  SensKINParams sens;
  KINMem kin_mem;

  sens = (SensKINParams) SENSKIN_smem;
  kin_mem = (KINMem) sens->kin_mem;

  binopt = *inopt;

  /* Since uuvec, uscalevec, and fscalevec aren't deallocated       */
  /* locally, deallocate existing copies before allocating new      */
  /* ones. If the kin_mem pointers are null, these are just no-ops. */

  N_VDISPOSE(kin_mem->kin_uu);
  N_VDISPOSE(kin_mem->kin_uscale);
  N_VDISPOSE(kin_mem->kin_fscale);

  N_VMAKE(uuvec, uu, *Neq);
  N_VMAKE(uscalevec, uscale, *Neq);
  N_VMAKE(fscalevec, fscale, *Neq);
  N_VMAKE(constraintsvec, constraints, *Neq);

  /***********************************************************************
   * Call SensKINSol:
   *
   *  SENSKIN_smem   is the pointer to the Sens_KINSol memory block
   *  Neq            is the total system size
   *
   *  uuvec          is the N_Vector containing the solution found
   *                    by KINSol
   *
   *  SensKINfunc    is the standard interface function to the user-
   *                    supplied KFUN
   *
   *  globalstragegy is an integer defining the global strategy choice
   *
   *  uscalevec      is the N_Vector containing the u scaling
   *
   *  fscalevec      is the N_Vector containing the f scaling
   *
   *  fnormtol       is the real value of fnormtol (stopping tolerance)
   *
   *  scsteptol      is the real value of scsteptol (scaled step tolerance)
   *
   *  constraintsvec is the N_Vector containing the constraints
   *
   *  inopt          is the flag indicating whether or not the ropt/iopt
   *                    arrays are to be used as input parameters
   *
   *  iopt           is the array of integer user input/output
   *
   *  ropt           is the array of real  user input/output
   *
   *  NULL           is the (unused) pointer to f_data
   *
   ***********************************************************************/

  *ier = SensKINSol (SENSKIN_smem, *Neq, uuvec, SensKINfunc,
            *globalstrategy, uscalevec, fscalevec, *fnormtol, *scsteptol,
            constraintsvec, binopt, iopt, ropt, NULL);

   N_VDISPOSE(constraintsvec);  /* Not used in sensitivity analysis */

}

/***************************************************************************/

void F_SENSKINFREE()
{
  SensKINParams sens;
  KINMem kin_mem;

  sens = (SensKINParams) SENSKIN_smem;
  kin_mem = (KINMem) sens->kin_mem;


  /* Deallocate N_Vectors allocated by F_SENSKINSOL and not freed locally */

  N_VDISPOSE(kin_mem->kin_uu);
  N_VDISPOSE(kin_mem->kin_uscale);
  N_VDISPOSE(kin_mem->kin_fscale);

  /* Call SensKINFree:
     SENSKIN_smem is the global pointer to the Sens_KINSOL memory block */

  SensKINFree (SENSKIN_smem);
}

/***************************************************************************/

void F_SENSKINLININIT (integer *ier)
{
   /* Initialize the sensitivity solution state */

   *ier = SensKINLinInit (SENSKIN_smem);
}

/***************************************************************************/

void F_SENSKINDIFF (integer *diffrhs, real *diffscl, integer *ier)
{
  /* Call SensKINDiff to reset sensitivity right-hand-side differencing */
  /* option and/or the difference increment scaling factor              */

  *ier = SensKINDiff (SENSKIN_smem, *diffrhs, *diffscl);
}

/***************************************************************************/

void F_SENSKINLINSOLVE (integer *ip, real *ww, integer *ier)
{
  SensKINParams sens;
  KINMem kin_mem;
  N_Vector wwvec;
  integer Neq;

  sens = (SensKINParams) SENSKIN_smem;
  kin_mem = (KINMem) sens->kin_mem;

  /* Get the size of the system that was specified to F_SENSKINSOL */
  Neq = kin_mem->kin_Neq;
  N_VMAKE(wwvec, ww, Neq);


  /* Call SensKINLINSolve to solve the sensitivity problem for    */
  /* the ip'th parameter. The value of ip in the call is adjusted */
  /* assuming that ip is based on a one-origin loop in FORTRAN,   */
  /* while the C routines take it to be of zero-origin.           */

  *ier = SensKINLinSolve (SENSKIN_smem, (*ip-1), wwvec);
  N_VDISPOSE(wwvec);
}

/***************************************************************************/

void F_SENSKININIT (integer *Neq, real *uu, integer *globalstrategy,
	       real *uscale , real *fscale, real *fnormtol, real *scsteptol,
          real *constraints, integer *inopt, long int *iopt, real *ropt,
          integer *ier)
{
  N_Vector uuvec, uscalevec, fscalevec, constraintsvec;
  boole binopt;
  SensKINParams sens;
  KINMem kin_mem;

  sens = (SensKINParams) SENSKIN_smem;
  kin_mem = (KINMem) sens->kin_mem;

  binopt = *inopt;

  /* Since uuvec, uscalevec, and fscalevec aren't deallocated       */
  /* locally, deallocate existing copies before allocating new      */
  /* ones. If the kin_mem pointers are null, these are just no-ops. */

  N_VDISPOSE(kin_mem->kin_uu);
  N_VDISPOSE(kin_mem->kin_uscale);
  N_VDISPOSE(kin_mem->kin_fscale);

  N_VMAKE(uuvec, uu, *Neq);
  N_VMAKE(uscalevec, uscale, *Neq);
  N_VMAKE(fscalevec, fscale, *Neq);
  N_VMAKE(constraintsvec, constraints, *Neq);

  /* Call SensKINInit to initialize KINSol memory without calling KINSOL
     when the non-linear solution is already known. See F_SENSKINSOL for
     a description of the arguments */


  *ier = SensKINInit (SENSKIN_smem, *Neq, uuvec, KINfunc, *globalstrategy,
		 uscalevec, fscalevec, *fnormtol, *scsteptol,
		 constraintsvec, binopt, iopt, ropt, NULL);

  N_VDISPOSE(constraintsvec);  /* Not used in sensitivity analysis */

}

/***************************************************************************/


/**************************************************************************
* C function SensKINfunc acts as an interface between Sens_KINSol and
* the FORTRAN user-supplied subroutine SENSKFUN. Addresses of Nlocal,
* uu and fval are passed to SENSKFUN, using the macros N_VLENGTH
* and N_VDATA from the NVECTOR module. The array of the sensitivity
* parameters p[] is retrieved from the sensitivity memory structure
* and passed to SENSKFUN. Auxiliary data is assumed to be communicated
* by common.
**************************************************************************/

void SensKINfunc(integer Neq, N_Vector uu, N_Vector fval, void *f_data)
{
  real *udata, *fdata, *p;
  integer Nlocal;
  SensKINParams sens;
  integer np;

  sens = (SensKINParams) SENSKIN_smem;
  np = sens->np;
  p  = sens->p;

  Nlocal = N_VLENGTH(uu);
  udata = N_VDATA(uu);
  fdata = N_VDATA(fval);

  /* NOTE: The C user data f_data is not passed to KFUN... it is NULL */

  K_FUN (&Nlocal, udata, fdata, &np, p);

}
