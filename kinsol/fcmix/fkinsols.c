/******************************************************************
 *                                                                *
 * File          : fkinsols.c                                     *
 * Programmers   : Allan G Taylor and Alan C. Hindmarsh @ LLNL    *
 * Version of    : 17 January 2001                                *
 *----------------------------------------------------------------*
 * This is the implementation file for the Fortran interface to   *
 * the KINSOL package. See fkinsol.h for usage.                   *
 * NOTE: some routines are necessarily stored elsewhere to avoid  *
 * linking problems. See also, therefore, fkinpreco.c, fkinpsol.c *
 * fkinuatimes.c, and the five fkinspgmr**.c where ** = 01, 10, 11*
 * 20, 21, for all the options available.                          *
 ******************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "llnltyps.h" /* definitions of types real and integer            */
#include "nvector.h"  /* definitions of type N_Vector and related macros  */
#include "kinsol.h"   /* KINSOL constants and prototypes                  */
#include "kinspgmr.h" /* prototypes of KINSPGMR interface routines        */
#include "fkinsol.h"  /* prototypes of interfaces, global variables       */

/**************************************************************************/


void F_SKINMALLOC(integer *neq, integer *ier)
{

  /* Call KINMalloc to initialize memory for KINSOL: 
     *neq    is the problem size
     NULL    is the pointer to the error message file
     KIN_machEnv is the pointer to the machine environment block

     A pointer to KINSOL problem memory is returned and stored in KIN_kmem. */

  KIN_kmem = KINMalloc(*neq,  NULL, NULL);

  *ier = (KIN_kmem == NULL) ? -1 : 0 ;
}


/***************************************************************************/

void F_KINSPGMR00( int *maxl, int *maxlrst, int *msbpre)
{
  /* Call KINSpgmr to specify the SPGMR linear solver:

     This is the base case: no preconditioning routines and no user ATimes rtn.

     KIN_kmem is the pointer to the KINSOL memory block
     *maxl       is the maximum Krylov dimension
     *maxlrst    is the max number of lin solver restarts
     *msbpre     is the max number of steps calling the preconditioning solver
                  w/o calling the preconditioner setup routine
     NULL        is a pointer to the preconditioner setup interface routine
     NULL        is a pointer to the preconditioner solve interface routine
     NULL        is a pointer to the user ATimes interface routine     
     NULL        is a pointer to the P_data memory structure  */

  KINSpgmr (KIN_kmem, *maxl, *maxlrst, *msbpre, NULL, NULL, NULL, NULL);
}

/***************************************************************************/

void F_KINSOL(int *Neq, real *uu, int *globalstrategy, 
	       real *uscale , real *fscale, real *fnormtol, real *scsteptol,
	       real *constraints, int *inopt, long int *iopt, real *ropt,
	       int *ier)

{ N_Vector uuvec, uscalevec, fscalevec, constraintsvec; integer i;
  boole binopt;

  binopt = *inopt;
  
  N_VMAKE(uuvec, uu, *Neq);
  N_VMAKE(uscalevec, uscale, *Neq);
  N_VMAKE(fscalevec, fscale, *Neq);
  N_VMAKE(constraintsvec, constraints, *Neq);

  /* Call KINSol:

     KIN_kmem is the pointer to the CVODE memory block
     Neq is the total system size
     uuvec is the N_Vector containing the solution found by KINSol
     KINfunc   is the standard interface function to the user-supplied KFUN
     globalstragegy is an integer defining the global strategy choice
     uscalevec is the N_Vector containing the u scaling 
     fscalevec is the N_Vector containing the f scaling
     fnormtol is the real value of fnormtol (stopping tolerance)
     scsteptol is the real value of scsteptol (scaled step tolerance)
     constraintsvec is the N_Vector containing the constraints
     inopt is the flag indicating whether or not the ropt/iopt arrays are
             to be used as input parameters
     iopt  is the array of integer user input/output
     ropt  is the array of real  user input/output
     NULL is the pointer to f_data
     
     
      */

  *ier = KINSol (KIN_kmem, *Neq, uuvec, KINfunc, *globalstrategy, 
		 uscalevec, fscalevec, *fnormtol, *scsteptol,
		 constraintsvec, binopt, iopt, ropt, NULL);

  N_VDISPOSE(uuvec);
  N_VDISPOSE(uscalevec);
  N_VDISPOSE(fscalevec);
  N_VDISPOSE(constraintsvec);

}

/***************************************************************************/

void F_KINFREE()
{
  /* Call KINFree:
     KIN_kmem is the pointer to the CVODE memory block */

  KINFree (KIN_kmem);
}

/***************************************************************************/


/* C function KINfunc acts as an interface between KINSol and the Fortran 
              user-supplied subroutine KFUN
   Addresses of Nlocal, uu and fval are passed to KFUN, using the
   macros N_VLENGTH and N_VDATA from the NVECTOR module.
   Auxiliary data is assumed to be communicated by Common. */

void KINfunc(integer Neq, N_Vector uu, N_Vector fval, void *f_data)
{
  real *udata, *fdata;
  integer Nlocal;

  /* NOTE: f_data is not passed to KFUN... it is NULL */

  Nlocal = N_VLENGTH(uu);
  udata = N_VDATA(uu);
  fdata = N_VDATA(fval);

  K_FUN (&Nlocal, udata, fdata);

}
