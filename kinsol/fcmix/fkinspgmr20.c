/******************************************************************
 *                                                                *
 * File          : fkinspgmr20.c                                  *
 * Programmers   : Allan G Taylor and Alan C. Hindmarsh @ LLNL    *
 * Version of    : 30 July 2002                                   *
 *----------------------------------------------------------------*
 * Routine used to interface between a Fortran main and the       *
 * various options available re preconditioning and user-supplied *
 * A-times routines (Jacobian A times vector v)                   *
 ******************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "kinspgmr.h" /* prototypes of KINSPGMR interface routines   */
#include "fkinsol.h"  /* prototypes of interfaces, global variables  */

/***************************************************************************/

void F_KINSPGMR20(int *maxl, int *maxlrst, int *msbpre, int *ier)
{
  /* Call KINSpgmr to specify the SPGMR linear solver:

     This is case 20: preconditioning setup and solve routines but 
          no user ATimes rtn.

     KIN_kmem is the pointer to the KINSOL memory block
     *maxl       is the maximum Krylov dimension
     *maxlrst    is the max number of lin solver restarts
     *msbpre     is the max number of steps calling the preconditioning solver
                  w/o calling the preconditioner setup routine
     KINPreco    is a pointer to the preconditioner setup interface routine
     KINPSol    is a pointer to the preconditioner solve interface routine
     NULL        is a pointer to the user ATimes interface routine     
     NULL        is a pointer to the P_data memory structure  */

  *ier = KINSpgmr (KIN_kmem, *maxl, *maxlrst, *msbpre, KINPreco, KINPSol,
                   NULL, NULL);
}
