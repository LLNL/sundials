/******************************************************************
 *                                                                *
 * File          : fkinspgmr10.c                                  *
 * Programmers   : Allan G Taylor and Alan C. Hindmarsh @ LLNL    *
 * Version of    : 17 January 2001                                *
 *                                                                *
 *----------------------------------------------------------------*
 * This file has the interface routine fkinspgmr10. See fkinsol.h *
 * for usage.                                                     *
 ******************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "kinspgmr.h" /* prototypes of KINSPGMR interface routines   */
#include "fkinsol.h"  /* prototypes of interfaces, global variables  */

/**************************************************************************/

void F_KINSPGMR10(int *maxl, int *maxlrst, int *msbpre)
{
  /* Call KINSpgmr to specify the SPGMR linear solver:

     This is case 10: preconditioning solve routines but no user ATimes rtn.

     KIN_kmem is the pointer to the KINSOL memory block
     *maxl       is the maximum Krylov dimension
     *maxlrst    is the max number of lin solver restarts
     *msbpre     is the max number of steps calling the preconditioning solver
                  w/o calling the preconditioner setup  routine
     NULL        is a pointer to the preconditioner setup interface routine
     KINPSol     is a pointer to the preconditioner solve interface routine
     NULL        is a pointer to the user ATimes interface routine     
     NULL        is a pointer to the P_data memory structure  */

  KINSpgmr (KIN_kmem, *maxl, *maxlrst, *msbpre, NULL, KINPSol, NULL, NULL);
}
