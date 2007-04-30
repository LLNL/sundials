/*
 * -----------------------------------------------------------------
 * $Revision: 1.4 $
 * $Date: 2007-04-30 19:29:01 $
 * -----------------------------------------------------------------
 * Programmer(s): Allan Taylor, Alan Hindmarsh, Radu Serban, and
 *                Aaron Collier @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * -----------------------------------------------------------------
 * This module contains the routines necessary to interface with
 * the KINBBDPRE module and user-supplied Fortran routines. Generic
 * names are used (e.g. FK_COMMFN). The routines here call the
 * generically named routines and provide a standard interface to
 * the C code of the KINBBDPRE package.
 * ----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include "fkinsol.h"               /* standard interfaces and global variables     */
#include "fkinbbd.h"               /* prototypes of interfaces to KINBBDPRE        */

#include <kinsol/kinsol_bbdpre.h>  /* prototypes of KINBBDPRE functions and macros */
#include <kinsol/kinsol_sptfqmr.h> /* prototypes of KINSPTFQMR interface routines  */
#include <kinsol/kinsol_spbcgs.h>  /* prototypes of KINSPBCG interface routines    */
#include <kinsol/kinsol_spgmr.h>   /* prototypes of KINSPGMR interface routines    */

/*
 * ----------------------------------------------------------------
 * private constants
 * ----------------------------------------------------------------
 */

#define ZERO RCONST(0.0)

/*
 * ----------------------------------------------------------------
 * prototypes of the user-supplied fortran routines
 * ----------------------------------------------------------------
 */

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

extern void FK_LOCFN(int*, realtype*, realtype*, int*);
extern void FK_COMMFN(int*, realtype*, int*);

#ifdef __cplusplus
}
#endif

/*
 * ----------------------------------------------------------------
 * Function : FKIN_BBDINIT
 * ----------------------------------------------------------------
 */

void FKIN_BBDINIT(int *nlocal, int *mudq, int *mldq,
		  int *mu, int *ml, int *ier)
{
  *ier = KINBBDPrecInit(KIN_kinmem, *nlocal, *mudq, *mldq,
                        *mu, *ml, ZERO, FKINgloc, FKINgcomm);

  return;
}

/*
 * ----------------------------------------------------------------
 * Function : FKINgloc
 * ----------------------------------------------------------------
 * C function FKINgloc is the interface between the KINBBDPRE
 * module and the Fortran subroutine FK_LOCFN.
 * ----------------------------------------------------------------
 */

int FKINgloc(int Nloc, N_Vector uu, N_Vector gval, void *user_data)
{
  realtype *uloc, *gloc;
  int ier;

  uloc = gloc = NULL;

  uloc = N_VGetArrayPointer(uu);
  gloc = N_VGetArrayPointer(gval);

  FK_LOCFN(&Nloc, uloc, gloc, &ier);

  N_VSetArrayPointer(gloc, gval);

  return(0);
}

/*
 * ----------------------------------------------------------------
 * Function : FKINgcomm
 * ----------------------------------------------------------------
 * C function FKINgcomm is the interface between the KINBBDPRE
 * module and the Fortran subroutine FK_COMMFN.
 * ----------------------------------------------------------------
 */

int FKINgcomm(int Nloc, N_Vector uu, void *user_data)
{
  realtype *uloc;
  int ier;

  uloc = NULL;

  uloc = N_VGetArrayPointer(uu);
  
  FK_COMMFN(&Nloc, uloc, &ier);

  return(0);
}

/*
 * ----------------------------------------------------------------
 * Function : FKIN_BBDOPT
 * ----------------------------------------------------------------
 * C function FKIN_BBDOPT is used to access optional outputs 
 * realated to the BBD preconditioner.
 * ----------------------------------------------------------------
 */

void FKIN_BBDOPT(long int *lenrpw, long int *lenipw, long int *nge)
{
  KINBBDPrecGetWorkSpace(KIN_kinmem, lenrpw, lenipw);
  KINBBDPrecGetNumGfnEvals(KIN_kinmem, nge);

  return;
}

