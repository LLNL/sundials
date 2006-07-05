/*
 * -----------------------------------------------------------------
 * $Revision: 1.1 $
 * $Date: 2006-07-05 15:32:37 $
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

extern void FK_LOCFN(long int*, realtype*, realtype*, int*);
extern void FK_COMMFN(long int*, realtype*, int*);

#ifdef __cplusplus
}
#endif

/*
 * ----------------------------------------------------------------
 * Function : FKIN_BBDINIT
 * ----------------------------------------------------------------
 */

void FKIN_BBDINIT(long int *nlocal, long int *mudq, long int *mldq,
		  long int *mu, long int *ml, int *ier)
{
  KBBD_Data = NULL;

  KBBD_Data = KINBBDPrecAlloc(KIN_kinmem, *nlocal, *mudq, *mldq,
			      *mu, *ml, ZERO, FKINgloc, FKINgcomm);
  if (KBBD_Data == NULL) *ier = -1;
  else *ier = 0;

  return;
}

/*
 * ----------------------------------------------------------------
 * Function : FKIN_BBDSPTFQMR
 * ----------------------------------------------------------------
 */

void FKIN_BBDSPTFQMR(int *maxl, int *ier)
{
  *ier = 0;

  *ier = KINBBDSptfqmr(KIN_kinmem, *maxl, KBBD_Data);

  return;
}

/*
 * ----------------------------------------------------------------
 * Function : FKIN_BBDSPBCG
 * ----------------------------------------------------------------
 */

void FKIN_BBDSPBCG(int *maxl, int *ier)
{
  *ier = 0;

  *ier = KINBBDSpbcg(KIN_kinmem, *maxl, KBBD_Data);

  return;
}

/*
 * ----------------------------------------------------------------
 * Function : FKIN_BBDSPGMR
 * ----------------------------------------------------------------
 */

void FKIN_BBDSPGMR(int *maxl, int *maxlrst, int *ier)
{
  *ier = 0;

  *ier = KINBBDSpgmr(KIN_kinmem, *maxl, KBBD_Data);
  if (*ier != KINSPILS_SUCCESS) return;

  *ier = KINSpilsSetMaxRestarts(KIN_kinmem, *maxlrst);

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

int FKINgloc(long int Nloc, N_Vector uu, N_Vector gval, void *f_data)
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

int FKINgcomm(long int Nloc, N_Vector uu, void *f_data)
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
 * C function FKIN_BBDOPT is used to access optional outputs stored
 * in the KBBD_Data data structure.
 * ----------------------------------------------------------------
 */

void FKIN_BBDOPT(long int *lenrpw, long int *lenipw, long int *nge)
{
  KINBBDPrecGetWorkSpace(KBBD_Data, lenrpw, lenipw);
  KINBBDPrecGetNumGfnEvals(KBBD_Data, nge);

  return;
}

/*
 * ----------------------------------------------------------------
 * Function : FKIN_BBDFREE
 * ----------------------------------------------------------------
 * C function FKIN_BBDFREE is used to interface with KINBBDPrecFree
 * in order to free memory created by KINBBDPrecAlloc.
 * ----------------------------------------------------------------
 */

void FKIN_BBDFREE(void)
{
  KINBBDPrecFree(&KBBD_Data);

  return;
}
