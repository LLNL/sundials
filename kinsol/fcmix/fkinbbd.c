/*
 * ----------------------------------------------------------------
 * $Revision: 1.14 $
 * $Date: 2004-06-18 21:36:28 $
 * ----------------------------------------------------------------
 * Programmer(s): Allan Taylor, Alan Hindmarsh and
 *                Radu Serban @ LLNL
 * ----------------------------------------------------------------
 * This module contains the routines necessary to interface with
 * the KINBBDPRE module and user-supplied Fortran routines. Generic
 * names are used (e.g. FK_COMMFN). The routines here call the
 * generically named routines and provide a standard interface to
 * the C code of the KINBBDPRE package.
 * ----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>
#include "sundialstypes.h"  /* definition of type realtype */
#include "nvector.h"        /* definition of type N_Vector */
#include "kinsol.h"         /* KINSOL constants and prototypes */
#include "fkinsol.h"        /* prototypes of standard interfaces and
                               global variables */
#include "fkinbbd.h"        /* prototypes of interfaces to KINBBDPRE */
#include "kinspgmr.h"       /* prototypes of KINSPGMR interface routines */
#include "kinbbdpre.h"      /* prototypes of KINBBDPRE functions and macros */

/*
 * ----------------------------------------------------------------
 * prototypes of the user-supplied fortran routines
 * ----------------------------------------------------------------
 */

void FK_LOCFN(long int*, realtype*, realtype*);
void FK_COMMFN(long int*, realtype*);

/*
 * ----------------------------------------------------------------
 * Function : FKIN_BBDINIT
 * ----------------------------------------------------------------
 */

void FKIN_BBDINIT(long int *nlocal, long int *mu, long int *ml, int *ier)
{
  KBBD_Data = KINBBDPrecAlloc(KIN_mem, *nlocal, *mu, *ml, 0.0, FKINgloc, FKINgcomm);
  if (KBBD_Data == NULL) { *ier = -1; return; }
}

/*
 * ----------------------------------------------------------------
 * Function : FKIN_BBDSPGMR
 * ----------------------------------------------------------------
 */

void FKIN_BBDSPGMR(int *maxl, int *maxlrst, int *ier)
{
  *ier = KINBBDSpgmr(KIN_mem, *maxl, KBBD_Data);
  if (*ier != 0) return;

  *ier = KINSpgmrSetMaxRestarts(KIN_mem, *maxlrst);
  if (*ier != 0) return;
}

/*
 * ----------------------------------------------------------------
 * Function : FKINgloc
 * ----------------------------------------------------------------
 * C function FKINgloc is the interface between the KINBBDPRE
 * module and the Fortran subroutine FK_LOCFN.
 * ----------------------------------------------------------------
 */

void FKINgloc(long int Nloc, N_Vector uu, N_Vector gval, void *f_data)
{
  realtype *uloc, *gloc;

  uloc = (realtype *) N_VGetData(uu);
  gloc = (realtype *) N_VGetData(gval);

  FK_LOCFN(&Nloc, uloc, gloc);

  N_VSetData((void *)gloc, gval);
}

/*
 * ----------------------------------------------------------------
 * Function : FKINgcomm
 * ----------------------------------------------------------------
 * C function FKINgcomm is the interface between the KINBBDPRE
 * module and the Fortran subroutine FK_COMMFN.
 * ----------------------------------------------------------------
 */

void FKINgcomm(long int Nloc, N_Vector uu, void *f_data)
{
  realtype *uloc;

  uloc = (realtype *) N_VGetData(uu);
  
  FK_COMMFN(&Nloc, uloc);
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
  KINBBDPrecGetIntWorkSpace(KBBD_Data, lenipw);
  KINBBDPrecGetRealWorkSpace(KBBD_Data, lenrpw);
  KINBBDPrecGetNumGfnEvals(KBBD_Data, nge);
}

/*
 * ----------------------------------------------------------------
 * Function : FKIN_BBDFREE
 * ----------------------------------------------------------------
 * C function FKIN_BBDFREE is used to interface with KINBBDPrecFree
 * in order to free memory created by KINBBDPrecAlloc.
 * ----------------------------------------------------------------
 */

void FKIN_BBDFREE()
{
  KINBBDPrecFree(KBBD_Data);
}
