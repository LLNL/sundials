/*
 * -----------------------------------------------------------------
 * $Revision: 1.1 $
 * $Date: 2005-05-11 23:10:54 $
 * ----------------------------------------------------------------- 
 * Programmer(s): Aaron Collier @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2005, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see sundials/cvode/LICENSE.
 * -----------------------------------------------------------------
 * This module contains the routines necessary to interface with the
 * IDABBDPRE module and user-supplied Fortran routines.
 * The routines here call the generically named routines and provide
 * a standard interface to the C code of the IDABBDPRE package.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include "idabbdpre.h"      /* prototypes of IDABBDPRE functions and macros */
#include "ida.h"            /* IDA constants and prototypes                 */
#include "idaspgmr.h"       /* prototypes of IDASPGMR interface routines    */
#include "idaspbcg.h"       /* prototypes of IDASPBCG interface routines    */
#include "fidabbd.h"        /* prototypes of interfaces to IDABBD           */
#include "fida.h"           /* actual function names, prototypes and global
			       variables                                    */
#include "nvector.h"        /* definition of type N_Vector                  */
#include "sundialstypes.h"  /* definition of type realtype                  */

/*************************************************/

/* private constant(s) */

#define ZERO RCONST(0.0)

/*************************************************/

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

extern void FIDA_GLOCFN(long int*, realtype*, realtype*, realtype*, realtype*, int*);
extern void FIDA_COMMFN(long int*, realtype*, realtype*, realtype*, int*);

#ifdef __cplusplus
}
#endif

/*************************************************/

void FIDA_BBDINIT(long int *Nloc, long int *mudq, long int *mldq,
		  long int *mu, long int *ml, realtype *dqrely, int *ier)
{
  IDABBD_Data = IDABBDPrecAlloc(IDA_idamem, *Nloc, *mudq, *mldq, *mu, *ml,
				*dqrely, (IDABBDLocalFn) FIDAgloc, (IDABBDCommFn) FIDAcfn);
  if (IDABBD_Data == NULL) *ier = -1;
  else *ier = 0;

  return;
}

/*************************************************/

void FIDA_BBDSPBCG(int *maxl, realtype *eplifac, realtype *dqincfac, int *ier)
{
  *ier = 0;

  *ier = IDABBDSpbcg(IDA_idamem, *maxl, IDABBD_Data);
  if (*ier != IDASPBCG_SUCCESS) return;

  *ier = IDASpbcgSetEpsLin(IDA_idamem, *eplifac);
  if (*ier != IDASPBCG_SUCCESS) return;

  if (*dqincfac != ZERO) {
    *ier = IDASpbcgSetIncrementFactor(IDA_idamem, *dqincfac);
    if (*ier != IDASPBCG_SUCCESS) return;
  }

  IDA_ls = IDA_SPBCG;

  return;
}

/*************************************************/

void FIDA_BBDSPGMR(int *maxl, int *gstype, int *maxrs,
		   realtype *eplifac, realtype *dqincfac, int *ier)
{
  *ier = 0;

  *ier = IDABBDSpgmr(IDA_idamem, *maxl, IDABBD_Data);
  if (*ier != IDASPGMR_SUCCESS) return;

  if (*gstype != 0) {
    *ier = IDASpgmrSetGSType(IDA_idamem, *gstype);
    if (*ier != IDASPGMR_SUCCESS) return;
  }

  if (*maxrs != 0) {
    *ier = IDASpgmrSetMaxRestarts(IDA_idamem, *maxrs);
    if (*ier != IDASPGMR_SUCCESS) return;
  }

  *ier = IDASpgmrSetEpsLin(IDA_idamem, *eplifac);
  if (*ier != IDASPGMR_SUCCESS) return;

  if (*dqincfac != ZERO) {
    *ier = IDASpgmrSetIncrementFactor(IDA_idamem, *dqincfac);
    if (*ier != IDASPGMR_SUCCESS) return;
  }

  IDA_ls = IDA_SPGMR;

  return;
}

/*************************************************/

void FIDA_BBDREINIT(long int *Nloc, long int *mudq, long int *mldq,
		    realtype *dqrely, int *ier)
{
  *ier = 0;

  *ier = IDABBDPrecReInit(IDABBD_Data, *mudq, *mldq,
			  *dqrely, (IDABBDLocalFn) FIDAgloc, (IDABBDCommFn) FIDAcfn);

  return;
}

/*************************************************/

int FIDAgloc(long int Nloc, realtype t, N_Vector yy, N_Vector yp,
	     N_Vector gval, void *res_data)
{
  realtype *yy_data, *yp_data, *gval_data;
  int ier;

  /* Initialize all pointers to NULL */
  yy_data = yp_data = gval_data = NULL;

  /* NOTE: The user-supplied routine should set ier to an
     appropriate value, but we preset the value to zero
     (meaning SUCCESS) so the user need only reset the
     value if an error occurred */
  ier = 0;

  /* Get pointers to vector data */
  yy_data   = N_VGetArrayPointer(yy);
  yp_data   = N_VGetArrayPointer(yp);
  gval_data = N_VGetArrayPointer(gval);

  /* Call user-supplied routine */
  FIDA_GLOCFN(&Nloc, &t, yy_data, yp_data, gval_data, &ier);

  return(ier);
}

/*************************************************/

int FIDAcfn(long int Nloc, realtype t, N_Vector yy, N_Vector yp,
	    void *res_data)
{
  realtype *yy_data, *yp_data;
  int ier;

  /* Initialize all pointers to NULL */
  yy_data = yp_data = NULL;

  /* NOTE: The user-supplied routine should set ier to an
     appropriate value, but we preset the value to zero
     (meaning SUCCESS) so the user need only reset the
     value if an error occurred */
  ier = 0;

  /* Get pointers to vector data */
  yy_data = N_VGetArrayPointer(yy);
  yp_data = N_VGetArrayPointer(yp);

  /* Call user-supplied routine */
  FIDA_COMMFN(&Nloc, &t, yy_data, yp_data, &ier);

  return(ier);
}

/*************************************************/

void FIDA_BBDOPT(long int *lenrpw, long int *lenipw, long int *nge)
{
  IDABBDPrecGetWorkSpace(IDABBD_Data, lenrpw, lenipw);
  IDABBDPrecGetNumGfnEvals(IDABBD_Data, nge);

  return;
}

/*************************************************/

void FIDA_BBDFREE(void)
{
  IDABBDPrecFree(IDABBD_Data);

  return;
}
