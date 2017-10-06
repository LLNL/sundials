/*
 * -----------------------------------------------------------------
 * $Revision: 4075 $
 * $Date: 2014-04-24 10:46:58 -0700 (Thu, 24 Apr 2014) $
 * ----------------------------------------------------------------- 
 * Programmer(s): Aaron Collier @ LLNL
 * -----------------------------------------------------------------
 * LLNS Copyright Start
 * Copyright (c) 2014, Lawrence Livermore National Security
 * This work was performed under the auspices of the U.S. Department 
 * of Energy by Lawrence Livermore National Laboratory in part under 
 * Contract W-7405-Eng-48 and in part under Contract DE-AC52-07NA27344.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * LLNS Copyright End
 * -----------------------------------------------------------------
 * This module contains the routines necessary to interface with the
 * IDABBDPRE module and user-supplied Fortran routines.
 * The routines here call the generically named routines and provide
 * a standard interface to the C code of the IDABBDPRE package.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include "fida.h"            /* function names, prototypes, global variables */
#include "fidabbd.h"         /* prototypes of interfaces to IDABBD           */

#include <ida/ida_bbdpre.h>  /* prototypes of IDABBDPRE functions and macros */
#include <ida/ida_spgmr.h>   /* prototypes of IDASPGMR interface routines    */
#include <ida/ida_spbcgs.h>  /* prototypes of IDASPBCG interface routines    */
#include <ida/ida_sptfqmr.h> /* prototypes of IDASPTFQMR interface routines  */

/*************************************************/

/* private constant(s) */

#define ZERO RCONST(0.0)

/*************************************************/

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

  extern void FIDA_GLOCFN(long int*, 
                          realtype*, realtype*, realtype*, realtype*, 
                          long int*, realtype*,
                          int*);
  extern void FIDA_COMMFN(long int*, 
                          realtype*, realtype*, realtype*, 
                          long int*, realtype*,
                          int*);

#ifdef __cplusplus
}
#endif

/*************************************************/

void FIDA_BBDINIT(long int *Nloc, long int *mudq, long int *mldq,
		  long int *mu, long int *ml, realtype *dqrely, int *ier)
{
  *ier = IDABBDPrecInit(IDA_idamem, *Nloc, *mudq, *mldq, *mu, *ml,
                        *dqrely, (IDABBDLocalFn) FIDAgloc, (IDABBDCommFn) FIDAcfn);

  return;
}

/*************************************************/

void FIDA_BBDREINIT(long int *Nloc, long int *mudq, long int *mldq,
		    realtype *dqrely, int *ier)
{
  *ier = 0;

  *ier = IDABBDPrecReInit(IDA_idamem, *mudq, *mldq, *dqrely);

  return;
}

/*************************************************/

int FIDAgloc(long int Nloc, realtype t, N_Vector yy, N_Vector yp,
	     N_Vector gval, void *user_data)
{
  realtype *yy_data, *yp_data, *gval_data;
  int ier;
  FIDAUserData IDA_userdata;

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

  IDA_userdata = (FIDAUserData) user_data;

  /* Call user-supplied routine */
  FIDA_GLOCFN(&Nloc, &t, yy_data, yp_data, gval_data, 
              IDA_userdata->ipar, IDA_userdata->rpar, &ier);

  return(ier);
}

/*************************************************/

int FIDAcfn(long int Nloc, realtype t, N_Vector yy, N_Vector yp,
	    void *user_data)
{
  realtype *yy_data, *yp_data;
  int ier;
  FIDAUserData IDA_userdata;

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

  IDA_userdata = (FIDAUserData) user_data;

  /* Call user-supplied routine */
  FIDA_COMMFN(&Nloc, &t, yy_data, yp_data, 
              IDA_userdata->ipar, IDA_userdata->rpar, &ier);

  return(ier);
}

/*************************************************/

void FIDA_BBDOPT(long int *lenrwbbd, long int *leniwbbd, long int *ngebbd)
{
  IDABBDPrecGetWorkSpace(IDA_idamem, lenrwbbd, leniwbbd);
  IDABBDPrecGetNumGfnEvals(IDA_idamem, ngebbd);

  return;
}
