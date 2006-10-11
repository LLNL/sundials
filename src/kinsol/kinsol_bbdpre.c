/*
 *-----------------------------------------------------------------
 * $Revision: 1.2 $
 * $Date: 2006-10-11 16:34:19 $
 *-----------------------------------------------------------------
 * Programmer(s): Allan Taylor, Alan Hindmarsh, Radu Serban, and
 *                Aaron Collier @ LLNL
 *-----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 *-----------------------------------------------------------------
 * This file contains implementations of routines for a
 * band-block-diagonal preconditioner, i.e. a block-diagonal
 * matrix with banded blocks, for use with KINSol, KINSp*
 * and the parallel implementation of NVECTOR.
 *
 * Note: With only one process, a banded matrix results
 * rather than a b-b-d matrix with banded blocks. Diagonal
 * blocking occurs at the process level.
 *-----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include <kinsol/kinsol_sptfqmr.h>
#include <kinsol/kinsol_spbcgs.h>
#include <kinsol/kinsol_spgmr.h>

#include "kinsol_impl.h"
#include "kinsol_bbdpre_impl.h"

#include <sundials/sundials_math.h>

/*
 *-----------------------------------------------------------------
 * private constants
 *-----------------------------------------------------------------
 */

#define ZERO RCONST(0.0)
#define ONE  RCONST(1.0)

/*
 *-----------------------------------------------------------------
 * prototype for difference quotient jacobian calculation routine
 *-----------------------------------------------------------------
 */

static int KBBDDQJac(KBBDPrecData pdata,
                     N_Vector uu, N_Vector uscale,
                     N_Vector gu, N_Vector gtemp, N_Vector utemp);

/*
 *-----------------------------------------------------------------
 * redability replacements
 *-----------------------------------------------------------------
 */

#define errfp    (kin_mem->kin_errfp)
#define uround   (kin_mem->kin_uround)
#define vec_tmpl (kin_mem->kin_vtemp1)

/*
 *-----------------------------------------------------------------
 * user-callable functions
 *-----------------------------------------------------------------
 */

/*
 *-----------------------------------------------------------------
 * Function : KINBBDPrecAlloc
 *-----------------------------------------------------------------
 */

void *KINBBDPrecAlloc(void *kinmem, long int Nlocal, 
		      long int mudq, long int mldq,
		      long int mukeep, long int mlkeep,
		      realtype dq_rel_uu, 
		      KINLocalFn gloc, KINCommFn gcomm)
{
  KBBDPrecData pdata;
  KINMem kin_mem;
  N_Vector vtemp3;
  long int muk, mlk, storage_mu;

  pdata = NULL;

  if (kinmem == NULL) {
    KINProcessError(NULL, 0, "KINBBDPRE", "KINBBDPrecAlloc", MSGBBD_KINMEM_NULL);
    return(NULL);
  }
  kin_mem = (KINMem) kinmem;

  /* test if the NVECTOR package is compatible with BLOCK BAND preconditioner */

  /* Note: do NOT need to check for N_VScale since it is required by KINSOL and
     so has already been checked for (see KINMalloc) */

  if (vec_tmpl->ops->nvgetarraypointer == NULL) {
    KINProcessError(kin_mem, 0, "KINBBDPRE", "KINBBDPrecAlloc", MSGBBD_BAD_NVECTOR);
    return(NULL);
  }

  pdata = NULL;
  pdata = (KBBDPrecData) malloc(sizeof *pdata);  /* allocate data memory */
  if (pdata == NULL) {
    KINProcessError(kin_mem, 0, "KINBBDPRE", "KINBBDPrecAlloc", MSGBBD_MEM_FAIL);
    return(NULL);
  }

  /* set pointers to gloc and gcomm and load half-bandwiths */

  pdata->kin_mem = kinmem;
  pdata->gloc = gloc;
  pdata->gcomm = gcomm;
  pdata->mudq = MIN(Nlocal-1, MAX(0, mudq));
  pdata->mldq = MIN(Nlocal-1, MAX(0, mldq));
  muk = MIN(Nlocal-1, MAX(0,mukeep));
  mlk = MIN(Nlocal-1, MAX(0,mlkeep));
  pdata->mukeep = muk;
  pdata->mlkeep = mlk;

  /* allocate memory for preconditioner matrix */

  storage_mu = MIN(Nlocal-1, muk+mlk);
  pdata->PP = NULL;
  pdata->PP = BandAllocMat(Nlocal, muk, mlk, storage_mu);
  if (pdata->PP == NULL) {
    free(pdata); pdata = NULL;
    return(NULL);
  }

  /* allocate memory for pivots */

  pdata->pivots = NULL;
  pdata->pivots = BandAllocPiv(Nlocal);
  if (pdata->pivots == NULL) {
    BandFreeMat(pdata->PP);
    free(pdata); pdata = NULL;
    return(NULL);
  }

  /* allocate vtemp3 for use by KBBDDQJac routine */

  vtemp3 = NULL;
  vtemp3 = N_VClone(kin_mem->kin_vtemp1);
  if (vtemp3 == NULL) {
    BandFreePiv(pdata->pivots);
    BandFreeMat(pdata->PP);
    free(pdata); pdata = NULL;
    return(NULL);
  }
  pdata->vtemp3 = vtemp3;

  /* set rel_uu based on input value dq_rel_uu */

  if (dq_rel_uu > ZERO) pdata->rel_uu = dq_rel_uu;
  else pdata->rel_uu = RSqrt(uround);  /* using dq_rel_uu = 0.0 means use default */

  /* store Nlocal to be used by the preconditioner routines */

  pdata->n_local = Nlocal;

  /* set work space sizes and initialize nge */

  pdata->rpwsize = Nlocal * (storage_mu*mlk + 1) + 1;
  pdata->ipwsize = Nlocal + 1;
  pdata->nge = 0;

  return((void *) pdata);
}

/*
 *-----------------------------------------------------------------
 * Function : KINBBDSptfqmr
 *-----------------------------------------------------------------
 */

int KINBBDSptfqmr(void *kinmem, int maxl, void *p_data)
{
  KINMem kin_mem;
  int flag;

  flag = KINSptfqmr(kinmem, maxl);
  if (flag != KINSPILS_SUCCESS) return(flag);

  kin_mem = (KINMem) kinmem;

  if (p_data == NULL) {
    KINProcessError(kin_mem, KINBBDPRE_PDATA_NULL, "KINBBDPRE", "KINBBDSptfqmr", MSGBBD_PDATA_NULL);
    return(KINBBDPRE_PDATA_NULL);
  }

  flag = KINSpilsSetPreconditioner(kinmem,
                                   KINBBDPrecSetup,
                                   KINBBDPrecSolve,
                                   p_data);
  if (flag != KINSPILS_SUCCESS) return(flag);

  return(KINSPILS_SUCCESS);
}

/*
 *-----------------------------------------------------------------
 * Function : KINBBDSpbcg
 *-----------------------------------------------------------------
 */

int KINBBDSpbcg(void *kinmem, int maxl, void *p_data)
{
  KINMem kin_mem;
  int flag;

  flag = KINSpbcg(kinmem, maxl);
  if (flag != KINSPILS_SUCCESS) return(flag);

  kin_mem = (KINMem) kinmem;

  if (p_data == NULL) {
    KINProcessError(kin_mem, KINBBDPRE_PDATA_NULL, "KINBBDPRE", "KINBBDSpbcg", MSGBBD_PDATA_NULL);
    return(KINBBDPRE_PDATA_NULL);
  }

  flag = KINSpilsSetPreconditioner(kinmem,
				   KINBBDPrecSetup,
				   KINBBDPrecSolve,
				   p_data);
  if (flag != KINSPILS_SUCCESS) return(flag);

  return(KINSPILS_SUCCESS);
}

/*
 *-----------------------------------------------------------------
 * Function : KINBBDSpgmr
 *-----------------------------------------------------------------
 */

int KINBBDSpgmr(void *kinmem, int maxl, void *p_data)
{
  KINMem kin_mem;
  int flag;

  flag = KINSpgmr(kinmem, maxl);
  if (flag != KINSPILS_SUCCESS) return(flag);

  kin_mem = (KINMem) kinmem;

  if (p_data == NULL) {
    KINProcessError(kin_mem, KINBBDPRE_PDATA_NULL, "KINBBDPRE", "KINBBDSpgmr", MSGBBD_PDATA_NULL);
    return(KINBBDPRE_PDATA_NULL);
  }

  flag = KINSpilsSetPreconditioner(kinmem,
				   KINBBDPrecSetup,
				   KINBBDPrecSolve,
				   p_data);
  if (flag != KINSPILS_SUCCESS) return(flag);

  return(KINSPILS_SUCCESS);
}

/*
 *-----------------------------------------------------------------
 * Function : KINBBDPrecFree
 *-----------------------------------------------------------------
 */

void KINBBDPrecFree(void **p_data)
{
  KBBDPrecData pdata;

  if (*p_data == NULL) return;

  pdata = (KBBDPrecData) (*p_data);
  N_VDestroy(pdata->vtemp3);
  BandFreeMat(pdata->PP);
  BandFreePiv(pdata->pivots);

  free(*p_data);
  *p_data = NULL;

}

/*
 *-----------------------------------------------------------------
 * Function : KINBBDPrecGetWorkSpace
 *-----------------------------------------------------------------
 */

int KINBBDPrecGetWorkSpace(void *p_data, long int *lenrwBBDP, long int *leniwBBDP)
{
  KBBDPrecData pdata;

  if (p_data == NULL) {
    KINProcessError(NULL, KINBBDPRE_PDATA_NULL, "KINBBDPRE", "KINBBDPrecGetWorkSpace", MSGBBD_PDATA_NULL);
    return(KINBBDPRE_PDATA_NULL);
  } 

  pdata = (KBBDPrecData) p_data;

  *lenrwBBDP = pdata->rpwsize;
  *leniwBBDP = pdata->ipwsize;

  return(KINBBDPRE_SUCCESS);
}

/*
 *-----------------------------------------------------------------
 * Function : KINBBDPrecGetNumGfnEvals
 *-----------------------------------------------------------------
 */

int KINBBDPrecGetNumGfnEvals(void *p_data, long int *ngevalsBBDP)
{
  KBBDPrecData pdata;

  if (p_data == NULL) {
    KINProcessError(NULL, KINBBDPRE_PDATA_NULL, "KINBBDPRE", "KINBBDPrecGetNumGfnEvals", MSGBBD_PDATA_NULL);
    return(KINBBDPRE_PDATA_NULL);
  } 

  pdata = (KBBDPrecData) p_data;

  *ngevalsBBDP = pdata->nge;

  return(KINBBDPRE_SUCCESS);
}

/*
 *-----------------------------------------------------------------
 * Function : KINBBDPrecGetReturnFlagName
 *-----------------------------------------------------------------
 */

char *KINBBDPrecGetReturnFlagName(int flag)
{
  char *name;

  name = (char *)malloc(30*sizeof(char));

  switch(flag) {
  case KINBBDPRE_SUCCESS:
    sprintf(name, "KINBBDPRE_SUCCESS");
    break;
  case KINBBDPRE_PDATA_NULL:
    sprintf(name, "KINBBDPRE_PDATA_NULL");
    break;
  default:
    sprintf(name, "NONE");
  }

  return(name);
}

/*
 *-----------------------------------------------------------------
 * preconditioner setup and solve functions
 *-----------------------------------------------------------------
 */
 
/*
 *-----------------------------------------------------------------
 * readability replacements
 *-----------------------------------------------------------------
 */

#define Nlocal (pdata->n_local)
#define mudq   (pdata->mudq)
#define mldq   (pdata->mldq)
#define mukeep (pdata->mukeep)
#define mlkeep (pdata->mlkeep)
#define gloc   (pdata->gloc)
#define gcomm  (pdata->gcomm)
#define pivots (pdata->pivots)
#define PP     (pdata->PP)
#define nge    (pdata->nge)
#define rel_uu (pdata->rel_uu)

/*
 *-----------------------------------------------------------------
 * Function : KINBBDPrecSetup
 *-----------------------------------------------------------------
 * KINBBDPrecSetup generates and factors a banded block of the
 * preconditioner matrix on each processor, via calls to the
 * user-supplied gloc and gcomm functions. It uses difference
 * quotient approximations to the Jacobian elements.
 *
 * KINBBDPrecSetup calculates a new Jacobian, stored in banded
 * matrix PP and does an LU factorization of P in place in PP.
 *
 * The parameters of KINBBDPrecSetup are as follows:
 *
 * uu      is the current value of the dependent variable vector,
 *         namely the solutin to func(uu)=0
 *
 * uscale  is the dependent variable scaling vector (i.e. uu)
 *
 * fval    is the vector f(u)
 *
 * fscale  is the function scaling vector
 *
 * p_data  is a pointer to user data - the same as the p_data
 *         parameter passed to KINSp*. For KINBBDPrecSetup,
 *         this should be of type KBBDData.
 *
 * vtemp1, vtemp2 are pointers to memory allocated for vectors of
 *                length N which are be used by KINBBDPrecSetup
 *                as temporary storage or work space. A third
 *                vector (vtemp3) required for KINBBDPrecSetup
 *                was previously allocated as pdata->vtemp3.
 *
 * Note: The value to be returned by the KINBBDPrecSetup function
 * is a flag indicating whether it was successful. This value is:
 *   0 if successful,
 *   > 0 for a recoverable error - step will be retried.
 *-----------------------------------------------------------------
 */

int KINBBDPrecSetup(N_Vector uu, N_Vector uscale,
		    N_Vector fval, N_Vector fscale, 
		    void *p_data,
		    N_Vector vtemp1, N_Vector vtemp2)
{
  long int ier;
  KBBDPrecData pdata;
  KINMem kin_mem;
  N_Vector vtemp3;
  int retval;

  pdata = (KBBDPrecData) p_data;

  kin_mem = (KINMem) pdata->kin_mem;

  vtemp3 = pdata->vtemp3;

  /* call KBBDDQJac for a new jacobian and store in PP */

  BandZero(PP);
  retval = KBBDDQJac(pdata, uu, uscale, vtemp1, vtemp2, vtemp3);
  if (retval != 0) {
    KINProcessError(kin_mem, KINBBDPRE_FUNC_UNRECVR, "KINBBDPRE", "KINBBDPrecSetup", MSGBBD_FUNC_FAILED);
    return(-1);
  }

  nge += (1 + MIN(mldq+mudq+1, Nlocal));

  /* do LU factorization of P in place (in PP) */

  ier = BandGBTRF(PP, pivots);

  /* return 0 if the LU was complete, else return 1 */

  if (ier > 0) return(1);
  else return(0);
}

/*
 *-----------------------------------------------------------------
 * Function : KINBBDPrecSolve
 *-----------------------------------------------------------------
 * KINBBDPrecSolve solves a linear system Pz = r, with the
 * banded blocked preconditioner matrix P generated and factored
 * by KINBBDPrecSetup. Here, r comes in as vtemp and z is
 * returned in vtemp as well.
 *
 * The parameters for KINBBDPrecSolve are as follows:
 *
 * uu     an N_Vector giving the current iterate for the system
 *
 * uscale an N_Vector giving the diagonal entries of the
 *        uu scaling matrix
 *
 * fval   an N_Vector giving the current function value
 *
 * fscale an N_Vector giving the diagonal entries of the
 *        function scaling matrix
 *
 * p_data is a pointer to user data - the same as the P_data
 *        parameter passed to KINSp*. For KINBBDPrecSolve,
 *        this should be of type KBBDData.
 *
 * vtemp  an N_Vector (temporary storage), usually the scratch
 *        vector vtemp from SPGMR/SPBCG/SPTFQMR (typical calling
 *        routine)
 *
 * Note: The value returned by the KINBBDPrecSolve function is a
 * flag indicating whether it was successful. Here this value is
 * always 0 which indicates success.
 *-----------------------------------------------------------------
 */

int KINBBDPrecSolve(N_Vector uu, N_Vector uscale,
		    N_Vector fval, N_Vector fscale, 
		    N_Vector vv, void *p_data,
		    N_Vector vtemp)
{
  KBBDPrecData pdata;
  realtype *vd;

  pdata = (KBBDPrecData) p_data;

  /* do the backsolve and return */

  vd = N_VGetArrayPointer(vv);
  BandGBTRS(PP, pivots, vd);

  return(0);
}

/*
 *-----------------------------------------------------------------
 * Function : KBBDDQJac
 *-----------------------------------------------------------------
 * This routine generates a banded difference quotient
 * approximation to the Jacobian of f(u). It assumes that a band
 * matrix of type BandMat is stored column-wise, and that elements
 * within each column are contiguous. All matrix elements are
 * generated as difference quotients, by way of calls to the user
 * routine gloc. By virtue of the band structure, the number of
 * these calls is bandwidth + 1, where bandwidth = ml + mu + 1.
 * This routine also assumes that the local elements of a vector
 * are stored contiguously.
 *-----------------------------------------------------------------
 */

#define f_data (kin_mem->kin_f_data)

static int KBBDDQJac(KBBDPrecData pdata,
                     N_Vector uu, N_Vector uscale,
                     N_Vector gu, N_Vector gtemp, N_Vector utemp)
{
  realtype inc, inc_inv;
  long int group, i, j, width, ngroups, i1, i2;
  KINMem kin_mem;
  realtype *udata, *uscdata, *gudata, *gtempdata, *utempdata, *col_j;
  int retval;

  kin_mem = (KINMem) pdata->kin_mem;

  /* set pointers to the data for all vectors */

  udata     = N_VGetArrayPointer(uu);
  uscdata   = N_VGetArrayPointer(uscale);
  gudata    = N_VGetArrayPointer(gu);
  gtempdata = N_VGetArrayPointer(gtemp);
  utempdata = N_VGetArrayPointer(utemp);

  /* load utemp with uu = predicted solution vector */

  N_VScale(ONE, uu, utemp);

  /* call gcomm and gloc to get base value of g(uu) */

  if (gcomm != NULL) {
    retval = gcomm(Nlocal, uu, f_data);
    if (retval != 0) return(retval);
  }

  retval = gloc(Nlocal, uu, gu, f_data);
  if (retval != 0) return(retval);

  /* set bandwidth and number of column groups for band differencing */

  width = mldq + mudq + 1;
  ngroups = MIN(width, Nlocal);

  /* loop over groups */
  
  for (group = 1; group <= ngroups; group++) {
  
    /* increment all u_j in group */

    for(j = group - 1; j < Nlocal; j += width) {
      inc = rel_uu * MAX(ABS(udata[j]), (ONE / uscdata[j]));
      utempdata[j] += inc;
    }
  
    /* evaluate g with incremented u */

    retval = gloc(Nlocal, utemp, gtemp, f_data);
    if (retval != 0) return(retval);

    /* restore utemp, then form and load difference quotients */

    for (j = group - 1; j < Nlocal; j += width) {
      utempdata[j] = udata[j];
      col_j = BAND_COL(PP,j);
      inc = rel_uu * MAX(ABS(udata[j]) , (ONE / uscdata[j]));
      inc_inv = ONE / inc;
      i1 = MAX(0, (j - mukeep));
      i2 = MIN((j + mlkeep), (Nlocal - 1));
      for (i = i1; i <= i2; i++)
	BAND_COL_ELEM(col_j, i, j) = inc_inv * (gtempdata[i] - gudata[i]);
    }
  }

  return(0);
}
