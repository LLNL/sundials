/*******************************************************************
 * File          : kinbbdpre.c                                     *
 * Programmers   : Allan G Taylor, Alan C Hindmarsh, and           *
 *                 Radu Serban @ LLNL                              *
 * Version of    : 07 February 2004                                *
 *-----------------------------------------------------------------*
 * Copyright (c) 2002, The Regents of the University of California * 
 * Produced at the Lawrence Livermore National Laboratory          *
 * All rights reserved                                             *
 * For details, see sundials/kinsol/LICENSE                        *
 *-----------------------------------------------------------------*
 * This file contains implementations of routines for a            *
 * band-block-diagonal preconditioner, i.e. a block-diagonal       *
 * matrix with banded blocks, for use with KINSol, KINSpgmr, and   *
 * the parallel implementation of NVECTOR.                         *
 * NOTE: with only one processor in use, a banded matrix results   *
 * rather than a b-b-d matrix with banded blocks. Diagonal         *
 * blocking occurs at the processor level.                         *
 *******************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "kinbbdpre.h"
#include "kinspgmr.h"
#include "sundialsmath.h"
#include "iterative.h"

#define ZERO         RCONST(0.0)
#define ONE          RCONST(1.0)

/* Error Messages */
#define KBBDALLOC      "KBBDAlloc-- "
#define MSG_KINMEM_NULL KBBDALLOC "KINSOL Memory is NULL.\n\n"
#define MSG_WRONG_NVEC  KBBDALLOC "Incompatible NVECTOR implementation.\n\n"
#define MSG_PDATA_NULL  "KBBDPrecGet*-- BBDPrecData is NULL. \n\n"

#define MSG_NO_PDATA   "KBBDSpgmr-- BBDPrecData is NULL. \n\n"

/* Prototype for difference quotient Jacobian calculation routine */

static void KBBDDQJac(KBBDPrecData pdata,
                      N_Vector uu, N_Vector uscale,
                      N_Vector gu, N_Vector gtemp, N_Vector utemp);

/* Redability replacements */
#define nvspec  (kin_mem->kin_nvspec)
#define errfp   (kin_mem->kin_errfp)
#define uround  (kin_mem->kin_uround)

/***************** User-Callable Functions: malloc and free ******************/

void *KINBBDPrecAlloc(void *kinmem, long int Nlocal, 
		      long int mu, long int ml,
		      realtype dq_rel_uu, 
		      KINLocalFn gloc, KINCommFn gcomm)
{
  KBBDPrecData pdata;
  KINMem kin_mem;
  N_Vector vtemp3;
  realtype rel_uu;

  if (kinmem == NULL) {
    fprintf(stdout, MSG_KINMEM_NULL);
    return(NULL);
  }
  kin_mem = (KINMem) kinmem;

  /* Test if the NVECTOR package is compatible with the BLOCK BAND preconditioner */
  if (nvspec->ops->nvgetdata == NULL || nvspec->ops->nvsetdata == NULL) {
    fprintf(errfp, MSG_WRONG_NVEC);
    return(NULL);
  }  

  pdata = (KBBDPrecData) malloc(sizeof *pdata);  /* Allocate data memory. */
  if (pdata == NULL) return(NULL);

  /* Set pointers to gloc and gcomm; load half-bandwiths. */
  pdata->kin_mem = kin_mem;
  pdata->ml = ml;
  pdata->mu = mu;
  pdata->gloc = gloc;
  pdata->gcomm = gcomm;
 
  /* Allocate memory for preconditioner matrix. */
  pdata->PP = BandAllocMat(Nlocal, mu, ml, mu+ml);
  if (pdata->PP == NULL) {
    free(pdata);
    return(NULL);
  }

  /* Allocate memory for pivots. */
  pdata->pivots = BandAllocPiv(Nlocal);
  if (pdata->PP == NULL) {
    BandFreeMat(pdata->PP);
    free(pdata);
    return(NULL);
  }

  /* Allocate vtemp3 for use by KBBDDQJac. */

  vtemp3 = N_VNew(nvspec);

  if (vtemp3 == NULL) {
    free(pdata);
    return(NULL);
  }
  pdata->vtemp3 = vtemp3;

  /* set rel_uu based on input value dq_rel_uu */

  if (dq_rel_uu > ZERO) rel_uu = dq_rel_uu;
  else rel_uu = RSqrt(kin_mem->kin_uround); /* A value dq_rel_uu = 0.0 received
                                               by this routine implies using
                                               the default instead.           */

  pdata->rel_uu = rel_uu;

  /* Store Nlocal to be used in KBBDPrecon. */
  pdata->n_local = Nlocal;

  /* Set work space sizes and initialize nge. */
  pdata->rpwsize = Nlocal*(2*mu + ml + 1);
  pdata->ipwsize = Nlocal;
  pdata->nge = 0;

  return((void *)pdata);
}

int KINBBDSpgmr(void *kinmem, int maxl, void *p_data)
{
  int flag;

  if ( p_data == NULL ) {
    fprintf(stdout, MSG_NO_PDATA);
    return(BBDP_NO_PDATA);
  }

  flag = KINSpgmr(kinmem, maxl);
  if(flag != SUCCESS) return(flag);

  flag = KINSpgmrSetPrecData(kinmem, p_data);
  if(flag != SUCCESS) return(flag);

  flag = KINSpgmrSetPrecSetupFn(kinmem, KINBBDPrecSetup);
  if(flag != SUCCESS) return(flag);

  flag = KINSpgmrSetPrecSolveFn(kinmem, KINBBDPrecSolve);
  if(flag != SUCCESS) return(flag);

  return(SUCCESS);
}

void KINBBDPrecFree(void *p_data)
{
  KBBDPrecData pdata;

  if (p_data != NULL) {
    pdata = (KBBDPrecData) p_data;
    BandFreeMat(pdata->PP);
    BandFreePiv(pdata->pivots);
    free(pdata);
  }
}

int KINBBDPrecGetIntWorkSpace(void *p_data, long int *leniwBBDP)
{
  KBBDPrecData pdata;

  if ( p_data == NULL ) {
    fprintf(stdout, MSG_PDATA_NULL);
    return(BBDP_NO_PDATA);
  } 

  pdata = (KBBDPrecData) p_data;

  *leniwBBDP = pdata->ipwsize;

  return(OKAY);
}

int KINBBDPrecGetRealWorkSpace(void *p_data, long int *lenrwBBDP)
{
  KBBDPrecData pdata;

  if ( p_data == NULL ) {
    fprintf(stdout, MSG_PDATA_NULL);
    return(BBDP_NO_PDATA);
  } 

  pdata = (KBBDPrecData) p_data;

  *lenrwBBDP = pdata->rpwsize;

  return(OKAY);
}

int KINBBDPrecGetNumGfnEvals(void *p_data, long int *ngevalsBBDP)
{
  KBBDPrecData pdata;

  if ( p_data == NULL ) {
    fprintf(stdout, MSG_PDATA_NULL);
    return(BBDP_NO_PDATA);
  } 

  pdata = (KBBDPrecData) p_data;

  *ngevalsBBDP = pdata->nge;

  return(OKAY);
}


/***************** Preconditioner setup and solve Functions ****************/
 

/* Readability Replacements */

#define Nlocal    (pdata->n_local)
#define mu        (pdata->mu)
#define ml        (pdata->ml)
#define gloc      (pdata->gloc)
#define gcomm     (pdata->gcomm)
#define pivots    (pdata->pivots)
#define PP        (pdata->PP)
#define nge       (pdata->nge)
#define rel_uu    (pdata->rel_uu)


/******************************************************************
 * Function : KINBBDPrecSetup                                     *
 *----------------------------------------------------------------*
 * KINBBDPrecSetup generates and factors a banded block of the    *
 * preconditioner matrix on each processor, via calls to the      *
 * user-supplied gloc and gcomm functions. It uses difference     *
 * quotient approximations to the Jacobian elements.              *
 *                                                                *
 * KINBBDPrecSetup calculates a new Jacobian, stored in banded    *
 * matrix PP and does an LU factorization of P in place  in PP.   *
 *                                                                *
 * The parameters of KINBBDPrecSetup are as follows:              *
 *                                                                *
 * uu      is the current value of the dependent variable vector, *
 *         namely the solutin to func(uu)=0.                      *
 *                                                                *
 * uscale  is the dependent variable scaling vector (i.e. uu)     *
 *                                                                *
 * fval    is the vector f(u).                                    *
 *                                                                *
 * fscale  is the function scaling vector.                        *
 *                                                                *
 * vtemp1, vtemp2 are pointers to memory allocated                *
 *           for vectors of length N which are be used by         *
 *           KINBBDPrecSetup as temporary storage or work space.  *
 *           A third vector vtemp3 required for KINBBDPrecSetup   *
 *           has been previously allocated as pdata->vtemp3       *
 *                                                                *
 * p_data  is a pointer to user data - the same as the p_data     *
 *           parameter passed to KINSpgmr. For KINBBDPrecSetup,   *
 *           this should be of type KBBDData.                     *
 *                                                                *
 * Return value:                                                  *
 * The value to be returned by the KINBBDPrecSetup function is a  *
 * flag indicating whether it was successful.  This value is      *
 *   0   if successful,                                           *
 *   > 0 for a recoverable error (step will be retried).          *
 ******************************************************************/

int KINBBDPrecSetup(N_Vector uu, N_Vector uscale,
		    N_Vector fval, N_Vector fscale, 
		    void *p_data,
		    N_Vector vtemp1, N_Vector vtemp2)
{
  long int ier;
  KBBDPrecData pdata;
  N_Vector vtemp3;

  pdata = (KBBDPrecData) p_data;
  vtemp3 = pdata->vtemp3;

  /* Call KBBDDQJac for a new Jacobian and store in PP. */
  BandZero(PP);
  KBBDDQJac(pdata, uu, uscale, vtemp1, vtemp2, vtemp3);
  nge += 1 + MIN(ml + mu + 1, Nlocal);

  /* Do LU factorization of P in place (in PP). */
  ier = BandFactor(PP, pivots);

  /* Return 0 if the LU was complete; otherwise return 1. */
  if (ier > 0) return(1);
  return(0);
}


/******************************************************************
 * Function : KINBBDPrecSolve                                     *
 *----------------------------------------------------------------*
 * KINBBDPrecSolve solves a linear system P z = r, with the       *
 * banded blocked preconditioner matrix P generated and factored  *
 * by KINBBDPrecSetup.  Here, r comes in as vtem and z is         *
 * returned in vtem as well.                                      *
 *                                                                *
 * The parameters of KINBBDPrecSolve are as follows:              *
 *                                                                *
 * uu     an N_Vector giving the current iterate for the system   *
 *                                                                *
 * uscale an N_Vector giving the diagonal entries                 *
 *          of the uu-scaling matrix                              *
 *                                                                *
 * fval   an N_Vector giving the current function value           *
 *                                                                *
 * fscale an N_Vector giving the diagonal entries of the          *
 *         function scaling matrix                                *
 *                                                                *
 * vtem   an N_Vector temporary, usually the scratch vector vtemp *
 *          from spgmr, the typical calling routine               *
 *                                                                *
 * p_data is a pointer to user data - the same as the P_data      *
 *          parameter passed to KINSpgmr. For KINBBDPrecSolve,    *
 *          this should be of type KBBDData.                      *
 *                                                                *
 * Return value:                                                  *
 * The value returned by the KINBBDPrecSolve function is a flag   *
 * indicating whether it was successful.  Here this value is      *
 * always 0, indicating success.                                  *
 ******************************************************************/

int KINBBDPrecSolve(N_Vector uu, N_Vector uscale,
		    N_Vector fval, N_Vector fscale, 
		    N_Vector vv, void *p_data,
		    N_Vector vtemp)
{
  KBBDPrecData pdata;
  realtype *vd;

  pdata = (KBBDPrecData) p_data;

  /* Do the backsolve and return. */
  vd = N_VGetData(vv);
  BandBacksolve(PP, pivots, vd);
  N_VSetData(vd, vv);

  return(0);
}

/*************** KBBDDQJac *****************************************

 This routine generates a banded difference quotient approximation to
 the Jacobian of f(u).  It assumes that a band matrix of type
 BandMat is stored column-wise, and that elements within each column
 are contiguous.  All matrix elements are generated as difference
 quotients, by way of calls to the user routine gloc.
 By virtue of the band structure, the number of these calls is
 bandwidth + 1, where bandwidth = ml + mu + 1.
 This routine also assumes that the local elements of a vector are
 stored contiguously.

**********************************************************************/

#define f_data (kin_mem->kin_f_data)

static void KBBDDQJac(KBBDPrecData pdata,
                      N_Vector uu, N_Vector uscale,
                      N_Vector gu, N_Vector gtemp, N_Vector utemp)
{
  KINMem kin_mem;
  realtype inc, inc_inv;
  long int group, i, j, width, ngroups, i1, i2;
  realtype *udata, *uscdata, *gudata, *gtempdata, *utempdata, *col_j;

  kin_mem = pdata->kin_mem;

  /* Set pointers to the data for all vectors. */
  udata     = N_VGetData(uu);
  uscdata   = N_VGetData(uscale);
  gudata    = N_VGetData(gu);
  gtempdata = N_VGetData(gtemp);
  utempdata = N_VGetData(utemp);

  /* Load utemp with uu = predicted solution vector. */
  N_VScale(ONE, uu, utemp);

  /* Call gcomm and gloc to get base value of g(uu). */
  gcomm (Nlocal, uu, f_data);
  gloc (Nlocal, uu, gu, f_data);

  /* Set bandwidth and number of column groups for band differencing. */
  width = ml + mu + 1;
  ngroups = MIN(width, Nlocal);

  /* Loop over groups. */  
  for (group=1; group <= ngroups; group++) {
  
    /* Increment all u_j in group. */
    for(j=group-1; j < Nlocal; j+=width) {
      inc = rel_uu * MAX(ABS(udata[j]), ONE/uscdata[j]);
      utempdata[j] += inc;
    }
  
    /* Evaluate g with incremented u. */
    N_VSetData(utempdata, utemp);
    gloc (Nlocal, utemp, gtemp, f_data);
    gtempdata = N_VGetData(gtemp);

    /* Restore utemp, then form and load difference quotients. */
    for (j=group-1; j < Nlocal; j+=width) {
      utempdata[j] = udata[j];
      col_j = BAND_COL(PP,j);
      inc =rel_uu * MAX( ABS(udata[j]) , ONE/uscdata[j] );
      inc_inv = ONE/inc;
      i1 = MAX(0, j-mu);
      i2 = MIN(j+ml, Nlocal-1);
      for (i=i1; i <= i2; i++)
        BAND_COL_ELEM(col_j,i,j) = inc_inv * (gtempdata[i] - gudata[i]);
    }
  }
}
