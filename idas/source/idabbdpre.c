/*
 * -----------------------------------------------------------------
 * $Revision: 1.11 $
 * $Date: 2004-10-26 20:17:12 $
 * ----------------------------------------------------------------- 
 * Programmers: Alan C. Hindmarsh, and Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California  
 * Produced at the Lawrence Livermore National Laboratory
 * All rights reserved
 * For details, see sundials/idas/LICENSE
 * -----------------------------------------------------------------
 * This file contains implementations of routines for a            
 * band-block-diagonal preconditioner, i.e. a block-diagonal       
 * matrix with banded blocks, for use with IDA and IDASpgmr.       
 * NOTE: with only one processor in use, a banded matrix results   
 * rather than a block-diagonal matrix with banded blocks.         
 * Diagonal blocking occurs at the processor level.                
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "idas_impl.h"
#include "idabbdpre_impl.h"
#include "idaspgmr_impl.h"
#include "sundialsmath.h"
#include "iterative.h"

#define ZERO         RCONST(0.0)
#define ONE          RCONST(1.0)
#define TWO          RCONST(2.0)

/* Error messages */

#define IDABBDALLOC      "IBBDPrecAlloc-- "
#define MSG_IDAMEM_NULL  IDABBDALLOC "Integrator Memory is NULL.\n\n"
#define MSG_BAD_NVECTOR  IDABBDALLOC "A required vector operation is not implemented.\n\n"
#define MSG_WRONG_NVEC   IDABBDALLOC "Incompatible NVECTOR implementation.\n\n"
#define MSG_PDATA_NULL   "IBBDPrecGet*-- BBDPrecData is NULL. \n\n"

#define MSG_NO_PDATA     "IBBDSpgmr-- BBDPrecData is NULL. \n\n"

/* Prototype for difference quotient Jacobian calculation routine */

static int IBBDDQJac(IBBDPrecData pdata, realtype tt, realtype cj,
                     N_Vector yy, N_Vector yp, N_Vector gref, 
                     N_Vector ytemp, N_Vector yptemp, N_Vector gtemp);

/**********************************************************************/

/* Readability Replacements */

#define errfp    (IDA_mem->ida_errfp)
#define uround   (IDA_mem->ida_uround)
#define vec_tmpl (IDA_mem->ida_tempv1)

/*
 * -----------------------------------------------------------------
 * User-Callable Functions: malloc, reinit, and free
 * -----------------------------------------------------------------
 */

void *IDABBDPrecAlloc(void *ida_mem, long int Nlocal, 
		      long int mudq, long int mldq, 
		      long int mukeep, long int mlkeep, 
		      realtype dq_rel_yy, 
		      IDABBDLocalFn Gres, IDABBDCommFn Gcomm)
{
  IDAMem IDA_mem;
  IBBDPrecData pdata;
  N_Vector tempv4;
  long int muk, mlk, storage_mu;

  if (ida_mem == NULL) {
    fprintf(stderr, MSG_IDAMEM_NULL);
    return(NULL);
  }

  IDA_mem = (IDAMem) ida_mem;

  /* Test if the NVECTOR package is compatible with BLOCK BAND preconditioner */
  if(vec_tmpl->ops->nvgetarraypointer == NULL) {
    if(errfp!=NULL) fprintf(errfp, MSG_BAD_NVECTOR);
    return(NULL);
  }

  /* Allocate data memory. */
  pdata = (IBBDPrecData) malloc(sizeof *pdata);
  if (pdata == NULL) return(NULL);

  /* Set pointers to glocal and gcomm; load half-bandwidths. */
  pdata->IDA_mem = IDA_mem;
  pdata->glocal = Gres;
  pdata->gcomm = Gcomm;
  pdata->mudq = MIN( Nlocal-1, MAX(0,mudq) );
  pdata->mldq = MIN( Nlocal-1, MAX(0,mldq) );
  muk = MIN( Nlocal-1, MAX(0,mukeep) );
  mlk = MIN( Nlocal-1, MAX(0,mlkeep) );
  pdata->mukeep = muk;
  pdata->mlkeep = mlk;

  /* Set extended upper half-bandwidth for PP (required for pivoting). */
  storage_mu = MIN(Nlocal-1, muk + mlk);

  /* Allocate memory for preconditioner matrix. */
  pdata->PP = BandAllocMat(Nlocal, muk, mlk, storage_mu);
  if (pdata->PP == NULL) { free(pdata); return(NULL); }

  /* Allocate memory for pivots. */
  pdata->pivots = BandAllocPiv(Nlocal);
  if (pdata->PP == NULL) {
    BandFreeMat(pdata->PP);
    free(pdata);
    return(NULL);
  }

  /* Allocate tempv4 for use by IBBDDQJac */
  tempv4 = N_VClone(vec_tmpl); 
  if (tempv4 == NULL){
    BandFreeMat(pdata->PP);
    BandFreePiv(pdata->pivots);
    free(pdata);
    return(NULL);
  }
  pdata->tempv4 = tempv4;
  
  /* Set rel_yy based on input value dq_rel_yy (0 implies default). */
  pdata->rel_yy = (dq_rel_yy > ZERO) ? dq_rel_yy : RSqrt(uround); 

  /* Store Nlocal to be used in IDABBDPrecSetup */
  pdata->n_local = Nlocal;
  
  /* Set work space sizes and initialize nge. */
  pdata->rpwsize = Nlocal*(mlk + storage_mu + 1);
  pdata->ipwsize = Nlocal;
  pdata->nge = 0;

  return((void *)pdata);
}

int IDABBDSpgmr(void *ida_mem, int maxl, void *bbd_data)
{
  int flag;

  if ( bbd_data == NULL ) {
    fprintf(stderr, MSG_NO_PDATA);
    return(IDA_PDATA_NULL);
  }

  flag = IDASpgmr(ida_mem, maxl);
  if(flag != IDASPGMR_SUCCESS) return(flag);

  flag = IDASpgmrSetPrecData(ida_mem, bbd_data);
  if(flag != IDASPGMR_SUCCESS) return(flag);

  flag = IDASpgmrSetPrecSetupFn(ida_mem, IDABBDPrecSetup);
  if(flag != IDASPGMR_SUCCESS) return(flag);

  flag = IDASpgmrSetPrecSolveFn(ida_mem, IDABBDPrecSolve);
  if(flag != IDASPGMR_SUCCESS) return(flag);

  return(IDASPGMR_SUCCESS);
}

int IDABBDPrecReInit(void *bbd_data,
		     long int mudq, long int mldq, 
		     realtype dq_rel_yy,  
		     IDABBDLocalFn Gres, IDABBDCommFn Gcomm)
{
  IBBDPrecData pdata;
  IDAMem IDA_mem;
  long int Nlocal;

  pdata =(IBBDPrecData) bbd_data;
  IDA_mem = pdata->IDA_mem;

  Nlocal = pdata->n_local;

  /* Set pointers to res_data, glocal, and gcomm; load half-bandwidths. */
  pdata->mudq = MIN( Nlocal-1, MAX(0,mudq) );
  pdata->mldq = MIN( Nlocal-1, MAX(0,mldq) );
  pdata->glocal = Gres;
  pdata->gcomm = Gcomm;

  /* Set rel_yy based on input value dq_rel_yy (0 implies default). */
  pdata->rel_yy = (dq_rel_yy > ZERO) ? dq_rel_yy : RSqrt(uround); 

  /* Re-initialize nge */
  pdata->nge = 0;

  return(0);
}


void IDABBDPrecFree(void *bbd_data)
{
  IBBDPrecData pdata;
  
  if ( bbd_data != NULL ) {
    pdata = (IBBDPrecData) bbd_data;
    BandFreeMat(pdata->PP);
    BandFreePiv(pdata->pivots);
    N_VDestroy(pdata->tempv4);
    free(pdata);
  }
}

int IDABBDPrecGetWorkSpace(void *bbd_data, long int *lenrwBBDP, long int *leniwBBDP)
{
  IBBDPrecData pdata;

  if ( bbd_data == NULL ) {
    fprintf(stderr, MSG_PDATA_NULL);
    return(IDA_PDATA_NULL);
  } 

  pdata = (IBBDPrecData) bbd_data;

  *lenrwBBDP = pdata->rpwsize;
  *leniwBBDP = pdata->ipwsize;

  return(IDA_SUCCESS);
}

int IDABBDPrecGetNumGfnEvals(void *bbd_data, long int *ngevalsBBDP)
{
  IBBDPrecData pdata;

  if ( bbd_data == NULL ) {
    fprintf(stderr, MSG_PDATA_NULL);
    return(IDA_PDATA_NULL);
  } 

  pdata = (IBBDPrecData) bbd_data;

  *ngevalsBBDP = pdata->nge;

  return(IDA_SUCCESS);
}

/* Readability Replacements */

#define Nlocal      (pdata->n_local)
#define mudq        (pdata->mudq)
#define mldq        (pdata->mldq)
#define mukeep      (pdata->mukeep)
#define mlkeep      (pdata->mlkeep)
#define glocal      (pdata->glocal)
#define gcomm       (pdata->gcomm)
#define pivots      (pdata->pivots)
#define PP          (pdata->PP)
#define nge         (pdata->nge)
#define rel_yy      (pdata->rel_yy)

/*
 * -----------------------------------------------------------------
 * Function : IDABBDPrecSetup                                     
 *----------------------------------------------------------------
 * IDABBDPrecSetup generates a band-block-diagonal preconditioner 
 * matrix, where the local block (on this processor) is a band    
 * matrix.  Each local block is computed by a difference quotient 
 * scheme via calls to the user-supplied routines glocal, gcomm.  
 * After generating the block in the band matrix PP, this routine 
 * does an LU factorization in place in PP.                       
 *                                                                
 * The IDABBDPrecSetup parameters used here are as follows:       
 *                                                                
 * tt is the current value of the independent variable t.        
 *                                                                
 * yy is the current value of the dependent variable vector, 
 *    namely the predicted value of y(t).                    
 *                                                                
 * yp is the current value of the derivative vector y',          
 *    namely the predicted value of y'(t).                    
 *                                                                
 * c_j is the scalar in the system Jacobian, proportional to 1/hh.
 *                                                                
 * prec_data is a pointer to user preconditioner data - the same as 
 *    the p_data parameter passed to IDASpgmr.                
 *                                                                
 * tmp1, tmp2, tmp3 are pointers to vectors of type         
 *    N_Vector, used for temporary storage or work space.            
 *                                                                  
 * The arguments Neq, rr, res, uround, and nrePtr are not used.   
 *                                                                  
 * Return value:                                                  
 * The value returned by this IDABBDPrecSetup function is a int   
 * flag indicating whether it was successful.  This value is      
 *    0    if successful,                                         
 *  > 0    for a recoverable error (step will be retried).        
 *  < 0    for a nonrecoverable error (step fails).               
 *
 * -----------------------------------------------------------------
 */

int IDABBDPrecSetup(realtype tt, 
		    N_Vector yy, N_Vector yp, N_Vector rr, 
		    realtype c_j, void *prec_data,
		    N_Vector tempv1, N_Vector tempv2, N_Vector tempv3)
{
  long int retfac;
  int retval;
  IBBDPrecData pdata;

  pdata =(IBBDPrecData) prec_data;

  /* Call IBBDDQJac for a new Jacobian calculation and store in PP. */
  BandZero(PP);
  retval = IBBDDQJac(pdata, tt, c_j, yy, yp,
                     tempv1, tempv2, tempv3, pdata->tempv4);
  if (retval < 0) return(-1);
  if (retval > 0) return(+1);
 
  /* Do LU factorization of preconditioner block in place (in PP). */
  retfac = BandFactor(PP, pivots);

  /* Return 0 if the LU was complete, or +1 otherwise. */
  if (retfac > 0) return(+1);
  return(0);
}


/*
 * -----------------------------------------------------------------
 * Function: IDABBDPrecSolve                                      
 *----------------------------------------------------------------
 * The function IDABBDPrecSolve computes a solution to the linear 
 * system P z = r, where P is the left preconditioner defined by  
 * the routine IDABBDPrecSetup.                                   
 *                                                                
 * The IDABBDPrecSolve parameters used here are as follows:       
 *                                                                
 * rvec is the input right-hand side vector r.                  
 *                                                                
 * zvec is the computed solution vector z.                      
 *                                                                
 * prec_data is a pointer to user preconditioner data - the same as  
 *     the p_data parameter passed to IDASpgmr.                
 *                                                                
 * The arguments tt, yy, yp, rr, c_j, delta, and tmp are NOT used.
 *                                                                
 * IDABBDPrecSolve always returns 0, indicating success.          
 *                                                                
 * -----------------------------------------------------------------
 */

int IDABBDPrecSolve(realtype tt, 
		    N_Vector yy, N_Vector yp, N_Vector rr, 
		    N_Vector rvec, N_Vector zvec,
		    realtype c_j, realtype delta, void *prec_data, 
                    N_Vector tmp)
{
  IBBDPrecData pdata;
  realtype *zd;

  pdata = (IBBDPrecData) prec_data;

  /* Copy rvec to zvec, do the backsolve, and return. */
  N_VScale(ONE, rvec, zvec);

  zd = N_VGetArrayPointer(zvec);

  BandBacksolve(PP, pivots, zd);

  return(0);
}

/*
 * -----------------------------------------------------------------
 * IBBDDQJac
 * -----------------------------------------------------------------
 * This routine generates a banded difference quotient approximation to
 * the local block of the Jacobian of G(t,y,y').  It assumes that a
 * band matrix of type BandMat is stored column-wise, and that elements
 * within each column are contiguous.
 *
 * All matrix elements are generated as difference quotients, by way
 * of calls to the user routine glocal.
 * By virtue of the band structure, the number of these calls is
 * bandwidth + 1, where bandwidth = mldq + mudq + 1.
 * But the band matrix kept has bandwidth = mlkeep + mukeep + 1.
 * This routine also assumes that the local elements of a vector are
 * stored contiguously.
 *
 * Return values are: 0 (success), > 0 (recoverable error),
 * or < 0 (nonrecoverable error).
 * -----------------------------------------------------------------
 */

#define ewt         (IDA_mem->ida_ewt)
#define res_data    (IDA_mem->ida_rdata)
#define hh          (IDA_mem->ida_hh)
#define constraints (IDA_mem->ida_constraints)

static int IBBDDQJac(IBBDPrecData pdata, realtype tt, realtype cj,
                     N_Vector yy, N_Vector yp, N_Vector gref, 
                     N_Vector ytemp, N_Vector yptemp, N_Vector gtemp)
{
  IDAMem IDA_mem;
  realtype inc, inc_inv;
  int  retval;
  long int group, i, j, width, ngroups, i1, i2;
  realtype *ydata, *ypdata, *ytempdata, *yptempdata, *grefdata, *gtempdata;
  realtype *cnsdata = NULL, *ewtdata;
  realtype *col_j, conj, yj, ypj, ewtj;

  IDA_mem = pdata->IDA_mem;

  /* Initialize ytemp and yptemp. */

  N_VScale(ONE, yy, ytemp);
  N_VScale(ONE, yp, yptemp);

  /* Obtain pointers as required to the data array of vectors. */

  ydata     = N_VGetArrayPointer(yy);
  ypdata    = N_VGetArrayPointer(yp);
  gtempdata = N_VGetArrayPointer(gtemp);
  ewtdata   = N_VGetArrayPointer(ewt);
  if (constraints != NULL) 
    cnsdata = N_VGetArrayPointer(constraints);
  ytempdata = N_VGetArrayPointer(ytemp);
  yptempdata= N_VGetArrayPointer(yptemp);
  grefdata = N_VGetArrayPointer(gref);

  /* Call gcomm and glocal to get base value of G(t,y,y'). */

  if (gcomm != NULL) {
    retval = gcomm(Nlocal, tt, yy, yp, res_data);
    if (retval != 0) return(retval);
  }

  retval = glocal(Nlocal, tt, yy, yp, gref, res_data); 
  nge++;
  if (retval != 0) return(retval);


  /* Set bandwidth and number of column groups for band differencing. */

  width = mldq + mudq + 1;
  ngroups = MIN(width, Nlocal);

  /* Loop over groups. */
  for(group = 1; group <= ngroups; group++) {
    
    /* Loop over the components in this group. */
    for(j = group-1; j < Nlocal; j += width) {
      yj = ydata[j];
      ypj = ypdata[j];
      ewtj = ewtdata[j];
      
      /* Set increment inc to yj based on rel_yy*abs(yj), with
         adjustments using ypj and ewtj if this is small, and a further
         adjustment to give it the same sign as hh*ypj. */
      inc = rel_yy*MAX(ABS(yj), MAX( ABS(hh*ypj), ONE/ewtj));
      if (hh*ypj < ZERO) inc = -inc;
      inc = (yj + inc) - yj;
      
      /* Adjust sign(inc) again if yj has an inequality constraint. */
      if (constraints != NULL) {
        conj = cnsdata[j];
        if (ABS(conj) == ONE)      {if ((yj+inc)*conj <  ZERO) inc = -inc;}
        else if (ABS(conj) == TWO) {if ((yj+inc)*conj <= ZERO) inc = -inc;}
      }

      /* Increment yj and ypj. */
      ytempdata[j] += inc;
      yptempdata[j] += cj*inc;
      
    }

    /* Evaluate G with incremented y and yp arguments. */

    retval = glocal(Nlocal, tt, ytemp, yptemp, gtemp, res_data); 
    nge++;
    if (retval != 0) return(retval);

    /* Loop over components of the group again; restore ytemp and yptemp. */
    for(j = group-1; j < Nlocal; j += width) {
      yj  = ytempdata[j]  = ydata[j];
      ypj = yptempdata[j] = ypdata[j];
      ewtj = ewtdata[j];

      /* Set increment inc as before .*/
      inc = rel_yy*MAX(ABS(yj), MAX( ABS(hh*ypj), ONE/ewtj));
      if (hh*ypj < ZERO) inc = -inc;
      inc = (yj + inc) - yj;
      if (constraints != NULL) {
        conj = cnsdata[j];
        if (ABS(conj) == ONE)      {if ((yj+inc)*conj <  ZERO) inc = -inc;}
        else if (ABS(conj) == TWO) {if ((yj+inc)*conj <= ZERO) inc = -inc;}
      }

      /* Form difference quotients and load into PP. */
      inc_inv = ONE/inc;
      col_j = BAND_COL(PP,j);
      i1 = MAX(0, j-mukeep);
      i2 = MIN(j+mlkeep, Nlocal-1);
      for(i = i1; i <= i2; i++) BAND_COL_ELEM(col_j,i,j) =
                                  inc_inv * (gtempdata[i] - grefdata[i]);
    }
  }
  
  return(0);
}
