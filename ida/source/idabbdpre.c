/******************************************************************
 *                                                                *
 * File          : idabbdpre.c                                    *
 * Programmers   : Allan G Taylor and Alan C Hindmarsh @ LLNL     *
 * Version of    : 2 March 2001                                   *
 *----------------------------------------------------------------*
 * This file contains implementations of routines for a           *
 * band-block-diagonal preconditioner, i.e. a block-diagonal      *
 * matrix with banded blocks, for use with IDA and IDASpgmr.      *
 * NOTE: with only one processor in use, a banded matrix results  *
 * rather than a block-diagonal matrix with banded blocks.        *
 * Diagonal blocking occurs at the processor level.               *
 ******************************************************************/

#include <stdlib.h>
#include "idabbdpre.h"
#include "ida.h"   
#include "llnltyps.h"
#include "nvector.h"
#include "llnlmath.h"
#include "iterativ.h"
#include "band.h"

#define ZERO         RCONST(0.0)
#define ONE          RCONST(1.0)
#define TWO          RCONST(2.0)


/* Prototype for difference quotient Jacobian calculation routine */

static int IBBDDQJac(integer Nlocal, integer mudq, integer mldq, 
                     integer mukeep, integer mlkeep, 
                     real cj, real hh, real rel_yy, real tt,
                     N_Vector ewt, N_Vector constraints,
                     IDALocalFn glocal, IDACommFn gcomm, BandMat JJ, 
                     integer *nginc, void *res_data, N_Vector yy,  
                     N_Vector yp, N_Vector gref, N_Vector gtemp, 
                     N_Vector ytemp, N_Vector yptemp);



/***************** User-Callable Functions: malloc and free ******************/

IBBDData IBBDAlloc(integer Nlocal, integer mudq, integer mldq, 
                   integer mukeep, integer mlkeep, real dq_rel_yy, 
                   IDALocalFn glocal, IDACommFn gcomm, 
		   void *idamem, void *res_data)
{
  IBBDData P_data;
  IDAMem ida_mem;
  N_Vector tempv4;
  real rel_yy;
  integer muk, mlk, storage_mu;

  ida_mem = (IDAMem)idamem;

  /* Allocate data memory. */
  P_data = (IBBDData) malloc(sizeof *P_data);
  if(P_data == NULL) return(NULL);

  /* Set pointers to res_data, glocal, and gcomm; load half-bandwidths. */
  P_data->res_data = res_data;
  P_data->mudq = MIN( Nlocal-1, MAX(0,mudq) );
  P_data->mldq = MIN( Nlocal-1, MAX(0,mldq) );
  muk = MIN( Nlocal-1, MAX(0,mukeep) );
  mlk = MIN( Nlocal-1, MAX(0,mlkeep) );
  P_data->mukeep = muk;
  P_data->mlkeep = mlk;
  P_data->glocal = glocal;
  P_data->gcomm = gcomm;

  /* Set extended upper half-bandwidth for PP (required for pivoting). */
  storage_mu = MIN(Nlocal-1, muk + mlk);

  /* Allocate memory for preconditioner matrix. */
  P_data->PP = BandAllocMat(Nlocal, muk, mlk, storage_mu);
  if(P_data->PP == NULL) { free(P_data); return(NULL); }

  /* Allocate memory for pivots. */
  P_data->pivots = BandAllocPiv(Nlocal);
  if(P_data->PP == NULL) {
    BandFreeMat(P_data->PP);
    free(P_data);
    return(NULL);
  }

  /* Allocate tempv4 for use by IBBDDQJac.  Note: Nlocal is a dummy,
     in that ida_machenv parameters are used to determine size. */
 tempv4 = N_VNew(Nlocal,ida_mem->ida_machenv); 
 if(tempv4 == NULL){
   BandFreeMat(P_data->PP);
   BandFreePiv(P_data->pivots);
   free(P_data);
   return(NULL);
 }
  P_data->tempv4 = tempv4;

  /* Set rel_yy based on input value dq_rel_yy (0 implies default). */
  if(dq_rel_yy > ZERO) rel_yy = dq_rel_yy;
  else               rel_yy = RSqrt(ida_mem->ida_uround); 
  P_data->rel_yy = rel_yy;

  /* Set work space sizes and initialize nge. */
  P_data->rpwsize = Nlocal*(mlk + storage_mu + 1);
  P_data->ipwsize = Nlocal;
  P_data->nge = 0;

  return(P_data);
}

void IBBDFree(IBBDData P_data)
{
  BandFreeMat(P_data->PP);
  BandFreePiv(P_data->pivots);
  N_VFree(P_data->tempv4);
  free(P_data);
}


/************* Preconditioner setup and solve functions **************/
 

/* Readability Replacements */

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


/******************************************************************
 * Function : IBBDPrecon                                          *
 *----------------------------------------------------------------*
 * IBBDPrecon generates a band-block-diagonal preconditioner      *
 * matrix, where the local block (on this processor) is a band    *
 * matrix.  Each local block is computed by a difference quotient *
 * scheme via calls to the user-supplied routines glocal, gcomm.  *
 * After generating the block in the band matrix PP, this routine *
 * does an LU factorization in place in PP.                       *
 *                                                                *
 * The IBBDPrecon parameters used here are as follows:            *
 *                                                                *
 * tt  is the current value of the independent variable t.        *
 *                                                                *
 * yy      is the current value of the dependent variable vector, *
 *         namely the predicted value of y(t).                    *
 *                                                                *
 * yp  is the current value of the derivative vector y',          *
 *        namely the predicted value of y'(t).                    *
 *                                                                *
 * cj  is the scalar in the system Jacobian, proportional to 1/hh.*
 *                                                                *
 * res_data  is a pointer to user data to be passed to res, the   *
 *        same as the res_data parameter passed to IDAMalloc.     *
 *                                                                *
 * P_data  is a pointer to user preconditioner data - the same as *
 *        the p_data parameter passed to IDASpgmr.                *
 *                                                                *
 * ewt    is the error weight vector.                             *
 *                                                                *
 * constraints  is the constraint vector.                         *
 *                                                                *
 * hh     is a tentative step size in t.                          *
 *                                                                *
 * tempv1, tempv2, tempv3 are pointers to vectors of type         *
 * N_Vector, used for temporary storage or work space.            *
 *                                                                *  
 * The arguments Neq, rr, res, uround, and nrePtr are not used.   *
 *                                                                *  
 * Return value:                                                  *
 * The value returned by this IBBDPrecon function is a int flag   *
 * indicating whether it was successful.  This value is           *
 *    0    if successful,                                         *
 *  > 0    for a recoverable error (step will be retried).        *
 *  < 0    for a nonrecoverable error (step fails).               *
 ******************************************************************/


int IBBDPrecon(integer Neq, real tt, N_Vector yy,
               N_Vector yp, N_Vector rr, real cj, ResFn res, 
               void *res_data, void *P_data, N_Vector ewt, 
               N_Vector constraints, real hh, real uround, long int *nrePtr,
               N_Vector tempv1, N_Vector tempv2, N_Vector tempv3)
{
  integer Nlocal, nginc, retfac;
  int retval;
  IBBDData pdata;
  N_Vector tempv4;

  pdata =(IBBDData)P_data;
  tempv4 = pdata->tempv4;
  Nlocal = N_VLOCLENGTH(yy);

  /* Call IBBDDQJac for a new Jacobian calculation and store in PP. */
  BandZero(PP);
  retval = IBBDDQJac(Nlocal, mudq, mldq, mukeep, mlkeep, cj, hh, rel_yy, tt,
                     ewt, constraints, glocal, gcomm, PP, &nginc, res_data, 
                     yy, yp, tempv1, tempv2, tempv3,tempv4);
  nge += nginc;
  if(retval<0)return(LSETUP_ERROR_NONRECVR);
  if(retval>0)return(LSETUP_ERROR_RECVR);
 
  /* Do LU factorization of preconditioner block in place (in PP). */
  retfac = BandFactor(PP, pivots);

  /* Return 0 if the LU was complete, or LSETUP_ERROR_RECVR otherwise. */
  if(retfac > 0) return(LSETUP_ERROR_RECVR);
  return(0);
}


/******************************************************************
 *                                                                *           
 * Function: IBBDPSol                                             *
 *----------------------------------------------------------------*
 * The function IBBDPSol computes a solution to the linear system *
 * P z = r, where P is the left preconditioner defined by the     *
 * routine IBBDPrecond.                                           *
 *                                                                *
 * The IBBDPSol parameters used here are as follows:              *
 *                                                                *
 * P_data is a pointer to user preconditioner data - the same as  *
 *        the p_data parameter passed to IDASpgmr.                *
 *                                                                *
 * rvec   is the input right-hand side vector r.                  *
 *                                                                *
 * zvec   is the computed solution vector z.                      *
 *                                                                *
 * The arguments Neq, tt, yy, yp, rr, cj, res, res_data, ewt,     *
 * delta, nrePtr, and tempv are NOT used.                         *
 *                                                                *
 * IBBDPSol always returns 0, indicating success.                 *
 *                                                                *
******************************************************************/

 int IBBDPSol(integer Neq, real tt, N_Vector yy, N_Vector yp, 
              N_Vector rr, real cj, ResFn res, void *res_data, 
              void *P_data, N_Vector ewt, real delta, N_Vector rvec, 
              N_Vector zvec, long int *nrePtr, N_Vector tempv)
{
  IBBDData pdata;

  pdata = (IBBDData)P_data;

  /* Copy rvec to zvec, do the backsolve, and return. */
  N_VScale(ONE, rvec, zvec);
  BandBacksolve(PP, pivots, zvec);
 
  return(0);
}

#undef mudq
#undef mukeep
#undef mldq
#undef mlkeep
#undef cj
#undef constraints
#undef ewt
#undef glocal
#undef gcomm
#undef pivots
#undef PP
#undef nge
#undef rel_yy


/*************** IBBDDQJac *****************************************

 This routine generates a banded difference quotient approximation to
 the local block of the Jacobian of G(t,y,y').  It assumes that a
 band matrix of type BandMat is stored column-wise, and that elements
 within each column are contiguous.

 All matrix elements are generated as difference quotients, by way
 of calls to the user routine glocal.
 By virtue of the band structure, the number of these calls is
 bandwidth + 1, where bandwidth = mldq + mudq + 1.
 But the band matrix kept has bandwidth = mlkeep + mukeep + 1.
 This routine also assumes that the local elements of a vector are
 stored contiguously.

 Return values are: 0 (success), > 0 (recoverable error),
 or < 0 (nonrecoverable error).
**********************************************************************/


/**********************************************************************/

static int IBBDDQJac(integer Nlocal, integer mudq, integer mldq, 
                     integer mukeep, integer mlkeep, 
                     real cj, real hh, real rel_yy, real tt,
                     N_Vector ewt, N_Vector constraints,
                     IDALocalFn glocal, IDACommFn gcomm, BandMat JJ, 
                     integer *nginc, void *res_data, N_Vector yy, N_Vector yp, 
                     N_Vector gref, N_Vector gtemp, 
                     N_Vector ytemp, N_Vector yptemp)
{
  real inc, inc_inv;
  int  retval;
  integer group, i, j, width, ngroups, i1, i2;
  real *ydata, *ypdata, *ytempdata, *yptempdata, *grefdata, *gtempdata;
  real *cnsdata, *ewtdata;
  real *col_j, conj, yj, ypj, ewtj;

  /* Obtain pointers as required to the data array of vectors. */
  ydata     = N_VDATA(yy);
  ypdata    = N_VDATA(yp);
  ytempdata = N_VDATA(ytemp);
  yptempdata= N_VDATA(yptemp);

  grefdata    = N_VDATA(gref);
  gtempdata = N_VDATA(gtemp);

  ewtdata = N_VDATA(ewt);
  if(constraints != NULL) cnsdata = N_VDATA(constraints);

  /* Initialize ytemp and yptemp. */

  N_VScale(ONE, yy, ytemp);
  N_VScale(ONE, yp, yptemp);

  /* Call gcomm and glocal to get base value of G(t,y,y'). */
  retval = gcomm(yy, yp, res_data);
  if(retval != 0) return(retval);

  retval = glocal(tt, yy, yp, gref, res_data); *nginc = 1;
  if(retval != 0) return(retval);

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
      if(hh*ypj < ZERO) inc = -inc;
      inc = (yj + inc) - yj;

      /* Adjust sign(inc) again if yj has an inequality constraint. */
      if(constraints != NULL) {
        conj = cnsdata[j];
        if(ABS(conj) == ONE)      {if((yj+inc)*conj <  ZERO) inc = -inc;}
        else if(ABS(conj) == TWO) {if((yj+inc)*conj <= ZERO) inc = -inc;}
      }

      /* Increment yj and ypj. */
      ytempdata[j] += inc;
      yptempdata[j] += cj*inc;
      
    }

    /* Evaluate G with incremented y and yp arguments. */
    retval = glocal(tt, ytemp, yptemp, gtemp, res_data); (*nginc)++;
    if(retval != 0) return(retval);

    /* Loop over components of the group again; restore ytemp and yptemp. */
    for(j = group-1; j < Nlocal; j += width) {
      yj  = ytempdata[j]  = ydata[j];
      ypj = yptempdata[j] = ypdata[j];
      ewtj = ewtdata[j];

      /* Set increment inc as before .*/
      inc = rel_yy*MAX(ABS(yj), MAX( ABS(hh*ypj), ONE/ewtj));
      if(hh*ypj < ZERO) inc = -inc;
      inc = (yj + inc) - yj;
      if(constraints != NULL) {
        conj = cnsdata[j];
        if(ABS(conj) == ONE)      {if((yj+inc)*conj <  ZERO) inc = -inc;}
        else if(ABS(conj) == TWO) {if((yj+inc)*conj <= ZERO) inc = -inc;}
      }

      /* Form difference quotients and load into JJ. */
      inc_inv = ONE/inc;
      col_j = BAND_COL(JJ,j);
      i1 = MAX(0, j-mukeep);
      i2 = MIN(j+mlkeep, Nlocal-1);
      for(i = i1; i <= i2; i++) BAND_COL_ELEM(col_j,i,j) =
	                inc_inv * (gtempdata[i] - grefdata[i]);       
    }
  }

 return(SUCCESS);
}
