/******************************************************************
 * File          : pvsbbdpre.c                                    *
 * Programmers   : Michael Wittman and Alan C. Hindmarsh @ LLNL   *
 * Version of    : 14 January 2002                                *
 *----------------------------------------------------------------*
 * This file contains implementations of routines for a           *
 * band-block-diagonal preconditioner, i.e. a block-diagonal      *
 * matrix with banded blocks, for use with PVODE and CVSpgmr.     *
 ******************************************************************/

#include "pvsbbdpre.h"
#include "cvodes.h"
#include "llnltyps.h"
#include "nvector.h"
#include "llnlmath.h"
#include "iterativ.h"
#include "band.h"

#define MIN_INC_MULT RCONST(1000.0)
#define ZERO         RCONST(0.0)
#define ONE          RCONST(1.0)


/* Prototype for difference quotient Jacobian calculation routine */

static void PVBBDDQJac(integer Nlocal, integer mudq, integer mldq, 
                       integer mukeep, integer mlkeep, real rely,
                       PVLocalFn gloc, PVCommFn cfn, BandMat J, void *f_data,
                       real t, N_Vector y, N_Vector ewt, real h, real uround,
                       N_Vector gy, N_Vector gtemp, N_Vector ytemp);


/***************** User-Callable Functions: malloc and free ******************/

PVBBDData PVBBDAlloc(integer Nlocal, integer mudq, integer mldq,
                     integer mukeep, integer mlkeep, real dqrely, 
                     PVLocalFn gloc, PVCommFn cfn, void *f_data)
{
  PVBBDData pdata;
  integer muk, mlk, storage_mu;

  pdata = (PVBBDData) malloc(sizeof *pdata);  /* Allocate data memory */
  if (pdata == NULL) return(NULL);

  /* Set pointers to f_data, gloc, and cfn; load half-bandwidths */
  pdata->f_data = f_data;
  pdata->gloc = gloc;
  pdata->cfn = cfn;
  pdata->mudq = MIN( Nlocal-1, MAX(0,mudq) );
  pdata->mldq = MIN( Nlocal-1, MAX(0,mldq) );
  muk = MIN( Nlocal-1, MAX(0,mukeep) );
  mlk = MIN( Nlocal-1, MAX(0,mlkeep) );
  pdata->mukeep = muk;
  pdata->mlkeep = mlk;

  /* Allocate memory for saved Jacobian */
  pdata->savedJ = BandAllocMat(Nlocal, muk, mlk, muk);
  if (pdata->savedJ == NULL) { free(pdata); return(NULL); }

  /* Allocate memory for preconditioner matrix */
  storage_mu = MIN(Nlocal-1, muk + mlk);
  pdata->savedP = BandAllocMat(Nlocal, muk, mlk, storage_mu);
  if (pdata->savedP == NULL) {
    BandFreeMat(pdata->savedJ);
    free(pdata);
    return(NULL);
  }
  /* Allocate memory for pivots */
  pdata->pivots = BandAllocPiv(Nlocal);
  if (pdata->savedJ == NULL) {
    BandFreeMat(pdata->savedP);
    BandFreeMat(pdata->savedJ);
    free(pdata);
    return(NULL);
  }

  /* Set pdata->dqrely based on input dqrely (0 implies default). */
  pdata->dqrely = (dqrely > ZERO) ? dqrely : RSqrt(UnitRoundoff());

  /* Set work space sizes and initialize nge */
  pdata->rpwsize = Nlocal*(muk + 2*mlk + storage_mu + 2);
  pdata->ipwsize = Nlocal;
  pdata->nge = 0;

  return(pdata);
}

void PVBBDFree(PVBBDData pdata)
{
  BandFreeMat(pdata->savedJ);
  BandFreeMat(pdata->savedP);
  BandFreePiv(pdata->pivots);
  free(pdata);
}


/***************** Preconditioner setup and solve Functions ****************/
 

/* Readability Replacements */

#define f_data    (pdata->f_data)
#define mudq      (pdata->mudq)
#define mldq      (pdata->mldq)
#define mukeep    (pdata->mukeep)
#define mlkeep    (pdata->mlkeep)
#define dqrely    (pdata->dqrely)
#define gloc      (pdata->gloc)
#define cfn       (pdata->cfn)
#define savedJ    (pdata->savedJ)
#define savedP    (pdata->savedP)
#define pivots    (pdata->pivots)
#define nge       (pdata->nge)


/******************************************************************
 * Function : PVBBDPrecon                                         *
 *----------------------------------------------------------------*
 * PVBBDPrecon generates and factors a banded block of the        *
 * preconditioner matrix on each processor, via calls to the      *
 * user-supplied gloc and cfn functions. It uses difference       *
 * quotient approximations to the Jacobian elements.              *
 *                                                                *
 * PVBBDPrecon calculates a new J, if necessary, then calculates  *
 * P = I - gamma*J, and does an LU factorization of P.            *
 *                                                                *
 * The parameters of PVBBDPrecon used here are as follows:        *
 *                                                                *
 * Neq     is the system size (global vector length).             *
 *                                                                *
 * t       is the current value of the independent variable.      *
 *                                                                *
 * y       is the current value of the dependent variable vector, *
 *         namely the predicted value of y(t).                    *
 *                                                                *
 * fy      is the vector f(t,y).                                  *
 *                                                                *
 * jok     is an input flag indicating whether Jacobian-related   *
 *         data needs to be recomputed, as follows:               *
 *           jok == FALSE means recompute Jacobian-related data   *
 *                  from scratch.                                 *
 *           jok == TRUE  means that Jacobian data from the       *
 *                  previous PVBBDPrecon call can be reused       *
 *                  (with the current value of gamma).            *
 *         A PVBBDPrecon call with jok == TRUE should only occur  *
 *         after a call with jok == FALSE.                        *
 *                                                                *
 * jcurPtr is a pointer to an output integer flag which is        *
 *         set by PVBBDPrecon as follows:                         *
 *           *jcurPtr = TRUE if Jacobian data was recomputed.     *
 *           *jcurPtr = FALSE if Jacobian data was not            * 
 *                      recomputed, but saved data was reused.    *
 *                                                                *
 * gamma   is the scalar appearing in the Newton matrix.          *
 *                                                                *
 * ewt     is the error weight vector.                            *
 *                                                                *
 * h       is a tentative step size in t.                         *
 *                                                                *
 * uround  is the machine unit roundoff.                          *
 *                                                                *
 * nfePtr  is a pointer to the memory location containing the     *
 *           CVODE counter nfe = number of calls to f. Not used.  *
 *                                                                *
 * P_data  is a pointer to user data - the same as the P_data     *
 *           parameter passed to CVSpgmr. For PVBBDPrecon, this   *
 *           should be of type PVBBDData.                         *
 *                                                                *
 * vtemp1, vtemp2, and vtemp3 are pointers to memory allocated    *
 *           for vectors of length N which are be used by         *
 *           PVBBDPrecon as temporary storage or work space.      *
 *                                                                *
 *                                                                *
 * Return value:                                                  *
 * The value returned by this PVBBDPrecon function is the int     *
 *   0  if successful,                                            *
 *   1  for a recoverable error (step will be retried).           *
 ******************************************************************/

int PVBBDPrecon(integer Neq, real t, N_Vector y, N_Vector fy,
                boole jok, boole *jcurPtr, real gamma, N_Vector ewt,
                real h, real uround, long int *nfePtr, void *P_data,
                N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3)
{
  integer Nlocal, ier;
  PVBBDData pdata;

  pdata = (PVBBDData)P_data;

  Nlocal = N_VLOCLENGTH(y);

  if (jok) {
    /* If jok = TRUE, use saved copy of J */
    *jcurPtr = FALSE;
    BandCopy(savedJ, savedP, mukeep, mlkeep);
  } else {
    /* Otherwise call PVBBDDQJac for new J value */
    *jcurPtr = TRUE;
    BandZero(savedJ);
    PVBBDDQJac(Nlocal, mudq, mldq, mukeep, mlkeep, dqrely, gloc, cfn, savedJ,
               f_data, t, y, ewt, h, uround, vtemp1, vtemp2, vtemp3);
    nge += 1 + MIN(mldq + mudq + 1, Nlocal);
    BandCopy(savedJ, savedP, mukeep, mlkeep);
  }
  
  /* Scale and add I to get P = I - gamma*J */
  BandScale(-gamma, savedP);
  BandAddI(savedP);
 
  /* Do LU factorization of P in place */
  ier = BandFactor(savedP, pivots);
 
  /* Return 0 if the LU was complete; otherwise return 1 */
  if (ier > 0) return(1);
  return(0);
}


/******************************************************************
 * Function : PVBBDPSol                                           *
 *----------------------------------------------------------------*
 * PVBBDPSol solves a linear system P z = r, with the band-block- *
 * diagonal preconditioner matrix P generated and factored by     *
 * PVBBDPrecon.                                                   *
 *                                                                *
 * The parameters of PVBBDPSol used here are as follows:          *
 *                                                                *
 * r      is the right-hand side vector of the linear system.     *
 *                                                                *
 * P_data is a pointer to the preconditioner data returned by     *
 *        PVBBDAlloc.                                             *
 *                                                                *
 * z      is the output vector computed by PVBBDPSol.             *
 *                                                                *
 * The value returned by the PVBBDPSol function is always 0,      *
 * indicating success.                                            *
 ******************************************************************/

int PVBBDPSol(integer N, real t, N_Vector y, N_Vector fy,
              N_Vector vtemp, real gamma, N_Vector ewt,
              real delta, long int *nfePtr, N_Vector r,
              int lr, void *P_data, N_Vector z)
{
  PVBBDData pdata;

  pdata = (PVBBDData)P_data;

  /* Copy r to z, then do backsolve and return */
  N_VScale(ONE, r, z);
  BandBacksolve(savedP, pivots, z);
 
  return(0);
}

#undef f_data
#undef mudq
#undef mldq
#undef mukeep
#undef mlkeep
#undef dqrely
#undef gloc
#undef cfn
#undef savedJ
#undef savedP
#undef pivots
#undef nge


/*************** PVBBDDQJac *****************************************

 This routine generates a banded difference quotient approximation to
 the local block of the Jacobian of g(t,y).  It assumes that a band 
 matrix of type BandMat is stored columnwise, and that elements within
 each column are contiguous.  All matrix elements are generated as
 difference quotients, by way of calls to the user routine gloc.
 By virtue of the band structure, the number of these calls is
 bandwidth + 1, where bandwidth = mldq + mudq + 1.
 But the band matrix kept has bandwidth = mlkeep + mukeep + 1.
 This routine also assumes that the local elements of a vector are
 stored contiguously.
**********************************************************************/

static void PVBBDDQJac(integer Nlocal, integer mudq, integer mldq, 
                       integer mukeep, integer mlkeep, real rely,
                       PVLocalFn gloc, PVCommFn cfn, BandMat J,
                       void *f_data, real t, N_Vector y,
                       N_Vector ewt, real h, real uround,
                       N_Vector gy, N_Vector gtemp, N_Vector ytemp)
{
  real    gnorm, minInc, inc, inc_inv;
  integer group, i, j, width, ngroups, i1, i2;
  real *y_data, *ewt_data, *gy_data, *gtemp_data, *ytemp_data, *col_j;

  /* Obtain pointers to the data for all vectors */
  y_data     = N_VDATA(y);
  ewt_data   = N_VDATA(ewt);
  gy_data    = N_VDATA(gy);
  gtemp_data = N_VDATA(gtemp);
  ytemp_data = N_VDATA(ytemp);

  /* Load ytemp with y = predicted solution vector */
  N_VScale(ONE, y, ytemp);

  /* Call cfn and gloc to get base value of g(t,y) */
  cfn (Nlocal, t, y, f_data);
  gloc(Nlocal, t, ytemp_data, gy_data, f_data);

  /* Set minimum increment based on uround and norm of g */
  gnorm = N_VWrmsNorm(gy, ewt);
  minInc = (gnorm != ZERO) ?
           (MIN_INC_MULT * ABS(h) * uround * Nlocal * gnorm) : ONE;

  /* Set bandwidth and number of column groups for band differencing */
  width = mldq + mudq + 1;
  ngroups = MIN(width, Nlocal);

  /* Loop over groups */  
  for (group=1; group <= ngroups; group++) {
    
    /* Increment all y_j in group */
    for(j=group-1; j < Nlocal; j+=width) {
      inc = MAX(rely*ABS(y_data[j]), minInc/ewt_data[j]);
      ytemp_data[j] += inc;
    }

    /* Evaluate g with incremented y */
    gloc(Nlocal, t, ytemp_data, gtemp_data, f_data);

    /* Restore ytemp, then form and load difference quotients */
    for (j=group-1; j < Nlocal; j+=width) {
      ytemp_data[j] = y_data[j];
      col_j = BAND_COL(J,j);
      inc = MAX(rely*ABS(y_data[j]), minInc/ewt_data[j]);
      inc_inv = ONE/inc;
      i1 = MAX(0, j-mukeep);
      i2 = MIN(j+mlkeep, Nlocal-1);
      for (i=i1; i <= i2; i++)
        BAND_COL_ELEM(col_j,i,j) =
          inc_inv * (gtemp_data[i] - gy_data[i]);
    }
  }
}
