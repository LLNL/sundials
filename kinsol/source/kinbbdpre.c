/******************************************************************
 *                                                                *
 * File          : kinbbdpre.c                                    *
 * Programmers   : Allan G Taylor and Alan C Hindmarsh @ LLNL     *
 * Version of    : 18 January 2001                                *
 *----------------------------------------------------------------*
 * This file contains implementations of routines for a           *
 * band-block-diagonal preconditioner, i.e. a block-diagonal      *
 * matrix with banded blocks, for use with KINSol and KINSpgmr.   *
 * NOTE: with only one processor in use, a banded matrix results  *
 * rather than a b-b matrix with banded blocks. Diagonal blocking *
 * occurs at the processor level.                                 *
 ******************************************************************/

#include "kinbbdpre.h"
#include "kinsol.h"
#include "llnltyps.h"
#include "nvector.h"
#include "llnlmath.h"
#include "iterativ.h"
#include "band.h"

#define ZERO         RCONST(0.0)
#define ONE          RCONST(1.0)


/* Prototype for difference quotient Jacobian calculation routine */

static void KBBDDQJac(integer Nlocal, BandMat J, void *P_data,
                      N_Vector uu, N_Vector uscale, N_Vector gu,
                      N_Vector gtemp, N_Vector utemp);


/***************** User-Callable Functions: malloc and free ******************/

KBBDData KBBDAlloc(integer Nlocal, integer mu, integer ml,
		   real dq_rel_uu, KINLocalFn gloc, KINCommFn gcomm, 
		   void *kinmem, void *f_data, machEnvType machEnv)
{
  KBBDData pdata;
  KINMem kin_mem;
  N_Vector vtemp3;
  real rel_uu;

  kin_mem = (KINMem)kinmem;

  pdata = (KBBDData) malloc(sizeof *pdata);  /* Allocate data memory */
  if (pdata == NULL)
    return(NULL);

  /* Set pointers to f_data, gloc, and gcomm; load half-bandwiths */
  pdata->f_data = f_data;
  pdata->ml = ml;
  pdata->mu = mu;
  pdata->gloc = gloc;
  pdata->gcomm = gcomm;
 
  /* Allocate memory for preconditioner matrix */
  pdata->PP = BandAllocMat(Nlocal, mu, ml, mu+ml);
  if (pdata->PP == NULL) {
    free(pdata);
    return(NULL);
  }

  /* Allocate memory for pivots */
  pdata->pivots = BandAllocPiv(Nlocal);
  if (pdata->PP == NULL) {
    BandFreeMat(pdata->PP);
    free(pdata);
    return(NULL);
  }

  /* allocate vtemp3 for use by KBBDDQJac  */

  vtemp3 = N_VNew(Nlocal, machEnv); /* Note: Nlocal here is a dummy in that 
                            machEnv parameters are used to determine size */
  if (vtemp3 == NULL) {
      free(pdata);
      return(NULL);
  }
  pdata->vtemp3 = vtemp3;

  /* set rel_uu based on input value dq_rel_uu */

  if (dq_rel_uu > ZERO) rel_uu=dq_rel_uu;
  else  rel_uu=RSqrt(kin_mem->kin_uround); /* a dq_rel_uu value of 0.0 received
				 	      by this routine implies using
					      the default instead */

  pdata->rel_uu = rel_uu;

  /* Set work space sizes and initialize nge */
  pdata->rpwsize = Nlocal*(2*mu + ml + 1);
  pdata->ipwsize = Nlocal;
  pdata->nge = 0;
  return(pdata);
}

void KBBDFree(KBBDData pdata)
{
  BandFreeMat(pdata->PP);
  BandFreePiv(pdata->pivots);
  free(pdata);
}


/***************** Preconditioner setup and solve Functions ****************/
 

/* Readability Replacements */

#define f_data    (pdata->f_data)
#define mu        (pdata->mu)
#define ml        (pdata->ml)
#define gloc      (pdata->gloc)
#define gcomm     (pdata->gcomm)
#define pivots    (pdata->pivots)
#define PP        (pdata->PP)
#define nge       (pdata->nge)
#define rel_uu    (pdata->rel_uu)

/******************************************************************
 * Function : KBBDPrecon                                          *
 *----------------------------------------------------------------*
 * KBBDPrecon generates and factors a banded block of the         *
 * preconditioner matrix on each processor, via calls to the      *
 * user-supplied gloc and gcomm functions. It uses difference     *
 * quotient approximations to the Jacobian elements.              *
 *                                                                *
 * KBBDPrecon calculates a new Jacobian, stored in banded         *
 * matrix PP and does an LU factorization of P in place  in PP.   *
 *                                                                *
 * The parameters of KBBDPrecon are as follows:                   *
 *                                                                *
 * Neq   is the system size and the global length of all vectors. *
 *           (not used)                                           *
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
 *           KBBDPrecon as temporary storage or work space.       *
 *           A third vector vtemp3 required for KBBDPrecon has    *
 *           been previously allocated as pdata->vtemp3           *
 *                                                                *
 * func     SysFn type function defining the system f(u)=0.       *
 *                                                                *
 * uround   real variable giving the unit roundoff  (not used)    *
 *                                                                *
 * nfePtr  is a pointer to the memory location containing the     *
 *           KINSOL counter nfe = number of calls to f. Not used. *
 *                                                                *
 * P_data  is a pointer to user data - the same as the P_data     *
 *           parameter passed to KINSpgmr. For KBBDPrecon, this   *
 *           should be of type KBBDData.                          *
 *                                                                *
 *                                                                *
 * Return value:                                                  *
 * The value to be returned by the KBBDPrecon function is a flag  *
 * indicating whether it was successful.  This value is           *
 *   0   if successful,                                           *
 *   > 0 for a recoverable error (step will be retried).          *
 ******************************************************************/

int KBBDPrecon(integer Neq, N_Vector uu, N_Vector uscale,
	        N_Vector fval, N_Vector fscale,
                N_Vector vtemp1, N_Vector vtemp2,
	        SysFn func, real uround,
                long int *nfePtr, void *P_data)

{
  integer Nlocal, ier;
  KBBDData pdata;
  N_Vector vtemp3;

  pdata = (KBBDData)P_data;
  vtemp3 = pdata->vtemp3;

  Nlocal = N_VLOCLENGTH(uu);

  /* Otherwise call KBBDDQJac for a new Jacobian calc and store in PP */
  BandZero(PP);
  KBBDDQJac(Nlocal, PP, P_data, uu, uscale, vtemp1, vtemp2, vtemp3);
  nge += 1 + MIN(ml + mu + 1, Nlocal);
 
  /* Do LU factorization of P in place  (in PP) */
  ier = BandFactor(PP, pivots);
 
  /* Return 0 if the LU was complete; otherwise return 1 */
  if (ier > 0) return(1);
  return(0);
}


/******************************************************************
 * Function : KBBDPSol                                            *
 *----------------------------------------------------------------*
 * KBBDPSol solves a linear system P z = r, with the banded       *
 * blocked preconditioner matrix P generated and factored by      *
 * KBBDPrecon.  Here, r comes in as vtem and z is returned in vtem*
 * as well.                                                       *
 *                                                                *
 * The parameters of KBBDPSol are as follows:                     *
 *                                                                *
 * Neq    is the global length of all vector arguments.           *
 *                  (not used)                                    *
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
 * func   the function func defines the system being solved:      *     
 *             func(uu) = 0    or  f(u)=0                         *
 *                                                                *
 * uround  is the machine unit roundoff  (not used)               *
 *                                                                *
 *                                                                *
 * nfePtr is a pointer to the memory location containing the      *
 *          KINSOL problem data nfe = number of calls to f. The   *
 *          KBBDPSol routine does not use this argument.          *
 *                                                                *
 * P_data is a pointer to user data - the same as the P_data      *
 *          parameter passed to KINSpgmr. For KBBDPSol, this      *
 *          should be of type KBBDData.                           *
 *                                                                *
 * Return value:                                                  *
 * The value returned by the KBBDPSol function is a flag          *
 * indicating whether it was successful.  Here this value is      *
 * always 0, indicating success.                                  *
 ******************************************************************/

int KBBDPSol(integer Nlocal, N_Vector uu, N_Vector uscale,
	      N_Vector fval, N_Vector fscale,
	      N_Vector vtem, N_Vector ftem,
	      SysFn func, real u_round,
	      long int *nfePtr, void *P_data)
{
  KBBDData pdata;

  pdata = (KBBDData)P_data;

  /* Do the backsolve and return */
  BandBacksolve(PP, pivots, vtem);
 
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

static void KBBDDQJac(integer Nlocal, BandMat J, void *P_data,
		       N_Vector uu, N_Vector uscale, N_Vector gu,
		       N_Vector gtemp, N_Vector utemp)
{
  real inc, inc_inv;
  integer group, i, j, width, ngroups, i1, i2;
  real *udata, *uscdata, *gudata, *gtempdata, *utempdata, *col_j;
  KBBDData pdata;

  pdata= (KBBDData)P_data;

  /* Set pointers to the data for all vectors */
  udata     = N_VDATA(uu);
  uscdata   = N_VDATA(uscale);
  gudata    = N_VDATA(gu);
  gtempdata = N_VDATA(gtemp);
  utempdata = N_VDATA(utemp);

  /* Load utemp with uu = predicted solution vector */
  N_VScale(ONE, uu, utemp);

  /* Call gcomm and gloc to get base value of g(uu) */
  gcomm (Nlocal, udata, f_data);
  gloc (Nlocal, uu, gu, f_data);


  /* Set bandwidth and number of column groups for band differencing */
  width = ml + mu + 1;
  ngroups = MIN(width, Nlocal);

  /* Loop over groups */  
  for (group=1; group <= ngroups; group++) {
    
    /* Increment all u_j in group */
    for(j=group-1; j < Nlocal; j+=width) {
      inc = rel_uu * MAX(ABS(udata[j]), ONE/uscdata[j]);
      utempdata[j] += inc;
    }

    /* Evaluate g with incremented u */
    gloc (Nlocal, utemp, gtemp, f_data);

    /* Restore utemp, then form and load difference quotients */
    for (j=group-1; j < Nlocal; j+=width) {
      utempdata[j] = udata[j];
      col_j = BAND_COL(J,j);
      inc =rel_uu * MAX( ABS(udata[j]) , ONE/uscdata[j] );
      inc_inv = ONE/inc;
      i1 = MAX(0, j-mu);
      i2 = MIN(j+ml, Nlocal-1);
      for (i=i1; i <= i2; i++)
	BAND_COL_ELEM(col_j,i,j) =
	  inc_inv * (gtempdata[i] - gudata[i]);
    }
  }
}
