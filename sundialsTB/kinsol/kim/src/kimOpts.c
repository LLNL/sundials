/*
 * -----------------------------------------------------------------
 * $Revision$
 * $Date$
 * -----------------------------------------------------------------
 * Programmer: Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2021, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * Option parsing functions for the KINSOL Matlab interface.
 * -----------------------------------------------------------------
 */

#include <string.h>
#include "kim.h"

/*
 * ---------------------------------------------------------------------------------
 * Global interface data variable (defined in kim.c)
 * ---------------------------------------------------------------------------------
 */

extern kimInterfaceData kimData;

/*
 * ---------------------------------------------------------------------------------
 * Redability replacements
 * ---------------------------------------------------------------------------------
 */

#define N            (kimData->n)
#define ls           (kimData->LS)
#define pm           (kimData->PM)

#define mtlb_data    (kimData->mtlb_data)

#define mtlb_JACfct  (kimData->JACfct)
#define mtlb_PSETfct (kimData->PSETfct)
#define mtlb_PSOLfct (kimData->PSOLfct)
#define mtlb_GLOCfct (kimData->GLOCfct)
#define mtlb_GCOMfct (kimData->GCOMfct)

/*
 * ---------------------------------------------------------------------------------
 * Option handling functions
 * ---------------------------------------------------------------------------------
 */

int get_SolverOptions(const mxArray *options,
                      booleantype *verbose, booleantype *errmsg,
                      int *mxiter, int *msbset, int *msbsetsub, 
                      int *etachoice, int *mxnbcf,
                      double *eta, double *egamma, double *ealpha, double *mxnewtstep, 
                      double *relfunc, double *fnormtol, double *scsteptol,
                      double **constraints,
                      booleantype *noInitSetup, booleantype *noMinEps)
{
  mxArray *opt;
  char *bufval;
  int buflen, status;
  sunindextype i, m, n;
  double *tmp;

  /* Set default values (pass 0 values. KINSOL does the rest) */

  *mxiter = 0;
  *msbset = 0;
  *msbsetsub = 0;
  *mxnbcf = 0;
  *etachoice = KIN_ETACHOICE1;

  *eta = 0.0;
  *egamma = 0.0;
  *ealpha = 0.0;
  *mxnewtstep = 0.0;
  *relfunc = 0.0;
  *fnormtol = 0.0;
  *scsteptol = 0.0;

  *noInitSetup = SUNFALSE;
  *noMinEps = SUNFALSE;

  *constraints = NULL;

  *verbose = SUNFALSE;
  *errmsg = SUNTRUE;

  /* Return now if options was empty */

  if (mxIsEmpty(options)) return(0);

  /* User data */

  opt = mxGetField(options,0,"UserData");
  if ( !mxIsEmpty(opt) ) {
    mxDestroyArray(mtlb_data);
    mtlb_data = mxDuplicateArray(opt);
  }

  /* Integer values */

  opt = mxGetField(options,0,"MaxNumIter");
  if ( !mxIsEmpty(opt) )
    *mxiter = (int)*mxGetPr(opt);
  
  opt = mxGetField(options,0,"MaxNumSetups");
  if ( !mxIsEmpty(opt) )
    *msbset = (int)*mxGetPr(opt);
  
  opt = mxGetField(options,0,"MaxNumSubSetups");
  if ( !mxIsEmpty(opt) )
    *msbsetsub = (int)*mxGetPr(opt);
  
  opt = mxGetField(options,0,"MaxNumBetaFails");
  if ( !mxIsEmpty(opt) )
    *mxnbcf = (int)*mxGetPr(opt);
  
  opt = mxGetField(options,0,"EtaForm");
  if ( !mxIsEmpty(opt) ) {
    buflen = mxGetM(opt) * mxGetN(opt) + 1;
    bufval = mxCalloc(buflen, sizeof(char));
    status = mxGetString(opt, bufval, buflen);
    if(status != 0) {
      kimErrHandler(-999, "KINSOL", "KINInit", "Cannot parse EtaForm.", NULL);
      return(-1);
    }
    if(!strcmp(bufval,"Type1"))         *etachoice = KIN_ETACHOICE1;
    else if(!strcmp(bufval,"Type2"))    *etachoice = KIN_ETACHOICE2;
    else if(!strcmp(bufval,"Constant")) *etachoice = KIN_ETACONSTANT;
    else {
      kimErrHandler(-999, "KINSOL", "KINInit", "EtaForm has an illegal value.", NULL);
      return(-1);
    }
  }
  
  /* Real values */

  opt = mxGetField(options,0,"Eta");
  if ( !mxIsEmpty(opt) )
    *eta = (double)*mxGetPr(opt);
  
  opt = mxGetField(options,0,"EtaAlpha");
  if ( !mxIsEmpty(opt) )
    *ealpha = (double)*mxGetPr(opt);
  
  opt = mxGetField(options,0,"EtaGamma");
  if ( !mxIsEmpty(opt) )
    *egamma = (double)*mxGetPr(opt);
  
  opt = mxGetField(options,0,"MaxNewtonStep");
  if ( !mxIsEmpty(opt) )
    *mxnewtstep = (double)*mxGetPr(opt);
  
  opt = mxGetField(options,0,"FuncRelErr");
  if ( !mxIsEmpty(opt) )
    *relfunc = (double)*mxGetPr(opt);

  opt = mxGetField(options,0,"FuncNormTol");
  if ( !mxIsEmpty(opt) )
    *fnormtol = (double)*mxGetPr(opt);

  opt = mxGetField(options,0,"ScaledStepTol");
  if ( !mxIsEmpty(opt) )
    *scsteptol = (double)*mxGetPr(opt);

  /* Boolean values */

  opt = mxGetField(options,0,"ErrorMessages");
  if ( !mxIsEmpty(opt) ) {
    if (!mxIsLogical(opt)) {
      kimErrHandler(-999, "KINSOL", "KINInit", "ErrorMessages is not a logical scalar.", NULL);
      return(-1);
    }
    if (mxIsLogicalScalarTrue(opt)) *errmsg = SUNTRUE;
    else                            *errmsg = SUNFALSE;
  }

  opt = mxGetField(options,0,"Verbose");
  if ( !mxIsEmpty(opt) ) {
    if (!mxIsLogical(opt)) {
      kimErrHandler(-999, "KINSOL", "KINInit", "Verbose is not a logical scalar.", NULL);
      return(-1);
    }
    if (mxIsLogicalScalarTrue(opt)) *verbose = SUNTRUE;
    else                            *verbose = SUNFALSE;
  }

  opt = mxGetField(options,0,"InitialSetup");
  if ( !mxIsEmpty(opt) ) {
    if (!mxIsLogical(opt)) {
      kimErrHandler(-999, "KINSOL", "KINInit", "InitialSetup is not a logical scalar.", NULL);
      return(-1);
    }
    if (mxIsLogicalScalarTrue(opt)) *noInitSetup = SUNFALSE;
    else                            *noInitSetup = SUNTRUE;
  }


  opt = mxGetField(options,0,"MinBoundEps");
  if ( !mxIsEmpty(opt) ) {
    if (!mxIsLogical(opt)) {
      kimErrHandler(-999, "KINSOL", "KINInit", "MinBoundEps is not a logical scalar.", NULL);
      return(-1);
    }
    if (mxIsLogicalScalarTrue(opt)) *noMinEps = SUNFALSE;
    else                            *noMinEps = SUNTRUE;
  }

  /* Constraints */

  opt = mxGetField(options,0,"Constraints");
  if ( !mxIsEmpty(opt) ) {
    m = mxGetM(opt);
    n = mxGetN(opt);
    if ( (n != 1) && (m != 1) ) {
      kimErrHandler(-999, "KINSOL", "KINInit", "constraints is not a vector.", NULL);
      return(-1);
    }
    if ( m > n ) n = m;
    if ( n != N ) {
      kimErrHandler(-999, "KINSOL", "KINInit", "constraints has wrong number of components.", NULL);
      return(-1);
    }
    tmp = mxGetPr(opt);
    *constraints = (double *) malloc(N*sizeof(double));
    for (i=0;i<N;i++) (*constraints)[i] = tmp[i];
  }

  /* We made it here without problems */

  return(0);
}


int get_LinSolvOptions(const mxArray *options,
                       sunindextype *mupper, sunindextype *mlower,
                       sunindextype *mudq, sunindextype *mldq, double *dqrely,
                       int *ptype, int *maxrs, int *maxl)
{
  mxArray *opt;
  char *bufval;
  int buflen, status;
  
  *mupper = 0;
  *mlower = 0;

  *mudq = 0;
  *mldq = 0;
  *dqrely = 0.0;

  *ptype = PREC_NONE;

  *maxl  = 0;
  *maxrs = 0;

  ls = LS_DENSE;
  pm = PM_NONE;

  /* Return now if options was empty */

  if (mxIsEmpty(options)) return(0);

  /* Linear solver type */

  opt = mxGetField(options,0,"LinearSolver");
  if ( !mxIsEmpty(opt) ) {
    buflen = mxGetM(opt) * mxGetN(opt) + 1;
    bufval = mxCalloc(buflen, sizeof(char));
    status = mxGetString(opt, bufval, buflen);
    if(status != 0) {
      kimErrHandler(-999, "KINSOL", "KINInit", "Cannot parse LinearSolver.", NULL);
      return(-1);
    }
    if(!strcmp(bufval,"Band"))          ls = LS_BAND;
    else if(!strcmp(bufval,"GMRES"))    ls = LS_SPGMR;
    else if(!strcmp(bufval,"BiCGStab")) ls = LS_SPBCG;
    else if(!strcmp(bufval,"TFQMR"))    ls = LS_SPTFQMR;
    else if(!strcmp(bufval,"Dense"))    ls = LS_DENSE;
    else {
      kimErrHandler(-999, "KINSOL", "KINInit", "LinearSolver has an illegal value.", NULL);
      return(-1);
    }
  }
  
  /* Jacobian function */

  opt = mxGetField(options,0,"JacobianFn");
  if ( !mxIsEmpty(opt) ) {
    mxDestroyArray(mtlb_JACfct);
    mtlb_JACfct  = mxDuplicateArray(opt);
  }

  /* Band linear solver */

  if (ls==LS_BAND) {

    opt = mxGetField(options,0,"UpperBwidth");
    if ( !mxIsEmpty(opt) )
      *mupper = (sunindextype)*mxGetPr(opt);
    
    opt = mxGetField(options,0,"LowerBwidth");
    if ( !mxIsEmpty(opt) )
      *mlower = (sunindextype)*mxGetPr(opt);

  }
  
  /*  SPGMR linear solver options */

  if (ls==LS_SPGMR) {

    opt = mxGetField(options,0,"MaxNumRestarts");
    if ( !mxIsEmpty(opt) )
      *maxrs = (int)*mxGetPr(opt);  

  }

  /* SPILS linear solver options */

  if ( (ls==LS_SPGMR) || (ls==LS_SPBCG) || (ls==LS_SPTFQMR) ) {

    /* Max. dimension of Krylov subspace */

    opt = mxGetField(options,0,"KrylovMaxDim");
    if ( !mxIsEmpty(opt) ) {
      *maxl = (int)*mxGetPr(opt);
      if (*maxl < 0) {
        kimErrHandler(-999, "KINSOL", "KINInit", "KrylovMaxDim is negative.", NULL);
        return(-1);
      }
    }

    /* Preconditioning type */

    opt = mxGetField(options,0,"PrecType");
    if ( !mxIsEmpty(opt) ) {
      buflen = mxGetM(opt) * mxGetN(opt) + 1;
      bufval = mxCalloc(buflen, sizeof(char));
      status = mxGetString(opt, bufval, buflen);
      if(status != 0) {
        kimErrHandler(-999, "KINSOL", "KINInit", "Cannot parse PrecType.", NULL);
        return(-1);
      }
      if(!strcmp(bufval,"Right"))     *ptype = PREC_RIGHT;
      else if(!strcmp(bufval,"None")) *ptype = PREC_NONE;
      else {
        kimErrHandler(-999, "KINSOL", "KINInit", "PrecType has an illegal value.", NULL);
        return(-1);
      }
    }

    /* User defined precoditioning */

    opt = mxGetField(options,0,"PrecSetupFn");
    if ( !mxIsEmpty(opt) ) {
      mxDestroyArray(mtlb_PSETfct);
      mtlb_PSETfct  = mxDuplicateArray(opt);
    }
  
    opt = mxGetField(options,0,"PrecSolveFn");
    if ( !mxIsEmpty(opt) ) {
      mxDestroyArray(mtlb_PSOLfct);
      mtlb_PSOLfct  = mxDuplicateArray(opt);
    }
    
    /* Preconditioner module */
  
    opt = mxGetField(options,0,"PrecModule");
    if ( !mxIsEmpty(opt) ) {
      buflen = mxGetM(opt) * mxGetN(opt) + 1;
      bufval = mxCalloc(buflen, sizeof(char));
      status = mxGetString(opt, bufval, buflen);
      if(status != 0) {
        kimErrHandler(-999, "KINSOL", "KINInit", "Cannot parse PrecModule.", NULL);
        return(-1);
      }
      if(!strcmp(bufval,"BBDPre"))           pm = PM_BBDPRE;
      else if(!strcmp(bufval,"UserDefined")) pm = PM_NONE;
      else {
        kimErrHandler(-999, "KINSOL", "KINInit", "PrecModule has an illegal value.", NULL);
        return(-1);
      }
    }

    
    if (pm == PM_BBDPRE) {

      opt = mxGetField(options,0,"UpperBwidth");
      if ( !mxIsEmpty(opt) )
        *mupper = (sunindextype)*mxGetPr(opt);
    
      opt = mxGetField(options,0,"LowerBwidth");
      if ( !mxIsEmpty(opt) )
        *mlower = (sunindextype)*mxGetPr(opt);

      opt = mxGetField(options,0,"UpperBwidthDQ");
      if ( !mxIsEmpty(opt) )
        *mudq = (sunindextype)*mxGetPr(opt);

      opt = mxGetField(options,0,"LowerBwidthDQ");
      if ( !mxIsEmpty(opt) )
        *mldq = (sunindextype)*mxGetPr(opt);
      
      opt = mxGetField(options,0,"GlocalFn");
      if ( !mxIsEmpty(opt) ) {
        mxDestroyArray(mtlb_GLOCfct);
        mtlb_GLOCfct  = mxDuplicateArray(opt);
      }
      else {
        kimErrHandler(-999, "KINSOL", "KINInit", "GlocalFn required for BBD preconditioner.", NULL);
        return(-1);
      }      

      opt = mxGetField(options,0,"GcommFn");
      if ( !mxIsEmpty(opt) ) {
        mxDestroyArray(mtlb_GCOMfct);
        mtlb_GCOMfct  = mxDuplicateArray(opt);
      }
      
    }

  }

  /* We made it here without problems */

  return(0);

}

