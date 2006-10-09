/*
 * -----------------------------------------------------------------
 * $Revision: 1.3 $
 * $Date: 2006-10-09 23:56:25 $
 * -----------------------------------------------------------------
 * Programmer: Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2005, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see sundials-x.y.z/src/kinsol/LICENSE.
 * -----------------------------------------------------------------
 * Option parsing functions for the KINSOL Matlab interface.
 * -----------------------------------------------------------------
 */

#include <string.h>
#include "kim.h"

extern kim_KINSOLdata kim_Kdata;  /* KINSOL data */
extern kim_MATLABdata kim_Mdata;  /* MATLAB data */

/* 
 * ---------------------------------------------------------------------------------
 * Private constants
 * ---------------------------------------------------------------------------------
 */

#define ONE RCONST(1.0)

/*
 * ---------------------------------------------------------------------------------
 * Redability replacements
 * ---------------------------------------------------------------------------------
 */

#define N  (kim_Kdata->N)
#define ls (kim_Kdata->ls)
#define pm (kim_Kdata->pm)

#define mx_JACfct  (kim_Mdata->mx_JACfct)
#define mx_PSETfct (kim_Mdata->mx_PSETfct)
#define mx_PSOLfct (kim_Mdata->mx_PSOLfct)
#define mx_GLOCfct (kim_Mdata->mx_GLOCfct)
#define mx_GCOMfct (kim_Mdata->mx_GCOMfct)

/*
 * ---------------------------------------------------------------------------------
 * Option handling functions
 * ---------------------------------------------------------------------------------
 */

int get_SolverOptions(const mxArray *options,
                      booleantype *verbose,
                      int *mxiter, int *msbset, int *msbsetsub, 
                      int *etachoice, int *mxnbcf,
                      double *eta, double *egamma, double *ealpha, double *mxnewtstep, 
                      double *relfunc, double *fnormtol, double *scsteptol,
                      double **constraints,
                      booleantype *noInitSetup, booleantype *noMinEps)
{
  mxArray *opt;
  char *bufval;
  int i, buflen, status, n;
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

  *noInitSetup = FALSE;
  *noMinEps = FALSE;
  *verbose = FALSE;

  *constraints = NULL;

  /* Return now if options was empty */

  if (mxIsEmpty(options)) return(0);

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
    if(status != 0) return(status);
    if(strcmp(bufval,"Type1"))         *etachoice = KIN_ETACHOICE1;
    else if(strcmp(bufval,"Type2"))    *etachoice = KIN_ETACHOICE2;
    else if(strcmp(bufval,"Constant")) *etachoice = KIN_ETACONSTANT;
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

  opt = mxGetField(options,0,"Verbose");
  if ( !mxIsEmpty(opt) )
    *verbose = mxIsLogicalScalarTrue(opt);

  opt = mxGetField(options,0,"InitialSetup");
  if ( !mxIsEmpty(opt) )
    *noInitSetup = !mxIsLogicalScalarTrue(opt);

  opt = mxGetField(options,0,"MinBoundEps");
  if ( !mxIsEmpty(opt) )
    *noMinEps = !mxIsLogicalScalarTrue(opt);


  /* Constraints */

  opt = mxGetField(options,0,"Constraints");
  if ( !mxIsEmpty(opt) ) {
    tmp = mxGetPr(opt);
    n = mxGetN(opt);
    if (n == N) {
      *constraints = (double *) malloc(N*sizeof(int));
      for (i=0;i<N;i++)
        (*constraints)[i] = tmp[i];
    } else {
      return(1); /* Error */
    }
  }

  /* We made it here without problems */

  return(0);
}




int get_LinSolvOptions(const mxArray *options,
                       int *mupper, int *mlower,
                       int *mudq, int *mldq, double *dqrely,
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
    if(status != 0) 
      mexErrMsgTxt("Cannot parse LinearSolver.");
    if(!strcmp(bufval,"Band")) ls = LS_BAND;
    else if(!strcmp(bufval,"GMRES")) ls = LS_SPGMR;
    else if(!strcmp(bufval,"BiCGStab")) ls = LS_SPBCG;
    else if(!strcmp(bufval,"TFQMR")) ls = LS_SPTFQMR;
    else if(!strcmp(bufval,"Dense")) ls = LS_DENSE;
    else mexErrMsgTxt("LinearSolver has an illegal value.");
  }
  
  /* Jacobian function */

  opt = mxGetField(options,0,"JacobianFn");
  if ( !mxIsEmpty(opt) ) {
    mxDestroyArray(mx_JACfct);
    mx_JACfct  = mxDuplicateArray(opt);
  }

  /* Band linear solver */

  if (ls==LS_BAND) {

    opt = mxGetField(options,0,"UpperBwidth");
    if ( !mxIsEmpty(opt) )
      *mupper = (int)*mxGetPr(opt);
    
    opt = mxGetField(options,0,"LowerBwidth");
    if ( !mxIsEmpty(opt) )
      *mlower = (int)*mxGetPr(opt);

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
      if (*maxl < 0) 
        mexErrMsgTxt("KrylovMaxDim is negative.");
    }

    /* Preconditioning type */

    opt = mxGetField(options,0,"PrecType");
    if ( !mxIsEmpty(opt) ) {
      buflen = mxGetM(opt) * mxGetN(opt) + 1;
      bufval = mxCalloc(buflen, sizeof(char));
      status = mxGetString(opt, bufval, buflen);
      if(status != 0)
        mexErrMsgTxt("Cannot parse PrecType.");
      if(!strcmp(bufval,"Right")) *ptype = PREC_RIGHT;
      else if(!strcmp(bufval,"None")) *ptype = PREC_NONE;
      else mexErrMsgTxt("PrecType has an illegal value.");
    }

    /* User defined precoditioning */

    opt = mxGetField(options,0,"PrecSetupFn");
    if ( !mxIsEmpty(opt) ) {
      mxDestroyArray(mx_PSETfct);
      mx_PSETfct  = mxDuplicateArray(opt);
    }
  
    opt = mxGetField(options,0,"PrecSolveFn");
    if ( !mxIsEmpty(opt) ) {
      mxDestroyArray(mx_PSOLfct);
      mx_PSOLfct  = mxDuplicateArray(opt);
    }
    
    /* Preconditioner module */
  
    opt = mxGetField(options,0,"PrecModule");
    if ( !mxIsEmpty(opt) ) {
      buflen = mxGetM(opt) * mxGetN(opt) + 1;
      bufval = mxCalloc(buflen, sizeof(char));
      status = mxGetString(opt, bufval, buflen);
      if(status != 0)
        mexErrMsgTxt("Cannot parse PrecModule.");
      if(!strcmp(bufval,"BBDPre")) pm = PM_BBDPRE;
      else if(!strcmp(bufval,"UserDefined")) pm = PM_NONE;
      else mexErrMsgTxt("PrecModule has an illegal value.");
      
    }

    
    if (pm == PM_BBDPRE) {

      opt = mxGetField(options,0,"UpperBwidth");
      if ( !mxIsEmpty(opt) )
        *mupper = (int)*mxGetPr(opt);
    
      opt = mxGetField(options,0,"LowerBwidth");
      if ( !mxIsEmpty(opt) )
        *mlower = (int)*mxGetPr(opt);

      opt = mxGetField(options,0,"UpperBwidthDQ");
      if ( !mxIsEmpty(opt) )
        *mudq = (int)*mxGetPr(opt);

      opt = mxGetField(options,0,"LowerBwidthDQ");
      if ( !mxIsEmpty(opt) )
        *mldq = (int)*mxGetPr(opt);
      
      opt = mxGetField(options,0,"GlocalFn");
      if ( !mxIsEmpty(opt) ) {
        mxDestroyArray(mx_GLOCfct);
        mx_GLOCfct  = mxDuplicateArray(opt);
      }
      else 
        mexErrMsgTxt("GlocalFn required for BBD preconditioner.");
      
      opt = mxGetField(options,0,"GcommFn");
      if ( !mxIsEmpty(opt) ) {
        mxDestroyArray(mx_GCOMfct);
        mx_GCOMfct  = mxDuplicateArray(opt);
      }
      
    }

  }

  /* We made it here without problems */

  return(0);

}

