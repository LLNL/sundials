/*
 * -----------------------------------------------------------------
 * $Revision: 1.1 $
 * $Date: 2006-01-06 19:00:23 $
 * -----------------------------------------------------------------
 * Programmer: Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2005, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see sundials/kinsol/LICENSE.
 * -----------------------------------------------------------------
 * Option parsing functions for the KINSOL Matlab interface.
 * -----------------------------------------------------------------
 */

#include <string.h>
#include "kim.h"

#define ONE RCONST(1.0)

/*
 * ---------------------------------------------------------------------------------
 * Option handling functions
 * ---------------------------------------------------------------------------------
 */

int get_SolverOptions(const mxArray *options,
                      int *mxiter, int *msbset, int *etachoice, int *mxnbcf,
                      double *eta, double *egamma, double *ealpha, double *mxnewtstep, 
                      double *relfunc, double *fnormtol, double *scsteptol,
                      double **constraints,
                      booleantype *noInitSetup, booleantype *noMinEps)
{
  mxArray *opt;
  char *bufval;
  int i, buflen, status, n;
  double *tmp;

  /* Set default values */

  *mxiter = 0;
  *msbset = 0;
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
                       int *ls_tmp,
                       int *maxl, int *maxrs, int *ptype, int *pm_tmp,
                       int *mudq, int *mldq, int *mukeep, int *mlkeep,
                       double *dq,
                       mxArray **mx_tmp_JACfct,
                       mxArray **mx_tmp_PSETfct, mxArray **mx_tmp_PSOLfct,
                       mxArray **mx_tmp_GLOCfct, mxArray **mx_tmp_GCOMfct)
{
  mxArray *opt, *empty;
  char *bufval;
  int buflen, status;
  

  *ls_tmp = LS_SPGMR;

  *mudq   = 0;
  *mldq   = 0;
  *mukeep = 0;
  *mlkeep = 0;
  *dq = 0.0;

  *maxl  = 0;
  *maxrs = 0;

  *ptype = PREC_NONE;
  *pm_tmp = PM_NONE;

  empty = mxCreateDoubleMatrix(0,0,mxREAL);
  *mx_tmp_JACfct  = mxDuplicateArray(empty); 
  *mx_tmp_PSETfct = mxDuplicateArray(empty); 
  *mx_tmp_PSOLfct = mxDuplicateArray(empty); 
  *mx_tmp_GLOCfct = mxDuplicateArray(empty); 
  *mx_tmp_GCOMfct = mxDuplicateArray(empty); 

  if (mxIsEmpty(options)) return(0);

  /* Linear solver type */

  opt = mxGetField(options,0,"LinearSolver");
  if ( !mxIsEmpty(opt) ) {
    buflen = mxGetM(opt) * mxGetN(opt) + 1;
    bufval = mxCalloc(buflen, sizeof(char));
    status = mxGetString(opt, bufval, buflen);
    if(status != 0) return(status);
    if(!strcmp(bufval,"Dense")) *ls_tmp = LS_DENSE;
    else if(!strcmp(bufval,"GMRES")) *ls_tmp = LS_SPGMR;
    else if(!strcmp(bufval,"BiCGStab")) *ls_tmp = LS_SPBCG;
    else mexErrMsgTxt("LinearSolver has an illegal value.");
  }
  
  opt = mxGetField(options,0,"JacobianFn");
  if ( !mxIsEmpty(opt) )
    *mx_tmp_JACfct = mxDuplicateArray(opt);
  
  opt = mxGetField(options,0,"KrylovMaxDim");
  if ( !mxIsEmpty(opt) )
    *maxl = (int)*mxGetPr(opt);

  opt = mxGetField(options,0,"MaxNumRestarts");
  if ( !mxIsEmpty(opt) )
    *maxrs = (int)*mxGetPr(opt);  
  
  /* Preconditioner? -- this property not used right now!! */

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

  /* Preconditioner module */
  
  opt = mxGetField(options,0,"PrecModule");
  if ( !mxIsEmpty(opt) ) {
    buflen = mxGetM(opt) * mxGetN(opt) + 1;
    bufval = mxCalloc(buflen, sizeof(char));
    status = mxGetString(opt, bufval, buflen);
    if(status != 0)
      mexErrMsgTxt("Cannot parse PrecModule.");
    if(!strcmp(bufval,"BBDPre")) *pm_tmp = PM_BBDPRE;
    else if(!strcmp(bufval,"UserDefined")) *pm_tmp = PM_NONE;
    else mexErrMsgTxt("PrecModule has an illegal value.");
  }

  /* BBDPre options */

  if (*pm_tmp == PM_BBDPRE) {

    opt = mxGetField(options,0,"UpperBwidth");
    if ( !mxIsEmpty(opt) )
      *mukeep = (int)*mxGetPr(opt);
    
    opt = mxGetField(options,0,"LowerBwidth");
    if ( !mxIsEmpty(opt) )
      *mlkeep = (int)*mxGetPr(opt);

    opt = mxGetField(options,0,"UpperBwidthDQ");
    if ( !mxIsEmpty(opt) )
      *mudq = (int)*mxGetPr(opt);

    opt = mxGetField(options,0,"LowerBwidthDQ");
    if ( !mxIsEmpty(opt) )
      *mldq = (int)*mxGetPr(opt);
      
    opt = mxGetField(options,0,"GlocalFn");
    if ( !mxIsEmpty(opt) ) 
      *mx_tmp_GLOCfct = mxDuplicateArray(opt);
    else 
      mexErrMsgTxt("GlocalFn required for BBD preconditioner.");
    
    opt = mxGetField(options,0,"GcommFn");
    if ( !mxIsEmpty(opt) ) 
      *mx_tmp_GCOMfct = mxDuplicateArray(opt);
    
  }
  
  /* User provided preconditioner */

  if (*pm_tmp == PM_NONE) {
  
    opt = mxGetField(options,0,"PrecSetupFn");
    if ( !mxIsEmpty(opt) ) 
      *mx_tmp_PSETfct = mxDuplicateArray(opt);
  
    opt = mxGetField(options,0,"PrecSolveFn");
    if ( !mxIsEmpty(opt) )
      *mx_tmp_PSOLfct = mxDuplicateArray(opt);

  }

  /* We made it here without problems */

  mxDestroyArray(empty);

  return(0);

}

