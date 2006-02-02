/*
 * -----------------------------------------------------------------
 * $Revision: 1.2 $
 * $Date: 2006-02-02 00:39:04 $
 * -----------------------------------------------------------------
 * Programmer: Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2005, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see sundials/cvodes/LICENSE.
 * -----------------------------------------------------------------
 * Option parsing functions for the CVODES Matlab interface.
 * -----------------------------------------------------------------
 */

#include <string.h>
#include "cvm.h"

#define ONE RCONST(1.0)

/*
 * ---------------------------------------------------------------------------------
 * Option handling functions
 * ---------------------------------------------------------------------------------
 */

int get_IntgrOptions(const mxArray *options,
                     int *lmm, int *iter, int *maxord, booleantype *sld,
                     long int *mxsteps,
                     int *itol, realtype *reltol, double *Sabstol, double **Vabstol,
                     double *hin, double *hmax, double *hmin, double *tstop, booleantype *tstopSet,
                     int *Ng_tmp, mxArray **mx_tmp_Gfct,
                     booleantype *quad, booleantype *fsa, booleantype *asa,
                     booleantype *monitor, mxArray **mx_tmp_MONfct, mxArray **mx_tmp_MONdata)
{
  mxArray *opt, *empty;
  char *bufval;
  int i, buflen, status, q, m, n;
  double *tmp;

  /* Set default values */

  *lmm = CV_BDF;
  *iter = CV_NEWTON;
  *maxord = 5;

  *mxsteps = 0;

  *itol = CV_SS;
  *reltol = 1.0e-3;
  *Sabstol = 1.0e-6;
  *Vabstol = NULL;

  *Ng_tmp = 0;
  
  *hin = 0.0;
  *hmax = 0.0;
  *hmin = 0.0;
  *tstopSet = FALSE;
  *sld = FALSE;

  *quad = FALSE;
  *fsa  = FALSE;
  *asa  = FALSE;

  empty = mxCreateDoubleMatrix(0,0,mxREAL);
  *mx_tmp_Gfct = mxDuplicateArray(empty);
  *mx_tmp_MONfct = mxDuplicateArray(empty);
  *mx_tmp_MONdata = mxDuplicateArray(empty);

  /* Return now if options was empty */

  if (mxIsEmpty(options)) return(0);

  /* Tolerances */

  opt = mxGetField(options,0,"RelTol");
  if ( !mxIsEmpty(opt) ) {
    *reltol = *mxGetPr(opt);
    if (*reltol < 0.0 )
      mexErrMsgTxt("RelTol is negative.");
  }

  opt = mxGetField(options,0,"AbsTol");
  if ( !mxIsEmpty(opt) ) {
    m = mxGetM(opt);
    n = mxGetN(opt);
    if ( (n != 1) && (m != 1) )
      mexErrMsgTxt("AbsTol is not a scalar or a vector.");
    if ( m > n ) n = m;
    tmp = mxGetPr(opt);
    if (n == 1) {
      *itol = CV_SS;
      *Sabstol = *tmp;
      if (*Sabstol < 0.0)
        mexErrMsgTxt("AbsTol is negative.");
    } else if (n == N) {
      *itol = CV_SV;
      *Vabstol = (double *) malloc(N*sizeof(double));
      for(i=0;i<N;i++) {
        (*Vabstol)[i] = tmp[i];
        if (tmp[i] < 0.0)
          mexErrMsgTxt("AbsTol has a negative component.");
      }
    } else {
      mexErrMsgTxt("AbsTol does not contain N elements.");
    }
  }

  /* LMM */

  opt = mxGetField(options,0,"Adams");
  if ( !mxIsEmpty(opt) ) {
    buflen = mxGetM(opt) * mxGetN(opt) + 1;
    bufval = mxCalloc(buflen, sizeof(char));
    status = mxGetString(opt, bufval, buflen);
    if(status != 0)
      mexErrMsgTxt("Cannot parse Adams.");
    if(!strcmp(bufval,"on")) {*lmm = CV_ADAMS; *maxord = 12;}
    else if(!strcmp(bufval,"off")) {*lmm = CV_BDF; *maxord = 5;}
    else mexErrMsgTxt("Adams has an illegal value.");
  }

  /* ITER */
  
  opt = mxGetField(options,0,"NonlinearSolver");
  if ( !mxIsEmpty(opt) ) {
    buflen = mxGetM(opt) * mxGetN(opt) + 1;
    bufval = mxCalloc(buflen, sizeof(char));
    status = mxGetString(opt, bufval, buflen);
    if(status != 0)
      mexErrMsgTxt("Cannot parse NonlinearSolver.");
    if(!strcmp(bufval,"Functional")) *iter = CV_FUNCTIONAL;
    else if(!strcmp(bufval,"Newton")) *iter = CV_NEWTON;
    else mexErrMsgTxt("NonlinearSolver has an illegal value.");

  }
    
  /* Number of root functions */

  opt = mxGetField(options,0,"NumRoots");
  if ( !mxIsEmpty(opt) ) {
    *Ng_tmp = (int)*mxGetPr(opt);
    if (*Ng_tmp < 0)
      mexErrMsgTxt("NumRoots is negative.");
    if (*Ng_tmp > 0) {
      /* Roots function */
      opt = mxGetField(options,0,"RootsFn");
      if ( !mxIsEmpty(opt) ) {
        *mx_tmp_Gfct = mxDuplicateArray(opt);
      } else {
        mexErrMsgTxt("RootsFn required for NumRoots > 0");
      }
    }
  }

  /* Maximum number of steps */
  opt = mxGetField(options,0,"MaxNumSteps");
  if ( !mxIsEmpty(opt) ) {
    *mxsteps = (int)*mxGetPr(opt);
    if (*mxsteps < 0)
      mexErrMsgTxt("MaxNumSteps is negative.");
  }

  /* Maximum order */
  
  opt = mxGetField(options,0,"MaxOrder");
  if ( !mxIsEmpty(opt) ) {
    q = (int)*mxGetPr(opt);
    if (q <= 0)
      mexErrMsgTxt("MaxOrder must be positive.");
    if (q > *maxord)
      mexErrMsgTxt("MaxOrder is too large for the Method specified.");
    *maxord = q;
  }

  /* Initial step size */

  opt = mxGetField(options,0,"InitialStep");
  if ( !mxIsEmpty(opt) ) {
    *hin = *mxGetPr(opt);
  }

  /* Maximum step size */

  opt = mxGetField(options,0,"MaxStep");
  if ( !mxIsEmpty(opt) ) {
    tmp = mxGetPr(opt);
    if (*tmp < 0.0)
      mexErrMsgTxt("MaxStep is negative.");
    if ( mxIsInf(*tmp) )
      *hmax = 0.0;
    else
      *hmax = *tmp;
  }

  /* Minimum step size */

  opt = mxGetField(options,0,"MinStep");
  if ( !mxIsEmpty(opt) ) {
    *hmin = *mxGetPr(opt);
    if (*hmin < 0.0)
      mexErrMsgTxt("MinStep is negative.");
  }

  /* Stopping time */

  opt = mxGetField(options,0,"StopTime");
  if ( !mxIsEmpty(opt) ) {
    *tstop = *mxGetPr(opt);
    *tstopSet = TRUE;
  }

  /* Stability Limit Detection */

  opt = mxGetField(options,0,"StabilityLimDet");
  if ( !mxIsEmpty(opt) ) {
    buflen = mxGetM(opt) * mxGetN(opt) + 1;
    bufval = mxCalloc(buflen, sizeof(char));
    status = mxGetString(opt, bufval, buflen);
    if(status != 0) 
      mexErrMsgTxt("Cannot parse StabilityLimDet.");
    if(!strcmp(bufval,"on")) *sld = TRUE;
    else if(!strcmp(bufval,"off")) *sld = FALSE;
    else mexErrMsgTxt("StabilityLimDet has an illegal value.");
  }

  /* Quadratures? */

  opt = mxGetField(options,0,"Quadratures");
  if ( !mxIsEmpty(opt) ) {
    buflen = mxGetM(opt) * mxGetN(opt) + 1;
    bufval = mxCalloc(buflen, sizeof(char));
    status = mxGetString(opt, bufval, buflen);
    if(status != 0)
      mexErrMsgTxt("Cannot parse Quadratures.");
    if(!strcmp(bufval,"on")) *quad = TRUE;
    else if(!strcmp(bufval,"off")) *quad = FALSE;
    else mexErrMsgTxt("Quadratures has an illegal value.");
  }

  /* Sensitivity analysis? */

  opt = mxGetField(options,0,"SensiAnalysis");
  if ( !mxIsEmpty(opt) ) {
    buflen = mxGetM(opt) * mxGetN(opt) + 1;
    bufval = mxCalloc(buflen, sizeof(char));
    status = mxGetString(opt, bufval, buflen);
    if(status != 0) 
      mexErrMsgTxt("Cannot parse SensiAnalysis.");
    if(!strcmp(bufval,"FSA")) {*fsa = TRUE; *asa = FALSE;}
    else if(!strcmp(bufval,"ASA")) {*fsa = FALSE; *asa = TRUE;}
    else if(!strcmp(bufval,"off")) {*fsa = FALSE; *asa = FALSE;}
    else mexErrMsgTxt("SensiAnalysis has an illegal value.");
  }  

  /* Monitor? */
  opt = mxGetField(options,0,"MonitorFn");
  if ( !mxIsEmpty(opt) ) {
    *monitor = TRUE;
    *mx_tmp_MONfct = mxDuplicateArray(opt);
    opt = mxGetField(options,0,"MonitorData");
    if ( !mxIsEmpty(opt) ) {
      *mx_tmp_MONdata = mxDuplicateArray(opt);
    }
  }

  /* We made it here without problems */

  mxDestroyArray(empty);

  return(0);
}


int get_LinSolvOptions(const mxArray *options,
                       int *ls_tmp,
                       int *mupper, int *mlower,
                       int *mudq, int *mldq, double *dqrely,
                       int *ptype, int *gstype, int *maxl, int *pm_tmp,
                       mxArray **mx_tmp_JACfct,
                       mxArray **mx_tmp_PSETfct, mxArray **mx_tmp_PSOLfct,
                       mxArray **mx_tmp_GLOCfct, mxArray **mx_tmp_GCOMfct)
{
  mxArray *opt, *empty;
  char *bufval;
  int buflen, status;
  
  *ls_tmp = LS_DENSE;

  *mupper = 0;
  *mlower = 0;

  *mudq = 0;
  *mldq = 0;
  *dqrely = 0.0;

  *ptype = PREC_NONE;
  *gstype = MODIFIED_GS;
  *maxl = 0;

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
    if(status != 0) 
      mexErrMsgTxt("Cannot parse LinearSolver.");
    if(!strcmp(bufval,"Diag")) *ls_tmp = LS_DIAG;
    else if(!strcmp(bufval,"Band")) *ls_tmp = LS_BAND;
    else if(!strcmp(bufval,"GMRES")) *ls_tmp = LS_SPGMR;
    else if(!strcmp(bufval,"BiCGStab")) *ls_tmp = LS_SPBCG;
    else if(!strcmp(bufval,"TFQMR")) *ls_tmp = LS_SPTFQMR;
    else if(!strcmp(bufval,"Dense")) *ls_tmp = LS_DENSE;
    else mexErrMsgTxt("LinearSolver has an illegal value.");
  }
  
  opt = mxGetField(options,0,"JacobianFn");
  if ( !mxIsEmpty(opt) )
    *mx_tmp_JACfct = mxDuplicateArray(opt);
  
  opt = mxGetField(options,0,"UpperBwidth");
  if ( !mxIsEmpty(opt) )
    *mupper = (int)*mxGetPr(opt);

  opt = mxGetField(options,0,"LowerBwidth");
  if ( !mxIsEmpty(opt) )
    *mlower = (int)*mxGetPr(opt);

  opt = mxGetField(options,0,"PrecType");
  if ( !mxIsEmpty(opt) ) {
    buflen = mxGetM(opt) * mxGetN(opt) + 1;
    bufval = mxCalloc(buflen, sizeof(char));
    status = mxGetString(opt, bufval, buflen);
    if(status != 0)
      mexErrMsgTxt("Cannot parse PrecType.");
    if(!strcmp(bufval,"Left")) *ptype = PREC_LEFT;
    else if(!strcmp(bufval,"Right")) *ptype = PREC_RIGHT;
    else if(!strcmp(bufval,"Both")) *ptype = PREC_BOTH;
    else if(!strcmp(bufval,"None")) *ptype = PREC_NONE;
    else mexErrMsgTxt("PrecType has an illegal value.");
  }

  opt = mxGetField(options,0,"GramSchmidtType");
  if ( !mxIsEmpty(opt) ) {
    buflen = mxGetM(opt) * mxGetN(opt) + 1;
    bufval = mxCalloc(buflen, sizeof(char));
    status = mxGetString(opt, bufval, buflen);
    if(status != 0)
      mexErrMsgTxt("Cannot parse GramSchmidtType.");
    if(!strcmp(bufval,"Classical")) *gstype = CLASSICAL_GS;
    else if(!strcmp(bufval,"Modified")) *gstype = MODIFIED_GS;
    else mexErrMsgTxt("GramSchmidtType has an illegal value.");
  }

  opt = mxGetField(options,0,"KrylovMaxDim");
  if ( !mxIsEmpty(opt) ) {
    *maxl = (int)*mxGetPr(opt);
    if (*maxl < 0) 
      mexErrMsgTxt("KrylovMaxDim is negative.");
  }
  
  opt = mxGetField(options,0,"PrecSetupFn");
  if ( !mxIsEmpty(opt) ) 
    *mx_tmp_PSETfct = mxDuplicateArray(opt);
  
  opt = mxGetField(options,0,"PrecSolveFn");
  if ( !mxIsEmpty(opt) )
    *mx_tmp_PSOLfct = mxDuplicateArray(opt);

  /* Preconditioner module */
  
  opt = mxGetField(options,0,"PrecModule");
  if ( !mxIsEmpty(opt) ) {
    buflen = mxGetM(opt) * mxGetN(opt) + 1;
    bufval = mxCalloc(buflen, sizeof(char));
    status = mxGetString(opt, bufval, buflen);
    if(status != 0)
      mexErrMsgTxt("Cannot parse PrecModule.");
    if(!strcmp(bufval,"BandPre")) *pm_tmp = PM_BANDPRE;
    else if(!strcmp(bufval,"BBDPre")) *pm_tmp = PM_BBDPRE;
    else if(!strcmp(bufval,"UserDefined")) *pm_tmp = PM_NONE;
    else mexErrMsgTxt("PrecModule has an illegal value.");

    opt = mxGetField(options,0,"UpperBwidthDQ");
    if ( !mxIsEmpty(opt) )
      *mudq = (int)*mxGetPr(opt);

    opt = mxGetField(options,0,"LowerBwidthDQ");
    if ( !mxIsEmpty(opt) )
      *mldq = (int)*mxGetPr(opt);

    if (*pm_tmp == PM_BBDPRE) {
      
      opt = mxGetField(options,0,"GlocalFn");
      if ( !mxIsEmpty(opt) ) 
        *mx_tmp_GLOCfct = mxDuplicateArray(opt);
      else 
        mexErrMsgTxt("GlocalFn required for BBD preconditioner.");
        
      opt = mxGetField(options,0,"GcommFn");
      if ( !mxIsEmpty(opt) ) 
        *mx_tmp_GCOMfct = mxDuplicateArray(opt);
      
    }

  }
  
  /* We made it here without problems */

  mxDestroyArray(empty);

  return(0);

}


int get_QuadOptions(const mxArray *options,
                    int *Nq_tmp, double **yQ0, mxArray **mx_tmp_QUADfct,
                    booleantype *errconQ,
                    int *itolQ, double *reltolQ, double *SabstolQ, double **VabstolQ)
{
  mxArray *opt, *empty;
  char *bufval;
  int i, buflen, status, m, n;
  double *tmp;

  *Nq_tmp = 0;
  *yQ0 = NULL;

  *errconQ = FALSE;
  *itolQ = CV_SS;
  *reltolQ = 1.0e-4;
  *SabstolQ = 1.0e-6;
  *VabstolQ = NULL;

  empty = mxCreateDoubleMatrix(0,0,mxREAL);
  *mx_tmp_QUADfct = mxDuplicateArray(empty);

  /* Note that options cannot be empty if we got here */

  /* Initial conditions for quadratures */

  opt = mxGetField(options,0,"QuadInitCond");
  if ( !mxIsEmpty(opt) ) {
    m = mxGetM(opt);
    n = mxGetN(opt);
    if ( (n != 1) && (m != 1) )
      mexErrMsgTxt("QuadInitCond is not a vector.");
    if ( n > m ) *Nq_tmp = n;
    else         *Nq_tmp = m;
    tmp = mxGetPr(opt);
    *yQ0 = (double *)malloc((*Nq_tmp)*sizeof(double));
    for(i=0;i<*Nq_tmp;i++) 
      (*yQ0)[i] = tmp[i];
  } else {
    mexErrMsgTxt("QuadInitCond required for quadrature integration.");
  }

  /* Quadrature function */

  opt = mxGetField(options,0,"QuadRhsFn");
  if ( !mxIsEmpty(opt) ) {
    *mx_tmp_QUADfct = mxDuplicateArray(opt);
  } else {
    mexErrMsgTxt("QuadRhsFn required for quadrature integration.");
  }

  /* Quadrature error control and tolerances */

  opt = mxGetField(options,0,"QuadErrControl");
  if ( !mxIsEmpty(opt) ) {
    buflen = mxGetM(opt) * mxGetN(opt) + 1;
    bufval = mxCalloc(buflen, sizeof(char));
    status = mxGetString(opt, bufval, buflen);
    if(status != 0)
      mexErrMsgTxt("Canot parse QuadErrControl.");
    if(!strcmp(bufval,"on")) {
      *errconQ = TRUE;
      opt = mxGetField(options,0,"QuadRelTol");
      if ( !mxIsEmpty(opt) ) {
        *reltolQ = *mxGetPr(opt);
        if (*reltolQ < 0.0)
          mexErrMsgTxt("QuadRelTol is negative.");
      } 
      opt = mxGetField(options,0,"QuadAbsTol");
      if ( !mxIsEmpty(opt) ) {
        m = mxGetN(opt);
        n = mxGetM(opt);
        if ( (n != 1) && (m != 1) )
          mexErrMsgTxt("QuadAbsTol is not a scalar or a vector.");
        if ( m > n ) n = m;
        tmp = mxGetPr(opt);
        if (n == 1) {
          *itolQ = CV_SS;
          *SabstolQ = *tmp;
          if (*SabstolQ < 0.0)
            mexErrMsgTxt("QuadAbsTol is negative.");
        } else if (n == Nq) {
          *itolQ = CV_SV;
          *VabstolQ = (double *)malloc((*Nq_tmp)*sizeof(double));
          for(i=0;i<Nq;i++) {
            (*VabstolQ)[i] = tmp[i];
            if (tmp[i] < 0.0)
              mexErrMsgTxt("QuadAbsTol has a negative component.");
          }
        } else {
          mexErrMsgTxt("QuadAbsTol does not contain Nq elements.");
        }
      }
    } else if(!strcmp(bufval,"off")) {
      *errconQ = FALSE;
    } else {
      mexErrMsgTxt("QuadErrControl has an illegal value.");
    }
  }

  /* We made it here without problems */

  mxDestroyArray(empty);

  return(0);
}

int get_ASAOptions(const mxArray *options, int *Nd_tmp, int *interp_tmp)
{
  mxArray *opt;
  char *bufval;
  int buflen, status;

  *Nd_tmp = 100;
  *interp_tmp = CV_HERMITE;

  opt = mxGetField(options,0,"ASANumPoints");
  if ( !mxIsEmpty(opt) ) {
    *Nd_tmp = (int)*mxGetPr(opt);
    if (*Nd_tmp <= 0)
      mexErrMsgTxt("ASaNumPoints must be positive");
  } else {
    mexErrMsgTxt("ASANumPoints required for ASA.");
  }

  opt = mxGetField(options,0,"ASAInterpType");
  if ( !mxIsEmpty(opt) ) {
    buflen = mxGetM(opt) * mxGetN(opt) + 1;
    bufval = mxCalloc(buflen, sizeof(char));
    status = mxGetString(opt, bufval, buflen);
    if(status != 0)
      mexErrMsgTxt("Cannot parse ASAInterpType");
    if(!strcmp(bufval,"Hermite")) *interp_tmp = CV_HERMITE;
    else if(!strcmp(bufval,"Polynomial")) *interp_tmp = CV_POLYNOMIAL;
    else mexErrMsgTxt("ASAInterpType has an illegal value.");
  }

  return(0);
}

int get_FSAOptions(const mxArray *options, 
                   int *Ns_tmp, double **yS0,
                   int *ism_tmp,
                   char **pfield_name, int **plist, double **pbar,
                   int *Srhs, mxArray **mx_tmp_SRHSfct, double *rho,
                   booleantype *errconS, int *itolS, double *reltolS, 
                   double **SabstolS, double **VabstolS)
{
  mxArray *opt, *empty;
  char *bufval;
  int i, is, m, n, buflen, status;
  double *tmp;

  /* Set default values */

  *Ns_tmp = 0;
  *yS0 = NULL;

  *ism_tmp = CV_STAGGERED;

  *errconS = TRUE;

  *rho = 0.0;

  *itolS = CV_EE;
  *SabstolS = NULL;
  *VabstolS = NULL;

  *pfield_name = NULL;
  *plist = NULL;
  *pbar  = NULL;

  empty = mxCreateDoubleMatrix(0,0,mxREAL);
  *mx_tmp_SRHSfct = mxDuplicateArray(empty);
  *Srhs = 0;

  /* Note that options cannot be empty if we got here */

  /* Initial conditions for sensitivities */

  opt = mxGetField(options,0,"FSAInitCond");

  if ( !mxIsEmpty(opt) ) {
    m = mxGetM(opt);
    n = mxGetN(opt);
    if ( m != N ) 
      mexErrMsgTxt("FSAInitCond does not have N rows.");
    *Ns_tmp = n;
    tmp = mxGetPr(opt);
    *yS0 = (double *)malloc((*Ns_tmp)*N*sizeof(double));
    for (i=0; i<(*Ns_tmp)*N; i++)
      (*yS0)[i] = tmp[i];
  } else {
    mexErrMsgTxt("FSAInitCond required for FSA.");
  }

  /* ISM */

  opt = mxGetField(options,0,"FSAMethod");
  if ( !mxIsEmpty(opt) ) {
    buflen = mxGetM(opt) * mxGetN(opt) + 1;
    bufval = mxCalloc(buflen, sizeof(char));
    status = mxGetString(opt, bufval, buflen);
    if(status != 0)
      mexErrMsgTxt("Cannot parse FSAMethod");
    if(!strcmp(bufval,"Simultaneous")) *ism_tmp = CV_SIMULTANEOUS;
    else if(!strcmp(bufval,"Staggered")) *ism_tmp = CV_STAGGERED;
    else if(!strcmp(bufval,"Staggered1")) *ism_tmp = CV_STAGGERED1;
    else mexErrMsgTxt("FSAMethod has an illegal value.");
  }

  /* Field name in data structure for params. */

  opt = mxGetField(options,0,"FSAParamField");
  if ( !mxIsEmpty(opt) ) {
    buflen = mxGetM(opt) * mxGetN(opt) + 1;
    bufval = mxCalloc(buflen, sizeof(char));
    status = mxGetString(opt, bufval, buflen);
    if(status != 0)
      mexErrMsgTxt("Could not parse FSAParamField.");
    *pfield_name = mxCalloc(buflen, sizeof(char));
    strcpy((*pfield_name), bufval);
  }  

  /* PLIST */

  opt = mxGetField(options,0,"FSAParamList");
  if ( !mxIsEmpty(opt) ) {
    tmp = mxGetPr(opt);
    m = mxGetM(opt);
    n = mxGetN(opt);
    if ( (n != 1) && (m != 1) )
      mexErrMsgTxt("FSAParamList is not a vector.");
    if (m > n) n = m;
    if ( n != *Ns_tmp)
      mexErrMsgTxt("FSAParamList does not contain Ns elements.");
    *plist = (int *) malloc((*Ns_tmp)*sizeof(int));
    for (is=0;is<*Ns_tmp;is++)
      (*plist)[is] = (int) tmp[is];
  }

  /* PBAR */

  opt = mxGetField(options,0,"FSAParamScales");
  if ( !mxIsEmpty(opt) ) {
    m = mxGetM(opt);
    n = mxGetN(opt);
    if ( (n != 1) && (m != 1) )
      mexErrMsgTxt("FSAParamScales is not a vector.");
    if ( m > n ) n = m;
    if ( n != *Ns_tmp)
      mexErrMsgTxt("FSAParamScales does not contain Ns elements.");
    tmp = mxGetPr(opt);
    *pbar = (double *) malloc((*Ns_tmp)*sizeof(double));
    for(i=0;i<*Ns_tmp;i++)
      (*pbar)[i] = tmp[i];
  }

  /* Rho */

  opt = mxGetField(options,0,"FSADQparam");
  if ( !mxIsEmpty(opt) )
    *rho = *mxGetPr(opt);

  /* Error control */

  opt = mxGetField(options,0,"FSAErrControl");
  if ( !mxIsEmpty(opt) ) {
    buflen = mxGetM(opt) * mxGetN(opt) + 1;
    bufval = mxCalloc(buflen, sizeof(char));
    status = mxGetString(opt, bufval, buflen);
    if(status != 0)
      mexErrMsgTxt("Canot parse FSAErrControl.");
    if(!strcmp(bufval,"off")) *errconS = FALSE;
    else if(!strcmp(bufval,"on")) *errconS = TRUE;
    else mexErrMsgTxt("FSAErrControl has an illegal value.");
  }

  /* Tolerances */
  
  opt = mxGetField(options,0,"FSARelTol");

  if ( !mxIsEmpty(opt) ) {
    *reltolS = *mxGetPr(opt);
    if (*reltolS < 0.0)
      mexErrMsgTxt("FSARelTol is negative.");
    opt = mxGetField(options,0,"FSAAbsTol");
    if ( !mxIsEmpty(opt) ) {
      m = mxGetM(opt);
      n = mxGetN(opt);
      if ( (m == 1) && (n == *Ns_tmp) ) {
        *itolS = CV_SS;
        tmp = mxGetPr(opt);
        *SabstolS = (double *) malloc((*Ns_tmp)*sizeof(double));
        for (is=0; is<*Ns_tmp; is++) {
          (*SabstolS)[is] = tmp[is];
          if ( tmp[is] < 0.0 )
            mexErrMsgTxt("FSAAbsTol has a negative component.");
        }
      } else if ( (m == N) && (n == *Ns_tmp) ) {
        *itolS = CV_SV;
        tmp = mxGetPr(opt);
        *VabstolS = (double *)malloc((*Ns_tmp)*N*sizeof(double));
        for (i=0; i<(*Ns_tmp)*N; i++) {
          (*VabstolS)[i] = tmp[i];
          if ( tmp[i] < 0.0 )
            mexErrMsgTxt("FSAAbsTol has a negative component.");
        }
      } else {
        mexErrMsgTxt("FSAAbsTol must be either a 1xNs vector or a NxNs matrix.");
      }
    } else {
      *itolS = CV_EE;
    }
  }

  /* Sensitivity RHS function type */

  opt = mxGetField(options,0,"FSARhsType");
  if ( !mxIsEmpty(opt) ) {
    buflen = mxGetM(opt) * mxGetN(opt) + 1;
    bufval = mxCalloc(buflen, sizeof(char));
    status = mxGetString(opt, bufval, buflen);
    if(status != 0)
      mexErrMsgTxt("Cannot parse FSARhsType.");
    if(!strcmp(bufval,"One")) *Srhs = 1;
    else if(!strcmp(bufval,"All")) *Srhs = 2;
    else if(!strcmp(bufval,"None")) *Srhs = 0;
    else mexErrMsgTxt("FSARhsType has an illegal value.");

    if (*Srhs != 0) {
      opt = mxGetField(options,0,"FSARhsFn");
      if ( !mxIsEmpty(opt) ) {
        *mx_tmp_SRHSfct = mxDuplicateArray(opt);
      } else {
        mexErrMsgTxt("FSARhsFn required for FSARhsType other than None");
      }
    }
  }

  /* We made it here without problems */

  mxDestroyArray(empty);

  return(0);

}
