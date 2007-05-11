/*
 * -----------------------------------------------------------------
 * $Revision: 1.7 $
 * $Date: 2007-05-11 18:51:32 $
 * -----------------------------------------------------------------
 * Programmer: Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2005, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see sundials-x.y.z/src/cvodes/LICENSE.
 * -----------------------------------------------------------------
 * Option parsing functions for the CVODES Matlab interface.
 * -----------------------------------------------------------------
 */

#include <string.h>
#include "cvm.h"

/*extern cvm_CVODESdata cvm_Cdata;*/

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


#define N           (thisPb->n) 
#define Ns          (thisPb->ns) 
#define Ng          (thisPb->ng) 
#define ls          (thisPb->LS) 
#define pm          (thisPb->PM) 

#define mx_JACfct   (thisPb->JACfct)
#define mx_PSETfct  (thisPb->PSETfct)
#define mx_PSOLfct  (thisPb->PSOLfct)
#define mx_GLOCfct  (thisPb->GLOCfct)
#define mx_GCOMfct  (thisPb->GCOMfct)
#define mx_Gfct     (thisPb->Gfct)

#define mon         (thisPb->Mon)
#define tstopSet    (thisPb->TstopSet)

#define mx_MONfct   (thisPb->MONfct)
#define mx_MONdata  (thisPb->MONdata)

/*
 * ---------------------------------------------------------------------------------
 * Option handling functions
 * ---------------------------------------------------------------------------------
 */

void get_IntgrOptions(const mxArray *options, cvmPbData thisPb, booleantype fwd,
                      int *lmm, int *iter, int *maxord, booleantype *sld,
                      long int *mxsteps,
                      int *itol, realtype *reltol, double *Sabstol, double **Vabstol,
                      double *hin, double *hmax, double *hmin, double *tstop)
{
  mxArray *opt;
  char *bufval;
  int i, buflen, status, q, m, n;
  double *tmp;
  char *fctName;
  char *fwd_fctName = "CVodeInit/CVodeReInit";
  char *bck_fctName = "CVodeInitB/CVodeReInitB";

  if (fwd) fctName = fwd_fctName;
  else     fctName = bck_fctName;

  /* Set default values */

  *lmm = CV_BDF;
  *iter = CV_NEWTON;
  *maxord = 5;
  *sld = FALSE;

  *mxsteps = 0;

  *itol = CV_SS;
  *reltol = 1.0e-3;
  *Sabstol = 1.0e-6;
  *Vabstol = NULL;

  *hin = 0.0;
  *hmax = 0.0;
  *hmin = 0.0;

  tstopSet = FALSE;
  mon = FALSE;

  /* Return now if options was empty */

  if (mxIsEmpty(options)) return;

  /* Tolerances */

  opt = mxGetField(options,0,"RelTol");
  if ( !mxIsEmpty(opt) ) {
    *reltol = *mxGetPr(opt);
    if (*reltol < 0.0 ) cvmErrHandler(-999, "CVODES", fctName,
                                      "RelTol is negative.", NULL);
  }

  opt = mxGetField(options,0,"AbsTol");
  if ( !mxIsEmpty(opt) ) {
    m = mxGetM(opt);
    n = mxGetN(opt);
    if ( (n != 1) && (m != 1) ) cvmErrHandler(-999, "CVODES", fctName,
                                              "AbsTol is not a scalar or a vector.", NULL);
    if ( m > n ) n = m;
    tmp = mxGetPr(opt);
    if (n == 1) {
      *itol = CV_SS;
      *Sabstol = *tmp;
      if (*Sabstol < 0.0) cvmErrHandler(-999, "CVODES", fctName,
                                        "AbsTol is negative.", NULL);
    } else if (n == N) {
      *itol = CV_SV;
      *Vabstol = (double *) malloc(N*sizeof(double));
      for(i=0;i<N;i++) {
        (*Vabstol)[i] = tmp[i];
        if (tmp[i] < 0.0) cvmErrHandler(-999, "CVODES", fctName,
                                        "AbsTol has a negative component.", NULL);
      }
    } else {
      cvmErrHandler(-999, "CVODES", fctName,
                    "AbsTol does not contain N elements.", NULL);
    }
  }

  /* LMM */

  opt = mxGetField(options,0,"LMM");
  if ( !mxIsEmpty(opt) ) {
    buflen = mxGetM(opt) * mxGetN(opt) + 1;
    bufval = mxCalloc(buflen, sizeof(char));
    status = mxGetString(opt, bufval, buflen);
    if(status != 0) cvmErrHandler(-999, "CVODES", fctName,
                                  "Cannot parse LMM.", NULL);
    if(!strcmp(bufval,"Adams")) {*lmm = CV_ADAMS; *maxord = 12;}
    else if(!strcmp(bufval,"BDF")) {*lmm = CV_BDF; *maxord = 5;}
    else cvmErrHandler(-999, "CVODES", fctName,
                       "LMM has an illegal value.", NULL);
  }

  /* ITER */
  
  opt = mxGetField(options,0,"NonlinearSolver");
  if ( !mxIsEmpty(opt) ) {
    buflen = mxGetM(opt) * mxGetN(opt) + 1;
    bufval = mxCalloc(buflen, sizeof(char));
    status = mxGetString(opt, bufval, buflen);
    if(status != 0) cvmErrHandler(-999, "CVODES", fctName,
                                  "Cannot parse NonlinearSolver.", NULL);
    if(!strcmp(bufval,"Functional")) *iter = CV_FUNCTIONAL;
    else if(!strcmp(bufval,"Newton")) *iter = CV_NEWTON;
    else cvmErrHandler(-999, "CVODES", fctName,
                       "NonlinearSolver has an illegal value.", NULL);
  }
    
  /* Maximum number of steps */

  opt = mxGetField(options,0,"MaxNumSteps");
  if ( !mxIsEmpty(opt) ) {
    *mxsteps = (int)*mxGetPr(opt);
    if (*mxsteps < 0) cvmErrHandler(-999, "CVODES", fctName,
                                    "MaxNumSteps is negative.", NULL);
  }

  /* Maximum order */
  
  opt = mxGetField(options,0,"MaxOrder");
  if ( !mxIsEmpty(opt) ) {
    q = (int)*mxGetPr(opt);
    if (q <= 0) cvmErrHandler(-999, "CVODES", fctName,
                              "MaxOrder must be positive.", NULL);
    if (q > *maxord) cvmErrHandler(-999, "CVODES", fctName,
                                   "MaxOrder is too large for the Method specified.", NULL);
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
    if (*tmp < 0.0) cvmErrHandler(-999, "CVODES", fctName,
                                  "MaxStep is negative.", NULL);
    if ( mxIsInf(*tmp) )
      *hmax = 0.0;
    else
      *hmax = *tmp;
  }

  /* Minimum step size */

  opt = mxGetField(options,0,"MinStep");
  if ( !mxIsEmpty(opt) ) {
    *hmin = *mxGetPr(opt);
    if (*hmin < 0.0) cvmErrHandler(-999, "CVODES", fctName,
                                   "MinStep is negative.", NULL);
  }

  /* Stopping time */

  opt = mxGetField(options,0,"StopTime");
  if ( !mxIsEmpty(opt) ) {
    *tstop = *mxGetPr(opt);
    tstopSet = TRUE;
  }

  /* Stability Limit Detection */

  opt = mxGetField(options,0,"StabilityLimDet");
  if ( !mxIsEmpty(opt) ) {
    buflen = mxGetM(opt) * mxGetN(opt) + 1;
    bufval = mxCalloc(buflen, sizeof(char));
    status = mxGetString(opt, bufval, buflen);
    if(status != 0) cvmErrHandler(-999, "CVODES", fctName,
                                  "Cannot parse StabilityLimDet.", NULL);
    if(!strcmp(bufval,"on")) *sld = TRUE;
    else if(!strcmp(bufval,"off")) *sld = FALSE;
    else cvmErrHandler(-999, "CVODES", fctName,
                       "StabilityLimDet has an illegal value.", NULL);
  }

  /* Monitor? */

  opt = mxGetField(options,0,"MonitorFn");
  if ( !mxIsEmpty(opt) ) {
    mon = TRUE;
    mxDestroyArray(mx_MONfct);
    mx_MONfct = mxDuplicateArray(opt);
    opt = mxGetField(options,0,"MonitorData");
    if ( !mxIsEmpty(opt) ) {
      mxDestroyArray(mx_MONdata);
      mx_MONdata  = mxDuplicateArray(opt);
    }
  }

  /* The remaining options are interpreted only for forward phase */

  if (!fwd) return;

  /* Number of root functions */
  opt = mxGetField(options,0,"NumRoots");
  if ( !mxIsEmpty(opt) ) {

    Ng = (int)*mxGetPr(opt);
    if (Ng < 0) cvmErrHandler(-999, "CVODES", fctName,
                              "NumRoots is negative.", NULL);
    if (Ng > 0) {
      /* Roots function */
      opt = mxGetField(options,0,"RootsFn");
      if ( !mxIsEmpty(opt) ) {
        mxDestroyArray(mx_Gfct);
        mx_Gfct = mxDuplicateArray(opt);
      } else {
        cvmErrHandler(-999, "CVODES", fctName,
                      "RootsFn required for NumRoots > 0", NULL);
      }
    }

  }

  /* We made it here without problems */

  return;
}


void get_LinSolvOptions(const mxArray *options, cvmPbData thisPb, booleantype fwd,
                        int *mupper, int *mlower,
                        int *mudq, int *mldq, double *dqrely,
                        int *ptype, int *gstype, int *maxl)
{
  mxArray *opt;
  char *bufval;
  int buflen, status, tmp_ls, tmp_pm;
  char *fctName;
  char *fwd_fctName = "CVodeInit/CVodeReInit";
  char *bck_fctName = "CVodeInitB/CVodeReInitB";

  if (fwd) fctName = fwd_fctName;
  else     fctName = bck_fctName;

  *mupper = 0;
  *mlower = 0;

  *mudq = 0;
  *mldq = 0;
  *dqrely = 0.0;

  *ptype = PREC_NONE;
  *gstype = MODIFIED_GS;
  *maxl = 0;

  /* Return now if options was empty */

  if (mxIsEmpty(options)) return;

  /* Linear solver type */

  opt = mxGetField(options,0,"LinearSolver");
  if ( !mxIsEmpty(opt) ) {
    buflen = mxGetM(opt) * mxGetN(opt) + 1;
    bufval = mxCalloc(buflen, sizeof(char));
    status = mxGetString(opt, bufval, buflen);
    if(status != 0) cvmErrHandler(-999, "CVODES", fctName,
                                  "Cannot parse LinearSolver.", NULL);
    if(!strcmp(bufval,"Diag"))          ls = LS_DIAG;
    else if(!strcmp(bufval,"Band"))     ls = LS_BAND;
    else if(!strcmp(bufval,"GMRES"))    ls = LS_SPGMR;
    else if(!strcmp(bufval,"BiCGStab")) ls = LS_SPBCG;
    else if(!strcmp(bufval,"TFQMR"))    ls = LS_SPTFQMR;
    else if(!strcmp(bufval,"Dense"))    ls = LS_DENSE;
    else cvmErrHandler(-999, "CVODES", fctName,
                       "LinearSolver has an illegal value.", NULL);
  }
  
  /* Jacobian function */

  opt = mxGetField(options,0,"JacobianFn");
  if ( !mxIsEmpty(opt) ) {
    mxDestroyArray(mx_JACfct);
    mx_JACfct  = mxDuplicateArray(opt);
  }

  /* Band linear solver */

  if (tmp_ls==LS_BAND) {

    opt = mxGetField(options,0,"UpperBwidth");
    if ( !mxIsEmpty(opt) )
      *mupper = (int)*mxGetPr(opt);
    
    opt = mxGetField(options,0,"LowerBwidth");
    if ( !mxIsEmpty(opt) )
      *mlower = (int)*mxGetPr(opt);

  }

  /* SPGMR linear solver options */
  
  if (tmp_ls==LS_SPGMR) {

    /* Type of Gram-Schmidt procedure */

    opt = mxGetField(options,0,"GramSchmidtType");
    if ( !mxIsEmpty(opt) ) {
      buflen = mxGetM(opt) * mxGetN(opt) + 1;
      bufval = mxCalloc(buflen, sizeof(char));
      status = mxGetString(opt, bufval, buflen);
      if(status != 0) cvmErrHandler(-999, "CVODES", fctName,
                                    "Cannot parse GramSchmidtType.", NULL);
      if(!strcmp(bufval,"Classical"))     *gstype = CLASSICAL_GS;
      else if(!strcmp(bufval,"Modified")) *gstype = MODIFIED_GS;
      else cvmErrHandler(-999, "CVODES", fctName,
                         "GramSchmidtType has an illegal value.", NULL);
    }

  }

  /* SPILS linear solver options */

  if ( (tmp_ls==LS_SPGMR) || (tmp_ls==LS_SPBCG) || (tmp_ls==LS_SPTFQMR) ) {

    /* Max. dimension of Krylov subspace */

    opt = mxGetField(options,0,"KrylovMaxDim");
    if ( !mxIsEmpty(opt) ) {
      *maxl = (int)*mxGetPr(opt);
      if (*maxl < 0) cvmErrHandler(-999, "CVODES", fctName,
                                   "KrylovMaxDim is negative.", NULL);
    }

    /* Preconditioning type */

    opt = mxGetField(options,0,"PrecType");
    if ( !mxIsEmpty(opt) ) {
      buflen = mxGetM(opt) * mxGetN(opt) + 1;
      bufval = mxCalloc(buflen, sizeof(char));
      status = mxGetString(opt, bufval, buflen);
      if(status != 0) cvmErrHandler(-999, "CVODES", fctName,
                                    "Cannot parse PrecType.", NULL);
      if(!strcmp(bufval,"Left")) *ptype = PREC_LEFT;
      else if(!strcmp(bufval,"Right")) *ptype = PREC_RIGHT;
      else if(!strcmp(bufval,"Both")) *ptype = PREC_BOTH;
      else if(!strcmp(bufval,"None")) *ptype = PREC_NONE;
      else cvmErrHandler(-999, "CVODES", fctName,
                         "PrecType has an illegal value.", NULL);
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
      if(status != 0) cvmErrHandler(-999, "CVODES", fctName,
                                    "Cannot parse PrecModule.", NULL);
      if(!strcmp(bufval,"BandPre"))          pm = PM_BANDPRE;
      else if(!strcmp(bufval,"BBDPre"))      pm = PM_BBDPRE;
      else if(!strcmp(bufval,"UserDefined")) tmp_pm = PM_NONE;
      else cvmErrHandler(-999, "CVODES", fctName,
                         "PrecModule has an illegal value.", NULL);
    }

    if (tmp_pm != PM_NONE) {
    
      opt = mxGetField(options,0,"UpperBwidth");
      if ( !mxIsEmpty(opt) )
        *mupper = (int)*mxGetPr(opt);
      
      opt = mxGetField(options,0,"LowerBwidth");
      if ( !mxIsEmpty(opt) )
        *mlower = (int)*mxGetPr(opt);
      
    }

    if (tmp_pm == PM_BBDPRE) {
      
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
      } else { 
        cvmErrHandler(-999, "CVODES", fctName,
                      "GlocalFn required for BBD preconditioner.", NULL);
      }      

      opt = mxGetField(options,0,"GcommFn");
      if ( !mxIsEmpty(opt) ) {
        mxDestroyArray(mx_GCOMfct);
        mx_GCOMfct  = mxDuplicateArray(opt);
      }

    }

  }

  
  /* We made it here without problems */

  return;

}


void get_QuadOptions(const mxArray *options, cvmPbData thisPb, booleantype fwd,
                     int Nq,
                     booleantype *errconQ,
                     int *itolQ, double *reltolQ, double *SabstolQ, double **VabstolQ)
{
  mxArray *opt;
  char *bufval;
  int i, buflen, status, m, n;
  double *tmp;
  char *fctName;
  char *fwd_fctName = "CVodeQuadInit/CVodeQuadReInit";
  char *bck_fctName = "CVodeQuadInitB/CVodeQuadReInitB";

  if (fwd) fctName = fwd_fctName;
  else     fctName = bck_fctName;

  *errconQ = FALSE;
  *itolQ = CV_SS;
  *reltolQ = 1.0e-4;
  *SabstolQ = 1.0e-6;
  *VabstolQ = NULL;

  /* Return now if options was empty */

  if (mxIsEmpty(options)) return;

  /* Quadrature error control and tolerances */

  opt = mxGetField(options,0,"ErrControl");
  if ( mxIsEmpty(opt) ) return;
  
  buflen = mxGetM(opt) * mxGetN(opt) + 1;
  bufval = mxCalloc(buflen, sizeof(char));
  status = mxGetString(opt, bufval, buflen);
  if(status != 0) cvmErrHandler(-999, "CVODES", fctName,
                                "Canot parse ErrControl.", NULL);

  if(strcmp(bufval,"off") == 0) return;

  if(strcmp(bufval,"on")  != 0) cvmErrHandler(-999, "CVODES", fctName,
                                              "ErrControl has an illegal value.", NULL);

  *errconQ = TRUE;

  opt = mxGetField(options,0,"RelTol");
  if ( !mxIsEmpty(opt) ) {
    *reltolQ = *mxGetPr(opt);
    if (*reltolQ < 0.0) cvmErrHandler(-999, "CVODES", fctName,
                                      "RelTol is negative.", NULL);
  } 

  opt = mxGetField(options,0,"AbsTol");
  if ( !mxIsEmpty(opt) ) {

    m = mxGetN(opt);
    n = mxGetM(opt);
    if ( (n != 1) && (m != 1) ) cvmErrHandler(-999, "CVODES", fctName,
                                              "AbsTol is not a scalar or a vector.", NULL);
    if ( m > n ) n = m;
    tmp = mxGetPr(opt);

    if (n == 1) {
      *itolQ = CV_SS;
      *SabstolQ = *tmp;
      if (*SabstolQ < 0.0) cvmErrHandler(-999, "CVODES", fctName,
                                         "AbsTol is negative.", NULL);
    } else if (n == Nq) {
      *itolQ = CV_SV;
      *VabstolQ = (double *)malloc(Nq*sizeof(double));
      for(i=0;i<Nq;i++) {
        (*VabstolQ)[i] = tmp[i];
        if (tmp[i] < 0.0) cvmErrHandler(-999, "CVODES", fctName,
                                        "AbsTol has a negative component.", NULL);
      }
    } else {
      cvmErrHandler(-999, "CVODES", fctName,
                    "AbsTol does not contain Nq elements.", NULL);
    }

  }

  /* We made it here without problems */

  return;
}

void get_FSAOptions(const mxArray *options, cvmPbData thisPb,
                    int *ism,
                    char **pfield_name, int **plist, double **pbar,
                    int *dqtype, double *rho,
                    booleantype *errconS, int *itolS, double *reltolS, 
                    double **SabstolS, double **VabstolS)
{
  mxArray *opt;
  char *bufval;
  int i, is, m, n, buflen, status, this_plist;
  double *tmp;

  /* Set default values */

  *ism = CV_STAGGERED;

  *dqtype = CV_CENTERED;
  *rho = 0.0;

  *errconS = TRUE;

  *itolS = CV_EE;
  *SabstolS = NULL;
  *VabstolS = NULL;

  *pfield_name = NULL;
  *plist = NULL;
  *pbar  = NULL;

  /* Return now if options was empty */

  if (mxIsEmpty(options)) return;

  /* Sensitivity method */

  opt = mxGetField(options,0,"method");
  if ( !mxIsEmpty(opt) ) {
  
    buflen = mxGetM(opt) * mxGetN(opt) + 1;
    bufval = mxCalloc(buflen, sizeof(char));
    status = mxGetString(opt, bufval, buflen);
    if(status != 0) cvmErrHandler(-999, "CVODES", "CVodeSensInit/CVodeSensReInit",
                                  "Could not parse method.", NULL);

    if(!strcmp(bufval,"Simultaneous"))   *ism = CV_SIMULTANEOUS;
    else if(!strcmp(bufval,"Staggered")) *ism = CV_STAGGERED;
    else cvmErrHandler(-999, "CVODES", "CVodeSensInit/CVodeSensReInit",
                       "method has an illegal value.", NULL);

  }


  /* Field name in data structure for params. */

  opt = mxGetField(options,0,"ParamField");
  if ( !mxIsEmpty(opt) ) {

    buflen = mxGetM(opt) * mxGetN(opt) + 1;
    bufval = mxCalloc(buflen, sizeof(char));
    status = mxGetString(opt, bufval, buflen);
    if(status != 0) cvmErrHandler(-999, "CVODES", "CVodeSensInit/CVodeSensReInit",
                                  "Could not parse ParamField.", NULL);

    *pfield_name = mxCalloc(buflen, sizeof(char));
    strcpy((*pfield_name), bufval);

  }  

  /* PLIST */

  opt = mxGetField(options,0,"ParamList");
  if ( !mxIsEmpty(opt) ) {

    tmp = mxGetPr(opt);
    m = mxGetM(opt);
    n = mxGetN(opt);
    if ( (n != 1) && (m != 1) ) cvmErrHandler(-999, "CVODES", "CVodeSensInit/CVodeSensReInit",
                                              "ParamList is not a vector.", NULL);
    if (m > n) n = m;
    if ( n != Ns) cvmErrHandler(-999, "CVODES", "CVodeSensInit/CVodeSensReInit",
                                "ParamList does not contain Ns elements.", NULL);
    *plist = (int *) malloc(Ns*sizeof(int));
    for (is=0;is<Ns;is++) {
      this_plist = (int) tmp[is];
      if (this_plist <= 0) cvmErrHandler(-999, "CVODES", "CVodeSensInit/CVodeSensReInit",
                                         "ParamList must contain only positive integers.", NULL);
      (*plist)[is] = this_plist - 1;
    }

  }

  /* PBAR */

  opt = mxGetField(options,0,"ParamScales");
  if ( !mxIsEmpty(opt) ) {

    m = mxGetM(opt);
    n = mxGetN(opt);
    if ( (n != 1) && (m != 1) ) cvmErrHandler(-999, "CVODES", "CVodeSensInit/CVodeSensReInit",
                                              "ParamScales is not a vector.", NULL);
    if ( m > n ) n = m;
    if ( n != Ns) cvmErrHandler(-999, "CVODES", "CVodeSensInit/CVodeSensReInit",
                                "ParamScales does not contain Ns elements.", NULL);
    tmp = mxGetPr(opt);
    *pbar = (double *) malloc(Ns*sizeof(double));
    for(i=0;i<Ns;i++)
      (*pbar)[i] = tmp[i];

  }

  /* DQ type */

  opt = mxGetField(options,0,"DQtype");
  if ( !mxIsEmpty(opt) ) {

    buflen = mxGetM(opt) * mxGetN(opt) + 1;
    bufval = mxCalloc(buflen, sizeof(char));
    status = mxGetString(opt, bufval, buflen);
    if(status != 0) cvmErrHandler(-999, "CVODES", "CVodeSensInit/CVodeSensReInit",
                                  "Cannot parse DQtype.", NULL);

    if(!strcmp(bufval,"Centered")) *dqtype = CV_CENTERED;
    else if(!strcmp(bufval,"Forward")) *dqtype = CV_FORWARD;
    else cvmErrHandler(-999, "CVODES", "CVodeSensInit/CVodeSensReInit",
                       "DQtype has an illegal value.", NULL);

  }
  
  /* DQ parameter */

  opt = mxGetField(options,0,"DQparam");
  if ( !mxIsEmpty(opt) )
    *rho = *mxGetPr(opt);

  /* Error control */

  opt = mxGetField(options,0,"ErrControl");
  if ( !mxIsEmpty(opt) ) {

    buflen = mxGetM(opt) * mxGetN(opt) + 1;
    bufval = mxCalloc(buflen, sizeof(char));
    status = mxGetString(opt, bufval, buflen);
    if(status != 0) cvmErrHandler(-999, "CVODES", "CVodeSensInit/CVodeSensReInit",
                                  "Canot parse ErrControl.", NULL);

    if(!strcmp(bufval,"off")) *errconS = FALSE;
    else if(!strcmp(bufval,"on")) *errconS = TRUE;
    else cvmErrHandler(-999, "CVODES", "CVodeSensInit/CVodeSensReInit",
                       "ErrControl has an illegal value.", NULL);

  }

  /* Tolerances */
  
  opt = mxGetField(options,0,"RelTol");
  if ( !mxIsEmpty(opt) ) {

    *reltolS = *mxGetPr(opt);
    if (*reltolS < 0.0) cvmErrHandler(-999, "CVODES", "CVodeSensInit/CVodeSensReInit",
                                      "RelTol is negative.", NULL);

    opt = mxGetField(options,0,"AbsTol");
    if ( !mxIsEmpty(opt) ) {

      m = mxGetM(opt);
      n = mxGetN(opt);
      if ( (m == 1) && (n == Ns) ) {
        *itolS = CV_SS;
        tmp = mxGetPr(opt);
        *SabstolS = (double *) malloc(Ns*sizeof(double));
        for (is=0; is<Ns; is++) {
          (*SabstolS)[is] = tmp[is];
          if ( tmp[is] < 0.0 ) cvmErrHandler(-999, "CVODES", "CVodeSensInit/CVodeSensReInit",
                                             "AbsTol has a negative component.", NULL);
        }
      } else if ( (m == N) && (n == Ns) ) {
        *itolS = CV_SV;
        tmp = mxGetPr(opt);
        *VabstolS = (double *)malloc(Ns*N*sizeof(double));
        for (i=0; i<Ns*N; i++) {
          (*VabstolS)[i] = tmp[i];
          if ( tmp[i] < 0.0 ) cvmErrHandler(-999, "CVODES", "CVodeSensInit/CVodeSensReInit",
                                            "AbsTol has a negative component.", NULL);
        }
      } else {
        cvmErrHandler(-999, "CVODES", "CVodeSensInit/CVodeSensReInit",
                      "AbsTol must be either a 1xNs vector or a NxNs matrix.", NULL);
      }

    } else {

      *itolS = CV_EE;

    }

  }

  /* We made it here without problems */

  return;

}

