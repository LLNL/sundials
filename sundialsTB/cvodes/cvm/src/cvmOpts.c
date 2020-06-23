/*
 * -----------------------------------------------------------------
 * $Revision$
 * $Date$
 * -----------------------------------------------------------------
 * Programmer: Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2020, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * Option parsing functions for the CVODES Matlab interface.
 * -----------------------------------------------------------------
 */

#include <string.h>
#include "cvm.h"

/*
 * ---------------------------------------------------------------------------------
 * Redability replacements
 * ---------------------------------------------------------------------------------
 */

#define N             (thisPb->n) 
#define Ns            (thisPb->ns) 
#define Ng            (thisPb->ng) 
#define ls            (thisPb->LS) 
#define pm            (thisPb->PM) 

#define mtlb_data     (thisPb->mtlb_data)

#define mtlb_JACfct   (thisPb->JACfct)
#define mtlb_PSETfct  (thisPb->PSETfct)
#define mtlb_PSOLfct  (thisPb->PSOLfct)
#define mtlb_GLOCfct  (thisPb->GLOCfct)
#define mtlb_GCOMfct  (thisPb->GCOMfct)
#define mtlb_Gfct     (thisPb->Gfct)

#define mon           (thisPb->Mon)
#define tstopSet      (thisPb->TstopSet)

#define mtlb_MONfct   (thisPb->MONfct)
#define mtlb_MONdata  (thisPb->MONdata)

/*
 * ---------------------------------------------------------------------------------
 * Option handling functions
 * ---------------------------------------------------------------------------------
 */

int get_IntgrOptions(const mxArray *options, cvmPbData thisPb, booleantype fwd, int lmm,
                     int *maxord, booleantype *sld, booleantype *errmsg,
                     sunindextype *mxsteps,
                     int *itol, realtype *reltol, double *Sabstol, double **Vabstol,
                     double *hin, double *hmax, double *hmin, double *tstop,
                     booleantype *rhs_s)
{
  mxArray *opt;
  int q;
  sunindextype i, m, n;
  double *tmp;
  char *fctName;
  char *fwd_fctName = "CVodeInit/CVodeReInit";
  char *bck_fctName = "CVodeInitB/CVodeReInitB";

  if (fwd) fctName = fwd_fctName;
  else     fctName = bck_fctName;
  
  /* Set default values */
  
  *maxord = (lmm == CV_ADAMS) ? 12 : 5;
  
  *sld = SUNFALSE;

  *mxsteps = 0;

  *itol = CV_SS;
  *reltol = 1.0e-3;
  *Sabstol = 1.0e-6;
  *Vabstol = NULL;

  *hin = 0.0;
  *hmax = 0.0;
  *hmin = 0.0;

  *rhs_s = SUNFALSE;

  Ng = 0;
  tstopSet = SUNFALSE;
  mon = SUNFALSE;

  *errmsg = SUNTRUE;


  /* Return now if options was empty */

  if (mxIsEmpty(options)) return(0);

  /* User data */

  opt = mxGetField(options,0,"UserData");
  if ( !mxIsEmpty(opt) ) {
    mxDestroyArray(mtlb_data);
    mtlb_data = mxDuplicateArray(opt);
  }
  
  /* Tolerances */

  opt = mxGetField(options,0,"RelTol");
  if ( !mxIsEmpty(opt) ) {
    *reltol = *mxGetPr(opt);
    if (*reltol < 0.0 ) {
      cvmErrHandler(-999, "CVODES", fctName, "RelTol is negative.", NULL);
      return(-1);
    }
  }
    
  opt = mxGetField(options,0,"AbsTol");
  if ( !mxIsEmpty(opt) ) {
    m = mxGetM(opt);
    n = mxGetN(opt);
    if ( (n != 1) && (m != 1) ) {
      cvmErrHandler(-999, "CVODES", fctName, "AbsTol is not a scalar or a vector.", NULL);
      return(-1);
    }
    if ( m > n ) n = m;
    tmp = mxGetPr(opt);
    if (n == 1) {
      *itol = CV_SS;
      *Sabstol = *tmp;
      if (*Sabstol < 0.0) {
        cvmErrHandler(-999, "CVODES", fctName, "AbsTol is negative.", NULL);
      return(-1);
      }
    } else if (n == N) {
      *itol = CV_SV;
      *Vabstol = (double *) malloc(N*sizeof(double));
      for(i=0;i<N;i++) {
        (*Vabstol)[i] = tmp[i];
        if (tmp[i] < 0.0) {
          cvmErrHandler(-999, "CVODES", fctName, "AbsTol has a negative component.", NULL);
          return(-1);
        }
      }
    } else {
      cvmErrHandler(-999, "CVODES", fctName, "AbsTol does not contain N elements.", NULL);
      return(-1);
    }
  }

  /* Maximum number of steps */

  opt = mxGetField(options,0,"MaxNumSteps");
  if ( !mxIsEmpty(opt) ) {
    *mxsteps = (int)*mxGetPr(opt);
    if (*mxsteps < 0) {
      cvmErrHandler(-999, "CVODES", fctName, "MaxNumSteps is negative.", NULL);
      return(-1);
    }
  }

  /* Maximum order */

  opt = mxGetField(options,0,"MaxOrder");
  if ( !mxIsEmpty(opt) ) {
    q = (int)*mxGetPr(opt);
    if (q <= 0) {
      cvmErrHandler(-999, "CVODES", fctName, "MaxOrder must be positive.", NULL);
      return(-1);
    }
    if (q > *maxord) {
      cvmErrHandler(-999, "CVODES", fctName, "MaxOrder is too large for the Method specified.", NULL);
      return(-1);
    }
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
    if (*tmp < 0.0) {
      cvmErrHandler(-999, "CVODES", fctName, "MaxStep is negative.", NULL);
      return(-1);
    }
    if ( mxIsInf(*tmp) ) *hmax = 0.0;
    else                 *hmax = *tmp;
  }

  /* Minimum step size */

  opt = mxGetField(options,0,"MinStep");
  if ( !mxIsEmpty(opt) ) {
    *hmin = *mxGetPr(opt);
    if (*hmin < 0.0) {
      cvmErrHandler(-999, "CVODES", fctName, "MinStep is negative.", NULL);
      return(-1);
    }
  }

  /* Stability Limit Detection */

  opt = mxGetField(options,0,"StabilityLimDet");
  if ( !mxIsEmpty(opt) ) {
    if (!mxIsLogical(opt)) {
      cvmErrHandler(-999, "CVODES", fctName, "StabilityLimDet is not a logical scalar.", NULL);
      return(-1);
    }
    if (mxIsLogicalScalarTrue(opt)) *sld = SUNTRUE;
    else                            *sld = SUNFALSE;
  }

  /* Monitor? */

  opt = mxGetField(options,0,"MonitorFn");
  if ( !mxIsEmpty(opt) ) {
    mon = SUNTRUE;
    mxDestroyArray(mtlb_MONfct);
    mtlb_MONfct = mxDuplicateArray(opt);
    opt = mxGetField(options,0,"MonitorData");
    if ( !mxIsEmpty(opt) ) {
      mxDestroyArray(mtlb_MONdata);
      mtlb_MONdata  = mxDuplicateArray(opt);
    }
  }

  /* The remaining options are interpreted either for 
   * forward problems only or backward problems only */

  if (fwd) {   /* FORWARD PROBLEM ONLY */

    /* Disable error/warning messages? */

    opt = mxGetField(options,0,"ErrorMessages");
    if ( !mxIsEmpty(opt) ) {
      if (!mxIsLogical(opt)) {
        cvmErrHandler(-999, "CVODES", fctName, "ErrorMessages is not a logical scalar.", NULL);
        return(-1);
      }
      if (mxIsLogicalScalarTrue(opt)) *errmsg = SUNTRUE;
      else                            *errmsg = SUNFALSE;
    }

    /* Stopping time */
    opt = mxGetField(options,0,"StopTime");
    if ( !mxIsEmpty(opt) ) {
      *tstop = *mxGetPr(opt);
      tstopSet = SUNTRUE;
    }

    /* Number of root functions */
    opt = mxGetField(options,0,"NumRoots");
    if ( !mxIsEmpty(opt) ) {

      Ng = (int)*mxGetPr(opt);
      if (Ng < 0) {
        cvmErrHandler(-999, "CVODES", fctName, "NumRoots is negative.", NULL);
        return(-1);
      }
      if (Ng > 0) {
        /* Roots function */
        opt = mxGetField(options,0,"RootsFn");
        if ( !mxIsEmpty(opt) ) {
          mxDestroyArray(mtlb_Gfct);
          mtlb_Gfct = mxDuplicateArray(opt);
        } else {
          cvmErrHandler(-999, "CVODES", fctName, "RootsFn required for NumRoots > 0", NULL);
          return(-1);
        }
      }
      
    }

  } else {   /* BACKWARD PROBLEM ONLY */

    /* Dependency on forward sensitivities */

    opt = mxGetField(options,0,"SensDependent");
    if ( !mxIsEmpty(opt) ) {
      if (!mxIsLogical(opt)) {
        cvmErrHandler(-999, "CVODES", fctName, "SensDependent is not a logical scalar.", NULL);
        return(-1);
      }
      if (mxIsLogicalScalarTrue(opt)) *rhs_s = SUNTRUE;
      else                            *rhs_s = SUNFALSE;
    }

  }

  /* We made it here without problems */

  return(0);
}


int get_LinSolvOptions(const mxArray *options, cvmPbData thisPb, booleantype fwd,
                       sunindextype *mupper, sunindextype *mlower,
                       sunindextype *mudq, sunindextype *mldq, double *dqrely,
                       int *ptype, int *gstype, int *maxl)
{
  mxArray *opt;
  char *bufval;
  int buflen, status;
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

  if (mxIsEmpty(options)) return(0);

  /* Linear solver type */

  opt = mxGetField(options,0,"LinearSolver");
  if ( !mxIsEmpty(opt) ) {
    buflen = mxGetM(opt) * mxGetN(opt) + 1;
    bufval = mxCalloc(buflen, sizeof(char));
    status = mxGetString(opt, bufval, buflen);
    if(status != 0) {
      cvmErrHandler(-999, "CVODES", fctName, "Cannot parse LinearSolver.", NULL);
      return(-1);
    }
    if(!strcmp(bufval,"Diag"))          ls = LS_DIAG;
    else if(!strcmp(bufval,"Band"))     ls = LS_BAND;
    else if(!strcmp(bufval,"GMRES"))    ls = LS_SPGMR;
    else if(!strcmp(bufval,"BiCGStab")) ls = LS_SPBCG;
    else if(!strcmp(bufval,"TFQMR"))    ls = LS_SPTFQMR;
    else if(!strcmp(bufval,"Dense"))    ls = LS_DENSE;
    else {
      cvmErrHandler(-999, "CVODES", fctName, "LinearSolver has an illegal value.", NULL);
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

  /* SPGMR linear solver options */
  
  if (ls==LS_SPGMR) {

    /* Type of Gram-Schmidt procedure */

    opt = mxGetField(options,0,"GramSchmidtType");
    if ( !mxIsEmpty(opt) ) {
      buflen = mxGetM(opt) * mxGetN(opt) + 1;
      bufval = mxCalloc(buflen, sizeof(char));
      status = mxGetString(opt, bufval, buflen);
      if(status != 0) {
        cvmErrHandler(-999, "CVODES", fctName, "Cannot parse GramSchmidtType.", NULL);
        return(-1);
      }
      if(!strcmp(bufval,"Classical"))     *gstype = CLASSICAL_GS;
      else if(!strcmp(bufval,"Modified")) *gstype = MODIFIED_GS;
      else {
        cvmErrHandler(-999, "CVODES", fctName, "GramSchmidtType has an illegal value.", NULL);
        return(-1);
      }
    }

  }

  /* SPILS linear solver options */

  if ( (ls==LS_SPGMR) || (ls==LS_SPBCG) || (ls==LS_SPTFQMR) ) {

    /* Max. dimension of Krylov subspace */

    opt = mxGetField(options,0,"KrylovMaxDim");
    if ( !mxIsEmpty(opt) ) {
      *maxl = (int)*mxGetPr(opt);
      if (*maxl < 0) {
        cvmErrHandler(-999, "CVODES", fctName, "KrylovMaxDim is negative.", NULL);
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
        cvmErrHandler(-999, "CVODES", fctName, "Cannot parse PrecType.", NULL);
        return(-1);
      }
      if(!strcmp(bufval,"Left")) *ptype = PREC_LEFT;
      else if(!strcmp(bufval,"Right")) *ptype = PREC_RIGHT;
      else if(!strcmp(bufval,"Both"))  *ptype = PREC_BOTH;
      else if(!strcmp(bufval,"None"))  *ptype = PREC_NONE;
      else {
        cvmErrHandler(-999, "CVODES", fctName, "PrecType has an illegal value.", NULL);
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
        cvmErrHandler(-999, "CVODES", fctName, "Cannot parse PrecModule.", NULL);
        return(-1);
      }
      if(!strcmp(bufval,"BandPre"))          pm = PM_BANDPRE;
      else if(!strcmp(bufval,"BBDPre"))      pm = PM_BBDPRE;
      else if(!strcmp(bufval,"UserDefined")) pm = PM_NONE;
      else {
        cvmErrHandler(-999, "CVODES", fctName, "PrecModule has an illegal value.", NULL);
        return(-1);
      }
    }

    if (pm != PM_NONE) {
    
      opt = mxGetField(options,0,"UpperBwidth");
      if ( !mxIsEmpty(opt) )
        *mupper = (sunindextype)*mxGetPr(opt);
      
      opt = mxGetField(options,0,"LowerBwidth");
      if ( !mxIsEmpty(opt) )
        *mlower = (sunindextype)*mxGetPr(opt);
      
    }

    if (pm == PM_BBDPRE) {
      
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
      } else { 
        cvmErrHandler(-999, "CVODES", fctName, "GlocalFn required for BBD preconditioner.", NULL);
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


int get_QuadOptions(const mxArray *options, cvmPbData thisPb, booleantype fwd,
                    sunindextype Nq, booleantype *rhs_s,
                    booleantype *errconQ,
                    int *itolQ, double *reltolQ, double *SabstolQ, double **VabstolQ)
{
  mxArray *opt;
  sunindextype i, m, n;
  double *tmp;
  char *fctName;
  char *fwd_fctName = "CVodeQuadInit/CVodeQuadReInit";
  char *bck_fctName = "CVodeQuadInitB/CVodeQuadReInitB";

  if (fwd) fctName = fwd_fctName;
  else     fctName = bck_fctName;

  *errconQ = SUNFALSE;
  *itolQ = CV_SS;
  *reltolQ = 1.0e-4;
  *SabstolQ = 1.0e-6;
  *VabstolQ = NULL;

  *rhs_s = SUNFALSE;

  /* Return now if options was empty */

  if (mxIsEmpty(options)) return(0);

  /* For backward problems only, check dependency on forward sensitivities */

  if (!fwd) {

    opt = mxGetField(options,0,"SensDependent");
    if ( !mxIsEmpty(opt) ) {
      if (!mxIsLogical(opt)) {
        cvmErrHandler(-999, "CVODES", fctName, "SensDependent is not a logical scalar.", NULL);
        return(-1);
      }
      if (mxIsLogicalScalarTrue(opt)) *rhs_s = SUNTRUE;
      else                            *rhs_s = SUNFALSE;
    }

  }

  /* Quadrature error control and tolerances */

  opt = mxGetField(options,0,"ErrControl");

  if ( mxIsEmpty(opt) ) return(0);

  if (!mxIsLogical(opt)) {
    cvmErrHandler(-999, "CVODES", fctName, "ErrControl is not a logical scalar.", NULL);
    return(-1);
  }

  if (!mxIsLogicalScalarTrue(opt)) return(0);
  
  /* the remining options are interpreted only if quadratures are included in error control */

  *errconQ = SUNTRUE;

  opt = mxGetField(options,0,"RelTol");
  if ( !mxIsEmpty(opt) ) {
    *reltolQ = *mxGetPr(opt);
    if (*reltolQ < 0.0) {
      cvmErrHandler(-999, "CVODES", fctName, "RelTol is negative.", NULL);
      return(-1);
    }
  } 

  opt = mxGetField(options,0,"AbsTol");
  if ( !mxIsEmpty(opt) ) {

    m = mxGetN(opt);
    n = mxGetM(opt);
    if ( (n != 1) && (m != 1) ) {
      cvmErrHandler(-999, "CVODES", fctName, "AbsTol is not a scalar or a vector.", NULL);
      return(-1);
    }
    if ( m > n ) n = m;
    tmp = mxGetPr(opt);

    if (n == 1) {
      *itolQ = CV_SS;
      *SabstolQ = *tmp;
      if (*SabstolQ < 0.0) {
        cvmErrHandler(-999, "CVODES", fctName, "AbsTol is negative.", NULL);
        return(-1);
      }
    } else if (n == Nq) {
      *itolQ = CV_SV;
      *VabstolQ = (double *)malloc(Nq*sizeof(double));
      for(i=0;i<Nq;i++) {
        (*VabstolQ)[i] = tmp[i];
        if (tmp[i] < 0.0) {
          cvmErrHandler(-999, "CVODES", fctName, "AbsTol has a negative component.", NULL);
          return(-1);
        }
      }
    } else {
      cvmErrHandler(-999, "CVODES", fctName, "AbsTol does not contain Nq elements.", NULL);
      return(-1);
    }

  }

  /* We made it here without problems */

  return(0);
}

int get_FSAOptions(const mxArray *options, cvmPbData thisPb,
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

  *errconS = SUNTRUE;

  *itolS = CV_EE;
  *SabstolS = NULL;
  *VabstolS = NULL;

  *pfield_name = NULL;
  *plist = NULL;
  *pbar  = NULL;

  /* Return now if options was empty */

  if (mxIsEmpty(options)) return(0);

  /* Sensitivity method */

  opt = mxGetField(options,0,"method");
  if ( !mxIsEmpty(opt) ) {
  
    buflen = mxGetM(opt) * mxGetN(opt) + 1;
    bufval = mxCalloc(buflen, sizeof(char));
    status = mxGetString(opt, bufval, buflen);
    if(status != 0) {
      cvmErrHandler(-999, "CVODES", "CVodeSensInit/CVodeSensReInit", "Could not parse method.", NULL);
      return(-1);
    }
    if(!strcmp(bufval,"Simultaneous"))   *ism = CV_SIMULTANEOUS;
    else if(!strcmp(bufval,"Staggered")) *ism = CV_STAGGERED;
    else {
      cvmErrHandler(-999, "CVODES", "CVodeSensInit/CVodeSensReInit", "method has an illegal value.", NULL);
      return(-1);
    }
  }


  /* Field name in data structure for params. */

  opt = mxGetField(options,0,"ParamField");
  if ( !mxIsEmpty(opt) ) {

    buflen = mxGetM(opt) * mxGetN(opt) + 1;
    bufval = mxCalloc(buflen, sizeof(char));
    status = mxGetString(opt, bufval, buflen);
    if(status != 0) {
      cvmErrHandler(-999, "CVODES", "CVodeSensInit/CVodeSensReInit", "Could not parse ParamField.", NULL);
      return(-1);
    }
    *pfield_name = mxCalloc(buflen, sizeof(char));
    strcpy((*pfield_name), bufval);

  }  

  /* PLIST */

  opt = mxGetField(options,0,"ParamList");
  if ( !mxIsEmpty(opt) ) {

    tmp = mxGetPr(opt);
    m = mxGetM(opt);
    n = mxGetN(opt);
    if ( (n != 1) && (m != 1) ) {
      cvmErrHandler(-999, "CVODES", "CVodeSensInit/CVodeSensReInit", "ParamList is not a vector.", NULL);
      return(-1);
    }
    if (m > n) n = m;
    if ( n != Ns) {
      cvmErrHandler(-999, "CVODES", "CVodeSensInit/CVodeSensReInit", "ParamList does not contain Ns elements.", NULL);
      return(-1);
    }
    *plist = (int *) malloc(Ns*sizeof(int));
    for (is=0;is<Ns;is++) {
      this_plist = (int) tmp[is];
      if (this_plist <= 0) {
        cvmErrHandler(-999, "CVODES", "CVodeSensInit/CVodeSensReInit", "ParamList must contain only positive integers.", NULL);
        return(-1);
      }
      (*plist)[is] = this_plist - 1;
    }

  }

  /* PBAR */

  opt = mxGetField(options,0,"ParamScales");
  if ( !mxIsEmpty(opt) ) {

    m = mxGetM(opt);
    n = mxGetN(opt);
    if ( (n != 1) && (m != 1) ) {
      cvmErrHandler(-999, "CVODES", "CVodeSensInit/CVodeSensReInit", "ParamScales is not a vector.", NULL);
      return(-1);
    }
    if ( m > n ) n = m;
    if ( n != Ns) {
      cvmErrHandler(-999, "CVODES", "CVodeSensInit/CVodeSensReInit", "ParamScales does not contain Ns elements.", NULL);
      return(-1);
    }
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
    if(status != 0) {
      cvmErrHandler(-999, "CVODES", "CVodeSensInit/CVodeSensReInit", "Cannot parse DQtype.", NULL);
      return(-1);
    }
    if(!strcmp(bufval,"Centered")) *dqtype = CV_CENTERED;
    else if(!strcmp(bufval,"Forward")) *dqtype = CV_FORWARD;
    else {
      cvmErrHandler(-999, "CVODES", "CVodeSensInit/CVodeSensReInit", "DQtype has an illegal value.", NULL);
      return(-1);
    }
  }
  
  /* DQ parameter */

  opt = mxGetField(options,0,"DQparam");
  if ( !mxIsEmpty(opt) )
    *rho = *mxGetPr(opt);

  /* Error control */

  opt = mxGetField(options,0,"ErrControl");
  if ( !mxIsEmpty(opt) ) {
    if (!mxIsLogical(opt)) {
      cvmErrHandler(-999, "CVODES", "CVodeSensInit/CVodeSensReInit", "ErrControl is not a logical scalar.", NULL);
      return(-1);
    }
    if (mxIsLogicalScalarTrue(opt)) *errconS = SUNTRUE;
    else                            *errconS = SUNFALSE;
  }

  /* Tolerances */
  
  opt = mxGetField(options,0,"RelTol");
  if ( !mxIsEmpty(opt) ) {

    *reltolS = *mxGetPr(opt);
    if (*reltolS < 0.0) {
      cvmErrHandler(-999, "CVODES", "CVodeSensInit/CVodeSensReInit", "RelTol is negative.", NULL);
      return(-1);
    }
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
          if ( tmp[is] < 0.0 ) {
            cvmErrHandler(-999, "CVODES", "CVodeSensInit/CVodeSensReInit", "AbsTol has a negative component.", NULL);
            return(-1);
          }
        }
      } else if ( (m == N) && (n == Ns) ) {
        *itolS = CV_SV;
        tmp = mxGetPr(opt);
        *VabstolS = (double *)malloc(Ns*N*sizeof(double));
        for (i=0; i<Ns*N; i++) {
          (*VabstolS)[i] = tmp[i];
          if ( tmp[i] < 0.0 ) {
            cvmErrHandler(-999, "CVODES", "CVodeSensInit/CVodeSensReInit", "AbsTol has a negative component.", NULL);
            return(-1);
          }
        }
      } else {
        cvmErrHandler(-999, "CVODES", "CVodeSensInit/CVodeSensReInit", "AbsTol must be either a 1xNs vector or an NxNs matrix.", NULL);
        return(-1);
      }

    } else {

      *itolS = CV_EE;

    }

  }

  /* We made it here without problems */

  return(0);

}
