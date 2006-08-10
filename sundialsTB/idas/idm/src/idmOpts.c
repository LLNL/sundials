/*
 * -----------------------------------------------------------------
 * $Revision: 1.6 $
 * $Date: 2006-08-10 18:01:04 $
 * -----------------------------------------------------------------
 * Programmer: Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2005, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see sundials/idas/LICENSE.
 * -----------------------------------------------------------------
 * Option parsing functions for the IDAS Matlab interface.
 * -----------------------------------------------------------------
 */

#include <string.h>
#include "idm.h"

extern idm_IDASdata idm_Cdata;   /* IDAS data */
extern booleantype idm_quad;     /* Quadratures? */
extern booleantype idm_asa;      /* Adjoint sensitivity? */
extern booleantype idm_fsa;      /* forward sensitivity? */

extern idm_MATLABdata idm_Mdata; /* MATLAB data */

/* 
 * ---------------------------------------------------------------------------------
 * Private constants
 * ---------------------------------------------------------------------------------
 */

#define ONE RCONST(1.0)







#define IDA_CENTERED 1
#define IDA_FORWARD  2




/*
 * ---------------------------------------------------------------------------------
 * Redability replacements
 * ---------------------------------------------------------------------------------
 */

#define N           (idm_Cdata->N) 
#define Ns          (idm_Cdata->Ns) 
#define Nq          (idm_Cdata->Nq) 
#define NqB         (idm_Cdata->NqB) 
#define Ng          (idm_Cdata->Ng) 
#define ls          (idm_Cdata->ls) 
#define pm          (idm_Cdata->pm) 
#define lsB         (idm_Cdata->lsB) 
#define pmB         (idm_Cdata->pmB) 

#define mx_QUADfct  (idm_Mdata->mx_QUADfct)
#define mx_JACfct   (idm_Mdata->mx_JACfct)
#define mx_PSETfct  (idm_Mdata->mx_PSETfct)
#define mx_PSOLfct  (idm_Mdata->mx_PSOLfct)
#define mx_GLOCfct  (idm_Mdata->mx_GLOCfct)
#define mx_GCOMfct  (idm_Mdata->mx_GCOMfct)
#define mx_Gfct     (idm_Mdata->mx_Gfct)
#define mx_SRESfct  (idm_Mdata->mx_SRESfct)

#define mx_QUADfctB (idm_Mdata->mx_QUADfctB)
#define mx_JACfctB  (idm_Mdata->mx_JACfctB)
#define mx_PSETfctB (idm_Mdata->mx_PSETfctB)
#define mx_PSOLfctB (idm_Mdata->mx_PSOLfctB)
#define mx_GLOCfctB (idm_Mdata->mx_GLOCfctB)
#define mx_GCOMfctB (idm_Mdata->mx_GCOMfctB)

#define mx_MONfct   (idm_Mdata->mx_MONfct)
#define mx_MONdata  (idm_Mdata->mx_MONdata)

#define mx_MONfctB  (idm_Mdata->mx_MONfctB)
#define mx_MONdataB (idm_Mdata->mx_MONdataB)

/*
 * ---------------------------------------------------------------------------------
 * Option handling functions
 * ---------------------------------------------------------------------------------
 */


/* MUST ADD: suppressALG, ID, constraints */

int get_IntgrOptions(const mxArray *options, booleantype fwd,
                     int *maxord, long int *mxsteps,
                     int *itol, realtype *reltol, double *Sabstol, double **Vabstol,
                     double *hin, double *hmax, double *tstop, booleantype *tstopSet,
                     booleantype *suppress,
                     double **id, double **cnstr)
{
  mxArray *opt;
  char *bufval;
  int i, buflen, status, q, m, n;
  double *tmp;
  booleantype tmp_quad;

  /* Set default values */

  *maxord = 5;

  *mxsteps = 0;

  *itol = IDA_SS;
  *reltol = 1.0e-3;
  *Sabstol = 1.0e-6;
  *Vabstol = NULL;

  *hin = 0.0;
  *hmax = 0.0;
  *tstopSet = FALSE;

  *suppress = FALSE;

  *id = NULL;
  *cnstr = NULL;

  tmp_quad = FALSE;

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
      *itol = IDA_SS;
      *Sabstol = *tmp;
      if (*Sabstol < 0.0)
        mexErrMsgTxt("AbsTol is negative.");
    } else if (n == N) {
      *itol = IDA_SV;
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

  /* Stopping time */

  opt = mxGetField(options,0,"StopTime");
  if ( !mxIsEmpty(opt) ) {
    *tstop = *mxGetPr(opt);
    *tstopSet = TRUE;
  }

  /* ID vector */

  opt = mxGetField(options,0,"VariableTypes");
  if ( !mxIsEmpty(opt) ) {
    m = mxGetM(opt);
    n = mxGetN(opt);
    if ( (n != 1) && (m != 1) )
      mexErrMsgTxt("VariableTypes is not a vector.");
    if ( m > n ) n = m;
    if (n == N) {
      tmp = mxGetPr(opt);
      *id = (double *)malloc(N*sizeof(double));
      for(i=0;i<N;i++) 
        (*id)[i] = tmp[i];
    } else {
      mexErrMsgTxt("VariableTypes has wrong dimension.");
    }
  }

  /* Constraints vector */

  opt = mxGetField(options,0,"ConstraintTypes");
  if ( !mxIsEmpty(opt) ) {
    m = mxGetM(opt);
    n = mxGetN(opt);
    if ( (n != 1) && (m != 1) )
      mexErrMsgTxt("ConstraintTypes is not a vector.");
    if ( m > n ) n = m;
    if (n == N) {
      tmp = mxGetPr(opt);
      *cnstr = (double *)malloc(N*sizeof(double));
      for(i=0;i<N;i++) 
        (*cnstr)[i] = tmp[i];
    } else {
      mexErrMsgTxt("ConstraintTypes has wrong dimension.");
    }
  }

  /* Suppress algebraic variables? */
  opt = mxGetField(options,0,"SuppressAlgVars");
  if ( !mxIsEmpty(opt) ) {
    buflen = mxGetM(opt) * mxGetN(opt) + 1;
    bufval = mxCalloc(buflen, sizeof(char));
    status = mxGetString(opt, bufval, buflen);
    if(status != 0)
      mexErrMsgTxt("Canot parse SuppressAlgVars.");
    if(!strcmp(bufval,"on")) {
      *suppress = TRUE;
    } else if(!strcmp(bufval,"off")) {
      *suppress = FALSE;
    } else {
      mexErrMsgTxt("SuppressAlgVars has an illegal value.");
    }
  }

  /* Options interpreted only for forward phase */

  if (fwd) {

    /* Number of root functions */
    opt = mxGetField(options,0,"NumRoots");
    if ( !mxIsEmpty(opt) ) {
      Ng = (int)*mxGetPr(opt);
      if (Ng < 0)
        mexErrMsgTxt("NumRoots is negative.");
      if (Ng > 0) {
        /* Roots function */
        opt = mxGetField(options,0,"RootsFn");
        if ( !mxIsEmpty(opt) ) {
          mxDestroyArray(mx_Gfct);
          mx_Gfct = mxDuplicateArray(opt);
        } else {
          mexErrMsgTxt("RootsFn required for NumRoots > 0");
        }
      }
    }

  }

  /* Quadratures? */

  opt = mxGetField(options,0,"Quadratures");
  if ( !mxIsEmpty(opt) ) {
    buflen = mxGetM(opt) * mxGetN(opt) + 1;
    bufval = mxCalloc(buflen, sizeof(char));
    status = mxGetString(opt, bufval, buflen);
    if(status != 0)
      mexErrMsgTxt("Cannot parse Quadratures.");
    if(!strcmp(bufval,"on")) tmp_quad = TRUE;
    else if(!strcmp(bufval,"off")) tmp_quad = FALSE;
    else mexErrMsgTxt("Quadratures has an illegal value.");
  }

  if (fwd) idm_quad  = tmp_quad;
  else     idm_quadB = tmp_quad;

  /* Monitor? */

  opt = mxGetField(options,0,"MonitorFn");
  if ( !mxIsEmpty(opt) ) {
    if (fwd) {
      idm_mon = TRUE;
      mxDestroyArray(mx_MONfct);
      mx_MONfct = mxDuplicateArray(opt);
    } else {
      idm_monB = TRUE;
      mxDestroyArray(mx_MONfctB);
      mx_MONfctB = mxDuplicateArray(opt);
    }
    opt = mxGetField(options,0,"MonitorData");
    if ( !mxIsEmpty(opt) ) {
      if (fwd) {
        mxDestroyArray(mx_MONdata);
        mx_MONdata  = mxDuplicateArray(opt);
      } else {
        mxDestroyArray(mx_MONdataB);
        mx_MONdataB = mxDuplicateArray(opt);
      }
    }
  }

  /* We made it here without problems */

  return(0);
}


/* MUST ADD: maxrs, epslin?, dqincfac? */

int get_LinSolvOptions(const mxArray *options, booleantype fwd,
                       int *mupper, int *mlower,
                       int *mudq, int *mldq,
                       int *gstype, int *maxl)
{
  mxArray *opt;
  char *bufval;
  int buflen, status, tmp_ls, tmp_pm;
  
  *mupper = 0;
  *mlower = 0;

  *mudq = 0;
  *mldq = 0;

  *gstype = MODIFIED_GS;
  *maxl = 0;

  tmp_ls = LS_DENSE;
  tmp_pm = PM_NONE;

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
    if(!strcmp(bufval,"Dense")) tmp_ls = LS_DENSE;
    else if(!strcmp(bufval,"Band")) tmp_ls = LS_BAND;
    else if(!strcmp(bufval,"GMRES")) tmp_ls = LS_SPGMR;
    else if(!strcmp(bufval,"BiCGStab")) tmp_ls = LS_SPBCG;
    else if(!strcmp(bufval,"TFQMR")) tmp_ls = LS_SPTFQMR;
    else mexErrMsgTxt("LinearSolver has an illegal value.");

    if (fwd) ls  = tmp_ls;
    else     lsB = tmp_ls;
  }
  
  /* Jacobian function */

  opt = mxGetField(options,0,"JacobianFn");
  if ( !mxIsEmpty(opt) ) {
    if (fwd) {
      mxDestroyArray(mx_JACfct);
      mx_JACfct  = mxDuplicateArray(opt);
    } else {
      mxDestroyArray(mx_JACfctB);
      mx_JACfctB = mxDuplicateArray(opt);
    }
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
      if(status != 0)
        mexErrMsgTxt("Cannot parse GramSchmidtType.");
      if(!strcmp(bufval,"Classical")) *gstype = CLASSICAL_GS;
      else if(!strcmp(bufval,"Modified")) *gstype = MODIFIED_GS;
      else mexErrMsgTxt("GramSchmidtType has an illegal value.");
    }

  }

  /* SPILS linear solver options */

  if ( (tmp_ls==LS_SPGMR) || (tmp_ls==LS_SPBCG) || (tmp_ls==LS_SPTFQMR) ) {

    /* Max. dimension of Krylov subspace */

    opt = mxGetField(options,0,"KrylovMaxDim");
    if ( !mxIsEmpty(opt) ) {
      *maxl = (int)*mxGetPr(opt);
      if (*maxl < 0) 
        mexErrMsgTxt("KrylovMaxDim is negative.");
    }

    /* User defined precoditioning */

    opt = mxGetField(options,0,"PrecSetupFn");
    if ( !mxIsEmpty(opt) ) {
      if (fwd) {
        mxDestroyArray(mx_PSETfct);
        mx_PSETfct  = mxDuplicateArray(opt);
      } else {
        mxDestroyArray(mx_PSETfctB);
        mx_PSETfctB = mxDuplicateArray(opt);
      }
    }

    opt = mxGetField(options,0,"PrecSolveFn");
    if ( !mxIsEmpty(opt) ) {
      if (fwd) {
        mxDestroyArray(mx_PSOLfct);
        mx_PSOLfct  = mxDuplicateArray(opt);
      } else {
        mxDestroyArray(mx_PSOLfctB);
        mx_PSOLfctB = mxDuplicateArray(opt);
      }
    }

    /* Preconditioner module */
  
    opt = mxGetField(options,0,"PrecModule");
    if ( !mxIsEmpty(opt) ) {
      buflen = mxGetM(opt) * mxGetN(opt) + 1;
      bufval = mxCalloc(buflen, sizeof(char));
      status = mxGetString(opt, bufval, buflen);
      if(status != 0)
        mexErrMsgTxt("Cannot parse PrecModule.");
      if(!strcmp(bufval,"BBDPre")) tmp_pm = PM_BBDPRE;
      else if(!strcmp(bufval,"UserDefined")) tmp_pm = PM_NONE;
      else mexErrMsgTxt("PrecModule has an illegal value.");
      
      if (fwd) pm  = tmp_pm;
      else     pmB = tmp_pm;
    }

    if (tmp_pm == PM_BBDPRE) {
      
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
        if (fwd) {
          mxDestroyArray(mx_GLOCfct);
          mx_GLOCfct  = mxDuplicateArray(opt);
        } else {
          mxDestroyArray(mx_GLOCfctB);
          mx_GLOCfctB = mxDuplicateArray(opt);
        }
      } else {
        mexErrMsgTxt("GlocalFn required for BBD preconditioner.");
      }      

      opt = mxGetField(options,0,"GcommFn");
      if ( !mxIsEmpty(opt) ) {
        if (fwd) {
          mxDestroyArray(mx_GCOMfct);
          mx_GCOMfct  = mxDuplicateArray(opt);
        } else {
          mxDestroyArray(mx_GCOMfctB);
          mx_GCOMfctB = mxDuplicateArray(opt);
        }
      }
    }

  }

  
  /* We made it here without problems */

  return(0);

}


int get_QuadOptions(const mxArray *options, booleantype fwd,
                    double **yQ0, booleantype *errconQ,
                    int *itolQ, double *reltolQ, double *SabstolQ, double **VabstolQ)
{
  mxArray *opt;
  char *bufval;
  int tmp_Nq, i, buflen, status, m, n;
  double *tmp;

  tmp_Nq = 0;
  *yQ0 = NULL;

  *errconQ = FALSE;
  *itolQ = IDA_SS;
  *reltolQ = 1.0e-4;
  *SabstolQ = 1.0e-6;
  *VabstolQ = NULL;

  /* Return now if options was empty */

  if (mxIsEmpty(options)) return(0);

  /* Initial conditions for quadratures */

  opt = mxGetField(options,0,"QuadInitCond");
  if ( !mxIsEmpty(opt) ) {
    m = mxGetM(opt);
    n = mxGetN(opt);
    if ( (n != 1) && (m != 1) )
      mexErrMsgTxt("QuadInitCond is not a vector.");
    if ( n > m ) tmp_Nq = n;
    else         tmp_Nq = m;
    tmp = mxGetPr(opt);
    *yQ0 = (double *)malloc((tmp_Nq)*sizeof(double));
    for(i=0;i<tmp_Nq;i++) 
      (*yQ0)[i] = tmp[i];
  } else {
    mexErrMsgTxt("QuadInitCond required for quadrature integration.");
  }

  if (fwd) Nq  = tmp_Nq;
  else     NqB = tmp_Nq;

  /* Quadrature function */

  opt = mxGetField(options,0,"QuadRhsFn");
  if ( !mxIsEmpty(opt) ) {
    if (fwd) {
      mxDestroyArray(mx_QUADfct);
      mx_QUADfct  = mxDuplicateArray(opt);
    } else {
      mxDestroyArray(mx_QUADfctB);
      mx_QUADfctB = mxDuplicateArray(opt);
    }
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
          *itolQ = IDA_SS;
          *SabstolQ = *tmp;
          if (*SabstolQ < 0.0)
            mexErrMsgTxt("QuadAbsTol is negative.");
        } else if (n == tmp_Nq) {
          *itolQ = IDA_SV;
          *VabstolQ = (double *)malloc(tmp_Nq*sizeof(double));
          for(i=0;i<tmp_Nq;i++) {
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

  return(0);
}

int get_FSAOptions(const mxArray *options, 
                   char **pfield_name, int **plist, double **pbar,
                   booleantype *userSRES, int *dqtype, double *rho,
                   booleantype *errconS, int *itolS, double *reltolS, 
                   double **SabstolS, double **VabstolS)
{
  mxArray *opt;
  char *bufval;
  int i, is, m, n, buflen, status;
  double *tmp;

  /* Set default values */

  *userSRES = FALSE;
  *errconS = TRUE;

  *dqtype = IDA_CENTERED;
  *rho = 0.0;

  /*  *itolS = IDA_EE;*/
  *SabstolS = NULL;
  *VabstolS = NULL;

  *pfield_name = NULL;
  *plist = NULL;
  *pbar  = NULL;

  /* Return now if options was empty */

  if (mxIsEmpty(options)) return(0);

  /* Field name in data structure for params. */

  opt = mxGetField(options,0,"ParamField");
  if ( !mxIsEmpty(opt) ) {
    buflen = mxGetM(opt) * mxGetN(opt) + 1;
    bufval = mxCalloc(buflen, sizeof(char));
    status = mxGetString(opt, bufval, buflen);
    if(status != 0)
      mexErrMsgTxt("Could not parse ParamField.");
    *pfield_name = mxCalloc(buflen, sizeof(char));
    strcpy((*pfield_name), bufval);
  }  

  /* PLIST */

  opt = mxGetField(options,0,"ParamList");
  if ( !mxIsEmpty(opt) ) {
    tmp = mxGetPr(opt);
    m = mxGetM(opt);
    n = mxGetN(opt);
    if ( (n != 1) && (m != 1) )
      mexErrMsgTxt("ParamList is not a vector.");
    if (m > n) n = m;
    if ( n != Ns)
      mexErrMsgTxt("ParamList does not contain Ns elements.");
    *plist = (int *) malloc(Ns*sizeof(int));
    for (is=0;is<Ns;is++)
      (*plist)[is] = (int) tmp[is];
  }

  /* PBAR */

  opt = mxGetField(options,0,"ParamScales");
  if ( !mxIsEmpty(opt) ) {
    m = mxGetM(opt);
    n = mxGetN(opt);
    if ( (n != 1) && (m != 1) )
      mexErrMsgTxt("ParamScales is not a vector.");
    if ( m > n ) n = m;
    if ( n != Ns)
      mexErrMsgTxt("ParamScales does not contain Ns elements.");
    tmp = mxGetPr(opt);
    *pbar = (double *) malloc(Ns*sizeof(double));
    for(i=0;i<Ns;i++)
      (*pbar)[i] = tmp[i];
  }

  /* DQtype and DQparam */

  opt = mxGetField(options,0,"SensDQtype");
  if ( !mxIsEmpty(opt) ) {
    buflen = mxGetM(opt) * mxGetN(opt) + 1;
    bufval = mxCalloc(buflen, sizeof(char));
    status = mxGetString(opt, bufval, buflen);
    if(status != 0)
      mexErrMsgTxt("Cannot parse SensDQtype.");
    if(!strcmp(bufval,"Centered")) *dqtype = IDA_CENTERED;
    else if(!strcmp(bufval,"Forward")) *dqtype = IDA_FORWARD;
    else mexErrMsgTxt("SensDQtype has an illegal value.");
  }


  opt = mxGetField(options,0,"SensDQparam");
  if ( !mxIsEmpty(opt) )
    *rho = *mxGetPr(opt);

  /* Error control */

  opt = mxGetField(options,0,"SensErrControl");
  if ( !mxIsEmpty(opt) ) {
    buflen = mxGetM(opt) * mxGetN(opt) + 1;
    bufval = mxCalloc(buflen, sizeof(char));
    status = mxGetString(opt, bufval, buflen);
    if(status != 0)
      mexErrMsgTxt("Canot parse SensErrControl.");
    if(!strcmp(bufval,"off")) *errconS = FALSE;
    else if(!strcmp(bufval,"on")) *errconS = TRUE;
    else mexErrMsgTxt("SensErrControl has an illegal value.");
  }

  /* Tolerances */
  
  opt = mxGetField(options,0,"SensRelTol");

  if ( !mxIsEmpty(opt) ) {
    *reltolS = *mxGetPr(opt);
    if (*reltolS < 0.0)
      mexErrMsgTxt("SensRelTol is negative.");
    opt = mxGetField(options,0,"SensAbsTol");
    if ( !mxIsEmpty(opt) ) {
      m = mxGetM(opt);
      n = mxGetN(opt);
      if ( (m == 1) && (n == Ns) ) {
        *itolS = IDA_SS;
        tmp = mxGetPr(opt);
        *SabstolS = (double *) malloc(Ns*sizeof(double));
        for (is=0; is<Ns; is++) {
          (*SabstolS)[is] = tmp[is];
          if ( tmp[is] < 0.0 )
            mexErrMsgTxt("SensAbsTol has a negative component.");
        }
      } else if ( (m == N) && (n == Ns) ) {
        *itolS = IDA_SV;
        tmp = mxGetPr(opt);
        *VabstolS = (double *)malloc(Ns*N*sizeof(double));
        for (i=0; i<Ns*N; i++) {
          (*VabstolS)[i] = tmp[i];
          if ( tmp[i] < 0.0 )
            mexErrMsgTxt("SensAbsTol has a negative component.");
        }
      } else {
        mexErrMsgTxt("SensAbsTol must be either a 1xNs vector or a NxNs matrix.");
      }
    } else {
      /**itolS = IDA_EE;*/
    }
  }

  /* Sensitivity RES function type */

  opt = mxGetField(options,0,"SensResFn");
  if ( !mxIsEmpty(opt) ) {
    *userSRES = TRUE;
    mxDestroyArray(mx_SRESfct);
    mx_SRESfct = mxDuplicateArray(opt);
  } 

  /* We made it here without problems */

  return(0);

}

