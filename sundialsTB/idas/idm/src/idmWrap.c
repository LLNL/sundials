/*
 * -----------------------------------------------------------------
 * $Revision: 1.4 $
 * $Date: 2006-10-09 23:56:24 $
 * -----------------------------------------------------------------
 * Programmer: Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2005, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see sundials-x.y.z/src/ida/LICENSE.
 * -----------------------------------------------------------------
 * IDA wrapper functions.
 * -----------------------------------------------------------------
 */

#include "idm.h"
#include "nvm.h"

static void UpdateUserData(mxArray *mx_data);
static void UpdateMonitorData(mxArray *mx_data);
static void UpdateMonitorDataB(mxArray *mx_data);

booleantype idm_quad;      /* Forward quadratures? */
booleantype idm_quadB;     /* Backward quadratures? */
booleantype idm_asa;       /* Adjoint sensitivity? */
booleantype idm_fsa;       /* Forward sensitivity? */

/*
 * ---------------------------------------------------------------------------------
 * Redability replacements
 * ---------------------------------------------------------------------------------
 */

#define N           (idm_Cdata->N) 
#define Nq          (idm_Cdata->Nq) 
#define Ng          (idm_Cdata->Ng) 
#define Ns          (idm_Cdata->Ns) 
#define Nd          (idm_Cdata->Nd) 
#define Nc          (idm_Cdata->Nc) 
#define ls          (idm_Cdata->ls) 
#define pm          (idm_Cdata->pm) 
#define ism         (idm_Cdata->ism) 
#define NB          (idm_Cdata->NB) 
#define NqB         (idm_Cdata->NqB) 
#define lsB         (idm_Cdata->lsB) 
#define pmB         (idm_Cdata->pmB) 

#define mx_data     (idm_Mdata->mx_data)

#define mx_RESfct   (idm_Mdata->mx_RESfct)
#define mx_QUADfct  (idm_Mdata->mx_QUADfct)
#define mx_JACfct   (idm_Mdata->mx_JACfct)
#define mx_PSETfct  (idm_Mdata->mx_PSETfct)
#define mx_PSOLfct  (idm_Mdata->mx_PSOLfct)
#define mx_GLOCfct  (idm_Mdata->mx_GLOCfct)
#define mx_GCOMfct  (idm_Mdata->mx_GCOMfct)
#define mx_Gfct     (idm_Mdata->mx_Gfct)
#define mx_SRESfct  (idm_Mdata->mx_SRESfct)

#define mx_RESfctB  (idm_Mdata->mx_RESfctB)
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
 * Error handler
 * ---------------------------------------------------------------------------------
 */

void mtlb_IdaErrHandler(int error_code, 
                        const char *module, const char *function, 
                        char *msg, void *eh_data)
{
  char err_type[10];

  if (error_code == IDA_WARNING)
    sprintf(err_type,"WARNING");
  else
    sprintf(err_type,"ERROR");

  mexPrintf("\n[%s %s]  %s\n",module,err_type,function);
  mexPrintf("  %s\n\n",msg);

}

/*
 * ---------------------------------------------------------------------------------
 * Wrapper functions
 * ---------------------------------------------------------------------------------
 */

int mtlb_IdaRes(realtype tt, N_Vector yy, N_Vector yp,
                N_Vector rr, void *res_data)
{
  mxArray *mx_in[6], *mx_out[3];
  int ret;

  /* Inputs to the Matlab function */
  mx_in[0] = mxCreateScalarDouble(1.0);        /* type=1: forward ODE */
  mx_in[1] = mxCreateScalarDouble(tt);         /* current t */
  mx_in[2] = mxCreateDoubleMatrix(N,1,mxREAL); /* current yy */
  mx_in[3] = mxCreateDoubleMatrix(N,1,mxREAL); /* current yp */
  mx_in[4] = mx_RESfct;                        /* matlab function handle */ 
  mx_in[5] = mx_data;                          /* matlab user data */

  /* Call matlab wrapper */
  GetData(yy, mxGetPr(mx_in[2]), N);
  GetData(yp, mxGetPr(mx_in[3]), N);

  mexCallMATLAB(3,mx_out,6,mx_in,"idm_res");

  PutData(rr, mxGetPr(mx_out[0]), N);
  ret = (int)*mxGetPr(mx_out[1]);

  if (!mxIsEmpty(mx_out[2])) {
    UpdateUserData(mx_out[2]);
  }

  /* Free temporary space */
  mxDestroyArray(mx_in[0]);
  mxDestroyArray(mx_in[1]);
  mxDestroyArray(mx_in[2]);
  mxDestroyArray(mx_in[3]);
  mxDestroyArray(mx_out[0]);
  mxDestroyArray(mx_out[1]);
  mxDestroyArray(mx_out[2]);

  return(ret);
}

int mtlb_IdaQuadFct(realtype tres, N_Vector yy, N_Vector yp, N_Vector ypQ,
                    void *rdataQ)
{
  mxArray *mx_in[6], *mx_out[3];
  int ret;

  /* Inputs to the Matlab function */
  mx_in[0] = mxCreateScalarDouble(1.0);        /* type=1: forward ODE */
  mx_in[1] = mxCreateScalarDouble(tres);       /* current t */
  mx_in[2] = mxCreateDoubleMatrix(N,1,mxREAL); /* current yy */
  mx_in[3] = mxCreateDoubleMatrix(N,1,mxREAL); /* current yp */
  mx_in[4] = mx_QUADfct;                       /* matlab function handle */ 
  mx_in[5] = mx_data;                          /* matlab user data */

  /* Call matlab wrapper */
  GetData(yy, mxGetPr(mx_in[2]), N);
  GetData(yp, mxGetPr(mx_in[3]), N);

  mexCallMATLAB(3,mx_out,6,mx_in,"idm_rhsQ");

  PutData(ypQ, mxGetPr(mx_out[0]), Nq);
  ret = (int)*mxGetPr(mx_out[1]);

  if (!mxIsEmpty(mx_out[2])) {
    UpdateUserData(mx_out[2]);
  }

  /* Free temporary space */
  mxDestroyArray(mx_in[0]);
  mxDestroyArray(mx_in[1]);
  mxDestroyArray(mx_in[2]);
  mxDestroyArray(mx_in[3]);
  mxDestroyArray(mx_out[0]);
  mxDestroyArray(mx_out[1]);
  mxDestroyArray(mx_out[2]);

  return(ret);
}

int mtlb_IdaGfct(realtype t, N_Vector yy, N_Vector yp,
                 realtype *gout, void *g_data)
{
  double *gdata;
  int i, ret;
  mxArray *mx_in[5], *mx_out[3];
  
  /* Inputs to the Matlab function */
  mx_in[0] = mxCreateScalarDouble(t);          /* current t */
  mx_in[1] = mxCreateDoubleMatrix(N,1,mxREAL); /* current yy */
  mx_in[2] = mxCreateDoubleMatrix(N,1,mxREAL); /* current yp */
  mx_in[3] = mx_Gfct;                          /* matlab function handle */
  mx_in[4] = mx_data;                          /* matlab user data */
  
  /* Call matlab wrapper */
  GetData(yy, mxGetPr(mx_in[1]), N);
  GetData(yp, mxGetPr(mx_in[2]), N);

  mexCallMATLAB(3,mx_out,5,mx_in,"idm_root");

  gdata = mxGetPr(mx_out[0]);
  for (i=0;i<Ng;i++) gout[i] = gdata[i];

  ret = (int)*mxGetPr(mx_out[1]);

  if (!mxIsEmpty(mx_out[2])) {
    UpdateUserData(mx_out[2]);
  }

  /* Free temporary space */
  mxDestroyArray(mx_in[0]);
  mxDestroyArray(mx_in[1]);
  mxDestroyArray(mx_in[2]);
  mxDestroyArray(mx_out[0]);
  mxDestroyArray(mx_out[1]);
  mxDestroyArray(mx_out[2]);

  return(ret);
}

int mtlb_IdaDenseJac(long int Neq, realtype tt, 
                     N_Vector yy, N_Vector yp, N_Vector rr,
                     realtype c_j, void *jac_data, 
                     DenseMat Jac, 
                     N_Vector tmp1, N_Vector tmp2, 
                     N_Vector tmp3)
{
  double *J_data;
  long int i;
  int ret;
  mxArray *mx_in[8], *mx_out[3];

  /* Inputs to the Matlab function */
  mx_in[0] = mxCreateScalarDouble(1.0);         /* type=1: forward ODE */
  mx_in[1] = mxCreateScalarDouble(tt);          /* current t */  
  mx_in[2] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current yy */
  mx_in[3] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current yp */
  mx_in[4] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current rr */
  mx_in[5] = mxCreateScalarDouble(c_j);         /* current c_j */  
  mx_in[6] = mx_JACfct;                         /* matlab function handle */
  mx_in[7] = mx_data;                           /* matlab user data */
  
  /* Call matlab wrapper */
  GetData(yy, mxGetPr(mx_in[2]), N);
  GetData(yp, mxGetPr(mx_in[3]), N);
  GetData(rr, mxGetPr(mx_in[4]), N);

  mexCallMATLAB(3,mx_out,8,mx_in,"idm_djac");

  /* Extract data */
  J_data = mxGetPr(mx_out[0]);
  for (i=0;i<N;i++)
    memcpy(DENSE_COL(Jac,i), J_data + i*N, N*sizeof(double));

  ret = (int)*mxGetPr(mx_out[1]);

  if (!mxIsEmpty(mx_out[2])) {
    UpdateUserData(mx_out[2]);
  }

  /* Free temporary space */
  mxDestroyArray(mx_in[0]);
  mxDestroyArray(mx_in[1]);
  mxDestroyArray(mx_in[2]);
  mxDestroyArray(mx_in[3]);
  mxDestroyArray(mx_in[4]);
  mxDestroyArray(mx_in[5]);
  mxDestroyArray(mx_out[0]);
  mxDestroyArray(mx_out[1]);
  mxDestroyArray(mx_out[2]);

  return(ret);
}

int mtlb_IdaBandJac(long int Neq, long int mupper, 
                    long int mlower, realtype tt, 
                    N_Vector yy, N_Vector yp, N_Vector rr, 
                    realtype c_j, void *jac_data, BandMat Jac, 
                    N_Vector tmp1, N_Vector tmp2, 
                    N_Vector tmp3)
{
  double *J_data;
  long int eband, i;
  int ret;
  mxArray *mx_in[8], *mx_out[3];

  /* Inputs to the Matlab function */
  mx_in[0] = mxCreateScalarDouble(1.0);         /* type=1: forward ODE */
  mx_in[1] = mxCreateScalarDouble(tt);          /* current t */
  mx_in[2] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current yy */
  mx_in[3] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current yp */
  mx_in[4] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current rr */
  mx_in[5] = mxCreateScalarDouble(c_j);         /* current c_j */
  mx_in[6] = mx_JACfct;                         /* matlab function handle */
  mx_in[7] = mx_data;                           /* matlab user data */
  
  /* Call matlab wrapper */
  GetData(yy, mxGetPr(mx_in[2]), N);
  GetData(yp, mxGetPr(mx_in[3]), N);
  GetData(rr, mxGetPr(mx_in[4]), N);

  mexCallMATLAB(3,mx_out,8,mx_in,"idm_bjac");

  /* Extract data */
  eband =  mupper + mlower + 1;
  J_data = mxGetPr(mx_out[0]);
  for (i=0;i<N;i++)
    memcpy(BAND_COL(Jac,i) - mupper, J_data + i*eband, eband*sizeof(double));

  ret = (int)*mxGetPr(mx_out[1]);
  
  if (!mxIsEmpty(mx_out[2])) {
    UpdateUserData(mx_out[2]);
  }

  /* Free temporary space */
  mxDestroyArray(mx_in[0]);
  mxDestroyArray(mx_in[1]);
  mxDestroyArray(mx_in[2]);
  mxDestroyArray(mx_in[3]);
  mxDestroyArray(mx_in[4]);
  mxDestroyArray(mx_in[5]);
  mxDestroyArray(mx_out[0]);
  mxDestroyArray(mx_out[1]);
  mxDestroyArray(mx_out[2]);

  return(ret);
}

int mtlb_IdaSpilsJac(realtype tt,
                     N_Vector yy, N_Vector yp, N_Vector rr,
                     N_Vector v, N_Vector Jv,
                     realtype c_j, void *jac_data,
                     N_Vector tmp1, N_Vector tmp2)
{
  mxArray *mx_in[9], *mx_out[3];
  int ret;

  /* Inputs to the Matlab function */
  mx_in[0] = mxCreateScalarDouble(1.0);         /* type=1: forward ODE */
  mx_in[1] = mxCreateScalarDouble(tt);          /* current t */ 
  mx_in[2] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current yy */
  mx_in[3] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current yp */
  mx_in[4] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current rr */
  mx_in[5] = mxCreateDoubleMatrix(N,1,mxREAL);  /* vector v */
  mx_in[6] = mxCreateScalarDouble(c_j);         /* current c_j */ 
  mx_in[7] = mx_JACfct;                         /* matlab function handle */
  mx_in[8] = mx_data;                           /* matlab user data */
  
  /* Call matlab wrapper */
  GetData(yy, mxGetPr(mx_in[2]), N);
  GetData(yp, mxGetPr(mx_in[3]), N);
  GetData(rr, mxGetPr(mx_in[4]), N);
  GetData(v, mxGetPr(mx_in[5]), N);

  mexCallMATLAB(3,mx_out,9,mx_in,"idm_jtv");

  PutData(Jv, mxGetPr(mx_out[0]), N);
  ret = (int)*mxGetPr(mx_out[1]);

  if (!mxIsEmpty(mx_out[2])) {
    UpdateUserData(mx_out[2]);
  }

  /* Free temporary space */
  mxDestroyArray(mx_in[0]);
  mxDestroyArray(mx_in[1]);
  mxDestroyArray(mx_in[2]);
  mxDestroyArray(mx_in[3]);
  mxDestroyArray(mx_in[4]);
  mxDestroyArray(mx_in[5]);
  mxDestroyArray(mx_in[6]);
  mxDestroyArray(mx_out[0]);
  mxDestroyArray(mx_out[1]);
  mxDestroyArray(mx_out[2]);

  return(ret);
}

int mtlb_IdaSpilsPset(realtype tt,
                      N_Vector yy, N_Vector yp, N_Vector rr,
                      realtype c_j, void *prec_data,
                      N_Vector tmp1, N_Vector tmp2,
                      N_Vector tmp3)
{
  mxArray *mx_in[8], *mx_out[2];
  int ret;

  /* Inputs to the Matlab function */
  mx_in[0] = mxCreateScalarDouble(1.0);         /* type=1: forward ODE */
  mx_in[1] = mxCreateScalarDouble(tt);          /* current t */
  mx_in[2] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current yy */
  mx_in[3] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current yp */
  mx_in[4] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current rr */
  mx_in[5] = mxCreateLogicalScalar(c_j);        /* current c_j */
  mx_in[6] = mx_PSETfct;                        /* matlab function handle */
  mx_in[7] = mx_data;                           /* matlab user data */
  
  /* Call matlab wrapper */
  GetData(yy, mxGetPr(mx_in[2]), N);
  GetData(yp, mxGetPr(mx_in[3]), N);
  GetData(rr, mxGetPr(mx_in[4]), N);

  mexCallMATLAB(2,mx_out,8,mx_in,"idm_pset");

  ret = (int)*mxGetPr(mx_out[0]);

  if (!mxIsEmpty(mx_out[1])) {
    UpdateUserData(mx_out[1]);
  }

  /* Free temporary space */
  mxDestroyArray(mx_in[0]);
  mxDestroyArray(mx_in[1]);
  mxDestroyArray(mx_in[2]);
  mxDestroyArray(mx_in[3]);
  mxDestroyArray(mx_in[4]);
  mxDestroyArray(mx_in[5]);
  mxDestroyArray(mx_out[0]);
  mxDestroyArray(mx_out[1]);

  return(ret);
}

int mtlb_IdaSpilsPsol(realtype tt,
                      N_Vector yy, N_Vector yp, N_Vector rr,
                      N_Vector rvec, N_Vector zvec,
                      realtype c_j, realtype delta, void *prec_data,
                      N_Vector tmp)
{
  mxArray *mx_in[9], *mx_out[3];
  int ret;

  /* Inputs to the Matlab function */
  mx_in[0] = mxCreateScalarDouble(1.0);        /* type=1: forward ODE */
  mx_in[1] = mxCreateScalarDouble(tt);         /* current t */   
  mx_in[2] = mxCreateDoubleMatrix(N,1,mxREAL); /* current yy */
  mx_in[3] = mxCreateDoubleMatrix(N,1,mxREAL); /* current yp */
  mx_in[4] = mxCreateDoubleMatrix(N,1,mxREAL); /* current rr */
  mx_in[5] = mxCreateDoubleMatrix(N,1,mxREAL); /* right hand side r */
  mx_in[6] = mxCreateScalarDouble(c_j);        /* current c_j */   
  mx_in[7] = mx_PSOLfct;                       /* matlab function handle */
  mx_in[8] = mx_data;                          /* matlab user data */
  
  /* Call matlab wrapper */
  GetData(yy, mxGetPr(mx_in[2]), N);
  GetData(yp, mxGetPr(mx_in[3]), N);
  GetData(rr, mxGetPr(mx_in[4]), N);
  GetData(rvec, mxGetPr(mx_in[5]), N);

  mexCallMATLAB(3,mx_out,9,mx_in,"idm_psol");

  PutData(zvec, mxGetPr(mx_out[0]), N);
  ret = (int)*mxGetPr(mx_out[1]);

  if (!mxIsEmpty(mx_out[2])) {
    UpdateUserData(mx_out[2]);
  }

  /* Free temporary space */
  mxDestroyArray(mx_in[0]);
  mxDestroyArray(mx_in[1]);
  mxDestroyArray(mx_in[2]);
  mxDestroyArray(mx_in[3]);
  mxDestroyArray(mx_in[4]);
  mxDestroyArray(mx_in[5]);
  mxDestroyArray(mx_in[6]);
  mxDestroyArray(mx_out[0]);
  mxDestroyArray(mx_out[1]);
  mxDestroyArray(mx_out[2]);

  return(ret);
}

/*
 * ----------------------------
 * BBD PRECONDITONER FUNCTIONS
 * ----------------------------
 */

int mtlb_IdaBBDgloc(long int Nlocal, realtype tt,
                    N_Vector yy, N_Vector yp, N_Vector gval,
                    void *res_data)
{
  mxArray *mx_in[6], *mx_out[3];
  int ret;

  /* Inputs to the Matlab function */
  mx_in[0] = mxCreateScalarDouble(1.0);         /* type=1: forward ODE */
  mx_in[1] = mxCreateScalarDouble(tt);          /* current t */
  mx_in[2] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current yy */
  mx_in[3] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current yp */
  mx_in[4] = mx_GLOCfct;                        /* matlab function handle */
  mx_in[5] = mx_data;                           /* matlab user data */
  
  /* Call matlab wrapper */
  GetData(yy, mxGetPr(mx_in[2]), N);
  GetData(yp, mxGetPr(mx_in[3]), N);

  mexCallMATLAB(3,mx_out,6,mx_in,"idm_gloc");

  PutData(gval, mxGetPr(mx_out[0]), N);
  ret = (int)*mxGetPr(mx_out[1]);

  if (!mxIsEmpty(mx_out[2])) {
    UpdateUserData(mx_out[2]);
  }

  /* Free temporary space */
  mxDestroyArray(mx_in[0]);
  mxDestroyArray(mx_in[1]);
  mxDestroyArray(mx_in[2]);
  mxDestroyArray(mx_in[3]);
  mxDestroyArray(mx_out[0]);
  mxDestroyArray(mx_out[1]);
  mxDestroyArray(mx_out[2]);

  return(ret);
}

int mtlb_IdaBBDgcom(long int Nlocal, realtype tt,
                    N_Vector yy, N_Vector yp,
                    void *res_data)
{
  mxArray *mx_in[6], *mx_out[2];
  int ret;

  /* Inputs to the Matlab function */
  mx_in[0] = mxCreateScalarDouble(1.0);         /* type=1: forward ODE */
  mx_in[1] = mxCreateScalarDouble(tt);          /* current t */
  mx_in[2] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current yy */
  mx_in[3] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current yp */
  mx_in[4] = mx_GCOMfct;                        /* matlab function handle */
  mx_in[5] = mx_data;                           /* matlab user data */
  
  /* Call matlab wrapper */
  GetData(yy, mxGetPr(mx_in[2]), N);
  GetData(yp, mxGetPr(mx_in[3]), N);

  mexCallMATLAB(2,mx_out,6,mx_in,"idm_gcom");

  ret = (int)*mxGetPr(mx_out[0]);

  if (!mxIsEmpty(mx_out[1])) {
    UpdateUserData(mx_out[1]);
  }

  /* Free temporary space */
  mxDestroyArray(mx_in[0]);
  mxDestroyArray(mx_in[1]);
  mxDestroyArray(mx_in[2]);
  mxDestroyArray(mx_in[3]);
  mxDestroyArray(mx_out[0]);
  mxDestroyArray(mx_out[1]);

  return(ret);
}

/*
 * ----------------------------
 * FORWARD SENSITVITY FUNCTIONS
 * ----------------------------
 */

int mtlb_IdaSensRes(int Nsens, realtype tres, 
                    N_Vector yy, N_Vector yp, N_Vector rr,
                    N_Vector *yyS, N_Vector *ypS, N_Vector *rrS,
                    void *rdataS,
                    N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  mxArray *mx_in[9], *mx_out[3];
  int is, ret;
  double *tmp_yyS, *tmp_ypS, *tmp_rrS;

  /* Inputs to the Matlab function */
  mx_in[0] = mxCreateScalarDouble(tres);          /* current t */
  mx_in[1] = mxCreateDoubleMatrix(N,1,mxREAL);    /* current yy */
  mx_in[2] = mxCreateDoubleMatrix(N,1,mxREAL);    /* current yp */
  mx_in[3] = mxCreateDoubleMatrix(N,1,mxREAL);    /* current rr */
  mx_in[4] = mxCreateScalarDouble(Ns);            /* number of sensitivities */
  mx_in[5] = mxCreateDoubleMatrix(N*Ns,1,mxREAL); /* current yyS */
  mx_in[6] = mxCreateDoubleMatrix(N*Ns,1,mxREAL); /* current ypS */
  mx_in[7] = mx_SRESfct;                          /* matlab function handle */      
  mx_in[8] = mx_data;                             /* matlab user data */
  
  /* Call matlab wrapper */
  GetData(yy, mxGetPr(mx_in[1]), N);
  GetData(yp, mxGetPr(mx_in[2]), N);
  GetData(rr, mxGetPr(mx_in[3]), N);
  tmp_yyS = mxGetPr(mx_in[5]);
  tmp_ypS = mxGetPr(mx_in[6]);
  for (is=0; is<Ns; is++) {
    GetData(yyS[is], &tmp_yyS[is*N], N);
    GetData(ypS[is], &tmp_ypS[is*N], N);
  }

  mexCallMATLAB(3,mx_out,9,mx_in,"idm_resS");
  
  tmp_rrS = mxGetPr(mx_out[0]);

  for(is=0;is<Ns;is++)
    PutData(rrS[is], &tmp_rrS[is*N], N);

  ret = (int)*mxGetPr(mx_out[1]);

  if (!mxIsEmpty(mx_out[2])) {
    UpdateUserData(mx_out[2]);
  }

  /* Free temporary space */

  mxDestroyArray(mx_in[0]);
  mxDestroyArray(mx_in[1]);
  mxDestroyArray(mx_in[2]);
  mxDestroyArray(mx_in[3]);
  mxDestroyArray(mx_in[4]);
  mxDestroyArray(mx_in[5]);
  mxDestroyArray(mx_in[6]);
  mxDestroyArray(mx_out[0]);
  mxDestroyArray(mx_out[1]);
  mxDestroyArray(mx_out[2]);

  return(ret);
}


/*
 * ----------------------------
 * ADJOINT SENSITVITY FUNCTIONS
 * ----------------------------
 */

int mtlb_IdaResB(realtype tt, 
                 N_Vector yy, N_Vector yp,
                 N_Vector yyB, N_Vector ypB, N_Vector rrB,
                 void *rdataB)
{
  mxArray *mx_in[8], *mx_out[3];
  int ret;

  /* Inputs to the Matlab function */
  mx_in[0] = mxCreateScalarDouble(-1.0);        /* type=-1: backward ODE */
  mx_in[1] = mxCreateScalarDouble(tt);          /* current t */
  mx_in[2] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current yy */
  mx_in[3] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current yp */
  mx_in[4] = mxCreateDoubleMatrix(NB,1,mxREAL); /* current yyB */
  mx_in[5] = mxCreateDoubleMatrix(NB,1,mxREAL); /* current ypB */
  mx_in[6] = mx_RESfctB;                        /* matlab function handle */ 
  mx_in[7] = mx_data;                           /* matlab user data */
  
  /* Call matlab wrapper */
  GetData(yy, mxGetPr(mx_in[2]), N);
  GetData(yp, mxGetPr(mx_in[3]), N);
  GetData(yyB, mxGetPr(mx_in[4]), NB);
  GetData(ypB, mxGetPr(mx_in[5]), NB);

  mexCallMATLAB(3,mx_out,8,mx_in,"idm_res");

  PutData(rrB, mxGetPr(mx_out[0]), NB);

  ret = (int)*mxGetPr(mx_out[1]);

  if (!mxIsEmpty(mx_out[2])) {
    UpdateUserData(mx_out[2]);
  }

  /* Free temporary space */
  mxDestroyArray(mx_in[0]);
  mxDestroyArray(mx_in[1]);
  mxDestroyArray(mx_in[2]);
  mxDestroyArray(mx_in[3]);
  mxDestroyArray(mx_in[4]);
  mxDestroyArray(mx_in[5]);
  mxDestroyArray(mx_out[0]);
  mxDestroyArray(mx_out[1]);
  mxDestroyArray(mx_out[2]);

  return(ret);
}

int mtlb_IdaQuadFctB(realtype tt, 
                     N_Vector yy, N_Vector yp, 
                     N_Vector yyB, N_Vector ypB,
                     N_Vector ypQB, void *rdataQB)
{
  mxArray *mx_in[8], *mx_out[3];
  int ret;

  /* Inputs to the Matlab function */
  mx_in[0] = mxCreateScalarDouble(-1.0);        /* type=-1: backward ODE */
  mx_in[1] = mxCreateScalarDouble(tt);          /* current t */
  mx_in[2] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current yy */
  mx_in[3] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current yp */
  mx_in[4] = mxCreateDoubleMatrix(NB,1,mxREAL); /* current yyB */
  mx_in[5] = mxCreateDoubleMatrix(NB,1,mxREAL); /* current ypB */
  mx_in[6] = mx_QUADfctB;                       /* matlab function handle */ 
  mx_in[7] = mx_data;                           /* matlab user data */

  /* Call matlab wrapper */
  GetData(yy, mxGetPr(mx_in[2]), N);
  GetData(yp, mxGetPr(mx_in[3]), N);
  GetData(yyB, mxGetPr(mx_in[4]), NB);
  GetData(ypB, mxGetPr(mx_in[5]), NB);

  mexCallMATLAB(3,mx_out,8,mx_in,"idm_rhsQ");

  PutData(ypQB, mxGetPr(mx_out[0]), NqB);

  ret = (int)*mxGetPr(mx_out[1]);

  if (!mxIsEmpty(mx_out[2])) {
    UpdateUserData(mx_out[2]);
  }

  /* Free temporary space */
  mxDestroyArray(mx_in[0]);
  mxDestroyArray(mx_in[1]);
  mxDestroyArray(mx_in[2]);
  mxDestroyArray(mx_in[3]);
  mxDestroyArray(mx_in[4]);
  mxDestroyArray(mx_in[5]);
  mxDestroyArray(mx_out[0]);
  mxDestroyArray(mx_out[1]);
  mxDestroyArray(mx_out[2]);

  return(ret);
}

int mtlb_IdaDenseJacB(long int NeqB, realtype tt, 
                      N_Vector yy, N_Vector yp,
                      N_Vector yyB, N_Vector ypB, N_Vector rrB,
                      realtype c_jB, void *jac_dataB, 
                      DenseMat JacB, 
                      N_Vector tmp1B, N_Vector tmp2B, 
                      N_Vector tmp3B)
{
  double *JB_data;
  long int i;
  mxArray *mx_in[10], *mx_out[3];
  int ret;

  /* Inputs to the Matlab function */
  mx_in[0] = mxCreateScalarDouble(-1.0);        /* type=-1: backward ODE */
  mx_in[1] = mxCreateScalarDouble(tt);          /* current t */  
  mx_in[2] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current yy */
  mx_in[3] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current yp */
  mx_in[4] = mxCreateDoubleMatrix(NB,1,mxREAL); /* current yyB */
  mx_in[5] = mxCreateDoubleMatrix(NB,1,mxREAL); /* current ypB */
  mx_in[6] = mxCreateDoubleMatrix(NB,1,mxREAL); /* current rrB */
  mx_in[7] = mxCreateScalarDouble(c_jB);        /* current c_jB */  
  mx_in[8] = mx_JACfctB;                        /* matlab function handle */
  mx_in[9] = mx_data;                           /* matlab user data */
  
  /* Call matlab wrapper */
  GetData(yy, mxGetPr(mx_in[2]), N);
  GetData(yp, mxGetPr(mx_in[3]), N);
  GetData(yyB, mxGetPr(mx_in[4]), NB);
  GetData(ypB, mxGetPr(mx_in[5]), NB);
  GetData(rrB, mxGetPr(mx_in[6]), NB);

  mexCallMATLAB(3,mx_out,10,mx_in,"idm_djac");

  JB_data = mxGetPr(mx_out[0]);
  for (i=0;i<NB;i++)
    memcpy(DENSE_COL(JacB,i), JB_data + i*NB, NB*sizeof(double));

  ret = (int)*mxGetPr(mx_out[1]);

  if (!mxIsEmpty(mx_out[2])) {
    UpdateUserData(mx_out[2]);
  }

  /* Free temporary space */
  mxDestroyArray(mx_in[0]);
  mxDestroyArray(mx_in[1]);
  mxDestroyArray(mx_in[2]);
  mxDestroyArray(mx_in[3]);
  mxDestroyArray(mx_in[4]);
  mxDestroyArray(mx_in[5]);
  mxDestroyArray(mx_in[6]);
  mxDestroyArray(mx_in[7]);
  mxDestroyArray(mx_out[0]);
  mxDestroyArray(mx_out[1]);
  mxDestroyArray(mx_out[2]);

  return(ret);
}

int mtlb_IdaBandJacB(long int NeqB, 
                     long int mupperB, long int mlowerB, 
                     realtype tt, 
                     N_Vector yy, N_Vector yp,
                     N_Vector yyB, N_Vector ypB, N_Vector rrB,
                     realtype c_jB, void *jac_dataB,
                     BandMat JacB, 
                     N_Vector tmp1B, N_Vector tmp2B, 
                     N_Vector tmp3B)
{
  double *JB_data;
  long int ebandB, i;
  mxArray *mx_in[10], *mx_out[3];
  int ret;

  /* Inputs to the Matlab function */
  mx_in[0] = mxCreateScalarDouble(-1.0);        /* type=-1: backward ODE */
  mx_in[1] = mxCreateScalarDouble(tt);          /* current t */
  mx_in[2] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current yy */
  mx_in[3] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current yp */
  mx_in[4] = mxCreateDoubleMatrix(NB,1,mxREAL); /* current yyB */
  mx_in[5] = mxCreateDoubleMatrix(NB,1,mxREAL); /* current ypB */
  mx_in[6] = mxCreateDoubleMatrix(NB,1,mxREAL); /* current rrB */
  mx_in[7] = mxCreateScalarDouble(c_jB);        /* current c_jB */
  mx_in[8] = mx_JACfctB;                        /* matlab function handle */
  mx_in[9] = mx_data;                           /* matlab user data */
  
  /* Call matlab wrapper */
  GetData(yy, mxGetPr(mx_in[2]), N);
  GetData(yp, mxGetPr(mx_in[3]), N);
  GetData(yyB, mxGetPr(mx_in[4]), NB);
  GetData(ypB, mxGetPr(mx_in[5]), NB);
  GetData(rrB, mxGetPr(mx_in[6]), NB);

  mexCallMATLAB(3,mx_out,10,mx_in,"idm_bjac");

  ebandB =  mupperB + mlowerB + 1;
  JB_data = mxGetPr(mx_out[0]);
  for (i=0;i<NB;i++)
    memcpy(BAND_COL(JacB,i) - mupperB, JB_data + i*ebandB, ebandB*sizeof(double));
    
  ret = (int)*mxGetPr(mx_out[1]);

  if (!mxIsEmpty(mx_out[2])) {
    UpdateUserData(mx_out[2]);
  }

  /* Free temporary space */
  mxDestroyArray(mx_in[0]);
  mxDestroyArray(mx_in[1]);
  mxDestroyArray(mx_in[2]);
  mxDestroyArray(mx_in[3]);
  mxDestroyArray(mx_in[4]);
  mxDestroyArray(mx_in[5]);
  mxDestroyArray(mx_in[6]);
  mxDestroyArray(mx_in[7]);
  mxDestroyArray(mx_out[0]);
  mxDestroyArray(mx_out[1]);
  mxDestroyArray(mx_out[2]);

  return(ret);
}

int mtlb_IdaSpilsJacB(realtype tt,
                      N_Vector yy, N_Vector yp,
                      N_Vector yyB, N_Vector ypB, N_Vector rrB,
                      N_Vector vB, N_Vector JvB, 
                      realtype c_jB, void *jac_dataB, 
                      N_Vector tmp1B, N_Vector tmp2B)
{
  mxArray *mx_in[11], *mx_out[3];
  int ret;

  /* Inputs to the Matlab function */
  mx_in[0] = mxCreateScalarDouble(-1.0);        /* type=-1: backward ODE */
  mx_in[1] = mxCreateScalarDouble(tt);          /* current t */ 
  mx_in[2] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current yy */
  mx_in[3] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current yp */
  mx_in[4] = mxCreateDoubleMatrix(NB,1,mxREAL); /* current yyB */
  mx_in[5] = mxCreateDoubleMatrix(NB,1,mxREAL); /* current ypB */
  mx_in[6] = mxCreateDoubleMatrix(NB,1,mxREAL); /* current rrB */
  mx_in[7] = mxCreateDoubleMatrix(NB,1,mxREAL); /* vector vB */
  mx_in[8] = mxCreateScalarDouble(c_jB);        /* current c_jB */ 
  mx_in[9] = mx_JACfctB;                        /* matlab function handle */
  mx_in[10] = mx_data;                          /* matlab user data */
  
  /* Call matlab wrapper */
  GetData(yy, mxGetPr(mx_in[2]), N);
  GetData(yp, mxGetPr(mx_in[3]), N);
  GetData(yyB, mxGetPr(mx_in[4]), NB);
  GetData(ypB, mxGetPr(mx_in[5]), NB);
  GetData(rrB, mxGetPr(mx_in[6]), NB);
  GetData(vB, mxGetPr(mx_in[7]), NB);

  mexCallMATLAB(3,mx_out,11,mx_in,"idm_jtv");

  PutData(JvB, mxGetPr(mx_out[0]), NB);
  ret = (int)*mxGetPr(mx_out[1]);
  
  if (!mxIsEmpty(mx_out[2])) {
    UpdateUserData(mx_out[2]);
  }

  /* Free temporary space */
  mxDestroyArray(mx_in[0]);
  mxDestroyArray(mx_in[1]);
  mxDestroyArray(mx_in[2]);
  mxDestroyArray(mx_in[3]);
  mxDestroyArray(mx_in[4]);
  mxDestroyArray(mx_in[5]);
  mxDestroyArray(mx_in[6]);
  mxDestroyArray(mx_in[7]);
  mxDestroyArray(mx_in[8]);
  mxDestroyArray(mx_out[0]);
  mxDestroyArray(mx_out[1]);
  mxDestroyArray(mx_out[2]);

  return(ret);

}

int mtlb_IdaSpilsPsetB(realtype tt, 
                       N_Vector yy, N_Vector yp,
                       N_Vector yyB, N_Vector ypB, N_Vector rrB, 
                       realtype c_jB, void *prec_dataB,
                       N_Vector tmp1B, N_Vector tmp2B, 
                       N_Vector tmp3B)
{
  mxArray *mx_in[10], *mx_out[2];
  int ret;

  /* Inputs to the Matlab function */
  mx_in[0] = mxCreateScalarDouble(-1.0);        /* type=-1: backward ODE */
  mx_in[1] = mxCreateScalarDouble(tt);          /* current t */
  mx_in[2] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current yy */
  mx_in[3] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current yp */
  mx_in[4] = mxCreateDoubleMatrix(NB,1,mxREAL); /* current yyB */
  mx_in[5] = mxCreateDoubleMatrix(NB,1,mxREAL); /* current ypB */
  mx_in[6] = mxCreateDoubleMatrix(NB,1,mxREAL); /* current rrB */
  mx_in[7] = mxCreateScalarDouble(c_jB);        /* current c_jB */
  mx_in[8] = mx_PSETfctB;                       /* matlab function handle */
  mx_in[9] = mx_data;                           /* matlab user data */
  
  /* Call matlab wrapper */
  GetData(yy, mxGetPr(mx_in[2]), N);
  GetData(yp, mxGetPr(mx_in[3]), N);
  GetData(yyB, mxGetPr(mx_in[4]), NB);
  GetData(ypB, mxGetPr(mx_in[5]), NB);
  GetData(rrB, mxGetPr(mx_in[6]), NB);

  mexCallMATLAB(2,mx_out,10,mx_in,"idm_pset");

  ret = (int)*mxGetPr(mx_out[0]);

  if (!mxIsEmpty(mx_out[1])) {
    UpdateUserData(mx_out[1]);
  }

  /* Free temporary space */
  mxDestroyArray(mx_in[0]);
  mxDestroyArray(mx_in[1]);
  mxDestroyArray(mx_in[2]);
  mxDestroyArray(mx_in[3]);
  mxDestroyArray(mx_in[4]);
  mxDestroyArray(mx_in[5]);
  mxDestroyArray(mx_in[6]);
  mxDestroyArray(mx_in[7]);
  mxDestroyArray(mx_out[0]);
  mxDestroyArray(mx_out[1]);

  return(ret);

}

int mtlb_IdaSpilsPsolB(realtype tt, 
                       N_Vector yy, N_Vector yp,
                       N_Vector yyB, N_Vector ypB, N_Vector rrB, 
                       N_Vector rvecB, N_Vector zvecB,
                       realtype c_jB, realtype deltaB,
                       void *prec_dataB, N_Vector tmpB)
{
  mxArray *mx_in[11], *mx_out[3];
  int ret;

  /* Inputs to the Matlab function */
  mx_in[0] = mxCreateScalarDouble(-1.0);        /* type=-1: backward ODE */
  mx_in[1] = mxCreateScalarDouble(tt);          /* current t */   
  mx_in[2] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current yy */
  mx_in[3] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current yp */
  mx_in[4] = mxCreateDoubleMatrix(NB,1,mxREAL); /* current yyB */
  mx_in[5] = mxCreateDoubleMatrix(NB,1,mxREAL); /* current ypB */
  mx_in[6] = mxCreateDoubleMatrix(NB,1,mxREAL); /* current rrB */
  mx_in[7] = mxCreateDoubleMatrix(NB,1,mxREAL); /* right hand side rB */
  mx_in[8] = mxCreateScalarDouble(c_jB);        /* current c_jB */   
  mx_in[9] = mx_PSOLfctB;                       /* matlab function handle */
  mx_in[10] = mx_data;                          /* matlab user data */
  
  /* Call matlab wrapper */
  GetData(yy, mxGetPr(mx_in[2]), N);
  GetData(yp, mxGetPr(mx_in[3]), N);
  GetData(yyB, mxGetPr(mx_in[4]), NB);
  GetData(ypB, mxGetPr(mx_in[5]), NB);
  GetData(rrB, mxGetPr(mx_in[6]), NB);
  GetData(rvecB, mxGetPr(mx_in[7]), NB);

  mexCallMATLAB(3,mx_out,11,mx_in,"idm_psol");

  PutData(zvecB, mxGetPr(mx_out[0]), NB);
  ret = (int)*mxGetPr(mx_out[1]);

  if (!mxIsEmpty(mx_out[2])) {
    UpdateUserData(mx_out[2]);
  }

  /* Free temporary space */
  mxDestroyArray(mx_in[0]);
  mxDestroyArray(mx_in[1]);
  mxDestroyArray(mx_in[2]);
  mxDestroyArray(mx_in[3]);
  mxDestroyArray(mx_in[4]);
  mxDestroyArray(mx_in[5]);
  mxDestroyArray(mx_in[6]);
  mxDestroyArray(mx_in[7]);
  mxDestroyArray(mx_in[8]);
  mxDestroyArray(mx_out[0]);
  mxDestroyArray(mx_out[1]);
  mxDestroyArray(mx_out[2]);

  return(ret);

}

int mtlb_IdaBBDglocB(long int NlocalB, realtype tt,
                     N_Vector yy, N_Vector yp, 
                     N_Vector yyB, N_Vector ypB, N_Vector gvalB,
                     void *res_dataB)
{
  mxArray *mx_in[8], *mx_out[3];
  int ret;

  /* Inputs to the Matlab function */
  mx_in[0] = mxCreateScalarDouble(-1.0);        /* type=-1: backward ODE */
  mx_in[1] = mxCreateScalarDouble(tt);          /* current t */
  mx_in[2] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current yy */
  mx_in[3] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current yp */
  mx_in[4] = mxCreateDoubleMatrix(NB,1,mxREAL); /* current yyB */
  mx_in[5] = mxCreateDoubleMatrix(NB,1,mxREAL); /* current ypB */
  mx_in[6] = mx_GLOCfctB;                       /* matlab function handle */
  mx_in[7] = mx_data;                           /* matlab user data */
  
  /* Call matlab wrapper */
  GetData(yy, mxGetPr(mx_in[2]), N);
  GetData(yp, mxGetPr(mx_in[3]), N);
  GetData(yyB, mxGetPr(mx_in[4]), NB);
  GetData(ypB, mxGetPr(mx_in[5]), NB);

  mexCallMATLAB(3,mx_out,8,mx_in,"idm_gloc");

  PutData(gvalB, mxGetPr(mx_out[0]), NB);

  ret = (int)*mxGetPr(mx_out[1]);

  if (!mxIsEmpty(mx_out[2])) {
    UpdateUserData(mx_out[2]);
  }

  /* Free temporary space */
  mxDestroyArray(mx_in[0]);
  mxDestroyArray(mx_in[1]);
  mxDestroyArray(mx_in[2]);
  mxDestroyArray(mx_in[3]);
  mxDestroyArray(mx_in[4]);
  mxDestroyArray(mx_in[5]);
  mxDestroyArray(mx_out[0]);
  mxDestroyArray(mx_out[1]);
  mxDestroyArray(mx_out[2]);

  return(ret);
}

int mtlb_IdaBBDgcomB(long int NlocalB, realtype tt,
                     N_Vector yy, N_Vector yp,
                     N_Vector yyB, N_Vector ypB,
                     void *res_dataB)
{
  mxArray *mx_in[8], *mx_out[2];
  int ret;

  /* Inputs to the Matlab function */
  mx_in[0] = mxCreateScalarDouble(-1.0);        /* type=-1: backward ODE */
  mx_in[1] = mxCreateScalarDouble(tt);          /* current t */
  mx_in[2] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current yy */
  mx_in[3] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current yp */
  mx_in[4] = mxCreateDoubleMatrix(NB,1,mxREAL); /* current yyB */
  mx_in[5] = mxCreateDoubleMatrix(NB,1,mxREAL); /* current ypB */
  mx_in[6] = mx_GCOMfctB;                       /* matlab function handle */
  mx_in[7] = mx_data;                           /* matlab user data */
  
  /* Call matlab wrapper */
  GetData(yy, mxGetPr(mx_in[2]), N);
  GetData(yp, mxGetPr(mx_in[3]), N);
  GetData(yyB, mxGetPr(mx_in[4]), NB);
  GetData(ypB, mxGetPr(mx_in[5]), NB);

  mexCallMATLAB(2,mx_out,8,mx_in,"idm_gcom");

  ret = (int)*mxGetPr(mx_out[0]);

  if (!mxIsEmpty(mx_out[1])) {
    UpdateUserData(mx_out[1]);
  }

  /* Free temporary space */
  mxDestroyArray(mx_in[0]);
  mxDestroyArray(mx_in[1]);
  mxDestroyArray(mx_in[2]);
  mxDestroyArray(mx_in[3]);
  mxDestroyArray(mx_in[4]);
  mxDestroyArray(mx_in[5]);
  mxDestroyArray(mx_out[0]);
  mxDestroyArray(mx_out[1]);

  return(ret);
}

/*
 * ---------------------------------------------------------------------------------
 * Wrapper around matlab monitor function
 * ---------------------------------------------------------------------------------
 */

void mtlb_IdaMonitor(int call, double t, 
                     N_Vector yy, N_Vector yp, 
                     N_Vector yQ, 
                     N_Vector *yyS, N_Vector *ypS)
{
  mxArray *mx_in[10], *mx_out[1];
  double *tmp_yyS, *tmp_ypS;
  int is;

  mx_in[0] = mxCreateScalarDouble(call);            /* call type (0:first, 1:interm. 2:last) */
  mx_in[1] = mxCreateScalarDouble(t);               /* current time */
  mx_in[2] = mxCreateDoubleMatrix(N,1,mxREAL);      /* current yy */
  mx_in[3] = mxCreateDoubleMatrix(N,1,mxREAL);      /* current yp */
  if (idm_quad) {
    mx_in[4] = mxCreateDoubleMatrix(Nq,1,mxREAL);   /* current quadratures */
  } else {
    mx_in[4] = mxCreateDoubleMatrix(0,0,mxREAL);
  }
  mx_in[5] = mxCreateScalarDouble(Ns);              /* number of sensitivities */
  if (idm_fsa) {
    mx_in[6] = mxCreateDoubleMatrix(N*Ns,1,mxREAL); /* current yyS */
    mx_in[7] = mxCreateDoubleMatrix(N*Ns,1,mxREAL); /* current ypS */
  } else {
    mx_in[6] = mxCreateDoubleMatrix(0,0,mxREAL);
    mx_in[7] = mxCreateDoubleMatrix(0,0,mxREAL);
  }
  mx_in[8] = mx_MONfct;                             /* Matlab monitor function */
  mx_in[9] = mx_MONdata;                            /* data for monitor function */

  if (call == 1) {

    GetData(yy, mxGetPr(mx_in[2]), N);
    GetData(yp, mxGetPr(mx_in[3]), N);

    if (idm_quad) {
      GetData(yQ, mxGetPr(mx_in[4]), Nq);
    }

    if (idm_fsa) {
      tmp_yyS = mxGetPr(mx_in[6]);
      tmp_ypS = mxGetPr(mx_in[7]);
      for (is=0; is<Ns; is++) {
        GetData(yyS[is], &tmp_yyS[is*N], N);
        GetData(ypS[is], &tmp_ypS[is*N], N);
      }
    }

  }

  mexCallMATLAB(1,mx_out,10,mx_in,"idm_monitor");

  if (!mxIsEmpty(mx_out[0])) {
    UpdateMonitorData(mx_out[0]);
  }

  mxDestroyArray(mx_in[0]);
  mxDestroyArray(mx_in[1]);
  mxDestroyArray(mx_in[2]);
  mxDestroyArray(mx_in[3]);
  mxDestroyArray(mx_in[4]);
  mxDestroyArray(mx_in[5]);
  mxDestroyArray(mx_in[6]);
  mxDestroyArray(mx_in[7]);
  mxDestroyArray(mx_out[0]);
}

void mtlb_IdaMonitorB(int call, double tB, N_Vector yyB, N_Vector ypB, N_Vector yQB)
{
  mxArray *mx_in[10], *mx_out[1];

  mx_in[0] = mxCreateScalarDouble(call);            /* 0: first, 1: interm. 2: last */
  mx_in[1] = mxCreateScalarDouble(tB);              /* current time */
  mx_in[2] = mxCreateDoubleMatrix(NB,1,mxREAL);     /* current yyB */
  mx_in[3] = mxCreateDoubleMatrix(NB,1,mxREAL);     /* current ypB */
  if (idm_quadB) {
    mx_in[4] = mxCreateDoubleMatrix(NqB,1,mxREAL);  /* current quadratures */
  } else {
    mx_in[4] = mxCreateDoubleMatrix(0,0,mxREAL);
  }
  mx_in[5] = mxCreateScalarDouble(0.0);             /* Ns is always zero here */
  mx_in[6] = mxCreateDoubleMatrix(0,0,mxREAL);      /* yyS is always empty here */
  mx_in[7] = mxCreateDoubleMatrix(0,0,mxREAL);      /* ypS is always empty here */
  mx_in[8] = mx_MONfctB;
  mx_in[9] = mx_MONdataB;

  if (call == 1) {
    
    GetData(yyB, mxGetPr(mx_in[2]), NB);
    GetData(ypB, mxGetPr(mx_in[3]), NB);

    if (idm_quadB)
      GetData(yQB, mxGetPr(mx_in[4]), NqB);
  }

  mexCallMATLAB(1,mx_out,10,mx_in,"idm_monitor");

  if (!mxIsEmpty(mx_out[0])) {
    UpdateMonitorDataB(mx_out[0]);
  }

  mxDestroyArray(mx_in[0]);
  mxDestroyArray(mx_in[1]);
  mxDestroyArray(mx_in[2]);
  mxDestroyArray(mx_in[3]);
  mxDestroyArray(mx_in[4]);
  mxDestroyArray(mx_in[5]);
  mxDestroyArray(mx_in[6]);
  mxDestroyArray(mx_in[7]);
  mxDestroyArray(mx_out[0]);
}


/*
 * ---------------------------------------------------------------------------------
 * Private functions to update the user data structures
 * ---------------------------------------------------------------------------------
 */

static void UpdateUserData(mxArray *mx_newdata)
{
  mexUnlock();
  mxDestroyArray(mx_data);
  mx_data = mxDuplicateArray(mx_newdata);
  mexMakeArrayPersistent(mx_data);
  mexLock();
}

static void UpdateMonitorData(mxArray *mx_newdata)
{
  mexUnlock();
  mxDestroyArray(mx_MONdata);
  mx_MONdata = mxDuplicateArray(mx_newdata);
  mexMakeArrayPersistent(mx_MONdata);
  mexLock();
}

static void UpdateMonitorDataB(mxArray *mx_newdata)
{
  mexUnlock();
  mxDestroyArray(mx_MONdataB);
  mx_MONdataB = mxDuplicateArray(mx_newdata);
  mexMakeArrayPersistent(mx_MONdataB);
  mexLock();
}

