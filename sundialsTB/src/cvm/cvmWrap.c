/*
 * -----------------------------------------------------------------
 * $Revision: 1.4 $
 * $Date: 2006-03-07 01:20:04 $
 * -----------------------------------------------------------------
 * Programmer: Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2005, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see sundials/cvodes/LICENSE.
 * -----------------------------------------------------------------
 * CVODES wrapper functions.
 * -----------------------------------------------------------------
 */

#include "cvm.h"

static void UpdateUserData(mxArray *mx_data);
static void UpdateMonitorData(mxArray *mx_data);
static void UpdateMonitorDataB(mxArray *mx_data);

booleantype cvm_quad;      /* Forward quadratures? */
booleantype cvm_quadB;     /* Backward quadratures? */
booleantype cvm_asa;       /* Adjoint sensitivity? */
booleantype cvm_fsa;       /* Forward sensitivity? */

/*
 * ---------------------------------------------------------------------------------
 * Redability replacements
 * ---------------------------------------------------------------------------------
 */

#define N           (cvm_Cdata->N) 
#define Nq          (cvm_Cdata->Nq) 
#define Ng          (cvm_Cdata->Ng) 
#define Ns          (cvm_Cdata->Ns) 
#define Nd          (cvm_Cdata->Nd) 
#define Nc          (cvm_Cdata->Nc) 
#define ls          (cvm_Cdata->ls) 
#define pm          (cvm_Cdata->pm) 
#define ism         (cvm_Cdata->ism) 
#define NB          (cvm_Cdata->NB) 
#define NqB         (cvm_Cdata->NqB) 
#define lsB         (cvm_Cdata->lsB) 
#define pmB         (cvm_Cdata->pmB) 

#define mx_data     (cvm_Mdata->mx_data)

#define mx_RHSfct   (cvm_Mdata->mx_RHSfct)
#define mx_QUADfct  (cvm_Mdata->mx_QUADfct)
#define mx_JACfct   (cvm_Mdata->mx_JACfct)
#define mx_PSETfct  (cvm_Mdata->mx_PSETfct)
#define mx_PSOLfct  (cvm_Mdata->mx_PSOLfct)
#define mx_GLOCfct  (cvm_Mdata->mx_GLOCfct)
#define mx_GCOMfct  (cvm_Mdata->mx_GCOMfct)
#define mx_Gfct     (cvm_Mdata->mx_Gfct)
#define mx_SRHSfct  (cvm_Mdata->mx_SRHSfct)

#define mx_RHSfctB  (cvm_Mdata->mx_RHSfctB)
#define mx_QUADfctB (cvm_Mdata->mx_QUADfctB)
#define mx_JACfctB  (cvm_Mdata->mx_JACfctB)
#define mx_PSETfctB (cvm_Mdata->mx_PSETfctB)
#define mx_PSOLfctB (cvm_Mdata->mx_PSOLfctB)
#define mx_GLOCfctB (cvm_Mdata->mx_GLOCfctB)
#define mx_GCOMfctB (cvm_Mdata->mx_GCOMfctB)

#define mx_MONfct   (cvm_Mdata->mx_MONfct)
#define mx_MONdata  (cvm_Mdata->mx_MONdata)

#define mx_MONfctB  (cvm_Mdata->mx_MONfctB)
#define mx_MONdataB (cvm_Mdata->mx_MONdataB)


/*
 * ---------------------------------------------------------------------------------
 * Error handler
 * ---------------------------------------------------------------------------------
 */

void mtlb_CVodeErrHandler(int error_code, 
                          const char *module, const char *function, 
                          char *msg, void *eh_data)
{
  char err_type[10];

  if (error_code == CV_WARNING)
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

int mtlb_CVodeRhs(realtype t, N_Vector y, N_Vector yd, void *f_data)
{
  mxArray *mx_in[5], *mx_out[3];
  int ret;

  /* Inputs to the Matlab function */
  mx_in[0] = mxCreateScalarDouble(1.0);        /* type=1: forward ODE */
  mx_in[1] = mxCreateScalarDouble(t);          /* current t */
  mx_in[2] = mxCreateDoubleMatrix(N,1,mxREAL); /* current y */
  mx_in[3] = mx_RHSfct;                   /* matlab function handle */ 
  mx_in[4] = mx_data;                     /* matlab user data */

  /* Call matlab wrapper */
  GetData(y, mxGetPr(mx_in[2]), N);

  mexCallMATLAB(3,mx_out,5,mx_in,"cvm_rhs");

  PutData(yd, mxGetPr(mx_out[0]), N);
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

int mtlb_CVodeQUADfct(realtype t, N_Vector y, N_Vector yQd, void *fQ_data)
{
  mxArray *mx_in[5], *mx_out[3];
  int ret;

  /* Inputs to the Matlab function */
  mx_in[0] = mxCreateScalarDouble(1.0);        /* type=1: forward ODE */
  mx_in[1] = mxCreateScalarDouble(t);          /* current t */
  mx_in[2] = mxCreateDoubleMatrix(N,1,mxREAL); /* current y */
  mx_in[3] = mx_QUADfct;                  /* matlab function handle */ 
  mx_in[4] = mx_data;                     /* matlab user data */

  /* Call matlab wrapper */
  GetData(y, mxGetPr(mx_in[2]), N);

  mexCallMATLAB(3,mx_out,5,mx_in,"cvm_rhsQ");

  PutData(yQd, mxGetPr(mx_out[0]), Nq);
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

int mtlb_CVodeGfct(realtype t, N_Vector y, double *g, void *g_data)
{
  double *gdata;
  int i, ret;
  mxArray *mx_in[4], *mx_out[3];
  
  /* Inputs to the Matlab function */
  mx_in[0] = mxCreateScalarDouble(t);          /* current t */
  mx_in[1] = mxCreateDoubleMatrix(N,1,mxREAL); /* current y */
  mx_in[2] = mx_Gfct;                     /* matlab function handle */
  mx_in[3] = mx_data;                     /* matlab user data */
  
  /* Call matlab wrapper */
  GetData(y, mxGetPr(mx_in[1]), N);

  mexCallMATLAB(3,mx_out,4,mx_in,"cvm_root");

  gdata = mxGetPr(mx_out[0]);
  for (i=0;i<Ng;i++) g[i] = gdata[i];

  ret = (int)*mxGetPr(mx_out[1]);

  if (!mxIsEmpty(mx_out[2])) {
    UpdateUserData(mx_out[2]);
  }

  /* Free temporary space */
  mxDestroyArray(mx_in[0]);
  mxDestroyArray(mx_in[1]);
  mxDestroyArray(mx_out[0]);
  mxDestroyArray(mx_out[1]);
  mxDestroyArray(mx_out[2]);

  return(ret);
}


int mtlb_CVodeDenseJac(long int Neq, DenseMat J, realtype t,
                       N_Vector y, N_Vector fy, void *jac_data,
                       N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  double *J_data;
  long int i;
  int ret;
  mxArray *mx_in[6], *mx_out[3];

  /* Inputs to the Matlab function */
  mx_in[0] = mxCreateScalarDouble(1.0);         /* type=1: forward ODE */
  mx_in[1] = mxCreateScalarDouble(t);           /* current t */  
  mx_in[2] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current y */
  mx_in[3] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current fy */
  mx_in[4] = mx_JACfct;                    /* matlab function handle */
  mx_in[5] = mx_data;                      /* matlab user data */
  
  /* Call matlab wrapper */
  GetData(y, mxGetPr(mx_in[2]), N);
  GetData(fy, mxGetPr(mx_in[3]), N);

  mexCallMATLAB(3,mx_out,6,mx_in,"cvm_djac");

  /* Extract data */
  J_data = mxGetPr(mx_out[0]);
  for (i=0;i<N;i++)
    memcpy(DENSE_COL(J,i), J_data + i*N, N*sizeof(double));

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

int mtlb_CVodeBandJac(long int Neq, long int mupper, long int mlower,
                      BandMat J, realtype t,
                      N_Vector y, N_Vector fy, void *jac_data,
                      N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  double *J_data;
  long int eband, i;
  int ret;
  mxArray *mx_in[6], *mx_out[3];

  /* Inputs to the Matlab function */
  mx_in[0] = mxCreateScalarDouble(1.0);         /* type=1: forward ODE */
  mx_in[1] = mxCreateScalarDouble(t);           /* current t */
  mx_in[2] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current y */
  mx_in[3] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current fy */
  mx_in[4] = mx_JACfct;                    /* matlab function handle */
  mx_in[5] = mx_data;                      /* matlab user data */
  
  /* Call matlab wrapper */
  GetData(y, mxGetPr(mx_in[2]), N);
  GetData(fy, mxGetPr(mx_in[3]), N);

  mexCallMATLAB(3,mx_out,6,mx_in,"cvm_bjac");

  /* Extract data */
  eband =  mupper + mlower + 1;
  J_data = mxGetPr(mx_out[0]);
  for (i=0;i<N;i++)
    memcpy(BAND_COL(J,i) - mupper, J_data + i*eband, eband*sizeof(double));

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


int mtlb_CVodeSpilsJac(N_Vector v, N_Vector Jv, realtype t,
                       N_Vector y, N_Vector fy,
                       void *jac_data, N_Vector tmp)
{
  mxArray *mx_in[7], *mx_out[3];
  int ret;

  /* Inputs to the Matlab function */
  mx_in[0] = mxCreateScalarDouble(1.0);         /* type=1: forward ODE */
  mx_in[1] = mxCreateScalarDouble(t);           /* current t */ 
  mx_in[2] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current y */
  mx_in[3] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current fy */
  mx_in[4] = mxCreateDoubleMatrix(N,1,mxREAL);  /* vector v */
  mx_in[5] = mx_JACfct;                    /* matlab function handle */
  mx_in[6] = mx_data;                      /* matlab user data */
  
  /* Call matlab wrapper */
  GetData(y, mxGetPr(mx_in[2]), N);
  GetData(fy, mxGetPr(mx_in[3]), N);
  GetData(v, mxGetPr(mx_in[4]), N);

  mexCallMATLAB(3,mx_out,7,mx_in,"cvm_jtv");

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
  mxDestroyArray(mx_out[0]);
  mxDestroyArray(mx_out[1]);
  mxDestroyArray(mx_out[2]);

  return(ret);
}


int mtlb_CVodeSpilsPset(realtype t, N_Vector y, N_Vector fy,
                        booleantype jok, booleantype *jcurPtr,
                        realtype gamma, void *P_data,
                        N_Vector tmp1, N_Vector tmp2,
                        N_Vector tmp3)
{
  mxArray *mx_in[8], *mx_out[3];
  int ret;

  /* Inputs to the Matlab function */
  mx_in[0] = mxCreateScalarDouble(1.0);         /* type=1: forward ODE */
  mx_in[1] = mxCreateScalarDouble(t);           /* current t */
  mx_in[2] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current y */
  mx_in[3] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current fy */
  mx_in[4] = mxCreateLogicalScalar(jok);        /* jok flag */
  mx_in[5] = mxCreateScalarDouble(gamma);       /* gamma value */
  mx_in[6] = mx_PSETfct;                   /* matlab function handle */
  mx_in[7] = mx_data;                      /* matlab user data */
  
  /* Call matlab wrapper */
  GetData(y, mxGetPr(mx_in[2]), N);
  GetData(fy, mxGetPr(mx_in[3]), N);

  mexCallMATLAB(3,mx_out,8,mx_in,"cvm_pset");

  *jcurPtr = mxIsLogicalScalarTrue(mx_out[0]);
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


int mtlb_CVodeSpilsPsol(realtype t, N_Vector y, N_Vector fy,
                        N_Vector r, N_Vector z,
                        realtype gamma, realtype delta,
                        int lr, void *P_data, N_Vector tmp)
{
  mxArray *mx_in[7], *mx_out[3];
  int ret;

  /* Inputs to the Matlab function */
  mx_in[0] = mxCreateScalarDouble(1.0);        /* type=1: forward ODE */
  mx_in[1] = mxCreateScalarDouble(t);          /* current t */   
  mx_in[2] = mxCreateDoubleMatrix(N,1,mxREAL); /* current y */
  mx_in[3] = mxCreateDoubleMatrix(N,1,mxREAL); /* current fy */
  mx_in[4] = mxCreateDoubleMatrix(N,1,mxREAL); /* right hand side r */
  mx_in[5] = mx_PSOLfct;                  /* matlab function handle */
  mx_in[6] = mx_data;                     /* matlab user data */
  
  /* Call matlab wrapper */
  GetData(y, mxGetPr(mx_in[2]), N);
  GetData(fy, mxGetPr(mx_in[3]), N);
  GetData(r, mxGetPr(mx_in[4]), N);

  mexCallMATLAB(3,mx_out,7,mx_in,"cvm_psol");

  PutData(z, mxGetPr(mx_out[0]), N);
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

int mtlb_CVodeBBDgloc(long int Nlocal, realtype t, N_Vector y,
                      N_Vector g, void *f_data)
{
  mxArray *mx_in[5], *mx_out[3];
  int ret;

  /* Inputs to the Matlab function */
  mx_in[0] = mxCreateScalarDouble(1.0);         /* type=1: forward ODE */
  mx_in[1] = mxCreateScalarDouble(t);           /* current t */
  mx_in[2] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current y */
  mx_in[3] = mx_GLOCfct;                   /* matlab function handle */
  mx_in[4] = mx_data;                      /* matlab user data */
  
  /* Call matlab wrapper */
  GetData(y, mxGetPr(mx_in[2]), N);

  mexCallMATLAB(3,mx_out,5,mx_in,"cvm_gloc");

  PutData(g, mxGetPr(mx_out[0]), N);
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

int mtlb_CVodeBBDgcom(long int Nlocal, realtype t, N_Vector y, void *f_data)
{
  mxArray *mx_in[5], *mx_out[2];
  int ret;

  /* Inputs to the Matlab function */
  mx_in[0] = mxCreateScalarDouble(1.0);         /* type=1: forward ODE */
  mx_in[1] = mxCreateScalarDouble(t);           /* current t */
  mx_in[2] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current y */
  mx_in[3] = mx_GCOMfct;                   /* matlab function handle */
  mx_in[4] = mx_data;                      /* matlab user data */
  
  /* Call matlab wrapper */
  GetData(y, mxGetPr(mx_in[3]), N);

  mexCallMATLAB(2,mx_out,5,mx_in,"cvm_gcom");

  ret = (int)*mxGetPr(mx_out[0]);

  if (!mxIsEmpty(mx_out[1])) {
    UpdateUserData(mx_out[1]);
  }

  /* Free temporary space */
  mxDestroyArray(mx_in[0]);
  mxDestroyArray(mx_in[1]);
  mxDestroyArray(mx_in[2]);
  mxDestroyArray(mx_out[0]);
  mxDestroyArray(mx_out[1]);

  return(ret);
}

/*
 * ----------------------------
 * FORWARD SENSITVITY FUNCTIONS
 * ----------------------------
 */

int mtlb_CVodeSensRhs(int Nsens, realtype t,
                      N_Vector y, N_Vector yd,
                      N_Vector *yS, N_Vector *ySd,
                      void *fS_data,
                      N_Vector tmp1, N_Vector tmp2)
{
  mxArray *mx_in[7], *mx_out[3];
  int is, ret;
  double *tmp;

  /* Inputs to the Matlab function */
  mx_in[0] = mxCreateScalarDouble(t);             /* current t */
  mx_in[1] = mxCreateDoubleMatrix(N,1,mxREAL);    /* current y */
  mx_in[2] = mxCreateDoubleMatrix(N,1,mxREAL);    /* current yd */
  mx_in[3] = mxCreateScalarDouble(Ns);            /* number of sensitivities */
  mx_in[4] = mxCreateDoubleMatrix(N*Ns,1,mxREAL); /* current yS */
  mx_in[5] = mx_SRHSfct;                          /* matlab function handle */      
  mx_in[6] = mx_data;                             /* matlab user data */
  
  /* Call matlab wrapper */
  GetData(y, mxGetPr(mx_in[1]), N);
  GetData(yd, mxGetPr(mx_in[2]), N);
  tmp = mxGetPr(mx_in[4]);
  for (is=0; is<Ns; is++)
    GetData(yS[is], &tmp[is*N], N);

  mexCallMATLAB(3,mx_out,7,mx_in,"cvm_rhsS");
  
  tmp = mxGetPr(mx_out[0]);

  for(is=0;is<Ns;is++)
    PutData(ySd[is], &tmp[is*N], N);

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


int mtlb_CVodeRhsB(realtype t, N_Vector y, N_Vector yB, N_Vector yBd, void *f_dataB)
{
  mxArray *mx_in[6], *mx_out[3];
  int ret;

  /* Inputs to the Matlab function */
  mx_in[0] = mxCreateScalarDouble(-1.0);        /* type=-1: backward ODE */
  mx_in[1] = mxCreateScalarDouble(t);           /* current t */
  mx_in[2] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current y */
  mx_in[3] = mxCreateDoubleMatrix(NB,1,mxREAL); /* current yB */
  mx_in[4] = mx_RHSfctB;                   /* matlab function handle */ 
  mx_in[5] = mx_data;                      /* matlab user data */
  
  /* Call matlab wrapper */
  GetData(y, mxGetPr(mx_in[2]), N);
  GetData(yB, mxGetPr(mx_in[3]), NB);

  mexCallMATLAB(3,mx_out,6,mx_in,"cvm_rhs");

  PutData(yBd, mxGetPr(mx_out[0]), NB);

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


int mtlb_CVodeQUADfctB(realtype t, N_Vector y, N_Vector yB, N_Vector yQBd, void *fQ_dataB)
{
  mxArray *mx_in[6], *mx_out[3];
  int ret;

  /* Inputs to the Matlab function */
  mx_in[0] = mxCreateScalarDouble(-1.0);        /* type=-1: backward ODE */
  mx_in[1] = mxCreateScalarDouble(t);           /* current t */
  mx_in[2] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current y */
  mx_in[3] = mxCreateDoubleMatrix(NB,1,mxREAL); /* current yB */
  mx_in[4] = mx_QUADfctB;                  /* matlab function handle */ 
  mx_in[5] = mx_data;                      /* matlab user data */

  /* Call matlab wrapper */
  GetData(y, mxGetPr(mx_in[2]), N);
  GetData(yB, mxGetPr(mx_in[3]), NB);

  mexCallMATLAB(3,mx_out,6,mx_in,"cvm_rhsQ");

  PutData(yQBd, mxGetPr(mx_out[0]), NqB);

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


int mtlb_CVodeDenseJacB(long int NeqB, DenseMat JB, realtype t,
                        N_Vector y, N_Vector yB, N_Vector fyB,
                        void *jac_dataB, N_Vector tmp1B,
                        N_Vector tmp2B, N_Vector tmp3B)
{
  double *JB_data;
  long int i;
  mxArray *mx_in[7], *mx_out[3];
  int ret;

  /* Inputs to the Matlab function */
  mx_in[0] = mxCreateScalarDouble(-1.0);        /* type=-1: backward ODE */
  mx_in[1] = mxCreateScalarDouble(t);           /* current t */  
  mx_in[2] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current y */
  mx_in[3] = mxCreateDoubleMatrix(NB,1,mxREAL); /* current yB */
  mx_in[4] = mxCreateDoubleMatrix(NB,1,mxREAL); /* current fyB */
  mx_in[5] = mx_JACfctB;                   /* matlab function handle */
  mx_in[6] = mx_data;                      /* matlab user data */
  
  /* Call matlab wrapper */
  GetData(y, mxGetPr(mx_in[2]), N);
  GetData(yB, mxGetPr(mx_in[3]), NB);
  GetData(fyB, mxGetPr(mx_in[4]), NB);

  mexCallMATLAB(3,mx_out,7,mx_in,"cvm_djac");

  JB_data = mxGetPr(mx_out[0]);
  for (i=0;i<NB;i++)
    memcpy(DENSE_COL(JB,i), JB_data + i*NB, NB*sizeof(double));

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
  mxDestroyArray(mx_out[0]);
  mxDestroyArray(mx_out[1]);
  mxDestroyArray(mx_out[2]);

  return(ret);
}


int mtlb_CVodeBandJacB(long int NeqB, long int mupperB,
                       long int mlowerB, BandMat JB,
                       realtype t, N_Vector y,
                       N_Vector yB, N_Vector fyB,
                       void *jac_dataB, N_Vector tmp1B,
                       N_Vector tmp2B, N_Vector tmp3B)
{
  double *JB_data;
  long int ebandB, i;
  mxArray *mx_in[7], *mx_out[3];
  int ret;

  /* Inputs to the Matlab function */
  mx_in[0] = mxCreateScalarDouble(-1.0);        /* type=-1: backward ODE */
  mx_in[1] = mxCreateScalarDouble(t);           /* current t */
  mx_in[2] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current y */
  mx_in[3] = mxCreateDoubleMatrix(NB,1,mxREAL); /* current yB */
  mx_in[4] = mxCreateDoubleMatrix(NB,1,mxREAL); /* current fyB */
  mx_in[5] = mx_JACfctB;                   /* matlab function handle */
  mx_in[6] = mx_data;                      /* matlab user data */
  
  /* Call matlab wrapper */
  GetData(y, mxGetPr(mx_in[2]), N);
  GetData(yB, mxGetPr(mx_in[3]), NB);
  GetData(fyB, mxGetPr(mx_in[4]), NB);

  mexCallMATLAB(3,mx_out,7,mx_in,"cvm_bjac");

  ebandB =  mupperB + mlowerB + 1;
  JB_data = mxGetPr(mx_out[0]);
  for (i=0;i<NB;i++)
    memcpy(BAND_COL(JB,i) - mupperB, JB_data + i*ebandB, ebandB*sizeof(double));
    
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
  mxDestroyArray(mx_out[0]);
  mxDestroyArray(mx_out[1]);
  mxDestroyArray(mx_out[2]);

  return(ret);
}


int mtlb_CVodeSpilsJacB(N_Vector vB, N_Vector JvB, realtype t,
                        N_Vector y, N_Vector yB, N_Vector fyB,
                        void *jac_dataB, N_Vector tmpB)
{
  mxArray *mx_in[8], *mx_out[3];
  int ret;

  /* Inputs to the Matlab function */
  mx_in[0] = mxCreateScalarDouble(-1.0);        /* type=-1: backward ODE */
  mx_in[1] = mxCreateScalarDouble(t);           /* current t */ 
  mx_in[2] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current y */
  mx_in[3] = mxCreateDoubleMatrix(NB,1,mxREAL); /* current yB */
  mx_in[4] = mxCreateDoubleMatrix(NB,1,mxREAL); /* current fyB */
  mx_in[5] = mxCreateDoubleMatrix(NB,1,mxREAL); /* vector vB */
  mx_in[6] = mx_JACfctB;                   /* matlab function handle */
  mx_in[7] = mx_data;                      /* matlab user data */
  
  /* Call matlab wrapper */
  GetData(y, mxGetPr(mx_in[2]), N);
  GetData(yB, mxGetPr(mx_in[3]), NB);
  GetData(fyB, mxGetPr(mx_in[4]), NB);
  GetData(vB, mxGetPr(mx_in[5]), NB);

  mexCallMATLAB(3,mx_out,8,mx_in,"cvm_jtv");

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
  mxDestroyArray(mx_out[0]);
  mxDestroyArray(mx_out[1]);
  mxDestroyArray(mx_out[2]);

  return(ret);

}


int mtlb_CVodeSpilsPsetB(realtype t, N_Vector y,
                         N_Vector yB, N_Vector fyB,
                         booleantype jokB,
                         booleantype *jcurPtrB, realtype gammaB,
                         void *P_dataB,
                         N_Vector tmp1B, N_Vector tmp2B,
                         N_Vector tmp3B)
{
  mxArray *mx_in[9], *mx_out[3];
  int ret;

  /* Inputs to the Matlab function */
  mx_in[0] = mxCreateScalarDouble(-1.0);        /* type=-1: backward ODE */
  mx_in[1] = mxCreateScalarDouble(t);           /* current t */
  mx_in[2] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current y */
  mx_in[3] = mxCreateDoubleMatrix(NB,1,mxREAL); /* current yB */
  mx_in[4] = mxCreateDoubleMatrix(NB,1,mxREAL); /* current fyB */
  mx_in[5] = mxCreateLogicalScalar(jokB);       /* jokB flag */
  mx_in[6] = mxCreateScalarDouble(gammaB);      /* gammaB value */
  mx_in[7] = mx_PSETfctB;                  /* matlab function handle */
  mx_in[8] = mx_data;                      /* matlab user data */
  
  /* Call matlab wrapper */
  GetData(y, mxGetPr(mx_in[2]), N);
  GetData(yB, mxGetPr(mx_in[3]), NB);
  GetData(fyB, mxGetPr(mx_in[4]), NB);

  mexCallMATLAB(3,mx_out,9,mx_in,"cvm_pset");

  *jcurPtrB = mxIsLogicalScalarTrue(mx_out[0]);
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


int mtlb_CVodeSpilsPsolB(realtype t, N_Vector y,
                         N_Vector yB, N_Vector fyB,
                         N_Vector rB, N_Vector zB,
                         realtype gammaB, realtype deltaB,
                         int lrB, void *P_dataB, N_Vector tmpB)
{
  mxArray *mx_in[8], *mx_out[3];
  int ret;

  /* Inputs to the Matlab function */
  mx_in[0] = mxCreateScalarDouble(-1.0);        /* type=-1: backward ODE */
  mx_in[1] = mxCreateScalarDouble(t);           /* current t */   
  mx_in[2] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current y */
  mx_in[3] = mxCreateDoubleMatrix(NB,1,mxREAL); /* current yB */
  mx_in[4] = mxCreateDoubleMatrix(NB,1,mxREAL); /* current fyB */
  mx_in[5] = mxCreateDoubleMatrix(NB,1,mxREAL); /* right hand side rB */
  mx_in[6] = mx_PSOLfctB;                  /* matlab function handle */
  mx_in[7] = mx_data;                      /* matlab user data */
  
  /* Call matlab wrapper */
  GetData(y, mxGetPr(mx_in[2]), N);
  GetData(yB, mxGetPr(mx_in[3]), NB);
  GetData(fyB, mxGetPr(mx_in[4]), NB);
  GetData(rB, mxGetPr(mx_in[5]), NB);

  mexCallMATLAB(3,mx_out,8,mx_in,"cvm_psol");

  PutData(zB, mxGetPr(mx_out[0]), NB);
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

int mtlb_CVodeBBDglocB(long int NlocalB, realtype t, N_Vector y,
                       N_Vector yB, N_Vector gB, void *f_dataB)
{
  mxArray *mx_in[6], *mx_out[3];
  int ret;

  /* Inputs to the Matlab function */
  mx_in[0] = mxCreateScalarDouble(-1.0);        /* type=-1: backward ODE */
  mx_in[1] = mxCreateScalarDouble(t);           /* current t */
  mx_in[2] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current y */
  mx_in[3] = mxCreateDoubleMatrix(NB,1,mxREAL); /* current yB */
  mx_in[4] = mx_GLOCfctB;                  /* matlab function handle */
  mx_in[5] = mx_data;                      /* matlab user data */
  
  /* Call matlab wrapper */
  GetData(y, mxGetPr(mx_in[2]), N);
  GetData(yB, mxGetPr(mx_in[3]), NB);

  mexCallMATLAB(3,mx_out,6,mx_in,"cvm_gloc");

  PutData(gB, mxGetPr(mx_out[0]), NB);

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

int mtlb_CVodeBBDgcomB(long int NlocalB, realtype t, N_Vector y, 
                       N_Vector yB, void *f_dataB)
{
  mxArray *mx_in[6], *mx_out[2];
  int ret;

  /* Inputs to the Matlab function */
  mx_in[0] = mxCreateScalarDouble(-1.0);        /* type=-1: backward ODE */
  mx_in[1] = mxCreateScalarDouble(t);           /* current t */
  mx_in[2] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current y */
  mx_in[3] = mxCreateDoubleMatrix(NB,1,mxREAL); /* current yB */
  mx_in[4] = mx_GCOMfctB;                       /* matlab function handle */
  mx_in[5] = mx_data;                           /* matlab user data */
  
  /* Call matlab wrapper */
  GetData(y, mxGetPr(mx_in[2]), N);
  GetData(yB, mxGetPr(mx_in[3]), NB);

  mexCallMATLAB(2,mx_out,6,mx_in,"cvm_gcom");

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
 * ---------------------------------------------------------------------------------
 * Wrapper around matlab monitor function
 * ---------------------------------------------------------------------------------
 */

void mtlb_CVodeMonitor(int call, double t, N_Vector y, N_Vector yQ, N_Vector *yS)
{
  mxArray *mx_in[8], *mx_out[1];
  double *tmp;
  int is;

  mx_in[0] = mxCreateScalarDouble(call);            /* call type (0:first, 1:interm. 2:last) */
  mx_in[1] = mxCreateScalarDouble(t);               /* current time */
  mx_in[2] = mxCreateDoubleMatrix(N,1,mxREAL);      /* current solution */
  if (cvm_quad)
    mx_in[3] = mxCreateDoubleMatrix(Nq,1,mxREAL);   /* current quadratures */
  else
    mx_in[3] = mxCreateDoubleMatrix(0,0,mxREAL);
  mx_in[4] = mxCreateScalarDouble(Ns);              /* number of sensitivities */
  if (cvm_fsa)
    mx_in[5] = mxCreateDoubleMatrix(N*Ns,1,mxREAL); /* current sensitivities */
  else
    mx_in[5] = mxCreateDoubleMatrix(0,0,mxREAL);
  mx_in[6] = mx_MONfct;                             /* Matlab monitor function */
  mx_in[7] = mx_MONdata;                            /* data for monitor function */

  if (call == 1) {

    GetData(y, mxGetPr(mx_in[2]), N);

    if (cvm_quad)
      GetData(yQ, mxGetPr(mx_in[3]), Nq);

    if (cvm_fsa) {
      tmp = mxGetPr(mx_in[5]);
      for (is=0; is<Ns; is++)
        GetData(yS[is], &tmp[is*N], N);
    }

  }

  mexCallMATLAB(1,mx_out,8,mx_in,"cvm_monitor");

  if (!mxIsEmpty(mx_out[0])) {
    UpdateMonitorData(mx_out[0]);
  }

  mxDestroyArray(mx_in[0]);
  mxDestroyArray(mx_in[1]);
  mxDestroyArray(mx_in[2]);
  mxDestroyArray(mx_in[3]);
  mxDestroyArray(mx_in[4]);
  mxDestroyArray(mx_in[5]);
  mxDestroyArray(mx_out[0]);
}

void mtlb_CVodeMonitorB(int call, double tB, N_Vector yB, N_Vector yQB)
{
  mxArray *mx_in[8], *mx_out[1];
  double *tmp;

  mx_in[0] = mxCreateScalarDouble(call);            /* 0: first, 1: interm. 2: last */
  mx_in[1] = mxCreateScalarDouble(tB);              /* current time */
  mx_in[2] = mxCreateDoubleMatrix(NB,1,mxREAL);     /* current solution */
  if (cvm_quadB)
    mx_in[3] = mxCreateDoubleMatrix(NqB,1,mxREAL);  /* current quadratures */
  else
    mx_in[3] = mxCreateDoubleMatrix(0,0,mxREAL);
  mx_in[4] = mxCreateScalarDouble(0.0);             /* Ns is always zero here */
  mx_in[5] = mxCreateDoubleMatrix(0,0,mxREAL);      /* yS is always empty here */
  mx_in[6] = mx_MONfctB;
  mx_in[7] = mx_MONdataB;

  if (call == 1) {
    
    GetData(yB, mxGetPr(mx_in[2]), NB);

    if (cvm_quadB)
      GetData(yQB, mxGetPr(mx_in[3]), NqB);
  }

  mexCallMATLAB(1,mx_out,8,mx_in,"cvm_monitor");

  if (!mxIsEmpty(mx_out[0])) {
    UpdateMonitorDataB(mx_out[0]);
  }

  mxDestroyArray(mx_in[0]);
  mxDestroyArray(mx_in[1]);
  mxDestroyArray(mx_in[2]);
  mxDestroyArray(mx_in[3]);
  mxDestroyArray(mx_in[4]);
  mxDestroyArray(mx_in[5]);
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

