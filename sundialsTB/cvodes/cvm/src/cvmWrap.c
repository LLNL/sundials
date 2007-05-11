/*
 * -----------------------------------------------------------------
 * $Revision: 1.8 $
 * $Date: 2007-05-11 21:42:53 $
 * -----------------------------------------------------------------
 * Programmer: Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2005, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see sundials-x.y.z/src/cvodes/LICENSE.
 * -----------------------------------------------------------------
 * CVODES wrapper functions.
 * -----------------------------------------------------------------
 */

#include "cvm.h"
#include "nvm.h"

static void UpdateUserData(mxArray *mx_data, cvmInterfaceData cvmData);
static void UpdateMonitorData(mxArray *mx_data, cvmInterfaceData cvmData);
static void UpdateMonitorDataB(mxArray *mx_data, cvmInterfaceData cvmData);

/*
 * ---------------------------------------------------------------------------------
 * Redability replacements
 * ---------------------------------------------------------------------------------
 */


#define mx_data     (cvmData->mx_data)


#define fsa         (fwdPb->Fsa)
#define quadr       (fwdPb->Quadr)
#define N           (fwdPb->n) 
#define Nq          (fwdPb->nq) 
#define Ng          (fwdPb->ng) 
#define Ns          (fwdPb->ns) 
#define ls          (fwdPb->LS) 
#define pm          (fwdPb->PM)

#define mx_RHSfct   (fwdPb->RHSfct)
#define mx_QUADfct  (fwdPb->QUADfct)
#define mx_JACfct   (fwdPb->JACfct)
#define mx_PSETfct  (fwdPb->PSETfct)
#define mx_PSOLfct  (fwdPb->PSOLfct)
#define mx_GLOCfct  (fwdPb->GLOCfct)
#define mx_GCOMfct  (fwdPb->GCOMfct)
#define mx_Gfct     (fwdPb->Gfct)
#define mx_SRHSfct  (fwdPb->SRHSfct)

#define mx_MONfct   (fwdPb->MONfct)
#define mx_MONdata  (fwdPb->MONdata)



#define quadrB      (bckPb->Quadr)
#define NB          (bckPb->n) 
#define NqB         (bckPb->nq) 
#define lsB         (bckPb->LS) 
#define pmB         (bckPb->PM) 

#define mx_RHSfctB  (bckPb->RHSfct)
#define mx_QUADfctB (bckPb->QUADfct)
#define mx_JACfctB  (bckPb->JACfct)
#define mx_PSETfctB (bckPb->PSETfct)
#define mx_PSOLfctB (bckPb->PSOLfct)
#define mx_GLOCfctB (bckPb->GLOCfct)
#define mx_GCOMfctB (bckPb->GCOMfct)

#define mx_MONfctB  (bckPb->MONfct)
#define mx_MONdataB (bckPb->MONdata)


/*
 * ---------------------------------------------------------------------------------
 * FORWARD PROBLEMS
 * ---------------------------------------------------------------------------------
 */

int mtlb_CVodeRhs(realtype t, N_Vector y, N_Vector yd, void *f_data)
{
  cvmInterfaceData cvmData;
  cvmPbData fwdPb;
  mxArray *mx_in[4], *mx_out[3];
  int ret;

  /* Extract global interface data from user-data */
  cvmData = (cvmInterfaceData) f_data;
  fwdPb = cvmData->fwdPb;

  /* Inputs to the Matlab function */
  mx_in[0] = mxCreateScalarDouble(t);          /* current t */
  mx_in[1] = mxCreateDoubleMatrix(N,1,mxREAL); /* current y */
  mx_in[2] = mx_RHSfct;                        /* matlab function handle */ 
  mx_in[3] = mx_data;                          /* matlab user data */

  /* Call matlab wrapper */
  GetData(y, mxGetPr(mx_in[1]), N);

  mexCallMATLAB(3,mx_out,4,mx_in,"cvm_rhs");

  PutData(yd, mxGetPr(mx_out[0]), N);
  ret = (int)*mxGetPr(mx_out[1]);

  if (!mxIsEmpty(mx_out[2])) {
    UpdateUserData(mx_out[2], cvmData);
  }

  /* Free temporary space */
  mxDestroyArray(mx_in[0]);
  mxDestroyArray(mx_in[1]);
  mxDestroyArray(mx_out[0]);
  mxDestroyArray(mx_out[1]);
  mxDestroyArray(mx_out[2]);

  return(ret);
}

int mtlb_CVodeQUADfct(realtype t, N_Vector y, N_Vector yQd, void *f_data)
{
  cvmInterfaceData cvmData;
  cvmPbData fwdPb;
  mxArray *mx_in[4], *mx_out[3];
  int ret;

  /* Extract global interface data from user-data */
  cvmData = (cvmInterfaceData) f_data;
  fwdPb = cvmData->fwdPb;

  /* Inputs to the Matlab function */
  mx_in[0] = mxCreateScalarDouble(t);          /* current t */
  mx_in[1] = mxCreateDoubleMatrix(N,1,mxREAL); /* current y */
  mx_in[2] = mx_QUADfct;                       /* matlab function handle */ 
  mx_in[3] = mx_data;                          /* matlab user data */

  /* Call matlab wrapper */
  GetData(y, mxGetPr(mx_in[1]), N);

  mexCallMATLAB(3,mx_out,4,mx_in,"cvm_rhsQ");

  PutData(yQd, mxGetPr(mx_out[0]), Nq);
  ret = (int)*mxGetPr(mx_out[1]);

  if (!mxIsEmpty(mx_out[2])) {
    UpdateUserData(mx_out[2], cvmData);
  }

  /* Free temporary space */
  mxDestroyArray(mx_in[0]);
  mxDestroyArray(mx_in[1]);
  mxDestroyArray(mx_out[0]);
  mxDestroyArray(mx_out[1]);
  mxDestroyArray(mx_out[2]);

  return(ret);
}

int mtlb_CVodeGfct(realtype t, N_Vector y, double *g, void *f_data)
{
  cvmInterfaceData cvmData;
  cvmPbData fwdPb;
  double *gdata;
  int i, ret;
  mxArray *mx_in[4], *mx_out[3];
  
  /* Extract global interface data from user-data */
  cvmData = (cvmInterfaceData) f_data;
  fwdPb = cvmData->fwdPb;

  /* Inputs to the Matlab function */
  mx_in[0] = mxCreateScalarDouble(t);          /* current t */
  mx_in[1] = mxCreateDoubleMatrix(N,1,mxREAL); /* current y */
  mx_in[2] = mx_Gfct;                          /* matlab function handle */
  mx_in[3] = mx_data;                          /* matlab user data */
  
  /* Call matlab wrapper */
  GetData(y, mxGetPr(mx_in[1]), N);

  mexCallMATLAB(3,mx_out,4,mx_in,"cvm_root");

  gdata = mxGetPr(mx_out[0]);
  for (i=0;i<Ng;i++) g[i] = gdata[i];

  ret = (int)*mxGetPr(mx_out[1]);

  if (!mxIsEmpty(mx_out[2])) {
    UpdateUserData(mx_out[2], cvmData);
  }

  /* Free temporary space */
  mxDestroyArray(mx_in[0]);
  mxDestroyArray(mx_in[1]);
  mxDestroyArray(mx_out[0]);
  mxDestroyArray(mx_out[1]);
  mxDestroyArray(mx_out[2]);

  return(ret);
}


int mtlb_CVodeDenseJac(int Neq, realtype t,
                       N_Vector y, N_Vector fy, 
                       DlsMat J, void *f_data,
                       N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  cvmInterfaceData cvmData;
  cvmPbData fwdPb;
  double *J_data;
  long int i;
  int ret;
  mxArray *mx_in[5], *mx_out[3];

  /* Extract global interface data from user-data */
  cvmData = (cvmInterfaceData) f_data;
  fwdPb = cvmData->fwdPb;

  /* Inputs to the Matlab function */
  mx_in[0] = mxCreateScalarDouble(t);           /* current t */  
  mx_in[1] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current y */
  mx_in[2] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current fy */
  mx_in[3] = mx_JACfct;                         /* matlab function handle */
  mx_in[4] = mx_data;                           /* matlab user data */
  
  /* Call matlab wrapper */
  GetData(y, mxGetPr(mx_in[1]), N);
  GetData(fy, mxGetPr(mx_in[2]), N);

  mexCallMATLAB(3,mx_out,5,mx_in,"cvm_djac");

  /* Extract data */
  J_data = mxGetPr(mx_out[0]);
  for (i=0;i<N;i++)
    memcpy(DENSE_COL(J,i), J_data + i*N, N*sizeof(double));

  ret = (int)*mxGetPr(mx_out[1]);

  if (!mxIsEmpty(mx_out[2])) {
    UpdateUserData(mx_out[2], cvmData);
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

int mtlb_CVodeBandJac(int Neq, int mupper, int mlower, realtype t,
                      N_Vector y, N_Vector fy, 
                      DlsMat J, void *f_data,
                      N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  cvmInterfaceData cvmData;
  cvmPbData fwdPb;
  double *J_data;
  int eband, i, ret;
  mxArray *mx_in[5], *mx_out[3];

  /* Extract global interface data from user-data */
  cvmData = (cvmInterfaceData) f_data;
  fwdPb = cvmData->fwdPb;

  /* Inputs to the Matlab function */
  mx_in[0] = mxCreateScalarDouble(t);           /* current t */
  mx_in[1] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current y */
  mx_in[2] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current fy */
  mx_in[3] = mx_JACfct;                         /* matlab function handle */
  mx_in[4] = mx_data;                           /* matlab user data */
  
  /* Call matlab wrapper */
  GetData(y, mxGetPr(mx_in[1]), N);
  GetData(fy, mxGetPr(mx_in[2]), N);

  mexCallMATLAB(3,mx_out,5,mx_in,"cvm_bjac");

  /* Extract data */
  eband =  mupper + mlower + 1;
  J_data = mxGetPr(mx_out[0]);
  for (i=0;i<N;i++)
    memcpy(BAND_COL(J,i) - mupper, J_data + i*eband, eband*sizeof(double));

  ret = (int)*mxGetPr(mx_out[1]);
  
  if (!mxIsEmpty(mx_out[2])) {
    UpdateUserData(mx_out[2], cvmData);
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


int mtlb_CVodeSpilsJac(N_Vector v, N_Vector Jv, realtype t,
                       N_Vector y, N_Vector fy,
                       void *f_data, N_Vector tmp)
{
  cvmInterfaceData cvmData;
  cvmPbData fwdPb;
  mxArray *mx_in[6], *mx_out[3];
  int ret;

  /* Extract global interface data from user-data */
  cvmData = (cvmInterfaceData) f_data;
  fwdPb = cvmData->fwdPb;

  /* Inputs to the Matlab function */
  mx_in[0] = mxCreateScalarDouble(t);           /* current t */ 
  mx_in[1] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current y */
  mx_in[2] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current fy */
  mx_in[3] = mxCreateDoubleMatrix(N,1,mxREAL);  /* vector v */
  mx_in[4] = mx_JACfct;                         /* matlab function handle */
  mx_in[5] = mx_data;                           /* matlab user data */
  
  /* Call matlab wrapper */
  GetData(y, mxGetPr(mx_in[1]), N);
  GetData(fy, mxGetPr(mx_in[2]), N);
  GetData(v, mxGetPr(mx_in[3]), N);

  mexCallMATLAB(3,mx_out,6,mx_in,"cvm_jtv");

  PutData(Jv, mxGetPr(mx_out[0]), N);
  ret = (int)*mxGetPr(mx_out[1]);

  if (!mxIsEmpty(mx_out[2])) {
    UpdateUserData(mx_out[2], cvmData);
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


int mtlb_CVodeSpilsPset(realtype t, N_Vector y, N_Vector fy,
                        booleantype jok, booleantype *jcurPtr,
                        realtype gamma, void *f_data,
                        N_Vector tmp1, N_Vector tmp2,
                        N_Vector tmp3)
{
  cvmInterfaceData cvmData;
  cvmPbData fwdPb;
  mxArray *mx_in[7], *mx_out[3];
  int ret;

  /* Extract global interface data from user-data */
  cvmData = (cvmInterfaceData) f_data;
  fwdPb = cvmData->fwdPb;

  /* Inputs to the Matlab function */
  mx_in[0] = mxCreateScalarDouble(t);           /* current t */
  mx_in[1] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current y */
  mx_in[2] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current fy */
  mx_in[3] = mxCreateLogicalScalar(jok);        /* jok flag */
  mx_in[4] = mxCreateScalarDouble(gamma);       /* gamma value */
  mx_in[5] = mx_PSETfct;                        /* matlab function handle */
  mx_in[6] = mx_data;                           /* matlab user data */
  
  /* Call matlab wrapper */
  GetData(y, mxGetPr(mx_in[1]), N);
  GetData(fy, mxGetPr(mx_in[2]), N);

  mexCallMATLAB(3,mx_out,7,mx_in,"cvm_pset");

  *jcurPtr = mxIsLogicalScalarTrue(mx_out[0]);
  ret = (int)*mxGetPr(mx_out[1]);

  if (!mxIsEmpty(mx_out[2])) {
    UpdateUserData(mx_out[2], cvmData);
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


int mtlb_CVodeSpilsPsol(realtype t, N_Vector y, N_Vector fy,
                        N_Vector r, N_Vector z,
                        realtype gamma, realtype delta,
                        int lr, void *f_data, N_Vector tmp)
{
  cvmInterfaceData cvmData;
  cvmPbData fwdPb;
  mxArray *mx_in[6], *mx_out[3];
  int ret;

  /* Extract global interface data from user-data */
  cvmData = (cvmInterfaceData) f_data;
  fwdPb = cvmData->fwdPb;

  /* Inputs to the Matlab function */
  mx_in[0] = mxCreateScalarDouble(t);          /* current t */   
  mx_in[1] = mxCreateDoubleMatrix(N,1,mxREAL); /* current y */
  mx_in[2] = mxCreateDoubleMatrix(N,1,mxREAL); /* current fy */
  mx_in[3] = mxCreateDoubleMatrix(N,1,mxREAL); /* right hand side r */
  mx_in[4] = mx_PSOLfct;                       /* matlab function handle */
  mx_in[5] = mx_data;                          /* matlab user data */
  
  /* Call matlab wrapper */
  GetData(y, mxGetPr(mx_in[1]), N);
  GetData(fy, mxGetPr(mx_in[2]), N);
  GetData(r, mxGetPr(mx_in[3]), N);

  mexCallMATLAB(3,mx_out,6,mx_in,"cvm_psol");

  PutData(z, mxGetPr(mx_out[0]), N);
  ret = (int)*mxGetPr(mx_out[1]);

  if (!mxIsEmpty(mx_out[2])) {
    UpdateUserData(mx_out[2], cvmData);
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

/*
 * ----------------------------
 * BBD PRECONDITONER FUNCTIONS
 * ----------------------------
 */

int mtlb_CVodeBBDgloc(int Nlocal, realtype t, N_Vector y,
                      N_Vector g, void *f_data)
{
  cvmInterfaceData cvmData;
  cvmPbData fwdPb;
  mxArray *mx_in[4], *mx_out[3];
  int ret;

  /* Extract global interface data from user-data */
  cvmData = (cvmInterfaceData) f_data;
  fwdPb = cvmData->fwdPb;

  /* Inputs to the Matlab function */
  mx_in[0] = mxCreateScalarDouble(t);           /* current t */
  mx_in[1] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current y */
  mx_in[2] = mx_GLOCfct;                        /* matlab function handle */
  mx_in[3] = mx_data;                           /* matlab user data */
  
  /* Call matlab wrapper */
  GetData(y, mxGetPr(mx_in[1]), N);

  mexCallMATLAB(3,mx_out,4,mx_in,"cvm_gloc");

  PutData(g, mxGetPr(mx_out[0]), N);
  ret = (int)*mxGetPr(mx_out[1]);

  if (!mxIsEmpty(mx_out[2])) {
    UpdateUserData(mx_out[2], cvmData);
  }

  /* Free temporary space */
  mxDestroyArray(mx_in[0]);
  mxDestroyArray(mx_in[1]);
  mxDestroyArray(mx_out[0]);
  mxDestroyArray(mx_out[1]);
  mxDestroyArray(mx_out[2]);

  return(ret);
}

int mtlb_CVodeBBDgcom(int Nlocal, realtype t, N_Vector y, void *f_data)
{
  cvmInterfaceData cvmData;
  cvmPbData fwdPb;
  mxArray *mx_in[4], *mx_out[2];
  int ret;

  /* Extract global interface data from user-data */
  cvmData = (cvmInterfaceData) f_data;
  fwdPb = cvmData->fwdPb;

  /* Inputs to the Matlab function */
  mx_in[0] = mxCreateScalarDouble(t);           /* current t */
  mx_in[1] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current y */
  mx_in[2] = mx_GCOMfct;                        /* matlab function handle */
  mx_in[3] = mx_data;                           /* matlab user data */
  
  /* Call matlab wrapper */
  GetData(y, mxGetPr(mx_in[1]), N);

  mexCallMATLAB(2,mx_out,4,mx_in,"cvm_gcom");

  ret = (int)*mxGetPr(mx_out[0]);

  if (!mxIsEmpty(mx_out[1])) {
    UpdateUserData(mx_out[1], cvmData);
  }

  /* Free temporary space */
  mxDestroyArray(mx_in[0]);
  mxDestroyArray(mx_in[1]);
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
                      void *f_data,
                      N_Vector tmp1, N_Vector tmp2)
{
  cvmInterfaceData cvmData;
  cvmPbData fwdPb;
  mxArray *mx_in[7], *mx_out[3];
  int is, ret;
  double *tmp;

  /* Extract global interface data from user-data */
  cvmData = (cvmInterfaceData) f_data;
  fwdPb = cvmData->fwdPb;

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
    UpdateUserData(mx_out[2], cvmData);
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
  cvmInterfaceData cvmData;
  cvmPbData fwdPb, bckPb;
  mxArray *mx_in[6], *mx_out[3];
  int ret;

  /* Extract global interface data from user-data */
  cvmData = (cvmInterfaceData) f_dataB;
  fwdPb = cvmData->fwdPb;
  bckPb = cvmData->bckPb;

  /* Inputs to the Matlab function */
  mx_in[0] = mxCreateScalarDouble(0.0);         /* type=0: not dependent on yS */
  mx_in[1] = mxCreateScalarDouble(t);           /* current t */
  mx_in[2] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current y */
  mx_in[3] = mxCreateDoubleMatrix(NB,1,mxREAL); /* current yB */
  mx_in[4] = mx_RHSfctB;                        /* matlab function handle */ 
  mx_in[5] = mx_data;                           /* matlab user data */
  
  /* Call matlab wrapper */
  GetData(y, mxGetPr(mx_in[2]), N);
  GetData(yB, mxGetPr(mx_in[3]), NB);

  mexCallMATLAB(3,mx_out,6,mx_in,"cvm_rhsB");

  PutData(yBd, mxGetPr(mx_out[0]), NB);

  ret = (int)*mxGetPr(mx_out[1]);

  if (!mxIsEmpty(mx_out[2])) {
    UpdateUserData(mx_out[2], cvmData);
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


int mtlb_CVodeRhsBS(realtype t, N_Vector y,  N_Vector *yS,
                    N_Vector yB, N_Vector yBd, void *f_dataB)
{
  cvmInterfaceData cvmData;
  cvmPbData fwdPb, bckPb;
  mxArray *mx_in[7], *mx_out[3];
  int is, ret;
  double *tmp;

  /* Extract global interface data from user-data */
  cvmData = (cvmInterfaceData) f_dataB;
  fwdPb = cvmData->fwdPb;
  bckPb = cvmData->bckPb;

  /* Inputs to the Matlab function */
  mx_in[0] = mxCreateScalarDouble(1.0);           /* type=1: dependent on yS */
  mx_in[1] = mxCreateScalarDouble(t);             /* current t */
  mx_in[2] = mxCreateDoubleMatrix(N,1,mxREAL);    /* current y */
  mx_in[3] = mxCreateDoubleMatrix(N*Ns,1,mxREAL); /* current yS */
  mx_in[4] = mxCreateDoubleMatrix(NB,1,mxREAL);   /* current yB */
  mx_in[5] = mx_RHSfctB;                          /* matlab function handle */ 
  mx_in[6] = mx_data;                             /* matlab user data */
  
  /* Call matlab wrapper */
  GetData(y, mxGetPr(mx_in[2]), N);

  tmp = mxGetPr(mx_in[4]);
  for (is=0; is<Ns; is++)
    GetData(yS[is], &tmp[is*N], N);

  GetData(yB, mxGetPr(mx_in[4]), NB);

  mexCallMATLAB(3,mx_out,7,mx_in,"cvm_rhs");

  PutData(yBd, mxGetPr(mx_out[0]), NB);

  ret = (int)*mxGetPr(mx_out[1]);

  if (!mxIsEmpty(mx_out[2])) {
    UpdateUserData(mx_out[2], cvmData);
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



int mtlb_CVodeQUADfctB(realtype t, N_Vector y, N_Vector yB, N_Vector yQBd, void *f_dataB)
{
  cvmInterfaceData cvmData;
  cvmPbData fwdPb, bckPb;
  mxArray *mx_in[6], *mx_out[3];
  int ret;

  /* Extract global interface data from user-data */
  cvmData = (cvmInterfaceData) f_dataB;
  fwdPb = cvmData->fwdPb;
  bckPb = cvmData->bckPb;

  /* Inputs to the Matlab function */
  mx_in[0] = mxCreateScalarDouble(0.0);         /* type=0: not dependent on yS */
  mx_in[1] = mxCreateScalarDouble(t);           /* current t */
  mx_in[2] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current y */
  mx_in[3] = mxCreateDoubleMatrix(NB,1,mxREAL); /* current yB */
  mx_in[4] = mx_QUADfctB;                       /* matlab function handle */ 
  mx_in[5] = mx_data;                           /* matlab user data */

  /* Call matlab wrapper */
  GetData(y, mxGetPr(mx_in[2]), N);
  GetData(yB, mxGetPr(mx_in[3]), NB);

  mexCallMATLAB(3,mx_out,6,mx_in,"cvm_rhsQB");

  PutData(yQBd, mxGetPr(mx_out[0]), NqB);

  ret = (int)*mxGetPr(mx_out[1]);

  if (!mxIsEmpty(mx_out[2])) {
    UpdateUserData(mx_out[2], cvmData);
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


int mtlb_CVodeQUADfctBS(realtype t, N_Vector y,  N_Vector *yS,
                        N_Vector yB, N_Vector yQBd, void *f_dataB)
{
  cvmInterfaceData cvmData;
  cvmPbData fwdPb, bckPb;
  mxArray *mx_in[7], *mx_out[3];
  int is, ret;
  double *tmp;

  /* Extract global interface data from user-data */
  cvmData = (cvmInterfaceData) f_dataB;
  fwdPb = cvmData->fwdPb;
  bckPb = cvmData->bckPb;

  /* Inputs to the Matlab function */
  mx_in[0] = mxCreateScalarDouble(1.0);           /* type=1: dependent on yS */
  mx_in[1] = mxCreateScalarDouble(t);             /* current t */
  mx_in[2] = mxCreateDoubleMatrix(N,1,mxREAL);    /* current y */
  mx_in[3] = mxCreateDoubleMatrix(N*Ns,1,mxREAL); /* current yS */
  mx_in[4] = mxCreateDoubleMatrix(NB,1,mxREAL);   /* current yB */
  mx_in[5] = mx_QUADfctB;                         /* matlab function handle */ 
  mx_in[6] = mx_data;                             /* matlab user data */

  /* Call matlab wrapper */

  GetData(y, mxGetPr(mx_in[2]), N);

  tmp = mxGetPr(mx_in[4]);
  for (is=0; is<Ns; is++)
    GetData(yS[is], &tmp[is*N], N);

  GetData(yB, mxGetPr(mx_in[4]), NB);

  mexCallMATLAB(3,mx_out,7,mx_in,"cvm_rhsQB");

  PutData(yQBd, mxGetPr(mx_out[0]), NqB);

  ret = (int)*mxGetPr(mx_out[1]);

  if (!mxIsEmpty(mx_out[2])) {
    UpdateUserData(mx_out[2], cvmData);
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




int mtlb_CVodeDenseJacB(int NeqB, realtype t,
                        N_Vector y, N_Vector yB, N_Vector fyB,
                        DlsMat JB, void *f_dataB, 
                        N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B)
{
  cvmInterfaceData cvmData;
  cvmPbData fwdPb, bckPb;
  double *JB_data;
  mxArray *mx_in[6], *mx_out[3];
  int i, ret;

  /* Extract global interface data from user-data */
  cvmData = (cvmInterfaceData) f_dataB;
  fwdPb = cvmData->fwdPb;
  bckPb = cvmData->bckPb;

  /* Inputs to the Matlab function */
  mx_in[0] = mxCreateScalarDouble(t);           /* current t */  
  mx_in[1] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current y */
  mx_in[2] = mxCreateDoubleMatrix(NB,1,mxREAL); /* current yB */
  mx_in[3] = mxCreateDoubleMatrix(NB,1,mxREAL); /* current fyB */
  mx_in[4] = mx_JACfctB;                        /* matlab function handle */
  mx_in[5] = mx_data;                           /* matlab user data */
  
  /* Call matlab wrapper */
  GetData(y, mxGetPr(mx_in[1]), N);
  GetData(yB, mxGetPr(mx_in[2]), NB);
  GetData(fyB, mxGetPr(mx_in[3]), NB);

  mexCallMATLAB(3,mx_out,6,mx_in,"cvm_djacB");

  JB_data = mxGetPr(mx_out[0]);
  for (i=0;i<NB;i++)
    memcpy(DENSE_COL(JB,i), JB_data + i*NB, NB*sizeof(double));

  ret = (int)*mxGetPr(mx_out[1]);

  if (!mxIsEmpty(mx_out[2])) {
    UpdateUserData(mx_out[2], cvmData);
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


int mtlb_CVodeBandJacB(int NeqB, int mupperB, int mlowerB, realtype t, 
                       N_Vector y, N_Vector yB, N_Vector fyB,
                       DlsMat JB, void *f_dataB, 
                       N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B)
{
  cvmInterfaceData cvmData;
  cvmPbData fwdPb, bckPb;
  double *JB_data;
  mxArray *mx_in[6], *mx_out[3];
  int ebandB, i, ret;

  /* Extract global interface data from user-data */
  cvmData = (cvmInterfaceData) f_dataB;
  fwdPb = cvmData->fwdPb;
  bckPb = cvmData->bckPb;

  /* Inputs to the Matlab function */
  mx_in[0] = mxCreateScalarDouble(t);           /* current t */
  mx_in[1] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current y */
  mx_in[2] = mxCreateDoubleMatrix(NB,1,mxREAL); /* current yB */
  mx_in[3] = mxCreateDoubleMatrix(NB,1,mxREAL); /* current fyB */
  mx_in[4] = mx_JACfctB;                        /* matlab function handle */
  mx_in[5] = mx_data;                           /* matlab user data */
  
  /* Call matlab wrapper */
  GetData(y, mxGetPr(mx_in[1]), N);
  GetData(yB, mxGetPr(mx_in[2]), NB);
  GetData(fyB, mxGetPr(mx_in[3]), NB);

  mexCallMATLAB(3,mx_out,6,mx_in,"cvm_bjacB");

  ebandB =  mupperB + mlowerB + 1;
  JB_data = mxGetPr(mx_out[0]);
  for (i=0;i<NB;i++)
    memcpy(BAND_COL(JB,i) - mupperB, JB_data + i*ebandB, ebandB*sizeof(double));
    
  ret = (int)*mxGetPr(mx_out[1]);

  if (!mxIsEmpty(mx_out[2])) {
    UpdateUserData(mx_out[2], cvmData);
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


int mtlb_CVodeSpilsJacB(N_Vector vB, N_Vector JvB, realtype t,
                        N_Vector y, N_Vector yB, N_Vector fyB,
                        void *f_dataB, N_Vector tmpB)
{
  cvmInterfaceData cvmData;
  cvmPbData fwdPb, bckPb;
  mxArray *mx_in[7], *mx_out[3];
  int ret;

  /* Extract global interface data from user-data */
  cvmData = (cvmInterfaceData) f_dataB;
  fwdPb = cvmData->fwdPb;
  bckPb = cvmData->bckPb;

  /* Inputs to the Matlab function */
  mx_in[0] = mxCreateScalarDouble(t);           /* current t */ 
  mx_in[1] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current y */
  mx_in[2] = mxCreateDoubleMatrix(NB,1,mxREAL); /* current yB */
  mx_in[3] = mxCreateDoubleMatrix(NB,1,mxREAL); /* current fyB */
  mx_in[4] = mxCreateDoubleMatrix(NB,1,mxREAL); /* vector vB */
  mx_in[5] = mx_JACfctB;                        /* matlab function handle */
  mx_in[6] = mx_data;                           /* matlab user data */
  
  /* Call matlab wrapper */
  GetData(y, mxGetPr(mx_in[1]), N);
  GetData(yB, mxGetPr(mx_in[2]), NB);
  GetData(fyB, mxGetPr(mx_in[3]), NB);
  GetData(vB, mxGetPr(mx_in[4]), NB);

  mexCallMATLAB(3,mx_out,7,mx_in,"cvm_jtvB");

  PutData(JvB, mxGetPr(mx_out[0]), NB);
  ret = (int)*mxGetPr(mx_out[1]);
  
  if (!mxIsEmpty(mx_out[2])) {
    UpdateUserData(mx_out[2], cvmData);
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


int mtlb_CVodeSpilsPsetB(realtype t, N_Vector y,
                         N_Vector yB, N_Vector fyB,
                         booleantype jokB,
                         booleantype *jcurPtrB, realtype gammaB,
                         void *f_dataB,
                         N_Vector tmp1B, N_Vector tmp2B,
                         N_Vector tmp3B)
{
  cvmInterfaceData cvmData;
  cvmPbData fwdPb, bckPb;
  mxArray *mx_in[8], *mx_out[3];
  int ret;

  /* Extract global interface data from user-data */
  cvmData = (cvmInterfaceData) f_dataB;
  fwdPb = cvmData->fwdPb;
  bckPb = cvmData->bckPb;

  /* Inputs to the Matlab function */
  mx_in[0] = mxCreateScalarDouble(t);           /* current t */
  mx_in[1] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current y */
  mx_in[2] = mxCreateDoubleMatrix(NB,1,mxREAL); /* current yB */
  mx_in[3] = mxCreateDoubleMatrix(NB,1,mxREAL); /* current fyB */
  mx_in[4] = mxCreateLogicalScalar(jokB);       /* jokB flag */
  mx_in[5] = mxCreateScalarDouble(gammaB);      /* gammaB value */
  mx_in[6] = mx_PSETfctB;                       /* matlab function handle */
  mx_in[7] = mx_data;                           /* matlab user data */
  
  /* Call matlab wrapper */
  GetData(y, mxGetPr(mx_in[1]), N);
  GetData(yB, mxGetPr(mx_in[2]), NB);
  GetData(fyB, mxGetPr(mx_in[3]), NB);

  mexCallMATLAB(3,mx_out,8,mx_in,"cvm_psetB");

  *jcurPtrB = mxIsLogicalScalarTrue(mx_out[0]);
  ret = (int)*mxGetPr(mx_out[1]);

  if (!mxIsEmpty(mx_out[2])) {
    UpdateUserData(mx_out[2], cvmData);
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


int mtlb_CVodeSpilsPsolB(realtype t, N_Vector y,
                         N_Vector yB, N_Vector fyB,
                         N_Vector rB, N_Vector zB,
                         realtype gammaB, realtype deltaB,
                         int lrB, void *f_dataB, N_Vector tmpB)
{
  cvmInterfaceData cvmData;
  cvmPbData fwdPb, bckPb;
  mxArray *mx_in[7], *mx_out[3];
  int ret;

  /* Extract global interface data from user-data */
  cvmData = (cvmInterfaceData) f_dataB;
  fwdPb = cvmData->fwdPb;
  bckPb = cvmData->bckPb;

  /* Inputs to the Matlab function */
  mx_in[0] = mxCreateScalarDouble(t);           /* current t */   
  mx_in[1] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current y */
  mx_in[2] = mxCreateDoubleMatrix(NB,1,mxREAL); /* current yB */
  mx_in[3] = mxCreateDoubleMatrix(NB,1,mxREAL); /* current fyB */
  mx_in[4] = mxCreateDoubleMatrix(NB,1,mxREAL); /* right hand side rB */
  mx_in[5] = mx_PSOLfctB;                       /* matlab function handle */
  mx_in[6] = mx_data;                           /* matlab user data */
  
  /* Call matlab wrapper */
  GetData(y, mxGetPr(mx_in[1]), N);
  GetData(yB, mxGetPr(mx_in[2]), NB);
  GetData(fyB, mxGetPr(mx_in[3]), NB);
  GetData(rB, mxGetPr(mx_in[4]), NB);

  mexCallMATLAB(3,mx_out,7,mx_in,"cvm_psolB");

  PutData(zB, mxGetPr(mx_out[0]), NB);
  ret = (int)*mxGetPr(mx_out[1]);

  if (!mxIsEmpty(mx_out[2])) {
    UpdateUserData(mx_out[2], cvmData);
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

int mtlb_CVodeBBDglocB(int NlocalB, realtype t, N_Vector y,
                       N_Vector yB, N_Vector gB, void *f_dataB)
{
  cvmInterfaceData cvmData;
  cvmPbData fwdPb, bckPb;
  mxArray *mx_in[5], *mx_out[3];
  int ret;

  /* Extract global interface data from user-data */
  cvmData = (cvmInterfaceData) f_dataB;
  fwdPb = cvmData->fwdPb;
  bckPb = cvmData->bckPb;

  /* Inputs to the Matlab function */
  mx_in[0] = mxCreateScalarDouble(t);           /* current t */
  mx_in[1] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current y */
  mx_in[2] = mxCreateDoubleMatrix(NB,1,mxREAL); /* current yB */
  mx_in[3] = mx_GLOCfctB;                       /* matlab function handle */
  mx_in[4] = mx_data;                           /* matlab user data */
  
  /* Call matlab wrapper */
  GetData(y, mxGetPr(mx_in[1]), N);
  GetData(yB, mxGetPr(mx_in[2]), NB);

  mexCallMATLAB(3,mx_out,5,mx_in,"cvm_glocB");

  PutData(gB, mxGetPr(mx_out[0]), NB);

  ret = (int)*mxGetPr(mx_out[1]);

  if (!mxIsEmpty(mx_out[2])) {
    UpdateUserData(mx_out[2], cvmData);
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

int mtlb_CVodeBBDgcomB(int NlocalB, realtype t, N_Vector y, 
                       N_Vector yB, void *f_dataB)
{
  cvmInterfaceData cvmData;
  cvmPbData fwdPb, bckPb;
  mxArray *mx_in[5], *mx_out[2];
  int ret;

  /* Extract global interface data from user-data */
  cvmData = (cvmInterfaceData) f_dataB;
  fwdPb = cvmData->fwdPb;
  bckPb = cvmData->bckPb;

  /* Inputs to the Matlab function */
  mx_in[0] = mxCreateScalarDouble(t);           /* current t */
  mx_in[1] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current y */
  mx_in[2] = mxCreateDoubleMatrix(NB,1,mxREAL); /* current yB */
  mx_in[3] = mx_GCOMfctB;                       /* matlab function handle */
  mx_in[4] = mx_data;                           /* matlab user data */
  
  /* Call matlab wrapper */
  GetData(y, mxGetPr(mx_in[1]), N);
  GetData(yB, mxGetPr(mx_in[2]), NB);

  mexCallMATLAB(2,mx_out,5,mx_in,"cvm_gcomB");

  ret = (int)*mxGetPr(mx_out[0]);

  if (!mxIsEmpty(mx_out[1])) {
    UpdateUserData(mx_out[1], cvmData);
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
 * ---------------------------------------------------------------------------------
 * Wrapper around matlab monitor function
 * ---------------------------------------------------------------------------------
 */

void mtlb_CVodeMonitor(int call, double t, N_Vector y, N_Vector yQ, N_Vector *yS, cvmInterfaceData cvmData)
{
  cvmPbData fwdPb;
  mxArray *mx_in[8], *mx_out[1];
  double *tmp;
  int is;

  fwdPb = cvmData->fwdPb;

  mx_in[0] = mxCreateScalarDouble(call);            /* call type (0:first, 1:interm. 2:last) */
  mx_in[1] = mxCreateScalarDouble(t);               /* current time */
  mx_in[2] = mxCreateDoubleMatrix(N,1,mxREAL);      /* current solution */
  if (quadr)
    mx_in[3] = mxCreateDoubleMatrix(Nq,1,mxREAL);   /* current quadratures */
  else
    mx_in[3] = mxCreateDoubleMatrix(0,0,mxREAL);
  mx_in[4] = mxCreateScalarDouble(Ns);              /* number of sensitivities */
  if (fsa)
    mx_in[5] = mxCreateDoubleMatrix(N*Ns,1,mxREAL); /* current sensitivities */
  else
    mx_in[5] = mxCreateDoubleMatrix(0,0,mxREAL);
  mx_in[6] = mx_MONfct;                             /* Matlab monitor function */
  mx_in[7] = mx_MONdata;                            /* data for monitor function */

  if (call == 1) {

    GetData(y, mxGetPr(mx_in[2]), N);

    if (quadr)
      GetData(yQ, mxGetPr(mx_in[3]), Nq);

    if (fsa) {
      tmp = mxGetPr(mx_in[5]);
      for (is=0; is<Ns; is++)
        GetData(yS[is], &tmp[is*N], N);
    }

  }

  mexCallMATLAB(1,mx_out,8,mx_in,"cvm_monitor");

  if (!mxIsEmpty(mx_out[0])) {
    UpdateMonitorData(mx_out[0], cvmData);
  }

  mxDestroyArray(mx_in[0]);
  mxDestroyArray(mx_in[1]);
  mxDestroyArray(mx_in[2]);
  mxDestroyArray(mx_in[3]);
  mxDestroyArray(mx_in[4]);
  mxDestroyArray(mx_in[5]);
  mxDestroyArray(mx_out[0]);
}

void mtlb_CVodeMonitorB(int call, int idxB, double tB, N_Vector yB, N_Vector yQB, cvmInterfaceData cvmData)
{
  cvmPbData bckPb;
  mxArray *mx_in[7], *mx_out[1];

  bckPb = cvmData->bckPb;

  mx_in[0] = mxCreateScalarDouble(call);            /* 0: first, 1: interm. 2: last */
  mx_in[1] = mxCreateScalarDouble(idxB);            /* index of current problem */
  mx_in[2] = mxCreateScalarDouble(tB);              /* current time */
  mx_in[3] = mxCreateDoubleMatrix(NB,1,mxREAL);     /* current solution */
  if (quadrB)
    mx_in[4] = mxCreateDoubleMatrix(NqB,1,mxREAL);  /* current quadratures */
  else
    mx_in[4] = mxCreateDoubleMatrix(0,0,mxREAL);
  mx_in[5] = mx_MONfctB;
  mx_in[6] = mx_MONdataB;

  if (call == 1) {
    
    GetData(yB, mxGetPr(mx_in[3]), NB);

    if (quadrB)
      GetData(yQB, mxGetPr(mx_in[4]), NqB);
  }

  mexCallMATLAB(1,mx_out,7,mx_in,"cvm_monitorB");

  if (!mxIsEmpty(mx_out[0])) {
    UpdateMonitorDataB(mx_out[0], cvmData);
  }

  mxDestroyArray(mx_in[0]);
  mxDestroyArray(mx_in[1]);
  mxDestroyArray(mx_in[2]);
  mxDestroyArray(mx_in[3]);
  mxDestroyArray(mx_in[4]);
  mxDestroyArray(mx_out[0]);
}


/*
 * ---------------------------------------------------------------------------------
 * Private functions to update the user data structures
 * ---------------------------------------------------------------------------------
 */

static void UpdateUserData(mxArray *mx_newdata, cvmInterfaceData cvmData)
{
  mexUnlock();
  mxDestroyArray(mx_data);
  mx_data = mxDuplicateArray(mx_newdata);
  mexMakeArrayPersistent(mx_data);
  mexLock();
}

static void UpdateMonitorData(mxArray *mx_newdata, cvmInterfaceData cvmData)
{
  cvmPbData fwdPb;

  fwdPb = cvmData->fwdPb;

  mexUnlock();
  mxDestroyArray(mx_MONdata);
  mx_MONdata = mxDuplicateArray(mx_newdata);
  mexMakeArrayPersistent(mx_MONdata);
  mexLock();
}

static void UpdateMonitorDataB(mxArray *mx_newdata, cvmInterfaceData cvmData)
{
  cvmPbData bckPb;

  bckPb = cvmData->bckPb;

  mexUnlock();
  mxDestroyArray(mx_MONdataB);
  mx_MONdataB = mxDuplicateArray(mx_newdata);
  mexMakeArrayPersistent(mx_MONdataB);
  mexLock();
}

