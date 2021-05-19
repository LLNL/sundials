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
 * CVODES wrapper functions.
 * -----------------------------------------------------------------
 */

#include "cvm.h"
#include "nvm.h"

static void UpdateUserData(mxArray *new_mtlb_data, cvmPbData pb);
static void UpdateMonitorData(mxArray *new_mtlb_data, cvmPbData pb);

/*
 * ---------------------------------------------------------------------------------
 * Redability replacements
 * ---------------------------------------------------------------------------------
 */

#define fsa         (fwdPb->Fsa)
#define quadr       (fwdPb->Quadr)
#define N           (fwdPb->n) 
#define Nq          (fwdPb->nq) 
#define Ng          (fwdPb->ng) 
#define Ns          (fwdPb->ns) 

#define quadrB      (bckPb->Quadr)
#define NB          (bckPb->n) 
#define NqB         (bckPb->nq) 

/*
 * ---------------------------------------------------------------------------------
 * FORWARD PROBLEMS
 * ---------------------------------------------------------------------------------
 */

int mxW_CVodeRhs(realtype t, N_Vector y, N_Vector yd, void *user_data)
{
  cvmPbData fwdPb;
  mxArray *mx_in[4], *mx_out[3];
  int ret;

  /* Extract global interface data from user-data */
  fwdPb = (cvmPbData) user_data;

  /* Inputs to the Matlab function */
  mx_in[0] = mxCreateDoubleScalar(t);          /* current t */
  mx_in[1] = mxCreateDoubleMatrix(N,1,mxREAL); /* current y */
  mx_in[2] = fwdPb->RHSfct;                    /* matlab function handle */ 
  mx_in[3] = fwdPb->mtlb_data;                 /* matlab user data */

  /* Call matlab wrapper */
  GetData(y, mxGetPr(mx_in[1]), N);

  mexCallMATLAB(3,mx_out,4,mx_in,"cvm_rhs");

  PutData(yd, mxGetPr(mx_out[0]), N);
  ret = (int)*mxGetPr(mx_out[1]);

  if (!mxIsEmpty(mx_out[2])) {
    UpdateUserData(mx_out[2], fwdPb);
  }

  /* Free temporary space */
  mxDestroyArray(mx_in[0]);
  mxDestroyArray(mx_in[1]);
  mxDestroyArray(mx_out[0]);
  mxDestroyArray(mx_out[1]);
  mxDestroyArray(mx_out[2]);

  return(ret);
}

int mxW_CVodeQUADfct(realtype t, N_Vector y, N_Vector yQd, void *user_data)
{
  cvmPbData fwdPb;
  mxArray *mx_in[4], *mx_out[3];
  int ret;

  /* Extract global interface data from user-data */
  fwdPb = (cvmPbData) user_data;

  /* Inputs to the Matlab function */
  mx_in[0] = mxCreateDoubleScalar(t);          /* current t */
  mx_in[1] = mxCreateDoubleMatrix(N,1,mxREAL); /* current y */
  mx_in[2] = fwdPb->QUADfct;                   /* matlab function handle */ 
  mx_in[3] = fwdPb->mtlb_data;                 /* matlab user data */

  /* Call matlab wrapper */
  GetData(y, mxGetPr(mx_in[1]), N);

  mexCallMATLAB(3,mx_out,4,mx_in,"cvm_rhsQ");

  PutData(yQd, mxGetPr(mx_out[0]), Nq);
  ret = (int)*mxGetPr(mx_out[1]);

  if (!mxIsEmpty(mx_out[2])) {
    UpdateUserData(mx_out[2], fwdPb);
  }

  /* Free temporary space */
  mxDestroyArray(mx_in[0]);
  mxDestroyArray(mx_in[1]);
  mxDestroyArray(mx_out[0]);
  mxDestroyArray(mx_out[1]);
  mxDestroyArray(mx_out[2]);

  return(ret);
}

int mxW_CVodeGfct(realtype t, N_Vector y, double *g, void *user_data)
{
  cvmPbData fwdPb;
  double *gdata;
  int i, ret;
  mxArray *mx_in[4], *mx_out[3];
  
  /* Extract global interface data from user-data */
  fwdPb = (cvmPbData) user_data;

  /* Inputs to the Matlab function */
  mx_in[0] = mxCreateDoubleScalar(t);          /* current t */
  mx_in[1] = mxCreateDoubleMatrix(N,1,mxREAL); /* current y */
  mx_in[2] = fwdPb->Gfct;                      /* matlab function handle */
  mx_in[3] = fwdPb->mtlb_data;                 /* matlab user data */
  
  /* Call matlab wrapper */
  GetData(y, mxGetPr(mx_in[1]), N);

  mexCallMATLAB(3,mx_out,4,mx_in,"cvm_root");

  gdata = mxGetPr(mx_out[0]);
  for (i=0;i<Ng;i++) g[i] = gdata[i];

  ret = (int)*mxGetPr(mx_out[1]);

  if (!mxIsEmpty(mx_out[2])) {
    UpdateUserData(mx_out[2], fwdPb);
  }

  /* Free temporary space */
  mxDestroyArray(mx_in[0]);
  mxDestroyArray(mx_in[1]);
  mxDestroyArray(mx_out[0]);
  mxDestroyArray(mx_out[1]);
  mxDestroyArray(mx_out[2]);

  return(ret);
}


int mxW_CVodeDenseJac(sunindextype Neq, realtype t,
                      N_Vector y, N_Vector fy, 
                      DlsMat J, void *user_data,
                      N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  cvmPbData fwdPb;
  double *J_data;
  sunindextype i;
  int ret;
  mxArray *mx_in[5], *mx_out[3];

  /* Extract global interface data from user-data */
  fwdPb = (cvmPbData) user_data;

  /* Inputs to the Matlab function */
  mx_in[0] = mxCreateDoubleScalar(t);           /* current t */  
  mx_in[1] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current y */
  mx_in[2] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current fy */
  mx_in[3] = fwdPb->JACfct;                     /* matlab function handle */
  mx_in[4] = fwdPb->mtlb_data;                  /* matlab user data */
  
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
    UpdateUserData(mx_out[2], fwdPb);
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

int mxW_CVodeBandJac(sunindextype Neq, sunindextype mupper, sunindextype mlower, realtype t,
                     N_Vector y, N_Vector fy, 
                     DlsMat J, void *user_data,
                     N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  cvmPbData fwdPb;
  double *J_data;
  sunindextype eband, i;
  int ret;
  mxArray *mx_in[5], *mx_out[3];

  /* Extract global interface data from user-data */
  fwdPb = (cvmPbData) user_data;

  /* Inputs to the Matlab function */
  mx_in[0] = mxCreateDoubleScalar(t);           /* current t */
  mx_in[1] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current y */
  mx_in[2] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current fy */
  mx_in[3] = fwdPb->JACfct;                     /* matlab function handle */
  mx_in[4] = fwdPb->mtlb_data;                  /* matlab user data */
  
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
    UpdateUserData(mx_out[2], fwdPb);
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


int mxW_CVodeSpilsJac(N_Vector v, N_Vector Jv, realtype t,
                      N_Vector y, N_Vector fy,
                      void *user_data, N_Vector tmp)
{
  cvmPbData fwdPb;
  mxArray *mx_in[6], *mx_out[3];
  int ret;

  /* Extract global interface data from user-data */
  fwdPb = (cvmPbData) user_data;

  /* Inputs to the Matlab function */
  mx_in[0] = mxCreateDoubleScalar(t);           /* current t */ 
  mx_in[1] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current y */
  mx_in[2] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current fy */
  mx_in[3] = mxCreateDoubleMatrix(N,1,mxREAL);  /* vector v */
  mx_in[4] = fwdPb->JACfct;                     /* matlab function handle */
  mx_in[5] = fwdPb->mtlb_data;                  /* matlab user data */
  
  /* Call matlab wrapper */
  GetData(y, mxGetPr(mx_in[1]), N);
  GetData(fy, mxGetPr(mx_in[2]), N);
  GetData(v, mxGetPr(mx_in[3]), N);

  mexCallMATLAB(3,mx_out,6,mx_in,"cvm_jtv");

  PutData(Jv, mxGetPr(mx_out[0]), N);
  ret = (int)*mxGetPr(mx_out[1]);

  if (!mxIsEmpty(mx_out[2])) {
    UpdateUserData(mx_out[2], fwdPb);
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


int mxW_CVodeSpilsPset(realtype t, N_Vector y, N_Vector fy,
                       booleantype jok, booleantype *jcurPtr,
                       realtype gamma, void *user_data,
                       N_Vector tmp1, N_Vector tmp2,
                       N_Vector tmp3)
{
  cvmPbData fwdPb;
  mxArray *mx_in[7], *mx_out[3];
  int ret;

  /* Extract global interface data from user-data */
  fwdPb = (cvmPbData) user_data;

  /* Inputs to the Matlab function */
  mx_in[0] = mxCreateDoubleScalar(t);           /* current t */
  mx_in[1] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current y */
  mx_in[2] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current fy */
  mx_in[3] = mxCreateLogicalScalar(jok);        /* jok flag */
  mx_in[4] = mxCreateDoubleScalar(gamma);       /* gamma value */
  mx_in[5] = fwdPb->PSETfct;                    /* matlab function handle */
  mx_in[6] = fwdPb->mtlb_data;                  /* matlab user data */
  
  /* Call matlab wrapper */
  GetData(y, mxGetPr(mx_in[1]), N);
  GetData(fy, mxGetPr(mx_in[2]), N);

  mexCallMATLAB(3,mx_out,7,mx_in,"cvm_pset");

  *jcurPtr = mxIsLogicalScalarTrue(mx_out[0]);
  ret = (int)*mxGetPr(mx_out[1]);

  if (!mxIsEmpty(mx_out[2])) {
    UpdateUserData(mx_out[2], fwdPb);
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


int mxW_CVodeSpilsPsol(realtype t, N_Vector y, N_Vector fy,
                       N_Vector r, N_Vector z,
                       realtype gamma, realtype delta,
                       int lr, void *user_data, N_Vector tmp)
{
  cvmPbData fwdPb;
  mxArray *mx_in[6], *mx_out[3];
  int ret;

  /* Extract global interface data from user-data */
  fwdPb = (cvmPbData) user_data;

  /* Inputs to the Matlab function */
  mx_in[0] = mxCreateDoubleScalar(t);          /* current t */   
  mx_in[1] = mxCreateDoubleMatrix(N,1,mxREAL); /* current y */
  mx_in[2] = mxCreateDoubleMatrix(N,1,mxREAL); /* current fy */
  mx_in[3] = mxCreateDoubleMatrix(N,1,mxREAL); /* right hand side r */
  mx_in[4] = fwdPb->PSOLfct;                   /* matlab function handle */
  mx_in[5] = fwdPb->mtlb_data;                 /* matlab user data */
  
  /* Call matlab wrapper */
  GetData(y, mxGetPr(mx_in[1]), N);
  GetData(fy, mxGetPr(mx_in[2]), N);
  GetData(r, mxGetPr(mx_in[3]), N);

  mexCallMATLAB(3,mx_out,6,mx_in,"cvm_psol");

  PutData(z, mxGetPr(mx_out[0]), N);
  ret = (int)*mxGetPr(mx_out[1]);

  if (!mxIsEmpty(mx_out[2])) {
    UpdateUserData(mx_out[2], fwdPb);
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

int mxW_CVodeBBDgloc(sunindextype Nlocal, realtype t, N_Vector y,
                     N_Vector g, void *user_data)
{
  cvmPbData fwdPb;
  mxArray *mx_in[4], *mx_out[3];
  int ret;

  /* Extract global interface data from user-data */
  fwdPb = (cvmPbData) user_data;

  /* Inputs to the Matlab function */
  mx_in[0] = mxCreateDoubleScalar(t);           /* current t */
  mx_in[1] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current y */
  mx_in[2] = fwdPb->GLOCfct;                    /* matlab function handle */
  mx_in[3] = fwdPb->mtlb_data;                  /* matlab user data */
  
  /* Call matlab wrapper */
  GetData(y, mxGetPr(mx_in[1]), N);

  mexCallMATLAB(3,mx_out,4,mx_in,"cvm_gloc");

  PutData(g, mxGetPr(mx_out[0]), N);
  ret = (int)*mxGetPr(mx_out[1]);

  if (!mxIsEmpty(mx_out[2])) {
    UpdateUserData(mx_out[2], fwdPb);
  }

  /* Free temporary space */
  mxDestroyArray(mx_in[0]);
  mxDestroyArray(mx_in[1]);
  mxDestroyArray(mx_out[0]);
  mxDestroyArray(mx_out[1]);
  mxDestroyArray(mx_out[2]);

  return(ret);
}

int mxW_CVodeBBDgcom(sunindextype Nlocal, realtype t, N_Vector y, void *user_data)
{
  cvmPbData fwdPb;
  mxArray *mx_in[4], *mx_out[2];
  int ret;

  /* Extract global interface data from user-data */
  fwdPb = (cvmPbData) user_data;

  /* Inputs to the Matlab function */
  mx_in[0] = mxCreateDoubleScalar(t);           /* current t */
  mx_in[1] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current y */
  mx_in[2] = fwdPb->GCOMfct;                    /* matlab function handle */
  mx_in[3] = fwdPb->mtlb_data;                  /* matlab user data */
  
  /* Call matlab wrapper */
  GetData(y, mxGetPr(mx_in[1]), N);

  mexCallMATLAB(2,mx_out,4,mx_in,"cvm_gcom");

  ret = (int)*mxGetPr(mx_out[0]);

  if (!mxIsEmpty(mx_out[1])) {
    UpdateUserData(mx_out[1], fwdPb);
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

int mxW_CVodeSensRhs(int Nsens, realtype t,
                     N_Vector y, N_Vector yd,
                     N_Vector *yS, N_Vector *ySd,
                     void *user_data,
                     N_Vector tmp1, N_Vector tmp2)
{
  cvmPbData fwdPb;
  mxArray *mx_in[7], *mx_out[3];
  int is, ret;
  double *tmp;

  /* Extract global interface data from user-data */
  fwdPb = (cvmPbData) user_data;

  /* Inputs to the Matlab function */
  mx_in[0] = mxCreateDoubleScalar(t);             /* current t */
  mx_in[1] = mxCreateDoubleMatrix(N,1,mxREAL);    /* current y */
  mx_in[2] = mxCreateDoubleMatrix(N,1,mxREAL);    /* current yd */
  mx_in[3] = mxCreateDoubleScalar(Ns);            /* number of sensitivities */
  mx_in[4] = mxCreateDoubleMatrix(N*Ns,1,mxREAL); /* current yS */
  mx_in[5] = fwdPb->SRHSfct;                      /* matlab function handle */      
  mx_in[6] = fwdPb->mtlb_data;                    /* matlab user data */
  
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
    UpdateUserData(mx_out[2], fwdPb);
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

int mxW_CVodeRhsB(realtype t, N_Vector y, N_Vector yB, N_Vector yBd, void *user_dataB)
{
  cvmPbData fwdPb, bckPb;
  mxArray *mx_in[6], *mx_out[3];
  int ret;

  /* Extract global interface data from user-data */
  bckPb = (cvmPbData) user_dataB;
  fwdPb = bckPb->fwd;

  /* Inputs to the Matlab function */
  mx_in[0] = mxCreateDoubleScalar(0.0);         /* type=0: not dependent on yS */
  mx_in[1] = mxCreateDoubleScalar(t);           /* current t */
  mx_in[2] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current y */
  mx_in[3] = mxCreateDoubleMatrix(NB,1,mxREAL); /* current yB */
  mx_in[4] = bckPb->RHSfct;                     /* matlab function handle */ 
  mx_in[5] = bckPb->mtlb_data;                  /* matlab user data */
  
  /* Call matlab wrapper */
  GetData(y, mxGetPr(mx_in[2]), N);
  GetData(yB, mxGetPr(mx_in[3]), NB);

  mexCallMATLAB(3,mx_out,6,mx_in,"cvm_rhsB");

  PutData(yBd, mxGetPr(mx_out[0]), NB);

  ret = (int)*mxGetPr(mx_out[1]);

  if (!mxIsEmpty(mx_out[2])) {
    UpdateUserData(mx_out[2], bckPb);
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


int mxW_CVodeRhsBS(realtype t, N_Vector y,  N_Vector *yS,
                   N_Vector yB, N_Vector yBd, void *user_dataB)
{
  cvmPbData fwdPb, bckPb;
  mxArray *mx_in[8], *mx_out[3];
  int is, ret;
  double *tmp;

  /* Extract global interface data from user-data */
  bckPb = (cvmPbData) user_dataB;
  fwdPb = bckPb->fwd;

  /* Inputs to the Matlab function */
  mx_in[0] = mxCreateDoubleScalar(1.0);           /* type=1: dependent on yS */
  mx_in[1] = mxCreateDoubleScalar(t);             /* current t */
  mx_in[2] = mxCreateDoubleMatrix(N,1,mxREAL);    /* current y */
  mx_in[3] = mxCreateDoubleScalar(Ns);            /* number of sensitivities */
  mx_in[4] = mxCreateDoubleMatrix(N*Ns,1,mxREAL); /* current yS */
  mx_in[5] = mxCreateDoubleMatrix(NB,1,mxREAL);   /* current yB */
  mx_in[6] = bckPb->RHSfct;                       /* matlab function handle */ 
  mx_in[7] = bckPb->mtlb_data;                    /* matlab user data */
  
  /* Call matlab wrapper */
  GetData(y, mxGetPr(mx_in[2]), N);

  tmp = mxGetPr(mx_in[4]);
  for (is=0; is<Ns; is++)
    GetData(yS[is], &tmp[is*N], N);

  GetData(yB, mxGetPr(mx_in[5]), NB);

  mexCallMATLAB(3,mx_out,8,mx_in,"cvm_rhsB");

  PutData(yBd, mxGetPr(mx_out[0]), NB);

  ret = (int)*mxGetPr(mx_out[1]);

  if (!mxIsEmpty(mx_out[2])) {
    UpdateUserData(mx_out[2], bckPb);
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



int mxW_CVodeQUADfctB(realtype t, N_Vector y, N_Vector yB, N_Vector yQBd, void *user_dataB)
{
  cvmPbData fwdPb, bckPb;
  mxArray *mx_in[6], *mx_out[3];
  int ret;

  /* Extract global interface data from user-data */
  bckPb = (cvmPbData) user_dataB;
  fwdPb = bckPb->fwd;

  /* Inputs to the Matlab function */
  mx_in[0] = mxCreateDoubleScalar(0.0);         /* type=0: not dependent on yS */
  mx_in[1] = mxCreateDoubleScalar(t);           /* current t */
  mx_in[2] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current y */
  mx_in[3] = mxCreateDoubleMatrix(NB,1,mxREAL); /* current yB */
  mx_in[4] = bckPb->QUADfct;                    /* matlab function handle */ 
  mx_in[5] = bckPb->mtlb_data;                  /* matlab user data */

  /* Call matlab wrapper */
  GetData(y, mxGetPr(mx_in[2]), N);
  GetData(yB, mxGetPr(mx_in[3]), NB);

  mexCallMATLAB(3,mx_out,6,mx_in,"cvm_rhsQB");

  PutData(yQBd, mxGetPr(mx_out[0]), NqB);

  ret = (int)*mxGetPr(mx_out[1]);

  if (!mxIsEmpty(mx_out[2])) {
    UpdateUserData(mx_out[2], bckPb);
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


int mxW_CVodeQUADfctBS(realtype t, N_Vector y,  N_Vector *yS,
                       N_Vector yB, N_Vector yQBd, void *user_dataB)
{
  cvmPbData fwdPb, bckPb;
  mxArray *mx_in[8], *mx_out[3];
  int is, ret;
  double *tmp;

  /* Extract global interface data from user-data */
  bckPb = (cvmPbData) user_dataB;
  fwdPb = bckPb->fwd;

  /* Inputs to the Matlab function */
  mx_in[0] = mxCreateDoubleScalar(1.0);           /* type=1: dependent on yS */
  mx_in[1] = mxCreateDoubleScalar(t);             /* current t */
  mx_in[2] = mxCreateDoubleMatrix(N,1,mxREAL);    /* current y */
  mx_in[3] = mxCreateDoubleScalar(Ns);            /* number of sensitivities */
  mx_in[4] = mxCreateDoubleMatrix(N*Ns,1,mxREAL); /* current yS */
  mx_in[5] = mxCreateDoubleMatrix(NB,1,mxREAL);   /* current yB */
  mx_in[6] = bckPb->QUADfct;                      /* matlab function handle */ 
  mx_in[7] = bckPb->mtlb_data;                    /* matlab user data */

  /* Call matlab wrapper */

  GetData(y, mxGetPr(mx_in[2]), N);

  tmp = mxGetPr(mx_in[4]);
  for (is=0; is<Ns; is++)
    GetData(yS[is], &tmp[is*N], N);

  GetData(yB, mxGetPr(mx_in[5]), NB);

  mexCallMATLAB(3,mx_out,8,mx_in,"cvm_rhsQB");

  PutData(yQBd, mxGetPr(mx_out[0]), NqB);

  ret = (int)*mxGetPr(mx_out[1]);

  if (!mxIsEmpty(mx_out[2])) {
    UpdateUserData(mx_out[2], bckPb);
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




int mxW_CVodeDenseJacB(sunindextype NeqB, realtype t,
                       N_Vector y, N_Vector yB, N_Vector fyB,
                       DlsMat JB, void *user_dataB, 
                       N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B)
{
  cvmPbData fwdPb, bckPb;
  double *JB_data;
  mxArray *mx_in[6], *mx_out[3];
  int i, ret;

  /* Extract global interface data from user-data */
  bckPb = (cvmPbData) user_dataB;
  fwdPb = bckPb->fwd;

  /* Inputs to the Matlab function */
  mx_in[0] = mxCreateDoubleScalar(t);           /* current t */  
  mx_in[1] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current y */
  mx_in[2] = mxCreateDoubleMatrix(NB,1,mxREAL); /* current yB */
  mx_in[3] = mxCreateDoubleMatrix(NB,1,mxREAL); /* current fyB */
  mx_in[4] = bckPb->JACfct;                     /* matlab function handle */
  mx_in[5] = bckPb->mtlb_data;                  /* matlab user data */
  
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
    UpdateUserData(mx_out[2], bckPb);
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


int mxW_CVodeBandJacB(sunindextype NeqB, sunindextype mupperB, sunindextype mlowerB, realtype t, 
                       N_Vector y, N_Vector yB, N_Vector fyB,
                       DlsMat JB, void *user_dataB, 
                       N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B)
{
  cvmPbData fwdPb, bckPb;
  double *JB_data;
  mxArray *mx_in[6], *mx_out[3];
  sunindextype ebandB, i;
  int ret;

  /* Extract global interface data from user-data */
  bckPb = (cvmPbData) user_dataB;
  fwdPb = bckPb->fwd;

  /* Inputs to the Matlab function */
  mx_in[0] = mxCreateDoubleScalar(t);           /* current t */
  mx_in[1] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current y */
  mx_in[2] = mxCreateDoubleMatrix(NB,1,mxREAL); /* current yB */
  mx_in[3] = mxCreateDoubleMatrix(NB,1,mxREAL); /* current fyB */
  mx_in[4] = bckPb->JACfct;                     /* matlab function handle */
  mx_in[5] = bckPb->mtlb_data;                  /* matlab user data */
  
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
    UpdateUserData(mx_out[2], bckPb);
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


int mxW_CVodeSpilsJacB(N_Vector vB, N_Vector JvB, realtype t,
                        N_Vector y, N_Vector yB, N_Vector fyB,
                        void *user_dataB, N_Vector tmpB)
{
  cvmPbData fwdPb, bckPb;
  mxArray *mx_in[7], *mx_out[3];
  int ret;

  /* Extract global interface data from user-data */
  bckPb = (cvmPbData) user_dataB;
  fwdPb = bckPb->fwd;

  /* Inputs to the Matlab function */
  mx_in[0] = mxCreateDoubleScalar(t);           /* current t */ 
  mx_in[1] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current y */
  mx_in[2] = mxCreateDoubleMatrix(NB,1,mxREAL); /* current yB */
  mx_in[3] = mxCreateDoubleMatrix(NB,1,mxREAL); /* current fyB */
  mx_in[4] = mxCreateDoubleMatrix(NB,1,mxREAL); /* vector vB */
  mx_in[5] = bckPb->JACfct;                     /* matlab function handle */
  mx_in[6] = bckPb->mtlb_data;                  /* matlab user data */
  
  /* Call matlab wrapper */
  GetData(y, mxGetPr(mx_in[1]), N);
  GetData(yB, mxGetPr(mx_in[2]), NB);
  GetData(fyB, mxGetPr(mx_in[3]), NB);
  GetData(vB, mxGetPr(mx_in[4]), NB);

  mexCallMATLAB(3,mx_out,7,mx_in,"cvm_jtvB");

  PutData(JvB, mxGetPr(mx_out[0]), NB);
  ret = (int)*mxGetPr(mx_out[1]);
  
  if (!mxIsEmpty(mx_out[2])) {
    UpdateUserData(mx_out[2], bckPb);
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


int mxW_CVodeSpilsPsetB(realtype t, N_Vector y,
                         N_Vector yB, N_Vector fyB,
                         booleantype jokB,
                         booleantype *jcurPtrB, realtype gammaB,
                         void *user_dataB,
                         N_Vector tmp1B, N_Vector tmp2B,
                         N_Vector tmp3B)
{
  cvmPbData fwdPb, bckPb;
  mxArray *mx_in[8], *mx_out[3];
  int ret;

  /* Extract global interface data from user-data */
  bckPb = (cvmPbData) user_dataB;
  fwdPb = bckPb->fwd;

  /* Inputs to the Matlab function */
  mx_in[0] = mxCreateDoubleScalar(t);           /* current t */
  mx_in[1] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current y */
  mx_in[2] = mxCreateDoubleMatrix(NB,1,mxREAL); /* current yB */
  mx_in[3] = mxCreateDoubleMatrix(NB,1,mxREAL); /* current fyB */
  mx_in[4] = mxCreateLogicalScalar(jokB);       /* jokB flag */
  mx_in[5] = mxCreateDoubleScalar(gammaB);      /* gammaB value */
  mx_in[6] = bckPb->PSETfct;                    /* matlab function handle */
  mx_in[7] = bckPb->mtlb_data;                  /* matlab user data */
  
  /* Call matlab wrapper */
  GetData(y, mxGetPr(mx_in[1]), N);
  GetData(yB, mxGetPr(mx_in[2]), NB);
  GetData(fyB, mxGetPr(mx_in[3]), NB);

  mexCallMATLAB(3,mx_out,8,mx_in,"cvm_psetB");

  *jcurPtrB = mxIsLogicalScalarTrue(mx_out[0]);
  ret = (int)*mxGetPr(mx_out[1]);

  if (!mxIsEmpty(mx_out[2])) {
    UpdateUserData(mx_out[2], bckPb);
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


int mxW_CVodeSpilsPsolB(realtype t, N_Vector y,
                         N_Vector yB, N_Vector fyB,
                         N_Vector rB, N_Vector zB,
                         realtype gammaB, realtype deltaB,
                         int lrB, void *user_dataB, N_Vector tmpB)
{
  cvmPbData fwdPb, bckPb;
  mxArray *mx_in[7], *mx_out[3];
  int ret;

  /* Extract global interface data from user-data */
  bckPb = (cvmPbData) user_dataB;
  fwdPb = bckPb->fwd;

  /* Inputs to the Matlab function */
  mx_in[0] = mxCreateDoubleScalar(t);           /* current t */   
  mx_in[1] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current y */
  mx_in[2] = mxCreateDoubleMatrix(NB,1,mxREAL); /* current yB */
  mx_in[3] = mxCreateDoubleMatrix(NB,1,mxREAL); /* current fyB */
  mx_in[4] = mxCreateDoubleMatrix(NB,1,mxREAL); /* right hand side rB */
  mx_in[5] = bckPb->PSOLfct;                    /* matlab function handle */
  mx_in[6] = bckPb->mtlb_data;                  /* matlab user data */
  
  /* Call matlab wrapper */
  GetData(y, mxGetPr(mx_in[1]), N);
  GetData(yB, mxGetPr(mx_in[2]), NB);
  GetData(fyB, mxGetPr(mx_in[3]), NB);
  GetData(rB, mxGetPr(mx_in[4]), NB);

  mexCallMATLAB(3,mx_out,7,mx_in,"cvm_psolB");

  PutData(zB, mxGetPr(mx_out[0]), NB);
  ret = (int)*mxGetPr(mx_out[1]);

  if (!mxIsEmpty(mx_out[2])) {
    UpdateUserData(mx_out[2], bckPb);
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

int mxW_CVodeBBDglocB(sunindextype NlocalB, realtype t, N_Vector y,
                       N_Vector yB, N_Vector gB, void *user_dataB)
{
  cvmPbData fwdPb, bckPb;
  mxArray *mx_in[5], *mx_out[3];
  int ret;

  /* Extract global interface data from user-data */
  bckPb = (cvmPbData) user_dataB;
  fwdPb = bckPb->fwd;

  /* Inputs to the Matlab function */
  mx_in[0] = mxCreateDoubleScalar(t);           /* current t */
  mx_in[1] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current y */
  mx_in[2] = mxCreateDoubleMatrix(NB,1,mxREAL); /* current yB */
  mx_in[3] = bckPb->GLOCfct;                    /* matlab function handle */
  mx_in[4] = bckPb->mtlb_data;                  /* matlab user data */
  
  /* Call matlab wrapper */
  GetData(y, mxGetPr(mx_in[1]), N);
  GetData(yB, mxGetPr(mx_in[2]), NB);

  mexCallMATLAB(3,mx_out,5,mx_in,"cvm_glocB");

  PutData(gB, mxGetPr(mx_out[0]), NB);

  ret = (int)*mxGetPr(mx_out[1]);

  if (!mxIsEmpty(mx_out[2])) {
    UpdateUserData(mx_out[2], bckPb);
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

int mxW_CVodeBBDgcomB(sunindextype NlocalB, realtype t, N_Vector y, 
                       N_Vector yB, void *user_dataB)
{
  cvmPbData fwdPb, bckPb;
  mxArray *mx_in[5], *mx_out[2];
  int ret;

  /* Extract global interface data from user-data */
  bckPb = (cvmPbData) user_dataB;
  fwdPb = bckPb->fwd;

  /* Inputs to the Matlab function */
  mx_in[0] = mxCreateDoubleScalar(t);           /* current t */
  mx_in[1] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current y */
  mx_in[2] = mxCreateDoubleMatrix(NB,1,mxREAL); /* current yB */
  mx_in[3] = bckPb->GCOMfct;                    /* matlab function handle */
  mx_in[4] = bckPb->mtlb_data;                  /* matlab user data */
  
  /* Call matlab wrapper */
  GetData(y, mxGetPr(mx_in[1]), N);
  GetData(yB, mxGetPr(mx_in[2]), NB);

  mexCallMATLAB(2,mx_out,5,mx_in,"cvm_gcomB");

  ret = (int)*mxGetPr(mx_out[0]);

  if (!mxIsEmpty(mx_out[1])) {
    UpdateUserData(mx_out[1], bckPb);
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

void mxW_CVodeMonitor(int call, double t, N_Vector y, N_Vector yQ, N_Vector *yS,
                       cvmPbData fwdPb)
{
  mxArray *mx_in[8], *mx_out[1];
  double *tmp;
  int is;

  mx_in[0] = mxCreateDoubleScalar(call);            /* call type (0:first, 1:interm. 2:last) */
  mx_in[1] = mxCreateDoubleScalar(t);               /* current time */
  mx_in[2] = mxCreateDoubleMatrix(N,1,mxREAL);      /* current solution */
  if (quadr)
    mx_in[3] = mxCreateDoubleMatrix(Nq,1,mxREAL);   /* current quadratures */
  else
    mx_in[3] = mxCreateDoubleMatrix(0,0,mxREAL);
  mx_in[4] = mxCreateDoubleScalar(Ns);              /* number of sensitivities */
  if (fsa)
    mx_in[5] = mxCreateDoubleMatrix(N*Ns,1,mxREAL); /* current sensitivities */
  else
    mx_in[5] = mxCreateDoubleMatrix(0,0,mxREAL);
  mx_in[6] = fwdPb->MONfct;                         /* Matlab monitor function */
  mx_in[7] = fwdPb->MONdata;                        /* data for monitor function */

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
    UpdateMonitorData(mx_out[0], fwdPb);
  }

  mxDestroyArray(mx_in[0]);
  mxDestroyArray(mx_in[1]);
  mxDestroyArray(mx_in[2]);
  mxDestroyArray(mx_in[3]);
  mxDestroyArray(mx_in[4]);
  mxDestroyArray(mx_in[5]);
  mxDestroyArray(mx_out[0]);
}

void mxW_CVodeMonitorB(int call, int idxB, double tB, N_Vector yB, N_Vector yQB,
                       cvmPbData bckPb)
{
  mxArray *mx_in[7], *mx_out[1];

  mx_in[0] = mxCreateDoubleScalar(call);            /* 0: first, 1: interm. 2: last */
  mx_in[1] = mxCreateDoubleScalar(idxB);            /* index of current problem */
  mx_in[2] = mxCreateDoubleScalar(tB);              /* current time */
  mx_in[3] = mxCreateDoubleMatrix(NB,1,mxREAL);     /* current solution */
  if (quadrB)
    mx_in[4] = mxCreateDoubleMatrix(NqB,1,mxREAL);  /* current quadratures */
  else
    mx_in[4] = mxCreateDoubleMatrix(0,0,mxREAL);
  mx_in[5] = bckPb->MONfct;
  mx_in[6] = bckPb->MONdata;

  if (call == 1) {
    
    GetData(yB, mxGetPr(mx_in[3]), NB);

    if (quadrB)
      GetData(yQB, mxGetPr(mx_in[4]), NqB);
  }

  mexCallMATLAB(1,mx_out,7,mx_in,"cvm_monitorB");

  if (!mxIsEmpty(mx_out[0])) {
    UpdateMonitorData(mx_out[0], bckPb);
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

static void UpdateUserData(mxArray *new_mtlb_data, cvmPbData pb)
{
  mexUnlock();
  mxDestroyArray(pb->mtlb_data);
  pb->mtlb_data = mxDuplicateArray(new_mtlb_data);
  mexMakeArrayPersistent(pb->mtlb_data);
  mexLock();
}

static void UpdateMonitorData(mxArray *new_mtlb_data, cvmPbData pb)
{
  mexUnlock();
  mxDestroyArray(pb->MONdata);
  pb->MONdata = mxDuplicateArray(new_mtlb_data);
  mexMakeArrayPersistent(pb->MONdata);
  mexLock();
}

