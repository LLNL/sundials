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
 * IDA wrapper functions.
 * -----------------------------------------------------------------
 */

#include "idm.h"
#include "nvm.h"

static void UpdateUserData(mxArray *new_mtlb_data, idmPbData pb);
static void UpdateMonitorData(mxArray *new_mtlb_data, idmPbData pb);

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
 * Wrapper functions
 * ---------------------------------------------------------------------------------
 */

int mxW_IDARes(realtype tt, N_Vector yy, N_Vector yp,
               N_Vector rr, void *user_data)
{
  idmPbData fwdPb;
  mxArray *mx_in[5], *mx_out[3];
  int ret;

  /* Extract global interface data from user-data */
  fwdPb = (idmPbData) user_data;

  /* Inputs to the Matlab function */
  mx_in[0] = mxCreateDoubleScalar(tt);         /* current t */
  mx_in[1] = mxCreateDoubleMatrix(N,1,mxREAL); /* current yy */
  mx_in[2] = mxCreateDoubleMatrix(N,1,mxREAL); /* current yp */
  mx_in[3] = fwdPb->RESfct;                    /* matlab function handle */ 
  mx_in[4] = fwdPb->mtlb_data;                 /* matlab user data */

  /* Call matlab wrapper */
  GetData(yy, mxGetPr(mx_in[1]), N);
  GetData(yp, mxGetPr(mx_in[2]), N);

  mexCallMATLAB(3,mx_out,5,mx_in,"idm_res");

  PutData(rr, mxGetPr(mx_out[0]), N);
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

int mxW_IDAQuadFct(realtype tres, N_Vector yy, N_Vector yp, N_Vector ypQ,
                   void *user_data)
{
  idmPbData fwdPb;
  mxArray *mx_in[5], *mx_out[3];
  int ret;

  /* Extract global interface data from user-data */
  fwdPb = (idmPbData) user_data;

  /* Inputs to the Matlab function */
  mx_in[0] = mxCreateDoubleScalar(tres);       /* current t */
  mx_in[1] = mxCreateDoubleMatrix(N,1,mxREAL); /* current yy */
  mx_in[2] = mxCreateDoubleMatrix(N,1,mxREAL); /* current yp */
  mx_in[3] = fwdPb->QUADfct;                   /* matlab function handle */ 
  mx_in[4] = fwdPb->mtlb_data;                 /* matlab user data */

  /* Call matlab wrapper */
  GetData(yy, mxGetPr(mx_in[1]), N);
  GetData(yp, mxGetPr(mx_in[2]), N);

  mexCallMATLAB(3,mx_out,5,mx_in,"idm_rhsQ");

  PutData(ypQ, mxGetPr(mx_out[0]), Nq);
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

int mxW_IDAGfct(realtype t, N_Vector yy, N_Vector yp,
                realtype *gout, void *user_data)
{
  idmPbData fwdPb;
  double *gdata;
  int i, ret;
  mxArray *mx_in[5], *mx_out[3];
  
  /* Extract global interface data from user-data */
  fwdPb = (idmPbData) user_data;

  /* Inputs to the Matlab function */
  mx_in[0] = mxCreateDoubleScalar(t);          /* current t */
  mx_in[1] = mxCreateDoubleMatrix(N,1,mxREAL); /* current yy */
  mx_in[2] = mxCreateDoubleMatrix(N,1,mxREAL); /* current yp */
  mx_in[3] = fwdPb->Gfct;                      /* matlab function handle */
  mx_in[4] = fwdPb->mtlb_data;                 /* matlab user data */
  
  /* Call matlab wrapper */
  GetData(yy, mxGetPr(mx_in[1]), N);
  GetData(yp, mxGetPr(mx_in[2]), N);

  mexCallMATLAB(3,mx_out,5,mx_in,"idm_root");

  gdata = mxGetPr(mx_out[0]);
  for (i=0;i<Ng;i++) gout[i] = gdata[i];

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

int mxW_IDADenseJac(sunindextype Neq, realtype tt, realtype c_j, 
                    N_Vector yy, N_Vector yp, N_Vector rr,
                    DlsMat Jac, void *user_data, 
                    N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  idmPbData fwdPb;
  double *J_data;
  int i, ret;
  mxArray *mx_in[7], *mx_out[3];

  /* Extract global interface data from user-data */
  fwdPb = (idmPbData) user_data;

  /* Inputs to the Matlab function */
  mx_in[0] = mxCreateDoubleScalar(tt);          /* current t */  
  mx_in[1] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current yy */
  mx_in[2] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current yp */
  mx_in[3] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current rr */
  mx_in[4] = mxCreateDoubleScalar(c_j);         /* current c_j */  
  mx_in[5] = fwdPb->JACfct;                     /* matlab function handle */
  mx_in[6] = fwdPb->mtlb_data;                  /* matlab user data */
  
  /* Call matlab wrapper */
  GetData(yy, mxGetPr(mx_in[1]), N);
  GetData(yp, mxGetPr(mx_in[2]), N);
  GetData(rr, mxGetPr(mx_in[3]), N);

  mexCallMATLAB(3,mx_out,7,mx_in,"idm_djac");

  /* Extract data */
  J_data = mxGetPr(mx_out[0]);
  for (i=0;i<N;i++)
    memcpy(DENSE_COL(Jac,i), J_data + i*N, N*sizeof(double));

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

int mxW_IDABandJac(sunindextype Neq, sunindextype mupper, sunindextype mlower, 
                   realtype tt, realtype c_j, 
                   N_Vector yy, N_Vector yp, N_Vector rr, 
                   DlsMat Jac, void *user_data,
                   N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  idmPbData fwdPb;
  double *J_data;
  sunindextype eband, i;
  int ret;
  mxArray *mx_in[7], *mx_out[3];

  /* Extract global interface data from user-data */
  fwdPb = (idmPbData) user_data;

  /* Inputs to the Matlab function */
  mx_in[0] = mxCreateDoubleScalar(tt);          /* current t */
  mx_in[1] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current yy */
  mx_in[2] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current yp */
  mx_in[3] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current rr */
  mx_in[4] = mxCreateDoubleScalar(c_j);         /* current c_j */
  mx_in[5] = fwdPb->JACfct;                     /* matlab function handle */
  mx_in[6] = fwdPb->mtlb_data;                  /* matlab user data */
  
  /* Call matlab wrapper */
  GetData(yy, mxGetPr(mx_in[1]), N);
  GetData(yp, mxGetPr(mx_in[2]), N);
  GetData(rr, mxGetPr(mx_in[3]), N);

  mexCallMATLAB(3,mx_out,7,mx_in,"idm_bjac");

  /* Extract data */
  eband =  mupper + mlower + 1;
  J_data = mxGetPr(mx_out[0]);
  for (i=0;i<N;i++)
    memcpy(BAND_COL(Jac,i) - mupper, J_data + i*eband, eband*sizeof(double));

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

int mxW_IDASpilsJac(realtype tt,
                    N_Vector yy, N_Vector yp, N_Vector rr,
                    N_Vector v, N_Vector Jv,
                    realtype c_j, void *user_data,
                    N_Vector tmp1, N_Vector tmp2)
{
  idmPbData fwdPb;
  mxArray *mx_in[8], *mx_out[3];
  int ret;

  /* Extract global interface data from user-data */
  fwdPb = (idmPbData) user_data;

  /* Inputs to the Matlab function */
  mx_in[0] = mxCreateDoubleScalar(tt);          /* current t */ 
  mx_in[1] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current yy */
  mx_in[2] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current yp */
  mx_in[3] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current rr */
  mx_in[4] = mxCreateDoubleMatrix(N,1,mxREAL);  /* vector v */
  mx_in[5] = mxCreateDoubleScalar(c_j);         /* current c_j */ 
  mx_in[6] = fwdPb->JACfct;                     /* matlab function handle */
  mx_in[7] = fwdPb->mtlb_data;                  /* matlab user data */
  
  /* Call matlab wrapper */
  GetData(yy, mxGetPr(mx_in[1]), N);
  GetData(yp, mxGetPr(mx_in[2]), N);
  GetData(rr, mxGetPr(mx_in[3]), N);
  GetData(v, mxGetPr(mx_in[4]), N);

  mexCallMATLAB(3,mx_out,8,mx_in,"idm_jtv");

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
  mxDestroyArray(mx_in[4]);
  mxDestroyArray(mx_in[5]);
  mxDestroyArray(mx_out[0]);
  mxDestroyArray(mx_out[1]);
  mxDestroyArray(mx_out[2]);

  return(ret);
}

int mxW_IDASpilsPset(realtype tt,
                     N_Vector yy, N_Vector yp, N_Vector rr,
                     realtype c_j, void *user_data,
                     N_Vector tmp1, N_Vector tmp2,
                     N_Vector tmp3)
{
  idmPbData fwdPb;
  mxArray *mx_in[7], *mx_out[2];
  int ret;

  /* Extract global interface data from user-data */
  fwdPb = (idmPbData) user_data;

  /* Inputs to the Matlab function */
  mx_in[0] = mxCreateDoubleScalar(tt);          /* current t */
  mx_in[1] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current yy */
  mx_in[2] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current yp */
  mx_in[3] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current rr */
  mx_in[4] = mxCreateLogicalScalar(c_j);        /* current c_j */
  mx_in[5] = fwdPb->PSETfct;                    /* matlab function handle */
  mx_in[6] = fwdPb->mtlb_data;                  /* matlab user data */
  
  /* Call matlab wrapper */
  GetData(yy, mxGetPr(mx_in[1]), N);
  GetData(yp, mxGetPr(mx_in[2]), N);
  GetData(rr, mxGetPr(mx_in[3]), N);

  mexCallMATLAB(2,mx_out,7,mx_in,"idm_pset");

  ret = (int)*mxGetPr(mx_out[0]);

  if (!mxIsEmpty(mx_out[1])) {
    UpdateUserData(mx_out[1], fwdPb);
  }

  /* Free temporary space */
  mxDestroyArray(mx_in[0]);
  mxDestroyArray(mx_in[1]);
  mxDestroyArray(mx_in[2]);
  mxDestroyArray(mx_in[3]);
  mxDestroyArray(mx_in[4]);
  mxDestroyArray(mx_out[0]);
  mxDestroyArray(mx_out[1]);

  return(ret);
}

int mxW_IDASpilsPsol(realtype tt,
                     N_Vector yy, N_Vector yp, N_Vector rr,
                     N_Vector rvec, N_Vector zvec,
                     realtype c_j, realtype delta, void *user_data,
                     N_Vector tmp)
{
  idmPbData fwdPb;
  mxArray *mx_in[8], *mx_out[3];
  int ret;

  /* Extract global interface data from user-data */
  fwdPb = (idmPbData) user_data;

  /* Inputs to the Matlab function */
  mx_in[0] = mxCreateDoubleScalar(tt);         /* current t */   
  mx_in[1] = mxCreateDoubleMatrix(N,1,mxREAL); /* current yy */
  mx_in[2] = mxCreateDoubleMatrix(N,1,mxREAL); /* current yp */
  mx_in[3] = mxCreateDoubleMatrix(N,1,mxREAL); /* current rr */
  mx_in[4] = mxCreateDoubleMatrix(N,1,mxREAL); /* right hand side r */
  mx_in[5] = mxCreateDoubleScalar(c_j);        /* current c_j */   
  mx_in[6] = fwdPb->PSOLfct;                   /* matlab function handle */
  mx_in[7] = fwdPb->mtlb_data;                 /* matlab user data */
  
  /* Call matlab wrapper */
  GetData(yy, mxGetPr(mx_in[1]), N);
  GetData(yp, mxGetPr(mx_in[2]), N);
  GetData(rr, mxGetPr(mx_in[3]), N);
  GetData(rvec, mxGetPr(mx_in[4]), N);

  mexCallMATLAB(3,mx_out,8,mx_in,"idm_psol");

  PutData(zvec, mxGetPr(mx_out[0]), N);
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
  mxDestroyArray(mx_in[5]);
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

int mxW_IDABBDgloc(sunindextype Nlocal, realtype tt,
                   N_Vector yy, N_Vector yp, N_Vector gval,
                   void *user_data)
{
  idmPbData fwdPb;
  mxArray *mx_in[5], *mx_out[3];
  int ret;

  /* Extract global interface data from user-data */
  fwdPb = (idmPbData) user_data;

  /* Inputs to the Matlab function */
  mx_in[0] = mxCreateDoubleScalar(tt);          /* current t */
  mx_in[1] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current yy */
  mx_in[2] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current yp */
  mx_in[3] = fwdPb->GLOCfct;                    /* matlab function handle */
  mx_in[4] = fwdPb->mtlb_data;                  /* matlab user data */
  
  /* Call matlab wrapper */
  GetData(yy, mxGetPr(mx_in[1]), N);
  GetData(yp, mxGetPr(mx_in[2]), N);

  mexCallMATLAB(3,mx_out,5,mx_in,"idm_gloc");

  PutData(gval, mxGetPr(mx_out[0]), N);
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

int mxW_IDABBDgcom(sunindextype Nlocal, realtype tt,
                   N_Vector yy, N_Vector yp,
                   void *user_data)
{
  idmPbData fwdPb;
  mxArray *mx_in[5], *mx_out[2];
  int ret;

  /* Extract global interface data from user-data */
  fwdPb = (idmPbData) user_data;

  /* Inputs to the Matlab function */
  mx_in[0] = mxCreateDoubleScalar(tt);          /* current t */
  mx_in[1] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current yy */
  mx_in[2] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current yp */
  mx_in[3] = fwdPb->GCOMfct;                    /* matlab function handle */
  mx_in[4] = fwdPb->mtlb_data;                  /* matlab user data */
  
  /* Call matlab wrapper */
  GetData(yy, mxGetPr(mx_in[1]), N);
  GetData(yp, mxGetPr(mx_in[2]), N);

  mexCallMATLAB(2,mx_out,5,mx_in,"idm_gcom");

  ret = (int)*mxGetPr(mx_out[0]);

  if (!mxIsEmpty(mx_out[1])) {
    UpdateUserData(mx_out[1], fwdPb);
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

int mxW_IDASensRes(int Nsens, realtype tres, 
                   N_Vector yy, N_Vector yp, N_Vector rr,
                   N_Vector *yyS, N_Vector *ypS, N_Vector *rrS,
                   void *user_data,
                   N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  idmPbData fwdPb;
  mxArray *mx_in[9], *mx_out[3];
  int is, ret;
  double *tmp_yyS, *tmp_ypS, *tmp_rrS;

  /* Extract global interface data from user-data */
  fwdPb = (idmPbData) user_data;

  /* Inputs to the Matlab function */
  mx_in[0] = mxCreateDoubleScalar(tres);          /* current t */
  mx_in[1] = mxCreateDoubleMatrix(N,1,mxREAL);    /* current yy */
  mx_in[2] = mxCreateDoubleMatrix(N,1,mxREAL);    /* current yp */
  mx_in[3] = mxCreateDoubleMatrix(N,1,mxREAL);    /* current rr */
  mx_in[4] = mxCreateDoubleScalar(Ns);            /* number of sensitivities */
  mx_in[5] = mxCreateDoubleMatrix(N*Ns,1,mxREAL); /* current yyS */
  mx_in[6] = mxCreateDoubleMatrix(N*Ns,1,mxREAL); /* current ypS */
  mx_in[7] = fwdPb->SRESfct;                      /* matlab function handle */      
  mx_in[8] = fwdPb->mtlb_data;                    /* matlab user data */
  
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
    UpdateUserData(mx_out[2], fwdPb);
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

int mxW_IDAResB(realtype tt, 
                N_Vector yy, N_Vector yp,
                N_Vector yyB, N_Vector ypB, N_Vector rrB,
                void *user_dataB)
{
  idmPbData fwdPb, bckPb;
  mxArray *mx_in[8], *mx_out[3];
  int ret;

  /* Extract global interface data from user-data */
  bckPb = (idmPbData) user_dataB;
  fwdPb = bckPb->fwd;

  /* Inputs to the Matlab function */
  mx_in[0] = mxCreateDoubleScalar(0.0);         /* type=0: not dependent on yS */
  mx_in[1] = mxCreateDoubleScalar(tt);          /* current t */
  mx_in[2] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current yy */
  mx_in[3] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current yp */
  mx_in[4] = mxCreateDoubleMatrix(NB,1,mxREAL); /* current yyB */
  mx_in[5] = mxCreateDoubleMatrix(NB,1,mxREAL); /* current ypB */
  mx_in[6] = bckPb->RESfct;                     /* matlab function handle */ 
  mx_in[7] = bckPb->mtlb_data;                  /* matlab user data */
  
  /* Call matlab wrapper */
  GetData(yy, mxGetPr(mx_in[2]), N);
  GetData(yp, mxGetPr(mx_in[3]), N);
  GetData(yyB, mxGetPr(mx_in[4]), NB);
  GetData(ypB, mxGetPr(mx_in[5]), NB);

  mexCallMATLAB(3,mx_out,8,mx_in,"idm_resB");

  PutData(rrB, mxGetPr(mx_out[0]), NB);

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

int mxW_IDAResBS(realtype tt, 
                 N_Vector yy, N_Vector yp,
                 N_Vector *yyS, N_Vector *ypS,
                 N_Vector yyB, N_Vector ypB, N_Vector rrB,
                 void *user_dataB)
{
  idmPbData fwdPb, bckPb;
  mxArray *mx_in[11], *mx_out[3];
  double *tmp_yyS, *tmp_ypS;
  int is, ret;

  /* Extract global interface data from user-data */
  bckPb = (idmPbData) user_dataB;
  fwdPb = bckPb->fwd;

  /* Inputs to the Matlab function */
  mx_in[0] = mxCreateDoubleScalar(1.0);           /* type=1: dependent on yS */
  mx_in[1] = mxCreateDoubleScalar(tt);            /* current t */
  mx_in[2] = mxCreateDoubleMatrix(N,1,mxREAL);    /* current yy */
  mx_in[3] = mxCreateDoubleMatrix(N,1,mxREAL);    /* current yp */
  mx_in[4] = mxCreateDoubleScalar(Ns);            /* number of sensitivities */
  mx_in[5] = mxCreateDoubleMatrix(N*Ns,1,mxREAL); /* current yyS */
  mx_in[6] = mxCreateDoubleMatrix(N*Ns,1,mxREAL); /* current ypS */
  mx_in[7] = mxCreateDoubleMatrix(NB,1,mxREAL);   /* current yyB */
  mx_in[8] = mxCreateDoubleMatrix(NB,1,mxREAL);   /* current ypB */
  mx_in[9] = bckPb->RESfct;                       /* matlab function handle */ 
  mx_in[10] = bckPb->mtlb_data;                   /* matlab user data */
  
  /* Call matlab wrapper */

  GetData(yy, mxGetPr(mx_in[2]), N);
  GetData(yp, mxGetPr(mx_in[3]), N);

  tmp_yyS = mxGetPr(mx_in[5]);
  tmp_ypS = mxGetPr(mx_in[6]);
  for (is=0; is<Ns; is++) {
    GetData(yyS[is], &tmp_yyS[is*N], N);
    GetData(ypS[is], &tmp_ypS[is*N], N);
  }

  GetData(yyB, mxGetPr(mx_in[7]), NB);
  GetData(ypB, mxGetPr(mx_in[8]), NB);

  mexCallMATLAB(3,mx_out,11,mx_in,"idm_resB");

  PutData(rrB, mxGetPr(mx_out[0]), NB);

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
  mxDestroyArray(mx_in[6]);
  mxDestroyArray(mx_in[7]);
  mxDestroyArray(mx_in[8]);
  mxDestroyArray(mx_out[0]);
  mxDestroyArray(mx_out[1]);
  mxDestroyArray(mx_out[2]);

  return(ret);
}

int mxW_IDAQuadFctB(realtype tt, 
                    N_Vector yy, N_Vector yp, 
                    N_Vector yyB, N_Vector ypB,
                    N_Vector ypQB, void *user_dataB)
{
  idmPbData fwdPb, bckPb;
  mxArray *mx_in[8], *mx_out[3];
  int ret;

  /* Extract global interface data from user-data */
  bckPb = (idmPbData) user_dataB;
  fwdPb = bckPb->fwd;

  /* Inputs to the Matlab function */
  mx_in[0] = mxCreateDoubleScalar(0.0);         /* type=0: not dependent on yS */
  mx_in[1] = mxCreateDoubleScalar(tt);          /* current t */
  mx_in[2] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current yy */
  mx_in[3] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current yp */
  mx_in[4] = mxCreateDoubleMatrix(NB,1,mxREAL); /* current yyB */
  mx_in[5] = mxCreateDoubleMatrix(NB,1,mxREAL); /* current ypB */
  mx_in[6] = bckPb->QUADfct;                    /* matlab function handle */ 
  mx_in[7] = bckPb->mtlb_data;                  /* matlab user data */

  /* Call matlab wrapper */
  GetData(yy, mxGetPr(mx_in[2]), N);
  GetData(yp, mxGetPr(mx_in[3]), N);
  GetData(yyB, mxGetPr(mx_in[4]), NB);
  GetData(ypB, mxGetPr(mx_in[5]), NB);

  mexCallMATLAB(3,mx_out,8,mx_in,"idm_rhsQB");

  PutData(ypQB, mxGetPr(mx_out[0]), NqB);

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

int mxW_IDAQuadFctBS(realtype tt, 
                     N_Vector yy, N_Vector yp, 
                     N_Vector *yyS, N_Vector *ypS,
                     N_Vector yyB, N_Vector ypB,
                     N_Vector ypQB, void *user_dataB)
{
  idmPbData fwdPb, bckPb;
  mxArray *mx_in[11], *mx_out[3];
  double *tmp_yyS, *tmp_ypS;
  int is, ret;

  /* Extract global interface data from user-data */
  bckPb = (idmPbData) user_dataB;
  fwdPb = bckPb->fwd;

  /* Inputs to the Matlab function */
  mx_in[0] = mxCreateDoubleScalar(1.0);           /* type=1: dependent on yS */
  mx_in[1] = mxCreateDoubleScalar(tt);            /* current t */
  mx_in[2] = mxCreateDoubleMatrix(N,1,mxREAL);    /* current yy */
  mx_in[3] = mxCreateDoubleMatrix(N,1,mxREAL);    /* current yp */
  mx_in[4] = mxCreateDoubleScalar(Ns);            /* number of sensitivities */
  mx_in[5] = mxCreateDoubleMatrix(N*Ns,1,mxREAL); /* current yyS */
  mx_in[6] = mxCreateDoubleMatrix(N*Ns,1,mxREAL); /* current ypS */
  mx_in[7] = mxCreateDoubleMatrix(NB,1,mxREAL);   /* current yyB */
  mx_in[8] = mxCreateDoubleMatrix(NB,1,mxREAL);   /* current ypB */
  mx_in[9] = bckPb->QUADfct;                      /* matlab function handle */ 
  mx_in[10] = bckPb->mtlb_data;                   /* matlab user data */

  /* Call matlab wrapper */

  GetData(yy, mxGetPr(mx_in[2]), N);
  GetData(yp, mxGetPr(mx_in[3]), N);

  tmp_yyS = mxGetPr(mx_in[5]);
  tmp_ypS = mxGetPr(mx_in[6]);
  for (is=0; is<Ns; is++) {
    GetData(yyS[is], &tmp_yyS[is*N], N);
    GetData(ypS[is], &tmp_ypS[is*N], N);
  }

  GetData(yyB, mxGetPr(mx_in[7]), NB);
  GetData(ypB, mxGetPr(mx_in[8]), NB);

  mexCallMATLAB(3,mx_out,11,mx_in,"idm_rhsQB");

  PutData(ypQB, mxGetPr(mx_out[0]), NqB);

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
  mxDestroyArray(mx_in[6]);
  mxDestroyArray(mx_in[7]);
  mxDestroyArray(mx_in[8]);
  mxDestroyArray(mx_out[0]);
  mxDestroyArray(mx_out[1]);
  mxDestroyArray(mx_out[2]);

  return(ret);
}

int mxW_IDADenseJacB(sunindextype NeqB,
                     realtype tt, realtype c_jB,
                     N_Vector yy, N_Vector yp,
                     N_Vector yyB, N_Vector ypB, N_Vector rrB,
                     DlsMat JacB, void *user_dataB, 
                     N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B)
{
  idmPbData fwdPb, bckPb;
  double *JB_data;
  mxArray *mx_in[9], *mx_out[3];
  int i, ret;

  /* Extract global interface data from user-data */
  bckPb = (idmPbData) user_dataB;
  fwdPb = bckPb->fwd;

  /* Inputs to the Matlab function */
  mx_in[0] = mxCreateDoubleScalar(tt);          /* current t */  
  mx_in[1] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current yy */
  mx_in[2] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current yp */
  mx_in[3] = mxCreateDoubleMatrix(NB,1,mxREAL); /* current yyB */
  mx_in[4] = mxCreateDoubleMatrix(NB,1,mxREAL); /* current ypB */
  mx_in[5] = mxCreateDoubleMatrix(NB,1,mxREAL); /* current rrB */
  mx_in[6] = mxCreateDoubleScalar(c_jB);        /* current c_jB */  
  mx_in[7] = bckPb->JACfct;                     /* matlab function handle */
  mx_in[8] = bckPb->mtlb_data;                  /* matlab user data */
  
  /* Call matlab wrapper */
  GetData(yy, mxGetPr(mx_in[1]), N);
  GetData(yp, mxGetPr(mx_in[2]), N);
  GetData(yyB, mxGetPr(mx_in[3]), NB);
  GetData(ypB, mxGetPr(mx_in[4]), NB);
  GetData(rrB, mxGetPr(mx_in[5]), NB);

  mexCallMATLAB(3,mx_out,9,mx_in,"idm_djacB");

  JB_data = mxGetPr(mx_out[0]);
  for (i=0;i<NB;i++)
    memcpy(DENSE_COL(JacB,i), JB_data + i*NB, NB*sizeof(double));

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
  mxDestroyArray(mx_in[6]);
  mxDestroyArray(mx_out[0]);
  mxDestroyArray(mx_out[1]);
  mxDestroyArray(mx_out[2]);

  return(ret);
}

int mxW_IDABandJacB(sunindextype NeqB, sunindextype mupperB, sunindextype mlowerB, 
                    realtype tt, realtype c_jB, 
                    N_Vector yy, N_Vector yp,
                    N_Vector yyB, N_Vector ypB, N_Vector rrB,
                    DlsMat JacB, void *user_dataB,
                    N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B)
{
  idmPbData fwdPb, bckPb;
  double *JB_data;
  mxArray *mx_in[9], *mx_out[3];
  sunindextype ebandB, i;
  int ret;

  /* Extract global interface data from user-data */
  bckPb = (idmPbData) user_dataB;
  fwdPb = bckPb->fwd;

  /* Inputs to the Matlab function */
  mx_in[0] = mxCreateDoubleScalar(tt);          /* current t */
  mx_in[1] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current yy */
  mx_in[2] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current yp */
  mx_in[3] = mxCreateDoubleMatrix(NB,1,mxREAL); /* current yyB */
  mx_in[4] = mxCreateDoubleMatrix(NB,1,mxREAL); /* current ypB */
  mx_in[5] = mxCreateDoubleMatrix(NB,1,mxREAL); /* current rrB */
  mx_in[6] = mxCreateDoubleScalar(c_jB);        /* current c_jB */
  mx_in[7] = bckPb->JACfct;                     /* matlab function handle */
  mx_in[8] = bckPb->mtlb_data;                  /* matlab user data */
  
  /* Call matlab wrapper */
  GetData(yy, mxGetPr(mx_in[1]), N);
  GetData(yp, mxGetPr(mx_in[2]), N);
  GetData(yyB, mxGetPr(mx_in[3]), NB);
  GetData(ypB, mxGetPr(mx_in[4]), NB);
  GetData(rrB, mxGetPr(mx_in[5]), NB);

  mexCallMATLAB(3,mx_out,9,mx_in,"idm_bjacB");

  ebandB =  mupperB + mlowerB + 1;
  JB_data = mxGetPr(mx_out[0]);
  for (i=0;i<NB;i++)
    memcpy(BAND_COL(JacB,i) - mupperB, JB_data + i*ebandB, ebandB*sizeof(double));
    
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
  mxDestroyArray(mx_in[6]);
  mxDestroyArray(mx_out[0]);
  mxDestroyArray(mx_out[1]);
  mxDestroyArray(mx_out[2]);

  return(ret);
}

int mxW_IDASpilsJacB(realtype tt,
                     N_Vector yy, N_Vector yp,
                     N_Vector yyB, N_Vector ypB, N_Vector rrB,
                     N_Vector vB, N_Vector JvB, 
                     realtype c_jB, void *user_dataB, 
                     N_Vector tmp1B, N_Vector tmp2B)
{
  idmPbData fwdPb, bckPb;
  mxArray *mx_in[10], *mx_out[3];
  int ret;

  /* Extract global interface data from user-data */
  bckPb = (idmPbData) user_dataB;
  fwdPb = bckPb->fwd;

  /* Inputs to the Matlab function */
  mx_in[0] = mxCreateDoubleScalar(tt);          /* current t */ 
  mx_in[1] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current yy */
  mx_in[2] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current yp */
  mx_in[3] = mxCreateDoubleMatrix(NB,1,mxREAL); /* current yyB */
  mx_in[4] = mxCreateDoubleMatrix(NB,1,mxREAL); /* current ypB */
  mx_in[5] = mxCreateDoubleMatrix(NB,1,mxREAL); /* current rrB */
  mx_in[6] = mxCreateDoubleMatrix(NB,1,mxREAL); /* vector vB */
  mx_in[7] = mxCreateDoubleScalar(c_jB);        /* current c_jB */ 
  mx_in[8] = bckPb->JACfct;                     /* matlab function handle */
  mx_in[9] = bckPb->mtlb_data;                  /* matlab user data */
  
  /* Call matlab wrapper */
  GetData(yy, mxGetPr(mx_in[1]), N);
  GetData(yp, mxGetPr(mx_in[2]), N);
  GetData(yyB, mxGetPr(mx_in[3]), NB);
  GetData(ypB, mxGetPr(mx_in[4]), NB);
  GetData(rrB, mxGetPr(mx_in[5]), NB);
  GetData(vB, mxGetPr(mx_in[6]), NB);

  mexCallMATLAB(3,mx_out,10,mx_in,"idm_jtvB");

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
  mxDestroyArray(mx_in[5]);
  mxDestroyArray(mx_in[6]);
  mxDestroyArray(mx_in[7]);
  mxDestroyArray(mx_out[0]);
  mxDestroyArray(mx_out[1]);
  mxDestroyArray(mx_out[2]);

  return(ret);

}

int mxW_IDASpilsPsetB(realtype tt, 
                      N_Vector yy, N_Vector yp,
                      N_Vector yyB, N_Vector ypB, N_Vector rrB, 
                      realtype c_jB, void *user_dataB,
                      N_Vector tmp1B, N_Vector tmp2B, 
                      N_Vector tmp3B)
{
  idmPbData fwdPb, bckPb;
  mxArray *mx_in[9], *mx_out[2];
  int ret;

  /* Extract global interface data from user-data */
  bckPb = (idmPbData) user_dataB;
  fwdPb = bckPb->fwd;

  /* Inputs to the Matlab function */
  mx_in[0] = mxCreateDoubleScalar(tt);          /* current t */
  mx_in[1] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current yy */
  mx_in[2] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current yp */
  mx_in[3] = mxCreateDoubleMatrix(NB,1,mxREAL); /* current yyB */
  mx_in[4] = mxCreateDoubleMatrix(NB,1,mxREAL); /* current ypB */
  mx_in[5] = mxCreateDoubleMatrix(NB,1,mxREAL); /* current rrB */
  mx_in[6] = mxCreateDoubleScalar(c_jB);        /* current c_jB */
  mx_in[7] = bckPb->PSETfct;                    /* matlab function handle */
  mx_in[8] = bckPb->mtlb_data;                  /* matlab user data */
  
  /* Call matlab wrapper */
  GetData(yy, mxGetPr(mx_in[1]), N);
  GetData(yp, mxGetPr(mx_in[2]), N);
  GetData(yyB, mxGetPr(mx_in[3]), NB);
  GetData(ypB, mxGetPr(mx_in[4]), NB);
  GetData(rrB, mxGetPr(mx_in[5]), NB);

  mexCallMATLAB(2,mx_out,9,mx_in,"idm_psetB");

  ret = (int)*mxGetPr(mx_out[0]);

  if (!mxIsEmpty(mx_out[1])) {
    UpdateUserData(mx_out[1], bckPb);
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

  return(ret);

}

int mxW_IDASpilsPsolB(realtype tt, 
                      N_Vector yy, N_Vector yp,
                      N_Vector yyB, N_Vector ypB, N_Vector rrB, 
                      N_Vector rvecB, N_Vector zvecB,
                      realtype c_jB, realtype deltaB,
                      void *user_dataB, N_Vector tmpB)
{
  idmPbData fwdPb, bckPb;
  mxArray *mx_in[10], *mx_out[3];
  int ret;

  /* Extract global interface data from user-data */
  bckPb = (idmPbData) user_dataB;
  fwdPb = bckPb->fwd;

  /* Inputs to the Matlab function */
  mx_in[0] = mxCreateDoubleScalar(tt);          /* current t */   
  mx_in[1] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current yy */
  mx_in[2] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current yp */
  mx_in[3] = mxCreateDoubleMatrix(NB,1,mxREAL); /* current yyB */
  mx_in[4] = mxCreateDoubleMatrix(NB,1,mxREAL); /* current ypB */
  mx_in[5] = mxCreateDoubleMatrix(NB,1,mxREAL); /* current rrB */
  mx_in[6] = mxCreateDoubleMatrix(NB,1,mxREAL); /* right hand side rB */
  mx_in[7] = mxCreateDoubleScalar(c_jB);        /* current c_jB */   
  mx_in[8] = bckPb->PSOLfct;                    /* matlab function handle */
  mx_in[9] = bckPb->mtlb_data;                  /* matlab user data */
  
  /* Call matlab wrapper */
  GetData(yy, mxGetPr(mx_in[1]), N);
  GetData(yp, mxGetPr(mx_in[2]), N);
  GetData(yyB, mxGetPr(mx_in[3]), NB);
  GetData(ypB, mxGetPr(mx_in[4]), NB);
  GetData(rrB, mxGetPr(mx_in[5]), NB);
  GetData(rvecB, mxGetPr(mx_in[6]), NB);

  mexCallMATLAB(3,mx_out,10,mx_in,"idm_psolB");

  PutData(zvecB, mxGetPr(mx_out[0]), NB);
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
  mxDestroyArray(mx_in[6]);
  mxDestroyArray(mx_in[7]);
  mxDestroyArray(mx_out[0]);
  mxDestroyArray(mx_out[1]);
  mxDestroyArray(mx_out[2]);

  return(ret);

}

int mxW_IDABBDglocB(sunindextype NlocalB, realtype tt,
                    N_Vector yy, N_Vector yp, 
                    N_Vector yyB, N_Vector ypB, N_Vector gvalB,
                    void *user_dataB)
{
  idmPbData fwdPb, bckPb;
  mxArray *mx_in[7], *mx_out[3];
  int ret;

  /* Extract global interface data from user-data */
  bckPb = (idmPbData) user_dataB;
  fwdPb = bckPb->fwd;

  /* Inputs to the Matlab function */
  mx_in[0] = mxCreateDoubleScalar(tt);          /* current t */
  mx_in[1] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current yy */
  mx_in[2] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current yp */
  mx_in[3] = mxCreateDoubleMatrix(NB,1,mxREAL); /* current yyB */
  mx_in[4] = mxCreateDoubleMatrix(NB,1,mxREAL); /* current ypB */
  mx_in[5] = bckPb->GLOCfct;                    /* matlab function handle */
  mx_in[6] = bckPb->mtlb_data;                  /* matlab user data */
  
  /* Call matlab wrapper */
  GetData(yy, mxGetPr(mx_in[1]), N);
  GetData(yp, mxGetPr(mx_in[2]), N);
  GetData(yyB, mxGetPr(mx_in[3]), NB);
  GetData(ypB, mxGetPr(mx_in[4]), NB);

  mexCallMATLAB(3,mx_out,7,mx_in,"idm_glocB");

  PutData(gvalB, mxGetPr(mx_out[0]), NB);

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

int mxW_IDABBDgcomB(sunindextype NlocalB, realtype tt,
                    N_Vector yy, N_Vector yp,
                    N_Vector yyB, N_Vector ypB,
                    void *user_dataB)
{
  idmPbData fwdPb, bckPb;
  mxArray *mx_in[7], *mx_out[2];
  int ret;

  /* Extract global interface data from user-data */
  bckPb = (idmPbData) user_dataB;
  fwdPb = bckPb->fwd;

  /* Inputs to the Matlab function */
  mx_in[0] = mxCreateDoubleScalar(tt);          /* current t */
  mx_in[1] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current yy */
  mx_in[2] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current yp */
  mx_in[3] = mxCreateDoubleMatrix(NB,1,mxREAL); /* current yyB */
  mx_in[4] = mxCreateDoubleMatrix(NB,1,mxREAL); /* current ypB */
  mx_in[5] = bckPb->GCOMfct;                    /* matlab function handle */
  mx_in[6] = bckPb->mtlb_data;                  /* matlab user data */
  
  /* Call matlab wrapper */
  GetData(yy, mxGetPr(mx_in[1]), N);
  GetData(yp, mxGetPr(mx_in[2]), N);
  GetData(yyB, mxGetPr(mx_in[3]), NB);
  GetData(ypB, mxGetPr(mx_in[4]), NB);

  mexCallMATLAB(2,mx_out,7,mx_in,"idm_gcomB");

  ret = (int)*mxGetPr(mx_out[0]);

  if (!mxIsEmpty(mx_out[1])) {
    UpdateUserData(mx_out[1], bckPb);
  }

  /* Free temporary space */
  mxDestroyArray(mx_in[0]);
  mxDestroyArray(mx_in[1]);
  mxDestroyArray(mx_in[2]);
  mxDestroyArray(mx_in[3]);
  mxDestroyArray(mx_in[4]);
  mxDestroyArray(mx_out[0]);
  mxDestroyArray(mx_out[1]);

  return(ret);
}

/*
 * ---------------------------------------------------------------------------------
 * Wrapper around matlab monitor function
 * ---------------------------------------------------------------------------------
 */

void mxW_IDAMonitor(int call, double t, 
                    N_Vector yy,
                    N_Vector yQ, 
                    N_Vector *yyS,
                    idmPbData fwdPb)
{
  mxArray *mx_in[8], *mx_out[1];
  double *tmp_yyS;
  int is;

  mx_in[0] = mxCreateDoubleScalar(call);            /* call type (0:first, 1:interm. 2:last) */
  mx_in[1] = mxCreateDoubleScalar(t);               /* current time */
  mx_in[2] = mxCreateDoubleMatrix(N,1,mxREAL);      /* current yy */
  if (quadr) {
    mx_in[3] = mxCreateDoubleMatrix(Nq,1,mxREAL);   /* current quadratures */
  } else {
    mx_in[3] = mxCreateDoubleMatrix(0,0,mxREAL);
  }
  mx_in[4] = mxCreateDoubleScalar(Ns);              /* number of sensitivities */
  if (fsa) {
    mx_in[5] = mxCreateDoubleMatrix(N*Ns,1,mxREAL); /* current yyS */
  } else {
    mx_in[5] = mxCreateDoubleMatrix(0,0,mxREAL);
  }
  mx_in[6] = fwdPb->MONfct;                         /* Matlab monitor function */
  mx_in[7] = fwdPb->MONdata;                        /* data for monitor function */

  if (call == 1) {

    GetData(yy, mxGetPr(mx_in[2]), N);

    if (quadr) {
      GetData(yQ, mxGetPr(mx_in[3]), Nq);
    }

    if (fsa) {
      tmp_yyS = mxGetPr(mx_in[5]);
      for (is=0; is<Ns; is++) {
        GetData(yyS[is], &tmp_yyS[is*N], N);
      }
    }

  }

  mexCallMATLAB(1,mx_out,8,mx_in,"idm_monitor");

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

void mxW_IDAMonitorB(int call, int idxB, double tB,
                     N_Vector yyB,
                     N_Vector yQB,
                     idmPbData bckPb)
{
  mxArray *mx_in[7], *mx_out[1];

  mx_in[0] = mxCreateDoubleScalar(call);            /* 0: first, 1: interm. 2: last */
  mx_in[1] = mxCreateDoubleScalar(idxB);            /* index of current problem */
  mx_in[2] = mxCreateDoubleScalar(tB);              /* current time */
  mx_in[3] = mxCreateDoubleMatrix(NB,1,mxREAL);     /* current yyB */
  if (quadrB) {
    mx_in[4] = mxCreateDoubleMatrix(NqB,1,mxREAL);  /* current quadratures */
  } else {
    mx_in[4] = mxCreateDoubleMatrix(0,0,mxREAL);
  }
  mx_in[5] = bckPb->MONfct;
  mx_in[6] = bckPb->MONdata;

  if (call == 1) {
    
    GetData(yyB, mxGetPr(mx_in[3]), NB);

    if (quadrB)
      GetData(yQB, mxGetPr(mx_in[4]), NqB);
  }

  mexCallMATLAB(1,mx_out,7,mx_in,"idm_monitorB");

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

static void UpdateUserData(mxArray *new_mtlb_data, idmPbData pb)
{
  mexUnlock();
  mxDestroyArray(pb->mtlb_data);
  pb->mtlb_data = mxDuplicateArray(new_mtlb_data);
  mexMakeArrayPersistent(pb->mtlb_data);
  mexLock();
}

static void UpdateMonitorData(mxArray *new_mtlb_data, idmPbData pb)
{
  mexUnlock();
  mxDestroyArray(pb->MONdata);
  pb->MONdata = mxDuplicateArray(new_mtlb_data);
  mexMakeArrayPersistent(pb->MONdata);
  mexLock();
}


