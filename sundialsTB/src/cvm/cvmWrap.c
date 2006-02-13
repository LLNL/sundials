/*
 * -----------------------------------------------------------------
 * $Revision: 1.3 $
 * $Date: 2006-02-13 23:01:29 $
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
  mx_in[3] = mx_mtlb_RHSfct;                   /* matlab function handle */ 
  mx_in[4] = mx_mtlb_data;                     /* matlab user data */

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
  mx_in[3] = mx_mtlb_QUADfct;                  /* matlab function handle */ 
  mx_in[4] = mx_mtlb_data;                     /* matlab user data */

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
  mx_in[2] = mx_mtlb_Gfct;                     /* matlab function handle */
  mx_in[3] = mx_mtlb_data;                     /* matlab user data */
  
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


int mtlb_CVodeDenseJac(long int N, DenseMat J, realtype t,
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
  mx_in[4] = mx_mtlb_JACfct;                    /* matlab function handle */
  mx_in[5] = mx_mtlb_data;                      /* matlab user data */
  
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

int mtlb_CVodeBandJac(long int N, long int mupper, long int mlower,
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
  mx_in[4] = mx_mtlb_JACfct;                    /* matlab function handle */
  mx_in[5] = mx_mtlb_data;                      /* matlab user data */
  
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
  mx_in[5] = mx_mtlb_JACfct;                    /* matlab function handle */
  mx_in[6] = mx_mtlb_data;                      /* matlab user data */
  
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
  mx_in[6] = mx_mtlb_PSETfct;                   /* matlab function handle */
  mx_in[7] = mx_mtlb_data;                      /* matlab user data */
  
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
  mx_in[5] = mx_mtlb_PSOLfct;                  /* matlab function handle */
  mx_in[6] = mx_mtlb_data;                     /* matlab user data */
  
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
  mx_in[3] = mx_mtlb_GLOCfct;                   /* matlab function handle */
  mx_in[4] = mx_mtlb_data;                      /* matlab user data */
  
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
  mx_in[3] = mx_mtlb_GCOMfct;                   /* matlab function handle */
  mx_in[4] = mx_mtlb_data;                      /* matlab user data */
  
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


int mtlb_CVodeSensRhs1(int Ns, realtype t,
                       N_Vector y, N_Vector yd,
                       int iS, N_Vector yS, N_Vector ySd,
                       void *fS_data,
                       N_Vector tmp1, N_Vector tmp2)
{
  double isd;
  mxArray *mx_in[7], *mx_out[3];
  int ret;

  isd = 1.0 + iS;

  /* Inputs to the Matlab function */
  mx_in[0] = mxCreateScalarDouble(isd);        /* iS */
  mx_in[1] = mxCreateScalarDouble(t);          /* current t */
  mx_in[2] = mxCreateDoubleMatrix(N,1,mxREAL); /* current y */
  mx_in[3] = mxCreateDoubleMatrix(N,1,mxREAL); /* current yd */
  mx_in[4] = mxCreateDoubleMatrix(N,1,mxREAL); /* current yS */
  mx_in[5] = mx_mtlb_SRHSfct;                  /* matlab function handle */      
  mx_in[6] = mx_mtlb_data;                     /* matlab user data */
  
  /* Call matlab wrapper */
  GetData(y, mxGetPr(mx_in[2]), N);
  GetData(yd, mxGetPr(mx_in[3]), N);
  GetData(yS, mxGetPr(mx_in[4]), N);

  mexCallMATLAB(3,mx_out,7,mx_in,"cvm_rhsS1");

  PutData(ySd, mxGetPr(mx_out[0]), N);

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


int mtlb_CVodeSensRhs(int Ns, realtype t,
                      N_Vector y, N_Vector yd,
                      N_Vector *yS, N_Vector *ySd,
                      void *fS_data,
                      N_Vector tmp1, N_Vector tmp2)
{
  mxArray *mx_in[6], *mx_out[3];
  int is, ret;
  double *tmp;

  /* Inputs to the Matlab function */
  mx_in[0] = mxCreateScalarDouble(t);          /* current t */
  mx_in[1] = mxCreateDoubleMatrix(N,1,mxREAL); /* current y */
  mx_in[2] = mxCreateDoubleMatrix(N,1,mxREAL); /* current yd */
  mx_in[3] = mxCreateDoubleMatrix(N*Ns,1,mxREAL); /* current yS */
  mx_in[4] = mx_mtlb_SRHSfct;                  /* matlab function handle */      
  mx_in[5] = mx_mtlb_data;                     /* matlab user data */
  
  /* Call matlab wrapper */
  GetData(y, mxGetPr(mx_in[1]), N);
  GetData(yd, mxGetPr(mx_in[2]), N);
  tmp = mxGetPr(mx_in[3]);
  for (is=0; is<Ns; is++)
    GetData(yS[is], &tmp[is*Ns], N);

  mexCallMATLAB(3,mx_out,6,mx_in,"cvm_rhsS");
  
  tmp = mxGetPr(mx_out[0]);
  PutData(ySd[is], &tmp[is*Ns], N);

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
  mx_in[4] = mx_mtlb_RHSfctB;                   /* matlab function handle */ 
  mx_in[5] = mx_mtlb_data;                      /* matlab user data */
  
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
  mx_in[4] = mx_mtlb_QUADfctB;                  /* matlab function handle */ 
  mx_in[5] = mx_mtlb_data;                      /* matlab user data */

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


int mtlb_CVodeDenseJacB(long int nB, DenseMat JB, realtype t,
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
  mx_in[5] = mx_mtlb_JACfctB;                   /* matlab function handle */
  mx_in[6] = mx_mtlb_data;                      /* matlab user data */
  
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


int mtlb_CVodeBandJacB(long int nB, long int mupperB,
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
  mx_in[5] = mx_mtlb_JACfctB;                   /* matlab function handle */
  mx_in[6] = mx_mtlb_data;                      /* matlab user data */
  
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
  mx_in[6] = mx_mtlb_JACfctB;                   /* matlab function handle */
  mx_in[7] = mx_mtlb_data;                      /* matlab user data */
  
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
  mx_in[7] = mx_mtlb_PSETfctB;                  /* matlab function handle */
  mx_in[8] = mx_mtlb_data;                      /* matlab user data */
  
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
  mx_in[6] = mx_mtlb_PSOLfctB;                  /* matlab function handle */
  mx_in[7] = mx_mtlb_data;                      /* matlab user data */
  
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
  mx_in[4] = mx_mtlb_GLOCfctB;                  /* matlab function handle */
  mx_in[5] = mx_mtlb_data;                      /* matlab user data */
  
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
  mx_in[4] = mx_mtlb_GCOMfctB;                  /* matlab function handle */
  mx_in[5] = mx_mtlb_data;                      /* matlab user data */
  
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
  mxArray *mx_in[7];
  double *tmp;
  int is;

  /* Inputs to the Matlab function */
  mx_in[0] = mxCreateScalarDouble(call);            /* 0: first, 1: interm. 2: last */
  mx_in[1] = mxCreateScalarDouble(t);               /* current t */
  mx_in[2] = mxCreateDoubleMatrix(N,1,mxREAL);      /* current y */  
  if (Nq > 0)
    mx_in[3] = mxCreateDoubleMatrix(Nq,1,mxREAL);   /* current yQ */
  else
    mx_in[3] = mxCreateDoubleMatrix(0,0,mxREAL);
  if (Ns > 0)
    mx_in[4] = mxCreateDoubleMatrix(N*Ns,1,mxREAL); /* current yS */
  else
    mx_in[4] = mxCreateDoubleMatrix(0,0,mxREAL);
  mx_in[5] = mx_mtlb_MONfct;                        /* matlab function handle */
  mx_in[6] = mx_mtlb_MONdata;                       /* matlab user data */  

  if (call !=2) {
    
    GetData(y, mxGetPr(mx_in[2]), N);

    if (Nq > 0)
      GetData(yQ, mxGetPr(mx_in[3]), Nq);

    if (Ns > 0) {
      tmp = mxGetPr(mx_in[4]);
      for (is=0; is<Ns; is++)
        GetData(yS[is], &tmp[is*Ns], N);
    }
  }

  /* Call matlab function */
  mexCallMATLAB(0,NULL,7,mx_in,"cvm_monitor");

  mxDestroyArray(mx_in[0]);
  mxDestroyArray(mx_in[1]);
  mxDestroyArray(mx_in[2]);
  mxDestroyArray(mx_in[3]);
  mxDestroyArray(mx_in[4]);
  
}

void mtlb_CVodeMonitorB(int call, double tB, N_Vector yB, N_Vector yQB)
{
  mxArray *mx_in[7];
  double *tmp;

  /* Inputs to the Matlab function */
  mx_in[0] = mxCreateScalarDouble(call);            /* 0: first, 1: interm. 2: last */
  mx_in[1] = mxCreateScalarDouble(tB);              /* current t */
  mx_in[2] = mxCreateDoubleMatrix(NB,1,mxREAL);     /* current yB */
  if (NqB > 0)
    mx_in[3] = mxCreateDoubleMatrix(NqB,1,mxREAL);  /* current yQB */
  else
    mx_in[3] = mxCreateDoubleMatrix(0,0,mxREAL);
  mx_in[4] = mxCreateDoubleMatrix(0,0,mxREAL);      /* no yS */
  mx_in[5] = mx_mtlb_MONfctB;                       /* matlab function handle */
  mx_in[6] = mx_mtlb_MONdataB;                      /* matlab user data */  

  if (call !=2) {
    
    GetData(yB, mxGetPr(mx_in[2]), NB);

    if (NqB > 0)
      GetData(yQB, mxGetPr(mx_in[3]), NqB);

  }

  /* Call matlab function */
  mexCallMATLAB(0,NULL,7,mx_in,"cvm_monitor");

  mxDestroyArray(mx_in[0]);
  mxDestroyArray(mx_in[1]);
  mxDestroyArray(mx_in[2]);
  mxDestroyArray(mx_in[3]);
  mxDestroyArray(mx_in[4]);

}



/*
 * ---------------------------------------------------------------------------------
 * Private function to update the user data structure
 * ---------------------------------------------------------------------------------
 */

static void UpdateUserData(mxArray *mx_data)
{
  mexUnlock();
  mxDestroyArray(mx_mtlb_data);
  mx_mtlb_data = mxDuplicateArray(mx_data);
  mexMakeArrayPersistent(mx_mtlb_data);
  mexLock();
}

