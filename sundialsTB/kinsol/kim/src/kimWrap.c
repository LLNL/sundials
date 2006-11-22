/*
 * -----------------------------------------------------------------
 * $Revision: 1.4 $
 * $Date: 2006-11-22 00:12:52 $
 * -----------------------------------------------------------------
 * Programmer: Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2005, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see sundials-x.y.z/src/kinsol/LICENSE.
 * -----------------------------------------------------------------
 * KINSOL wrapper functions.
 * -----------------------------------------------------------------
 */

#include "kim.h"
#include "nvm.h"

static void UpdateUserData(mxArray *mx_newdata);

/*
 * ---------------------------------------------------------------------------------
 * Redability replacements
 * ---------------------------------------------------------------------------------
 */

#define N  (kim_Kdata->N)
#define ls (kim_Kdata->ls)
#define pm (kim_Kdata->pm)

#define mx_data    (kim_Mdata->mx_data)

#define mx_SYSfct  (kim_Mdata->mx_SYSfct)
#define mx_JACfct  (kim_Mdata->mx_JACfct)
#define mx_PSETfct (kim_Mdata->mx_PSETfct)
#define mx_PSOLfct (kim_Mdata->mx_PSOLfct)
#define mx_GLOCfct (kim_Mdata->mx_GLOCfct)
#define mx_GCOMfct (kim_Mdata->mx_GCOMfct)

#define fig_handle (kim_Mdata->fig_handle)

/*
 * ---------------------------------------------------------------------------------
 * Error handler
 * ---------------------------------------------------------------------------------
 */

void mtlb_KINErrHandler(int error_code, 
                        const char *module, const char *function, 
                        char *msg, void *eh_data)
{
  char err_type[10];

  if (error_code == KIN_WARNING)
    sprintf(err_type,"WARNING");
  else
    sprintf(err_type,"ERROR");

  mexPrintf("\n[%s %s]  %s\n",module,err_type,function);
  mexPrintf("  %s\n\n",msg);

}


/*
 * ---------------------------------------------------------------------------------
 * Info handler
 * ---------------------------------------------------------------------------------
 */

void mtlb_KINInfoHandler(const char *module, const char *function, 
                         char *msg, void *ih_data)
{
  char my_msg[400];
  mxArray *mx_in[3];

  sprintf(my_msg,"[%s] %s\n  %s\n",module,function,msg);

  /* action=1 -> append */
  mx_in[0] = mxCreateScalarDouble(1);
  mx_in[1] = mxCreateScalarDouble((double)fig_handle);
  mx_in[2] = mxCreateString(my_msg);

  mexCallMATLAB(0,NULL,3,mx_in,"kim_info");

}


/*
 * ---------------------------------------------------------------------------------
 * Wrapper functions
 * ---------------------------------------------------------------------------------
 */

int mtlb_KINSys(N_Vector y, N_Vector fy, void *f_data )
{
  mxArray *mx_in[3], *mx_out[3];
  int ret;
  
  /* Inputs to the Matlab function */
  mx_in[0] = mxCreateDoubleMatrix(N,1,mxREAL); /* current y */
  mx_in[1] = mx_SYSfct;                        /* matlab function handle */ 
  mx_in[2] = mx_data;                          /* matlab user data */

  /* Call matlab wrapper */
  GetData(y, mxGetPr(mx_in[0]), N);
  mexCallMATLAB(3,mx_out,3,mx_in,"kim_sys");

  PutData(fy, mxGetPr(mx_out[0]), N);

  ret = (int)*mxGetPr(mx_out[1]);

  if (!mxIsEmpty(mx_out[2])) {
    UpdateUserData(mx_out[2]);
  }

  /* Free temporary space */
  mxDestroyArray(mx_in[0]);
  mxDestroyArray(mx_out[0]);
  mxDestroyArray(mx_out[1]);
  mxDestroyArray(mx_out[2]);

  return(ret);
}

int mtlb_KINDenseJac(int Neq,
                     N_Vector y, N_Vector fy, 
                     DlsMat J, void *jac_data,
                     N_Vector tmp1, N_Vector tmp2)
{
  double *J_data;
  mxArray *mx_in[4], *mx_out[3];
  int i, ret;

  /* Inputs to the Matlab function */
  mx_in[0] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current y */
  mx_in[1] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current fy */
  mx_in[2] = mx_JACfct;                         /* matlab function handle */
  mx_in[3] = mx_data;                           /* matlab user data */
  
  /* Call matlab wrapper */
  GetData(y, mxGetPr(mx_in[0]), N);
  GetData(fy, mxGetPr(mx_in[1]), N);
  mexCallMATLAB(3,mx_out,4,mx_in,"kim_djac");

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
  mxDestroyArray(mx_out[0]);
  mxDestroyArray(mx_out[1]);
  mxDestroyArray(mx_out[2]);

  return(ret);

}


int mtlb_KINBandJac(int Neq, int mupper, int mlower,
                    N_Vector y, N_Vector fy, 
                    DlsMat J, void *jac_data,
                    N_Vector tmp1, N_Vector tmp2)
{
  double *J_data;
  mxArray *mx_in[4], *mx_out[3];
  int eband, i, ret;

  /* Inputs to the Matlab function */
  mx_in[0] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current y */
  mx_in[1] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current fy */
  mx_in[2] = mx_JACfct;                         /* matlab function handle */
  mx_in[3] = mx_data;                           /* matlab user data */
  
  /* Call matlab wrapper */
  GetData(y, mxGetPr(mx_in[0]), N);
  GetData(fy, mxGetPr(mx_in[1]), N);

  mexCallMATLAB(3,mx_out,4,mx_in,"kim_bjac");

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
  mxDestroyArray(mx_out[0]);
  mxDestroyArray(mx_out[1]);
  mxDestroyArray(mx_out[2]);

  return(ret);

}

int mtlb_KINSpilsJac(N_Vector v, N_Vector Jv,
                     N_Vector y, booleantype *new_y, 
                     void *J_data)
{
  mxArray *mx_in[5], *mx_out[4];
  int ret;

  /* Inputs to the Matlab function */
  mx_in[0] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current y */
  mx_in[1] = mxCreateDoubleMatrix(N,1,mxREAL);  /* vector v */
  mx_in[2] = mxCreateLogicalScalar(*new_y);     /* */
  mx_in[3] = mx_JACfct;                         /* matlab function handle */
  mx_in[4] = mx_data;                           /* matlab user data */
 
  /* Call matlab wrapper */
  GetData(y, mxGetPr(mx_in[0]), N);
  GetData(v, mxGetPr(mx_in[1]), N);

  mexCallMATLAB(4,mx_out,5,mx_in,"kim_jtv");

  PutData(Jv, mxGetPr(mx_out[0]), N);
  *new_y = mxIsLogicalScalarTrue(mx_out[1]);
  ret = (int)*mxGetPr(mx_out[2]);

  if (!mxIsEmpty(mx_out[3])) {
    UpdateUserData(mx_out[3]);
  }

  /* Free temporary space */
  mxDestroyArray(mx_in[0]);
  mxDestroyArray(mx_in[1]);
  mxDestroyArray(mx_in[2]);
  mxDestroyArray(mx_out[0]);
  mxDestroyArray(mx_out[1]);
  mxDestroyArray(mx_out[2]);
  mxDestroyArray(mx_out[3]);

  return(ret);
}

int mtlb_KINSpilsPset(N_Vector y, N_Vector yscale,
                      N_Vector fy, N_Vector fscale,
                      void *P_data, N_Vector vtemp1,
                      N_Vector vtemp2)
{
  mxArray *mx_in[6], *mx_out[2];
  int ret;

  /* Inputs to the Matlab function */
  mx_in[0] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current y */
  mx_in[1] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current yscale */
  mx_in[2] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current fy */
  mx_in[3] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current fscale */
  mx_in[4] = mx_PSETfct;                        /* matlab function handle */
  mx_in[5] = mx_data;                           /* matlab user data */
  
  /* Call matlab wrapper */
  GetData(y,      mxGetPr(mx_in[0]), N);
  GetData(yscale, mxGetPr(mx_in[1]), N);
  GetData(fy,     mxGetPr(mx_in[2]), N);
  GetData(fscale, mxGetPr(mx_in[3]), N);

  mexCallMATLAB(2,mx_out,6,mx_in,"kim_pset");

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

int mtlb_KINSpilsPsol(N_Vector y, N_Vector yscale, 
                      N_Vector fy, N_Vector fscale, 
                      N_Vector v, void *P_data,
                      N_Vector vtemp)
{
  mxArray *mx_in[7], *mx_out[3];
  int ret;


  /* Inputs to the Matlab function */
  mx_in[0] = mxCreateDoubleMatrix(N,1,mxREAL); /* current y */
  mx_in[1] = mxCreateDoubleMatrix(N,1,mxREAL); /* current yscale */
  mx_in[2] = mxCreateDoubleMatrix(N,1,mxREAL); /* current fy */
  mx_in[3] = mxCreateDoubleMatrix(N,1,mxREAL); /* current fscale */
  mx_in[4] = mxCreateDoubleMatrix(N,1,mxREAL); /* right hand side */
  mx_in[5] = mx_PSOLfct;                       /* matlab function handle */
  mx_in[6] = mx_data;                          /* matlab user data */
  
  /* Call matlab wrapper */
  GetData(y,      mxGetPr(mx_in[0]), N);
  GetData(yscale, mxGetPr(mx_in[1]), N);
  GetData(fy,     mxGetPr(mx_in[2]), N);
  GetData(fscale, mxGetPr(mx_in[3]), N);
  GetData(v,      mxGetPr(mx_in[4]), N);

  mexCallMATLAB(3,mx_out,7,mx_in,"kim_psol");

  PutData(v, mxGetPr(mx_out[0]), N);
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

int mtlb_KINGloc(int Nlocal, N_Vector y, N_Vector gval, void *f_data)
{
  mxArray *mx_in[3], *mx_out[3];
  int ret;

  /* Inputs to the Matlab function */
  mx_in[0] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current y */
  mx_in[1] = mx_GLOCfct;                        /* matlab function handle */
  mx_in[2] = mx_data;                           /* matlab user data */
  
  /* Call matlab wrapper */
  GetData(y, mxGetPr(mx_in[0]), N);
  mexCallMATLAB(3,mx_out,3,mx_in,"kim_gloc");

  PutData(gval, mxGetPr(mx_out[0]), N);

  ret = (int)*mxGetPr(mx_out[1]);

  if (!mxIsEmpty(mx_out[2])) {
    UpdateUserData(mx_out[2]);
  }

  /* Free temporary space */
  mxDestroyArray(mx_in[0]);
  mxDestroyArray(mx_out[0]);
  mxDestroyArray(mx_out[1]);
  mxDestroyArray(mx_out[2]);

  return(ret);
}

int mtlb_KINGcom(int Nlocal, N_Vector y, void *f_data)
{
  mxArray *mx_in[5], *mx_out[2];
  int ret;

  /* Inputs to the Matlab function */
  mx_in[0] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current y */
  mx_in[1] = mx_GCOMfct;                        /* matlab function handle */
  mx_in[2] = mx_data;                           /* matlab user data */
  
  /* Call matlab wrapper */
  GetData(y, mxGetPr(mx_in[0]), N);
  mexCallMATLAB(2,mx_out,3,mx_in,"kim_gcom");

  ret = (int)*mxGetPr(mx_out[0]);

  if (!mxIsEmpty(mx_out[1])) {
    UpdateUserData(mx_out[1]);
  }

  /* Free temporary space */
  mxDestroyArray(mx_in[0]);
  mxDestroyArray(mx_out[0]);
  mxDestroyArray(mx_out[1]);

  return(ret);
}

/*
 * ---------------------------------------------------------------------------------
 * Private function to update the user data structure
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
