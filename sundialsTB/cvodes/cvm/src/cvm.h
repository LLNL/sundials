/*
 * -----------------------------------------------------------------
 * $Revision: 1.14 $
 * $Date: 2012-03-07 21:44:21 $
 * -----------------------------------------------------------------
 * Programmer: Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * LLNS Copyright Start
 * Copyright (c) 2014, Lawrence Livermore National Security
 * This work was performed under the auspices of the U.S. Department 
 * of Energy by Lawrence Livermore National Laboratory in part under 
 * Contract W-7405-Eng-48 and in part under Contract DE-AC52-07NA27344.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * LLNS Copyright End
 * -----------------------------------------------------------------
 * Header file for the CVODES Matlab interface.
 * -----------------------------------------------------------------
 */

#ifndef _CVM_H
#define _CVM_H

#include <string.h>
#include <stdlib.h>
#include "mex.h"
#include <cvodes/cvodes.h>
#include <cvodes/cvodes_dense.h>
#include <cvodes/cvodes_diag.h>
#include <cvodes/cvodes_band.h>
#include <cvodes/cvodes_spgmr.h>
#include <cvodes/cvodes_spbcgs.h>
#include <cvodes/cvodes_sptfqmr.h>
#include <cvodes/cvodes_bandpre.h>
#include <cvodes/cvodes_bbdpre.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/*
 * ---------------------------------------------------------------------------------
 * Constants
 * ---------------------------------------------------------------------------------
 */

/* Tolerance types */

enum {CV_SS, CV_SV, CV_EE};

/* Linear solver types */

enum {LS_NONE, LS_DENSE, LS_DIAG, LS_BAND, LS_SPGMR, LS_SPBCG, LS_SPTFQMR};

/* Preconditioner modules */

enum {PM_NONE, PM_BANDPRE, PM_BBDPRE};

/*
 * ---------------------------------------------------------------------------------
 * Types for global data structures
 * ---------------------------------------------------------------------------------
 */


typedef struct cvmPbData_ {


  long int n;        /* problem dimension */
  N_Vector Y;        /* solution vector */

  booleantype Quadr; /* integrate quadratures? */
  long int nq;       /* number of quadratures */
  N_Vector YQ;       /* quadratures vector */

  booleantype Fsa;   /* integrate sensitivities? */
  int ns;            /* number of sensitivities */
  N_Vector *YS;      /* sensitivity vectors */

  booleantype RootSet;  /* rootfinding active? */
  int ng;               /* number of root functions */

  booleantype TstopSet; /* tstop active? */

  int LS;            /* linear solver type */
  int PM;            /* preconditioner module */

  booleantype Mon;   /* monitoring? */

  /* Matlab functions and data associated with this problem */

  mxArray *RHSfct;
  mxArray *QUADfct;

  mxArray *JACfct;

  mxArray *PSETfct;
  mxArray *PSOLfct;

  mxArray *GLOCfct;
  mxArray *GCOMfct;
    
  mxArray *Gfct;
    
  mxArray *SRHSfct;

  mxArray *MONfct;
  mxArray *MONdata;

  /* Pointer to the global Matlab user data */

  mxArray *mtlb_data;

  /* Information for backward problems only */

  struct cvmPbData_ *fwd;
  int index;               /* index of this problem */
  struct cvmPbData_ *next; /* pointer to next problem in linked list */

} *cvmPbData;


typedef struct cvmInterfaceData_ {

  void *cvode_mem;       /* CVODES solver memory */

  booleantype asa;       /* Perform ASA? */
  int Nd;                /* number of data points */
  int Nc;                /* number of check points */

  struct cvmPbData_ *fwdPb;
  struct cvmPbData_ *bckPb;

  int NbckPb;            /* Number of backward problems in the linked list bckPb */

  booleantype errMsg;    /* post error/warning messages? */

} *cvmInterfaceData;


/*
 * ---------------------------------------------------------------------------------
 * Error handler function
 * ---------------------------------------------------------------------------------
 */

void cvmErrHandler(int error_code, 
                   const char *module, const char *function, 
                   char *msg, void *eh_data);

/*
 * ---------------------------------------------------------------------------------
 * Wrapper functions
 * ---------------------------------------------------------------------------------
 */

int mxW_CVodeRhs(realtype t, N_Vector y, N_Vector yd, void *user_data);

int mxW_CVodeGfct(realtype t, N_Vector y, double *g, void *user_data);

int mxW_CVodeQUADfct(realtype t, N_Vector y, N_Vector yQd, void *user_data);


int mxW_CVodeSensRhs1(int Ns, realtype t,
                      N_Vector y, N_Vector ydot,
                      int iS, N_Vector yS, N_Vector ySdot,
                      void *user_data,
                      N_Vector tmp1, N_Vector tmp2);
int mxW_CVodeSensRhs(int Ns, realtype t,
                     N_Vector y, N_Vector ydot,
                     N_Vector *yS, N_Vector *ySdot,
                     void *user_data,
                     N_Vector tmp1, N_Vector tmp2);


int mxW_CVodeDenseJac(long int N, realtype t,
                      N_Vector y, N_Vector fy, 
                      DlsMat J, void *user_data,
                      N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
int mxW_CVodeBandJac(long int N, long int mupper, long int mlower, realtype t,
                     N_Vector y, N_Vector fy, 
                     DlsMat J, void *user_data,
                     N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
int mxW_CVodeSpilsJac(N_Vector v, N_Vector Jv, realtype t,
                      N_Vector y, N_Vector fy,
                      void *user_data, N_Vector tmp);
int mxW_CVodeSpilsPset(realtype t, N_Vector y, N_Vector fy,
                       booleantype jok, booleantype *jcurPtr,
                       realtype gamma, void *user_data,
                       N_Vector tmp1, N_Vector tmp2,
                       N_Vector tmp3);
int mxW_CVodeSpilsPsol(realtype t, N_Vector y, N_Vector fy,
                       N_Vector r, N_Vector z,
                       realtype gamma, realtype delta,
                       int lr, void *user_data, N_Vector tmp);

  
int mxW_CVodeBBDgloc(long int Nlocal, realtype t, N_Vector y, N_Vector g, void *user_data);
int mxW_CVodeBBDgcom(long int Nlocal, realtype t, N_Vector y, void *user_data);

void mxW_CVodeMonitor(int call, double t, 
                      N_Vector y, N_Vector yQ, N_Vector *yS,
                      cvmPbData fwdPb);


int mxW_CVodeRhsB(realtype t, N_Vector y,
                  N_Vector yB, N_Vector yBdot, void *user_dataB);
int mxW_CVodeRhsBS(realtype t, N_Vector y,  N_Vector *yS,
                   N_Vector yB, N_Vector yBd, void *user_dataB);

int mxW_CVodeQUADfctB(realtype t, N_Vector y, 
                      N_Vector yB, N_Vector qBdot, void *user_dataB);
int mxW_CVodeQUADfctBS(realtype t, N_Vector y,  N_Vector *yS,
                       N_Vector yB, N_Vector yQBd, void *user_dataB);

int mxW_CVodeDenseJacB(long int nB, realtype t,
                       N_Vector y, N_Vector yB, N_Vector fyB,
                       DlsMat JB, void *user_dataB, 
                       N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B);
int mxW_CVodeBandJacB(long int nB, long int mupperB, long int mlowerB, realtype t,
                      N_Vector y, N_Vector yB, N_Vector fyB,
                      DlsMat JB, void *user_dataB, 
                      N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B);
int mxW_CVodeSpilsJacB(N_Vector vB, N_Vector JvB, realtype t,
                       N_Vector y, N_Vector yB, N_Vector fyB,
                       void *user_dataB, N_Vector tmpB);
int mxW_CVodeSpilsPsetB(realtype t, N_Vector y,
                        N_Vector yB, N_Vector fyB,
                        booleantype jokB,
                        booleantype *jcurPtrB, realtype gammaB,
                        void *user_dataB,
                        N_Vector tmp1B, N_Vector tmp2B,
                        N_Vector tmp3B);
int mxW_CVodeSpilsPsolB(realtype t, N_Vector y,
                        N_Vector yB, N_Vector fyB,
                        N_Vector rB, N_Vector zB,
                        realtype gammaB, realtype deltaB,
                        int lrB, void *user_dataB, N_Vector tmpB);
  
int mxW_CVodeBBDglocB(long int NlocalB, realtype t, N_Vector y, 
                      N_Vector yB, N_Vector gB, void *user_dataB);

int mxW_CVodeBBDgcomB(long int NlocalB, realtype t, N_Vector y, 
                      N_Vector yB, void *user_dataB);


void mxW_CVodeMonitorB(int call, int idxB, double tB,
                       N_Vector yB, N_Vector yQB,
                       cvmPbData bckPb);

/*
 * ---------------------------------------------------------------------------------
 * Option handling functions
 * ---------------------------------------------------------------------------------
 */

int get_IntgrOptions(const mxArray *options, cvmPbData thisPb, booleantype fwd, int lmm,
                     int *maxord, booleantype *sld, booleantype *errmsg,
                     long int *mxsteps,
                     int *itol, realtype *reltol, double *Sabstol, double **Vabstol,
                     double *hin, double *hmax, double *hmin, 
                     double *tstop, booleantype *rhs_s);

int get_LinSolvOptions(const mxArray *options, cvmPbData thisPb, booleantype fwd,
                       long int *mupper, long int *mlower,
                       long int *mudq, long int *mldq, double *dqrely,
                       int *ptype, int *gstype, int *maxl);

int get_QuadOptions(const mxArray *options, cvmPbData thisPb, booleantype fwd,
                    long int Nq, booleantype *rhs_s,
                    booleantype *errconQ,
                    int *itolQ, double *reltolQ, double *SabstolQ, double **VabstolQ);

int get_FSAOptions(const mxArray *options, cvmPbData thisPb,
                   int *ism,
                   char **pfield_name, int **plist, double **pbar,
                   int *dqtype, double *rho,
                   booleantype *errconS, int *itolS, double *reltolS, 
                   double **SabstolS, double **VabstolS);

#ifdef __cplusplus
}
#endif

#endif
