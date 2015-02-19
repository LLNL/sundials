/*
 * -----------------------------------------------------------------
 * $Revision: 1.10 $
 * $Date: 2012-03-07 21:49:18 $
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
 * Header file for the IDAS Matlab interface.
 * -----------------------------------------------------------------
 */

#ifndef _IDM_H
#define _IDM_H

#include <string.h>
#include <stdlib.h>
#include "mex.h"
#include <idas/idas.h>
#include <idas/idas_dense.h>
#include <idas/idas_band.h>
#include <idas/idas_spgmr.h>
#include <idas/idas_spbcgs.h>
#include <idas/idas_sptfqmr.h>
#include <idas/idas_bbdpre.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/*
 * ---------------------------------------------------------------------------------
 * Constants
 * ---------------------------------------------------------------------------------
 */

/* Tolerance types */

enum {IDA_SS, IDA_SV, IDA_EE};

/* Linear solver types */

enum {LS_DENSE, LS_BAND, LS_SPGMR, LS_SPBCG, LS_SPTFQMR};

/* Preconditioner modules */

enum {PM_NONE, PM_BBDPRE};

/*
 * ---------------------------------------------------------------------------------
 * Types for global data structures
 * ---------------------------------------------------------------------------------
 */


typedef struct idmPbData_ {


  long int n;        /* problem dimension */
  N_Vector YY;       /* solution vector */
  N_Vector YP;       /* derivative of solution vector */

  booleantype Quadr; /* integrate quadratures? */
  long int nq;       /* number of quadratures */
  N_Vector YQ;       /* quadratures vector */

  booleantype Fsa;   /* integrate sensitivities? */
  int ns;            /* number of sensitivities */
  N_Vector *YYS;     /* sensitivity vectors */
  N_Vector *YPS;     /* derivatives of sensitivity vectors */

  booleantype RootSet;  /* rootfinding active? */
  int ng;               /* number of root functions */

  booleantype TstopSet; /* tstop active? */

  int LS;            /* linear solver type */
  int PM;            /* preconditioner module */

  booleantype Mon;   /* monitoring? */

  /* Matlab functions and data associated with this problem */

  mxArray *RESfct;
  mxArray *QUADfct;

  mxArray *JACfct;

  mxArray *PSETfct;
  mxArray *PSOLfct;

  mxArray *GLOCfct;
  mxArray *GCOMfct;
    
  mxArray *Gfct;
    
  mxArray *SRESfct;

  mxArray *MONfct;
  mxArray *MONdata;

  /* Pointer to the global Matlab user data */

  mxArray *mtlb_data;

  /* Information for backward problems only */

  struct idmPbData_ *fwd;
  int index;               /* index of this problem */
  struct idmPbData_ *next; /* pointer to next problem in linked list */

} *idmPbData;


typedef struct idmInterfaceData_ {

  void *ida_mem;         /* IDAS solver memory */

  booleantype asa;       /* Perform ASA? */
  int Nd;                /* number of data points */
  int Nc;                /* number of check points */

  struct idmPbData_ *fwdPb;
  struct idmPbData_ *bckPb;

  int NbckPb;            /* Number of backward problems in the linked list bckPb */

  booleantype errMsg;    /* post error/warning messages? */

} *idmInterfaceData;



/*
 * ---------------------------------------------------------------------------------
 * Error handler function
 * ---------------------------------------------------------------------------------
 */

void idmErrHandler(int error_code, 
                   const char *module, const char *function, 
                   char *msg, void *eh_data);


/*
 * ---------------------------------------------------------------------------------
 * Wrapper functions
 * ---------------------------------------------------------------------------------
 */

int mxW_IDARes(realtype tt, N_Vector yy, N_Vector yp,
               N_Vector rr, void *user_data);

int mxW_IDAGfct(realtype t, N_Vector y, N_Vector yp,
                realtype *gout, void *user_data);

int mxW_IDAQuadFct(realtype tres, N_Vector yy, N_Vector yp,
                   N_Vector ypQ,
                   void *user_data);

int mxW_IDASensRes(int Nsens, realtype tres, 
                   N_Vector yy, N_Vector yp, N_Vector resval,
                   N_Vector *yyS, N_Vector *ypS, N_Vector *resvalS,
                   void *user_data,
                   N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

int mxW_IDADenseJac(long int Neq, 
                    realtype tt, realtype c_j, 
                    N_Vector yy, N_Vector yp, N_Vector rr,
                    DlsMat Jac, void *user_data, 
                    N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

int mxW_IDABandJac(long int Neq, long int mupper, long int mlower, 
                   realtype tt, realtype c_j, 
                   N_Vector yy, N_Vector yp, N_Vector rr, 
                   DlsMat Jac, void *user_data,
                   N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

int mxW_IDASpilsJac(realtype tt,
                    N_Vector yy, N_Vector yp, N_Vector rr,
                    N_Vector v, N_Vector Jv,
                    realtype c_j, void *user_data,
                    N_Vector tmp1, N_Vector tmp2);
int mxW_IDASpilsPset(realtype tt,
                     N_Vector yy, N_Vector yp, N_Vector rr,
                     realtype c_j, void *user_data,
                     N_Vector tmp1, N_Vector tmp2,
                     N_Vector tmp3);
int mxW_IDASpilsPsol(realtype tt,
                     N_Vector yy, N_Vector yp, N_Vector rr,
                     N_Vector rvec, N_Vector zvec,
                     realtype c_j, realtype delta, void *user_data,
                     N_Vector tmp);

int mxW_IDABBDgloc(long int Nlocal, realtype tt,
                   N_Vector yy, N_Vector yp, N_Vector gval,
                   void *user_data);
int mxW_IDABBDgcom(long int Nlocal, realtype tt,
                   N_Vector yy, N_Vector yp,
                   void *user_data);

void mxW_IDAMonitor(int call, double t, 
                    N_Vector yy,
                    N_Vector yQ, 
                    N_Vector *yyS,
                    idmPbData fwdPb);

int mxW_IDAResB(realtype tt, 
                N_Vector yy, N_Vector yp,
                N_Vector yyB, N_Vector ypB, N_Vector rrB,
                void *user_dataB);
int mxW_IDAResBS(realtype tt, 
                 N_Vector yy, N_Vector yp,
                 N_Vector *yyS, N_Vector *ypS,
                 N_Vector yyB, N_Vector ypB, N_Vector rrB,
                 void *user_dataB);

int mxW_IDAQuadFctB(realtype tt, 
                    N_Vector yy, N_Vector yp, 
                    N_Vector yyB, N_Vector ypB,
                    N_Vector ypQB,
                    void *user_dataB);
int mxW_IDAQuadFctBS(realtype t, 
                     N_Vector yy, N_Vector yp,
                     N_Vector *yyS, N_Vector *ypS,
                     N_Vector yyB, N_Vector ypB,
                     N_Vector ypQB,
                     void *user_dataB);

int mxW_IDADenseJacB(long int NeqB, 
                     realtype tt, realtype c_jB,
                     N_Vector yy, N_Vector yp,
                     N_Vector yyB, N_Vector ypB, N_Vector rrB,
                     DlsMat JacB, void *user_dataB, 
                     N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B);

int mxW_IDABandJacB(long int NeqB, long int mupperB, long int mlowerB, 
                    realtype tt, realtype c_jB, 
                    N_Vector yy, N_Vector yp,
                    N_Vector yyB, N_Vector ypB, N_Vector rrB,
                    DlsMat JacB, void *user_dataB,
                    N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B);

int mxW_IDASpilsJacB(realtype t,
                     N_Vector yy, N_Vector yp,
                     N_Vector yyB, N_Vector ypB, N_Vector rrB,
                     N_Vector vB, N_Vector JvB, 
                     realtype c_jB, void *user_dataB, 
                     N_Vector tmp1B, N_Vector tmp2B);
int mxW_IDASpilsPsetB(realtype tt, 
                      N_Vector yy, N_Vector yp,
                      N_Vector yyB, N_Vector ypB, N_Vector rrB, 
                      realtype c_jB, void *user_dataB,
                      N_Vector tmp1B, N_Vector tmp2B, 
                      N_Vector tmp3B);
int mxW_IDASpilsPsolB(realtype tt, 
                      N_Vector yy, N_Vector yp,
                      N_Vector yyB, N_Vector ypB, N_Vector rrB, 
                      N_Vector rvecB, N_Vector zvecB,
                      realtype c_jB, realtype deltaB,
                      void *user_dataB, N_Vector tmpB);

int mxW_IDABBDglocB(long int NlocalB, realtype tt,
                    N_Vector yy, N_Vector yp, 
                    N_Vector yyB, N_Vector ypB, N_Vector gvalB,
                    void *user_dataB);
int mxW_IDABBDgcomB(long int NlocalB, realtype tt,
                    N_Vector yy, N_Vector yp,
                    N_Vector yyB, N_Vector ypB,
                    void *user_dataB);

void mxW_IDAMonitorB(int call, int idxB, double tB, 
                     N_Vector yyB,
                     N_Vector yQB,
                     idmPbData bckPb);

/*
 * ---------------------------------------------------------------------------------
 * Option handling functions
 * ---------------------------------------------------------------------------------
 */

int get_IntgrOptions(const mxArray *options, idmPbData thisPb, booleantype fwd,
                     int *maxord,
                     long int *mxsteps,
                     int *itol, realtype *reltol, double *Sabstol, double **Vabstol,
                     double *hin, double *hmax,
                     double *tstop,
                     booleantype *suppress,
                     booleantype *errmsg,
                     double **id, double **cnstr,
                     booleantype *res_s);

int get_LinSolvOptions(const mxArray *options, idmPbData thisPb, booleantype fwd,
                       long int *mupper, long int *mlower,
                       long int *mudq, long int *mldq, double *dqrely,
                       int *gstype, int *maxl);

int get_QuadOptions(const mxArray *options, idmPbData thisPb, booleantype fwd,
                    long int Nq, booleantype *rhs_s,
                    booleantype *errconQ,
                    int *itolQ, double *reltolQ, double *SabstolQ, double **VabstolQ);

int get_FSAOptions(const mxArray *options, idmPbData thisPb,
                   int *ism,
                   char **pfield_name, int **plist, double **pbar,
                   int *dqtype, double *rho,
                   booleantype *errconS, int *itolS, double *reltolS, 
                   double **SabstolS, double **VabstolS);

#ifdef __cplusplus
}
#endif

#endif
